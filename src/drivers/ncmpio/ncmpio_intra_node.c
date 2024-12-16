/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * This file contains the implementation of intra-node aggregation feature,
 * which is designed for the I/O patterns that contain many noncontiguous
 * requests interleaved among processes, and spreading across a wide range of
 * file space. It is particularly useful when the number of MPI processes
 * allocated to a compute node is large.
 *
 * This feature is enabled by setting the PnetCDF hint 'nc_num_aggrs_per_node'
 * to a positive integral value indicating the desired number of processes per
 * compute node to be selected as the intra-node I/O aggregators. Each process
 * is assigned a unique aggregator. The non-aggregators send their requests to
 * the assigned aggregators, and then the aggregators make MPI-IO requests to
 * the file.
 *
 * Such strategy can effectively reduce communication congestion due to many
 * pending asynchronous messages produced in the collective write inside of
 * MPI-IO.
 *
 * The concept of intra-node request aggregation is based on the paper:
 * Q. Kang, S. Lee, K. Hou, R. Ross, A. Agrawal, A. Choudhary, and W. Liao.
 * Improving MPI Collective I/O for High Volume Non-Contiguous Requests With
 * Intra-Node Aggregation. IEEE Transactions on Parallel and Distributed
 * Systems (TPDS), 31(11):2682-2695, November 2020.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

int debug = 0;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* strcmp() strdup() */
#include <assert.h>
#include <errno.h>
#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

#ifdef HAVE_MPI_LARGE_COUNT
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    MPI_Count aint; \
    MPI_Count cint; \
    MPI_Count d0 = (x) - offsets; \
    MPI_Count d1 = (y) - offsets; \
    if (d0 != d1) { \
        cint = *(x) ; *(x) = *(y) ; *(y) = cint ; \
        cint = lengths[d0] ; lengths[d0] = lengths[d1] ; lengths[d1] = cint ; \
        aint = bufAddr[d0] ; bufAddr[d0] = bufAddr[d1] ; bufAddr[d1] = aint ; \
    } \
}
#else
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    int int4; \
    MPI_Offset aint; \
    MPI_Offset d0 = (x) - offsets; \
    MPI_Offset d1 = (y) - offsets; \
    if (d0 != d1) { \
        aint = *(x) ; *(x) = *(y) ; *(y) = aint ; \
        int4 = lengths[d0] ; lengths[d0] = lengths[d1] ; lengths[d1] = int4 ; \
        aint = bufAddr[d0] ; bufAddr[d0] = bufAddr[d1] ; bufAddr[d1] = aint ; \
    } \
}
#endif

#define MEDIAN(a,b,c) ((*(a) < *(b)) ? \
                      ((*(b) < *(c)) ? (b) : ((*(a) < *(c)) ? (c) : (a))) : \
                      ((*(b) > *(c)) ? (b) : ((*(a) < *(c)) ? (a) : (c))))

/*----< qsort_off_len_buf() >------------------------------------------------*/
/* Sort three arrays of offsets, lengths, and buffer addresses based on the
 * increasing order of offsets. This code is based on the qsort routine from
 * Bentley & McIlroy's "Engineering a Sort Function".
 */
static void
qsort_off_len_buf(MPI_Aint    num,
#ifdef HAVE_MPI_LARGE_COUNT
                  MPI_Count  *offsets,
                  MPI_Count  *lengths,
#else
                  MPI_Offset *offsets,
                  int        *lengths,
#endif
                  MPI_Aint   *bufAddr)
{
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *pa, *pb, *pc, *pd, *pl, *pm, *pn, cmp_result, swap_cnt, i, r;
#else
    MPI_Offset *pa, *pb, *pc, *pd, *pl, *pm, *pn, cmp_result, swap_cnt, i, r;
#endif

    while (1) {
        swap_cnt = 0;
        pm = offsets + (num / 2);
        if (num > 7) {
            pl = offsets;
            pn = offsets + (num - 1);
            if (num > 40) {
                size_t d = (num / 8);
                pl = MEDIAN(pl, pl + d, pl + 2 * d);
                pm = MEDIAN(pm - d, pm, pm + d);
                pn = MEDIAN(pn - 2 * d, pn - d, pn);
            }
            pm = MEDIAN(pl, pm, pn);
        }
        SWAP(offsets, lengths, bufAddr, offsets, pm);
        pa = pb = offsets;

        pc = pd = offsets + (num - 1);
        for (;;) {
            while (pb <= pc && (cmp_result = (*pb - *offsets)) <= 0) {
                if (cmp_result == 0) {
                    swap_cnt = 1;
                    SWAP(offsets, lengths, bufAddr, pa, pb);
                    pa++;
                }
                pb++;
            }
            while (pb <= pc && (cmp_result = (*pc - *offsets)) >= 0) {
                if (cmp_result == 0) {
                    swap_cnt = 1;
                    SWAP(offsets, lengths, bufAddr, pc, pd);
                    pd--;
                }
                pc--;
            }
            if (pb > pc)
                break;
            SWAP(offsets, lengths, bufAddr, pb, pc);
            swap_cnt = 1;
            pb++;
            pc--;
        }
        if (swap_cnt == 0) {  /* Switch to insertion sort */
            for (pm = offsets; pm < offsets + num; pm++)
                for (pl = pm; pl > offsets && (*(pl-1) > *pl); pl--)
                    SWAP(offsets, lengths, bufAddr, pl, pl-1);
            return;
        }

        pn = offsets + num;
        r = MIN(pa - offsets, pb - pa);
        for (i=0; i<r; i++) SWAP(offsets, lengths, bufAddr, offsets+i, pb-r+i)

        r = MIN(pd - pc, pn - pd - 1);
        for (i=0; i<r; i++) SWAP(offsets, lengths, bufAddr, pb+i, pn-r+i)

        if ((r = pb - pa) > 1)
            qsort_off_len_buf(r, offsets, lengths, bufAddr);
        if ((r = pd - pc) > 1) {
            /* Iterate rather than recursively call self to save stack space */
            lengths = lengths + (num - r);
            bufAddr = bufAddr + (num - r);
            offsets = pn - r;
            num = r;
        }
        else
            break;
    }
}

/* Heapify(a, i, heapsize); Algorithm from Cormen et al. pg. 143 modified for a
 * heap with smallest element at root. The recursion has been removed so that
 * there are no function calls. Function calls are too expensive.
 *
 * Requirement: all individual offsets lists must be already sorted !!!
 */
static
void heap_merge(int              nprocs,
                const MPI_Aint  *count,    /* [nprocs] */
                MPI_Aint         nelems,
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count       *offsets,  /* [nelems] */
                MPI_Count       *blklens,  /* [nelems] */
#else
                MPI_Offset      *offsets,  /* [nelems] */
                int             *blklens,  /* [nelems] */
#endif
                MPI_Aint        *bufAddr)  /* [nelems] */
{
    typedef struct {
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count  *off_list;
        MPI_Count  *len_list;
#else
        MPI_Offset *off_list;
        int        *len_list;
#endif
        MPI_Aint  *addr_list;
        MPI_Aint  count;
    } heap_struct;

    heap_struct *a, tmp;
    int i, j, heapsize, l, r, k, smallest;
    size_t sum;

    /* This heap_merge is not in-plance, taking too much memory footprint */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *srt_off = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
    MPI_Count *srt_len = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
#else
    MPI_Aint *srt_off = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nelems);
    int      *srt_len = (int*)      NCI_Malloc(sizeof(int)      * nelems);
#endif
    MPI_Aint *srt_addr = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nelems);

    a = (heap_struct *) NCI_Calloc(nprocs, sizeof(heap_struct));

    /* there are nprocs number of lists to be merged */
    j = 0;
    sum = 0;
    for (i = 0; i < nprocs; i++) {
        if (count[i]) {
            /* each of a[j].off_list is already sorted */
            a[j].off_list = offsets + sum;
            a[j].len_list = blklens + sum;
            a[j].addr_list = bufAddr + sum;
            sum += count[i];
            a[j].count = count[i];
            j++;
        }
    }

#define SWAP_HEAP(x, y, tmp) { tmp = x ; x = y ; y = tmp ; }

    heapsize = nprocs;

    /* Build a heap out of the first element from each list, with the smallest
     * element of the heap at the root. The first for loop is to find and move
     * the smallest a[*].off_list[0] to a[0].
     */
    for (i = heapsize / 2 - 1; i >= 0; i--) {
        k = i;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if (l < heapsize && a[l].off_list[0] < a[k].off_list[0])
                smallest = l;
            else
                smallest = k;

            if (r < heapsize && a[r].off_list[0] < a[smallest].off_list[0])
                smallest = r;

            if (smallest != k) {
                SWAP_HEAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }

    /* The heap keeps the smallest element in its first element, i.e.
     * a[0].off_list[0].
     */
    j = 0;
    for (i = 0; i < nelems; i++) {
        /* extract smallest element from heap, i.e. the root */
        srt_off[i] = a[0].off_list[0];
        srt_len[i] = a[0].len_list[0];
        srt_addr[i] = a[0].addr_list[0];
        a[0].count--;

        if (!a[0].count) {
            a[0] = a[heapsize - 1];
            heapsize--;
        } else {
            a[0].off_list++;
            a[0].len_list++;
            a[0].addr_list++;
        }

        /* Heapify(a, 0, heapsize); */
        k = 0;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if (l < heapsize && a[l].off_list[0] < a[k].off_list[0])
                smallest = l;
            else
                smallest = k;

            if (r < heapsize && a[r].off_list[0] < a[smallest].off_list[0])
                smallest = r;

            if (smallest != k) {
                SWAP_HEAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }

#ifdef HAVE_MPI_LARGE_COUNT
    memcpy(offsets, srt_off, sizeof(MPI_Count) * nelems);
    memcpy(blklens, srt_len, sizeof(MPI_Count) * nelems);
#else
    memcpy(offsets, srt_off, sizeof(MPI_Offset) * nelems);
    memcpy(blklens, srt_len, sizeof(int)        * nelems);
#endif
    memcpy(bufAddr, srt_addr, sizeof(MPI_Aint) * nelems);

    NCI_Free(a);
    NCI_Free(srt_addr);
    NCI_Free(srt_len);
    NCI_Free(srt_off);
}


/*----< ncmpio_intra_node_aggr_init() >--------------------------------------*/
/* When intra-node write aggregation is enabled, processes on the same node
 * will be divided into groups. The number of groups is the number of
 * aggregators on that node. The rank IDs of each group must be established.
 *
 * 1. Find the affinity of each MPI process to the compute node. This should
 *    have been done from a call to ncmpii_construct_node_list() during
 *    ncmpio_create() and ncmpio_open().
 * 2. Determine whether self process is an intra-node aggregator.
 * 3. For an aggregator, find the number of non-aggregators assigned to it and
 *    construct rank IDs of assigned non-aggregators.
 * 4. For a non-aggregator, find the rank ID of its assigned aggregator.
 *
 * This subroutine should be called only once per file at ncmpio_create() or
 * ncmpio_open(). It sets the following variables.
 *   ncp->my_aggr          rank ID of my aggregator
 *   ncp->num_nonaggrs     number of non-aggregators assigned to me
 *   ncp->nonaggr_ranks[]  rank IDs of assigned non-aggregators
 */
int
ncmpio_intra_node_aggr_init(NC *ncp)
{
    int i, naggrs_my_node, num_nonaggrs;
    int my_rank_index, *ranks_my_node, my_node_id, nprocs_my_node;

    /* initialize parameters of local-node aggregation */
    ncp->my_aggr = -1;         /* rank ID of my aggregator */
    ncp->num_nonaggrs = 0;     /* number of non-aggregators assigned */
    ncp->nonaggr_ranks = NULL; /* ranks of assigned non-aggregators */

    if (ncp->num_aggrs_per_node == 0 || ncp->num_aggrs_per_node == ncp->nprocs)
        /* disable intra-node aggregation */
        return NC_NOERR;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
    ncp->aggr_time[0] = ncp->aggr_time[1] = 0.0;
    ncp->aggr_time[2] = ncp->aggr_time[3] = 0.0;
    ncp->aggr_time[4] = ncp->aggr_time[5] = 0.0;
#endif

    /* allocate space for storing the rank IDs of non-aggregators assigned to
     * this rank. Note ncp->nonaggr_ranks[] will be freed when closing the
     * file, if allocated.
     */
    num_nonaggrs = ncp->nprocs / ncp->num_aggrs_per_node + 1;
    ncp->nonaggr_ranks = (int*) NCI_Malloc(sizeof(int) * num_nonaggrs);

    /* ncp->node_ids[] has been established in ncmpii_construct_node_list() */
    /* my_node_id is this rank's node ID */
    my_node_id = ncp->node_ids[ncp->rank];

    /* nprocs_my_node: the number of processes in my nodes
     * ranks_my_node[]: rank IDs of all processes in my node.
     * my_rank_index points to ranks_my_node[] where
     * ranks_my_node[my_rank_index] == ncp->rank
     */
    ranks_my_node = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
    my_rank_index = -1;
    nprocs_my_node = 0;
    for (i=0; i<ncp->nprocs; i++) {
        if (ncp->node_ids[i] == my_node_id) {
            if (i == ncp->rank)
                my_rank_index = nprocs_my_node;
            ranks_my_node[nprocs_my_node] = i;
            nprocs_my_node++;
        }
    }
    assert(my_rank_index >= 0);

    /* Now, ranks_my_node[my_rank_index] == ncp->rank */

    /* make sure number of aggregators in my node <= nprocs_my_node */
    naggrs_my_node = MIN(ncp->num_aggrs_per_node, nprocs_my_node);

    /* calculate the number of non-aggregators assigned to an aggregator.
     * Note num_nonaggrs includes self.
     */
    num_nonaggrs = nprocs_my_node / naggrs_my_node;
    if (nprocs_my_node % naggrs_my_node) num_nonaggrs++;

    int do_io = 0;
    if (num_nonaggrs == 1) {
        /* disable aggregation if the number of non-aggregators assigned to
         * this aggregator is 1. Note num_nonaggrs includes self. It is
         * possible for aggregation enabled or disabled on different nodes and
         * even different aggregation groups on the same node.
         *
         * Use whether ncp->my_aggr < 0 to tell if aggregation is disabled or
         * enabled.
         */
        ncp->my_aggr = -1;
        do_io = 1;
    }
    else {
        /* find the rank ID of aggregator assigned to this rank */
        ncp->my_aggr = ranks_my_node[my_rank_index - my_rank_index % num_nonaggrs];

        if (ncp->my_aggr == ncp->rank) { /* this rank is an aggregator */
            /* Set the number of non-aggregators assigned to this rank. For the
             * last group, make sure it does not go beyond nprocs_my_node.
             */
            do_io = 1;
            ncp->num_nonaggrs = MIN(num_nonaggrs, nprocs_my_node - my_rank_index);
            if (ncp->num_nonaggrs == 1)
                /* disable aggregation, as this aggregation group contains only
                 * self rank
                 */
                ncp->my_aggr = -1;
            else
                /* copy the rank IDs over to ncp->nonaggr_ranks[] */
                memcpy(ncp->nonaggr_ranks,
                       ranks_my_node + my_rank_index,
                       sizeof(int) * num_nonaggrs);
        }
    }
    NCI_Free(ranks_my_node);

    if (ncp->my_aggr < 0) {
        /* free ncp->nonaggr_ranks if aggregation is not enabled */
        NCI_Free(ncp->nonaggr_ranks);
        ncp->nonaggr_ranks = NULL;
    }

    /* Note when ncp->my_aggr < 0, it does not mean this rank will not write to
     * the file. It just means intra-node aggregation is disabled for this
     * rank. This rank may still have data to write to the file.
     */

#ifdef NOT_YET
/* NOTE: using a communicator for intra-node aggregators does NOT work, because
 * cb_nodes may not be the intra-node aggregator and communication in 2-phase
 * I/O needs all cb_nodes to participate.
 */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        /* create a new MPI communicator containing intra-node aggregators */
        int j, ina_nprocs, *ina_ranks;
        MPI_Group origin_group, ina_group;

        assert(ncp->adio_fh != NULL);

        /* construct an arry of ranks that access to file */
        ina_ranks = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
        MPI_Allgather(&do_io, 1, MPI_INT, ina_ranks, 1, MPI_INT, ncp->comm);

        /* calculate number of processes that will access to the file */
        for (ina_nprocs=0, i=0; i<ncp->nprocs; i++)
            if (ina_ranks[i]) ina_nprocs++;

// printf("ina_nprocs=%d ina_ranks=%d %d %d\n",ina_nprocs, ina_ranks[0],ina_ranks[1],ina_ranks[2]);

        /* construct an arry of rank IDs that access to file */
        for (j=0,i=0; i<ncp->nprocs; i++)
            if (ina_ranks[i]) ina_ranks[j++] = i;

// printf("ina_nprocs=%d ina_ranks=%d %d\n",ina_nprocs, ina_ranks[0],ina_ranks[1]);

        /* create a new communicator containing only ranks that access to file */
        MPI_Comm_group(ncp->comm, &origin_group);
        MPI_Group_incl(origin_group, ina_nprocs, ina_ranks, &ina_group);
        MPI_Comm_create(ncp->comm, ina_group, &ncp->adio_fh->ina_comm);
        MPI_Group_free(&ina_group);
        MPI_Group_free(&origin_group);
        NCI_Free(ina_ranks);
    }

if (debug && ncp->fstype != ADIO_FSTYPE_MPIIO) {
    int io_rank;
    printf("world rank %d ncp->adio_fh->ina_comm = %s\n", ncp->rank,
           (ncp->adio_fh->ina_comm == MPI_COMM_NULL) ? "MPI_COMM_NULL" : "NOT NULL");

    if (ncp->adio_fh->ina_comm != MPI_COMM_NULL)
        MPI_Comm_rank(ncp->adio_fh->ina_comm, &io_rank);
    else
        io_rank = -1;
    printf("World rank = %d I/O comm rank = %d\n",ncp->rank, io_rank);
}
#endif

    /* TODO: For automatically determine Whether to enable intra-node write
     * aggregation, this should be done right before each collective write
     * call.
     *   1. obtain hint cb_noddes, and striping_unit
     *   2. calculate aggregate access region
     * In each round of two-phase I/O, when the number of senders to each
     * cb_nodes is very large, then intra-node aggregation should be enabled.
     * Average of all nprocs_per_node may be a factor for determining whether
     * to enable intra-node aggregation. It indicates whether the high number
     * of processes are allocated on the same node.
     */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncp->aggr_time[0] = MPI_Wtime() - timing;
#endif

    return NC_NOERR;
}

/*----< flatten_subarray() >-------------------------------------------------*/
/* flatten a subarray request into a list of offset-length pairs */
static int
flatten_subarray(int                ndim,       /* number of dimensions */
                 int                el_size,    /* array element size */
                 MPI_Offset         var_begin,  /* starting file offset */
                 const MPI_Offset  *dimlen,     /* [ndim] dimension lengths */
                 const MPI_Offset  *start,      /* [ndim] starts of subarray */
                 const MPI_Offset  *count,      /* [ndim] counts of subarray */
                 const MPI_Offset  *stride,     /* [ndim] strides of subarray */
                 MPI_Aint          *npairs,     /* OUT: num of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
                 MPI_Count         *offsets,    /* OUT: array of offsets */
                 MPI_Count         *lengths     /* OUT: array of lengths */
#else
                 MPI_Offset        *offsets,    /* OUT: array of offsets */
                 int               *lengths     /* OUT: array of lengths */
#endif
                                     )
{
    int i, j;
    MPI_Offset length, nstride, array_len, off, subarray_len;
    size_t idx=0, idx0;

    *npairs = 0;
    if (ndim < 0) return NC_NOERR;

    if (ndim == 0) {  /* scalar record variable */
        *npairs = 1;
        offsets[0] = var_begin;
        lengths[0] = el_size;
        return NC_NOERR;
    }

    /* TODO: check if all stride[] >= 1
       Q: Is it legal if any stride[] <= 0 ? */

    /* calculate the number of offset-length pairs */
    *npairs = (stride[ndim-1] == 1) ? 1 : count[ndim-1];
    for (i=0; i<ndim-1; i++)
        *npairs *= count[i];
    if (*npairs == 0) /* not reachable, an error if count[] == 0 */
        return NC_NOERR;

    /* length of each row of the subarray are of the same size */
    length  = (stride[ndim-1] == 1) ? count[ndim-1] : 1;
    length *= el_size;
    nstride  = (stride[ndim-1] == 1) ? 1 : count[ndim-1];

    /* set the offset-length pairs for the lowest dimension */
    off = var_begin + start[ndim-1] * el_size;
    for (i=0; i<nstride; i++) {
        offsets[idx]  = off;
        lengths[idx]  = length;
        off          += stride[ndim-1] * el_size;
        idx++;
    }
    ndim--;

    subarray_len = nstride;
    array_len = 1;
    /* for higher dimensions */
    while (ndim > 0) {
        /* array_len is global array size from lowest up to ndim */
        array_len *= dimlen[ndim];

        /* off is the global array offset for this dimension, ndim-1 */
        off = start[ndim-1] * array_len * el_size;

        /* update all offsets from lowest up to dimension ndim-1 */
        idx0 = 0;
        for (j=0; j<subarray_len; j++) {
            offsets[idx0] += off;
            idx0++;
        }

        /* update each plan subarray of dimension ndim-1 */
        off = array_len * stride[ndim-1] * el_size;
        for (i=1; i<count[ndim-1]; i++) {
            idx0 = 0;
            for (j=0; j<subarray_len; j++) {
                offsets[idx] = offsets[idx0] + off;
                lengths[idx] = length;
                idx++;
                idx0++;
            }
            off += array_len * stride[ndim-1] * el_size;
        }
        ndim--;  /* move to next higher dimension */
        subarray_len *= count[ndim];
    }

    return NC_NOERR;
}

/*----< flatten_req() >-----------------------------------------------------*/
/* flatten one write request into offset-length pairs.
 * offsets and lengths are allocated here and need to be freed by the caller
 */
static int
flatten_req(NC                *ncp,
            NC_var            *varp,
            const MPI_Offset  *start,
            const MPI_Offset  *count,
            const MPI_Offset  *stride,
            int               *is_incr,   /* OUT: are offsets incrementing */
            MPI_Aint          *num_pairs, /* OUT: number of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count        **offsets,   /* OUT: array of flattened offsets */
            MPI_Count        **lengths    /* OUT: array of flattened lengths */
#else
            MPI_Offset       **offsets,   /* OUT: array of flattened offsets */
            int              **lengths    /* OUT: array of flattened lengths */
#endif
                                   )
{
    int i, j, err=NC_NOERR, ndims;
    MPI_Aint num, idx;
    MPI_Offset var_begin, *shape, count0, *ones=NULL;

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs, so we can malloc a contiguous memory
     * space for storing off-len pairs
     */
    if (varp->ndims == 0) { /* scalar variable */
#ifdef HAVE_MPI_LARGE_COUNT
        *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count));
        *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count));
#else
        *offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset));
        *lengths = (int*)       NCI_Malloc(sizeof(int));
#endif
        (*offsets)[0] = varp->begin;
        (*lengths)[0] = varp->xsz;
        *num_pairs = 1;
        return NC_NOERR;
    }
    else if (varp->ndims == 1 && IS_RECVAR(varp)) { /* scalar variable */
        num = count[0];
    }
    else {
        num = 1;
        if (stride != NULL && stride[varp->ndims-1] > 1)
            num = count[varp->ndims-1];  /* count of last dimension */
        for (i=0; i<varp->ndims-1; i++)
            num *= count[i];       /* all count[] except the last dimension */
    }
    *num_pairs = num;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count prev_offset;
    *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num);
    *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num);
#else
    MPI_Offset prev_offset;
    *offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * num);
    *lengths = (int*)       NCI_Malloc(sizeof(int)        * num);
#endif

    if (stride == NULL) { /* equivalent to {1, 1, ..., 1} */
        ones = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * varp->ndims);
        for (i=0; i<varp->ndims; i++) ones[i] = 1;
    }

    ndims = varp->ndims;
    var_begin = varp->begin;
    shape = varp->shape;
    if (IS_RECVAR(varp)) {
        count0 = count[0];
        var_begin += start[0] * ncp->recsize;
        ndims--;
        start++;
        count++;
        shape++;
        if (stride != NULL) stride++;
    }
    else
        count0 = 1;

    idx = 0;
    *is_incr = 1;
    prev_offset = -1;
    for (i=0; i<count0; i++) {
        /* flatten the request into a list of offset-length pairs */
        err = flatten_subarray(ndims, varp->xsz, var_begin, shape,
                               start, count, (stride == NULL) ? ones : stride,
                               &num,            /* OUT: num of off-len pairs */
                               *offsets + idx,  /* OUT: array of offsets */
                               *lengths + idx); /* OUT: array of lengths */

        if (num == 0) continue;

        /* check if (*offsets)[] are in an increasing order */
        for (j=0; j<num; j++) {
            if (prev_offset > (*offsets)[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_offset = (*offsets)[idx+j];
        }

        idx += num;
        assert(idx <= *num_pairs);

        if (IS_RECVAR(varp))
            var_begin += ncp->recsize;
    }
    if (ones != NULL)
        NCI_Free(ones);

    return err;
}

/*----< flatten_reqs() >-----------------------------------------------------*/
/* flatten all write requests into offset-length pairs.
 * offsets and lengths are allocated here and need to be freed by the caller
 */
static int
flatten_reqs(NC            *ncp,
             int            num_reqs,  /* IN: # requests */
             const NC_req  *reqs,      /* [num_reqs] requests */
             int           *is_incr,   /* OUT: are offsets incrementing */
             MPI_Aint      *num_pairs, /* OUT: total number of off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
             MPI_Count    **offsets,   /* OUT: array of flattened offsets */
             MPI_Count    **lengths    /* OUT: array of flattened lengths */
#else
             MPI_Offset   **offsets,   /* OUT: array of flattened offsets */
             int          **lengths    /* OUT: array of flattened lengths */
#endif
                                   )
{
    int i, j, status=NC_NOERR, ndims, max_ndims=0;
    MPI_Aint num, idx;
    MPI_Offset *start, *count, *shape, *stride, *ones;

    *num_pairs = 0;    /* total number of offset-length pairs */

    /* Count the number off-len pairs from reqs[], so we can malloc a
     * contiguous memory space for storing off-len pairs
     */
#if 1
    for (i=0; i<num_reqs; i++) {
        /* reqs[i].npairs is the number of offset-length pairs of this request,
         * calculated in ncmpio_igetput_varm() and igetput_varn()
         */
        *num_pairs += reqs[i].npairs;
        ndims = ncp->put_lead_list[reqs[i].lead_off].varp->ndims;
        max_ndims = MAX(max_ndims, ndims);
    }
#else
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;
        ndims = lead->varp->ndims;
        max_ndims = MAX(max_ndims, ndims);
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

        /* for record variable, each reqs[] is within a record */
        if (IS_RECVAR(lead->varp)) {
            ndims--;
            start++;
            count++;
            stride++;
        }
        if (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) stride = NULL;

        if (ndims < 0) continue;
        if (ndims == 0) {  /* 1D record variable */
            (*num_pairs)++;
            continue;
        }
        num = 1;
        if (stride != NULL && stride[ndims-1] > 1)
            num = count[ndims-1];  /* count of last dimension */
        for (j=0; j<ndims-1; j++)
            num *= count[j];  /* all count[] except the last dimension */

        (*num_pairs) += num;
    }

MPI_Aint nn = 0;
int wkl=0;
    for (i=0; i<num_reqs; i++) {
        nn += reqs[i].npairs;
// if (reqs[i].npairs > 1 && wkl == 0) { printf("reqs[%d].npairs=%ld > 1\n",i,reqs[i].npairs); wkl++; }
    }
if (nn != *num_pairs) printf("%s %d -------- num_reqs=%d num_pairs=%ld nn=%ld\n",__func__,__LINE__,num_reqs,(*num_pairs),nn);
// wkl=0;
#endif

    /* now we can allocate a contiguous memory space for the off-len pairs */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count prev_offset;
    *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs));
    *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs));
#else
    MPI_Offset prev_offset;
    *offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (*num_pairs));
    *lengths = (int*)       NCI_Malloc(sizeof(int)        * (*num_pairs));
#endif

    ones = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * max_ndims);
    for (i=0; i<max_ndims; i++) ones[i] = 1;

    idx = 0;
    prev_offset = -1;
    *is_incr = 1;

    /* now re-run the loop to fill in the off-len pairs */
    for (i=0; i<num_reqs; i++) {
        MPI_Offset var_begin;
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;

        /* The case when reqs[i] contains only one offset-length pair has been
         * specially handled to in ncmpio_igetput_varm() and igetput_varn()
         * which store the relative offset to the beginning of the variable,
         * varp->begin, to reqs[i].offset_start.  The length is stored in
         * reqs[i].offset_end.
         */
        if (reqs[i].npairs == 1) {
            (*offsets)[idx] = reqs[i].offset_start + lead->varp->begin;
            if (IS_RECVAR(lead->varp))
                (*offsets)[idx] += reqs[i].start[0] * ncp->recsize;
            (*lengths)[idx] = reqs[i].nelems * lead->varp->xsz;

            /* check if (*offsets)[] are in an increasing order */
            if (prev_offset > (*offsets)[idx])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_offset = (*offsets)[idx];
            idx++;
            continue;
        }

        ndims = lead->varp->ndims;
        if (ndims > 0) {
            start  = reqs[i].start;
            count  = start + ndims;
            stride = count + ndims;
        }
        else
            start = count = stride = NULL;

        shape = lead->varp->shape;

        /* find the starting file offset for this variable */
        var_begin = lead->varp->begin;

        /* for record variable, each reqs[] is within a record */
        if (IS_RECVAR(lead->varp)) {
            ndims--;
            start++;
            count++;
            stride++;
            shape++;
            /* find the starting file offset for this record */
            var_begin += reqs[i].start[0] * ncp->recsize;
        }

        if (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) stride = NULL;

        /* flatten each request into a list of offset-length pairs and
         * append to the end of offsets and lengths
         */
        flatten_subarray(ndims, lead->varp->xsz, var_begin, shape,
                         start, count, (stride == NULL) ? ones : stride,
                         &num,            /* OUT: number of off-len pairs */
                         *offsets + idx,  /* OUT: array of offsets */
                         *lengths + idx); /* OUT: array of lengths */
/*
if (num > 1 && wkl == 0) { printf("---------------------- num=%ld > 1\n",num); wkl++; }
if ((*offsets)[idx] != reqs[i].offset_start && wkl == 0) {printf("---------------------- offsets[idx=%ld] %lld != reqs[i=%d].offset_start %lld\n",idx, (*offsets)[idx], i, reqs[i].offset_start); wkl++; }
*/

        /* check if (*offsets)[] are in an increasing order */
        for (j=0; j<num; j++) {
            if (prev_offset > (*offsets)[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_offset = (*offsets)[idx+j];
        }
        idx += num;
    }
    NCI_Free(ones);

    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead = ncp->put_lead_list + reqs[i].lead_off;
        if (fIsSet(lead->flag, NC_REQ_TO_FREE)) {
            NCI_Free(lead->start);
            lead->start = NULL;
        }
    }

    return status;
}

/*----< construct_buf_type() >-----------------------------------------------*/
/* construct an MPI derived datatype for I/O buffers from the request list, by
 * concatenate all buffers.
 */
static int
construct_buf_type(const NC     *ncp,
                   int           num_reqs,  /* IN: # requests */
                   const NC_req *reqs,      /* [num_reqs] requests */
                   MPI_Aint     *bufLen,    /* OUT: buffer size in bytes */
                   MPI_Datatype *bufType)   /* OUT: buffer datatype */
{
    int i, err, mpireturn, status=NC_NOERR;
    NC_lead_req *lead;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *blocklens = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    MPI_Count *disps     = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
#else
    int       *blocklens = (int*)      NCI_Malloc(sizeof(int)       * num_reqs);
    MPI_Aint  *disps     = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint)  * num_reqs);
#endif

    *bufLen = 0;
    for (i=0; i<num_reqs; i++) {
        MPI_Aint addr;

        /* displacement uses MPI_BOTTOM */
        MPI_Get_address(reqs[i].xbuf, &addr);
        disps[i] = addr;

        /* blocklens[] in bytes */
        lead = ncp->put_lead_list + reqs[i].lead_off;
        blocklens[i] = reqs[i].nelems * lead->varp->xsz;

        *bufLen += blocklens[i];
    }

    /* construct buffer derived datatype */
#ifdef HAVE_MPI_LARGE_COUNT
    mpireturn = MPI_Type_create_hindexed_c(num_reqs, blocklens, disps,
                                           MPI_BYTE, bufType);
#else
    mpireturn = MPI_Type_create_hindexed(num_reqs, blocklens, disps,
                                         MPI_BYTE, bufType);
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_hindexed");
        /* return the first encountered error if there is any */
        if (status == NC_NOERR) status = err;

        *bufType = MPI_DATATYPE_NULL;
    }
    else {
        MPI_Type_commit(bufType);
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count typeSize;
        MPI_Type_size_c(*bufType, &typeSize);
#else
        int typeSize;
        MPI_Type_size(*bufType, &typeSize);
#endif
        assert(typeSize == *bufLen);
    }

    NCI_Free(blocklens);
    NCI_Free(disps);

    return status;
}

/*----< intra_node_aggregation() >-------------------------------------------*/
/* This is a collective call */
static int
intra_node_aggregation(NC           *ncp,
                       int           is_incr,  /* offsets are incrementing */
                       MPI_Aint      num_pairs,
#ifdef HAVE_MPI_LARGE_COUNT
                       MPI_Count    *offsets,
                       MPI_Count    *lengths,
#else
                       MPI_Offset   *offsets,
                       int          *lengths,
#endif
                       MPI_Offset    bufCount,
                       MPI_Datatype  bufType,
                       void         *buf)
{
    int i, j, err, mpireturn, status=NC_NOERR, nreqs, do_sort=0, indv_sorted=1;
    char *recv_buf=NULL, *wr_buf = NULL;
    MPI_Aint npairs=0, *msg, *count=NULL;
    MPI_Offset disp=0, buf_count=0;
    MPI_Datatype saved_fileType, fileType=MPI_BYTE;
    MPI_File fh;
    MPI_Request *req=NULL;
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double endT, startT = MPI_Wtime();
#endif
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufLen;
    MPI_Type_size_c(bufType, &bufLen);
#else
    int bufLen;
    MPI_Type_size(bufType, &bufLen);
#endif
    bufLen *= bufCount;

MPI_Offset maxm;
ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

    /* First, tell aggregator how much to receive by sending:
     * (num_pairs and bufLen). The message size to be sent by this rank
     * is num_pairs * 2 * sizeof(MPI_Offset) + bufLen
     */
    if (ncp->rank == ncp->my_aggr)
        msg = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs * 3);
    else
        msg = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 3);

    msg[0] = num_pairs;
    msg[1] = bufLen;
    msg[2] = is_incr;

    /* Aggregator collects each non-aggregator's num_pairs and bufLen */
    if (ncp->rank == ncp->my_aggr) {
        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
        nreqs = 0;
        for (i=1; i<ncp->num_nonaggrs; i++)
            MPI_Irecv(msg + i*3, 3, MPI_AINT, ncp->nonaggr_ranks[i], 0,
                      ncp->comm, &req[nreqs++]);

        mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    else /* non-aggregator */
        MPI_Send(msg, 3, MPI_AINT, ncp->my_aggr, 0, ncp->comm);

    /* Aggregator collects offset-length pairs from non-aggregators */
    if (ncp->rank == ncp->my_aggr) {
        MPI_Datatype recvTypes;

        /* calculate the total number of offset-length pairs */
        npairs = num_pairs;
        for (i=1; i<ncp->num_nonaggrs; i++) npairs += msg[i*3];

#ifdef HAVE_MPI_LARGE_COUNT
        if (npairs > num_pairs) {
            /* realloc to store all pairs in a contiguous buffer */
            offsets = (MPI_Count*) NCI_Realloc(offsets, sizeof(MPI_Count) * npairs);
            lengths = (MPI_Count*) NCI_Realloc(lengths, sizeof(MPI_Count) * npairs);
        }
#else
        if (npairs > num_pairs) {
            /* realloc to store all pairs in a contiguous buffer */
            offsets = (MPI_Offset*) NCI_Realloc(offsets, sizeof(MPI_Offset) * npairs);
            lengths = (int*) NCI_Realloc(lengths, sizeof(int) * npairs);
        }
#endif

        nreqs = 0;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        MPI_Get_address(offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Count) * msg[0]);
        MPI_Get_address(lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(MPI_Count) * msg[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (msg[i*3] == 0) continue;
            bklens[0] = msg[i*3] * sizeof(MPI_Count);
            bklens[1] = msg[i*3] * sizeof(MPI_Count);
            mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                                   &recvTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvTypes);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            MPI_Irecv_c(MPI_BOTTOM, 1, recvTypes, ncp->nonaggr_ranks[i],
                        0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvTypes);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#else
        int bklens[2];
        MPI_Aint aint, disps[2];

        MPI_Get_address(offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Aint) * msg[0]);
        MPI_Get_address(lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(int) * msg[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (msg[i*3] == 0) continue;
            bklens[0] = msg[i*3] * sizeof(MPI_Aint);
            bklens[1] = msg[i*3] * sizeof(int);
            mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                                 &recvTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvTypes);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            MPI_Irecv(MPI_BOTTOM, 1, recvTypes, ncp->nonaggr_ranks[i],
                      0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvTypes);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#endif
        mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    else if (num_pairs > 0) { /* non-aggregator */
        /* send offset-length pairs data to the aggregator */
        MPI_Datatype sendTypes;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        bklens[0] = msg[0] * sizeof(MPI_Count);
        bklens[1] = bklens[0];
        MPI_Get_address(offsets, &aint);
        disps[0] = aint;
        MPI_Get_address(lengths, &aint);
        disps[1] = aint;
        mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                               &sendTypes);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&sendTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        MPI_Send_c(MPI_BOTTOM, 1, sendTypes, ncp->my_aggr, 0, ncp->comm);
        MPI_Type_free(&sendTypes);
#else
        int bklens[2];
        MPI_Aint disps[2];

        bklens[0] = msg[0] * sizeof(MPI_Aint);
        bklens[1] = msg[0] * sizeof(int);
        MPI_Get_address(offsets, &disps[0]);
        MPI_Get_address(lengths, &disps[1]);
        mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                             &sendTypes);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&sendTypes);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        MPI_Send(MPI_BOTTOM, 1, sendTypes, ncp->my_aggr, 0, ncp->comm);
        MPI_Type_free(&sendTypes);
#endif

        /* Non-aggregators send write data to the aggregator */
        void *buf_ptr = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count num = (buf == NULL) ? 1 : bufCount;
        MPI_Send_c(buf_ptr, num, bufType, ncp->my_aggr, 0, ncp->comm);
#else
        int num = (buf == NULL) ? 1 : bufCount;
        MPI_Send(buf_ptr, num, bufType, ncp->my_aggr, 0, ncp->comm);
#endif
        NCI_Free(offsets);
        NCI_Free(lengths);
        offsets = NULL;
        lengths = NULL;
    }

ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

#ifdef PNETCDF_PROFILING
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[2] += endT - startT;
    startT = endT;
#endif

    if (ncp->rank == ncp->my_aggr && npairs > 0) {
        /* check if qsort is necessary. Find the first non-aggregator with
         * non-zero pairs
         */
        for (i=-1,j=0; j<ncp->num_nonaggrs; j++) {
            if (i == -1 && msg[j*3] > 0) /* find 1st one whose num_paris > 0 */
                i = j;
            if (msg[j*3+2] == 0) { /* individual j list is not sorted */
                indv_sorted = 0;
                do_sort = 1;
            }
        }
        if (i >= 0 && indv_sorted == 1) { /* further check if sort is needed */
            if (msg[i*3+2] == 0) /* first is_incr */
                do_sort = 1;
            else {
                MPI_Aint sum = msg[i*3];
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count prev_off = offsets[sum-1];
#else
                MPI_Offset prev_off = offsets[sum-1];
#endif
                for (++i; i<ncp->num_nonaggrs; i++) {
                    if (msg[i*3] == 0) continue;
                    if (msg[i*3+2] == 0 || prev_off > offsets[sum]) {
                        do_sort = 1; /* offsets are not incrementing */
                        break;
                    }
                    sum += msg[i*3];
                    prev_off = offsets[sum-1];
                }
            }
        }
        if (do_sort && indv_sorted) {
            /* count will be used in heap_merge() */
            count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs);
            for (i=0; i<ncp->num_nonaggrs; i++) count[i] = msg[i*3];
        }

        /* construct array of buffer addresses */
        MPI_Aint *bufAddr = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * npairs);
        bufAddr[0] = 0;
        for (i=1; i<npairs; i++)
            bufAddr[i] = bufAddr[i-1] + lengths[i-1];

        /* sort offsets, lengths, bufAddr altogether, based on offsets into
         * an increasing order
         */
        if (do_sort) {
            if (indv_sorted) {
                /* heap-merge of already sorted lists is much faster than qsort.
                 * Howevr, it has a much bigger memory footprint.
                 */
                heap_merge(ncp->num_nonaggrs, count, npairs, offsets, lengths,
                           bufAddr);
                NCI_Free(count);
            }
            else
                /* qsort is an in-place sorting */
                qsort_off_len_buf(npairs, offsets, lengths, bufAddr);
        }
ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

#ifdef PNETCDF_PROFILING
        endT = MPI_Wtime();
        if (ncp->rank == ncp->my_aggr) {
            ncp->aggr_time[3] += endT - startT;
            ncp->aggr_time[5] = MAX(ncp->aggr_time[5], npairs);
        }
        startT = endT;
#endif

        /* calculate the total write account */
        buf_count = bufLen;
        for (i=1; i<ncp->num_nonaggrs; i++) buf_count += msg[i*3 + 1];

        /* Allocate receive buffer, which will be sorted into an increasing
         * order based on the file offsets. Thus, after sorting pack recv_buf
         * to wr_buf to avoid creating another buffer datatype.
         */
        if (buf_count > 0) {
            recv_buf = (char*) NCI_Malloc(buf_count);
            wr_buf = (char*) NCI_Malloc(buf_count);
        }
ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

        /* First, pack self write data into front of the recv_buf */
        if (bufLen > 0) {
            if (bufType == MPI_BYTE)
                memcpy(recv_buf, buf, bufLen);
            else {
                void *inbuf = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count position=0;
                MPI_Count incount = (buf == NULL) ? 1 : bufCount;
                MPI_Pack_c(inbuf, incount, bufType, recv_buf, bufLen, &position,
                           MPI_COMM_SELF);
#else
                int position=0;
                int incount = (buf == NULL) ? 1 : bufCount;
                MPI_Pack(inbuf, incount, bufType, recv_buf, bufLen, &position,
                         MPI_COMM_SELF);
#endif
            }
        }

        /*
         * TODO, define a datatype to combine sends of offset-length pairs with
         * the write data into a single send call.
         */

        /* post requests to receive write data from non-aggregators */
        if (buf_count > 0) {
            char *ptr = recv_buf + bufLen;
            nreqs = 0;
            for (i=1; i<ncp->num_nonaggrs; i++) {
                if (msg[i*3 + 1] == 0) continue;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(ptr, msg[i*3 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                            0, ncp->comm, &req[nreqs++]);
#else
                MPI_Irecv(ptr, msg[i*3 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                          0, ncp->comm, &req[nreqs++]);
#endif
                ptr += msg[i*3 + 1];
            }
            mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }

        /* merge the overlapped buffer segments, skip the overlapped regions
         * for those with higher j indices (i.e. requests with lower j indices
         * win the writes to the overlapped regions)
         */
        for (i=0, j=1; j<npairs; j++) {
            if (offsets[i] + lengths[i] >= offsets[j] + lengths[j])
                /* segment i completely covers segment j, skip j */
                continue;

            MPI_Offset gap = offsets[i] + lengths[i] - offsets[j];
            if (gap >= 0) { /* segments i and j overlaps */
                if (bufAddr[i] + lengths[i] == bufAddr[j] + gap) {
                    /* buffers i and j are contiguous, merge j to i */
                    lengths[i] += lengths[j] - gap;
                }
                else { /* buffers are not contiguous, reduce j's len */
                    offsets[i+1] = offsets[j] + gap;
                    lengths[i+1] = lengths[j] - gap;
                    bufAddr[i+1] = bufAddr[j] + gap;
                    i++;
                }
            }
            else { /* i and j do not overlap */
                i++;
                if (i < j) {
                    offsets[i] = offsets[j];
                    lengths[i] = lengths[j];
                    bufAddr[i] = bufAddr[j];
                }
            }
        }
        /* update number of pairs, now all off-len pairs are not overlapped */
        npairs = i+1;

        /* pack recv_buf, data received from non-aggregators, into wr_buf, a
         * contiguous buffer, wr_buf, which will later be used in a call to
         * MPI_File_write_at_all()
         */
        char *ptr = wr_buf;
        buf_count = 0;
        if (npairs > 0) {
            memcpy(ptr, recv_buf + bufAddr[0], lengths[0]);
            ptr += lengths[0];
            buf_count = lengths[0];
        }
        for (i=0, j=1; j<npairs; j++) {
            memcpy(ptr, recv_buf + bufAddr[j], lengths[j]);
            ptr += lengths[j];
            /* overlap may be found, recalculate buf_count */
            buf_count += lengths[j];

            /* coalesce the offset-length pairs */
            if (offsets[i] + lengths[i] == offsets[j]) {
                /* coalesce j into i */
                lengths[i] += lengths[j];
            }
            else {
                i++;
                if (i < j) {
                    offsets[i] = offsets[j];
                    lengths[i] = lengths[j];
                }
            }
        }
        NCI_Free(bufAddr);
        if (recv_buf != NULL) NCI_Free(recv_buf);

        /* update number of pairs, now all off-len pairs are not overlapped */
        npairs = i+1;

        if (ncp->fstype == ADIO_FSTYPE_MPIIO) {

            if (npairs == 1) {
                /* Need not create fileType if writing to a contiguous space,
                 * as when fileType is MPI_BYTE, ncmpio_file_set_view() will
                 * make the entire file visible.
                 */
                disp = offsets[0];
            }
            else {
#ifdef HAVE_MPI_LARGE_COUNT
                /* construct fileview */
                mpireturn = MPI_Type_create_hindexed_c(npairs, lengths, offsets,
                                                       MPI_BYTE, &fileType);

#else
                assert(sizeof(*offsets) == sizeof(MPI_Aint));
                /* construct fileview */
                mpireturn = MPI_Type_create_hindexed(npairs, lengths,
                                                     (MPI_Aint*)offsets,
                                                     MPI_BYTE, &fileType);

#endif
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
                else {
                    mpireturn = MPI_Type_commit(&fileType);
                    if (mpireturn != MPI_SUCCESS) {
                        err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                        /* return the first encountered error if there is any */
                        if (status == NC_NOERR) status = err;
                    }
                }
            }
            NCI_Free(offsets);
            NCI_Free(lengths);
            offsets = NULL;
            lengths = NULL;
        }
    }

    NCI_Free(msg);
    if (ncp->rank == ncp->my_aggr)
        NCI_Free(req);

ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("==== %s line %d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[4] += endT - startT;
#endif

    if (ncp->rank != ncp->my_aggr) { /* non-aggregator writes nothing */
        buf_count = 0;
        npairs = 0;
    }

    /* Only aggregators writes non-zero sized of data to the file. The
     * non-aggregators participate the collective write call with zero-length
     * write requests.
     */
    fh = ncp->collective_fh;

    /* Preserve fd->filetype previously used */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        disp = 0;
        saved_fileType = ncp->adio_fh->filetype;

        /* When using internal Lustre driver, pass offsets and lengths to
         * ncmpio_file_set_view() and then pass to ADIO_File_set_view(), to
         * avoid creating a fileType and later be flattened. This can also save
         * memory space.
         *
         * Use MPI_DATATYPE_NULL to  indicate this setting of fileview comes
         * from intra-node aggregation.
         */
        fileType = MPI_DATATYPE_NULL;
    }

    /* set the MPI-IO fileview, this is a collective call */
    err = ncmpio_file_set_view(ncp, fh, &disp, fileType, npairs,
                               offsets, lengths);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        buf_count = 0;
    }
    if (ncp->fstype == ADIO_FSTYPE_MPIIO && fileType != MPI_BYTE)
        MPI_Type_free(&fileType);

ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

    /* call MPI_File_write_at_all, which returns disp that should point to the
     * very file offset to be written
     */
    err = ncmpio_read_write(ncp, NC_REQ_WR, NC_REQ_COLL, disp, buf_count,
                            MPI_BYTE, wr_buf, 1);
    if (status == NC_NOERR) status = err;

    if (wr_buf  != NULL) NCI_Free(wr_buf);
    if (offsets != NULL) NCI_Free(offsets);
    if (lengths != NULL) NCI_Free(lengths);

ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

    /* ncp->adio_fh->flat_file is allocated in ncmpio_file_set_view() */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        NCI_Free(ncp->adio_fh->flat_file);
        ncp->adio_fh->flat_file = NULL;

        /* restore the original filetype */
        ncp->adio_fh->filetype = saved_fileType;
    }

    return status;
}

/*----< ncmpio_intra_node_aggregation_nreqs() >------------------------------*/
/* This is a collective call */
int
ncmpio_intra_node_aggregation_nreqs(NC         *ncp,
                                    int         reqMode,
                                    int         num_reqs,
                                    NC_req     *put_list,
                                    MPI_Offset  newnumrecs)
{
    int err, status=NC_NOERR, is_incr=1;
    MPI_Aint bufLen, num_pairs;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *offsets=NULL, *lengths=NULL;
#else
    MPI_Offset *offsets=NULL;
    int *lengths=NULL;
#endif
    MPI_Datatype bufType=MPI_BYTE;
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif

MPI_Offset maxm;
ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

    /* currently supports write requests only */
    if (fIsSet(reqMode, NC_REQ_RD)) return NC_NOERR;

    assert(ncp->my_aggr >= 0);

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     */
    if (num_reqs > 0)
        flatten_reqs(ncp, num_reqs, put_list, &is_incr, &num_pairs, &offsets,
                     &lengths);
    else
        num_pairs = 0;

    /* construct write buffer datatype, bufType.
     * bufLen is the buffer size in bytes
     */
    if (num_reqs > 0) {
        construct_buf_type(ncp, num_reqs, put_list, &bufLen, &bufType);
        bufLen = 1;
    }
    else
        bufLen = 0;

    if (put_list != NULL)
        NCI_Free(put_list);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[1] += MPI_Wtime() - timing;
#endif

    err = intra_node_aggregation(ncp, is_incr, num_pairs, offsets, lengths,
                                 bufLen, bufType, NULL);
    if (status == NC_NOERR) status = err;

    /* free and reset bufType */
    if (bufType != MPI_BYTE && bufType != MPI_DATATYPE_NULL)
        MPI_Type_free(&bufType);

    /* Update the number of records if new records have been created.
     * For nonblocking APIs, there is no way for a process to know whether
     * others write to a record variable or not. Note newnumrecs has been
     * sync-ed and always >= ncp->numrecs.
     */
    if (newnumrecs > ncp->numrecs) {
        /* update new record number in file. Note newnumrecs is already
         * sync-ed among all processes and in collective mode
         * ncp->numrecs is always sync-ed in memory among processes,
         * thus no need another MPI_Allreduce to sync it. */
        err = ncmpio_write_numrecs(ncp, newnumrecs);
        if (status == NC_NOERR) status = err;
        /* retain the first error if there is any */
        if (ncp->numrecs < newnumrecs) ncp->numrecs = newnumrecs;
    }

    return status;
}

/*----< ncmpio_intra_node_aggregation() >------------------------------------*/
/* This is a collective call */
int
ncmpio_intra_node_aggregation(NC               *ncp,
                              int               reqMode,
                              NC_var           *varp,
                              const MPI_Offset *start,
                              const MPI_Offset *count,
                              const MPI_Offset *stride,
                              MPI_Offset        bufCount,
                              MPI_Datatype      bufType,
                              void             *buf)
{
    int err, status=NC_NOERR, is_incr=1;
    MPI_Aint num_pairs;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *offsets=NULL, *lengths=NULL;
#else
    MPI_Offset *offsets=NULL;
    int *lengths=NULL;
#endif
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif

MPI_Offset maxm;
ncmpi_inq_malloc_max_size(&maxm); if (debug && ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);

    /* currently supports write requests only */
    if (fIsSet(reqMode, NC_REQ_RD)) return NC_NOERR;

    if (buf == NULL) /* zero-length request */
        return intra_node_aggregation(ncp, 1, 0, NULL, NULL, 0, MPI_BYTE, NULL);

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     */
    err = flatten_req(ncp, varp, start, count, stride, &is_incr, &num_pairs,
                      &offsets, &lengths);
    if (err != NC_NOERR) {
        num_pairs = 0;
        if (offsets != NULL)
            NCI_Free(offsets);
        offsets = NULL;
    }
    status = err;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[1] += MPI_Wtime() - timing;
#endif

    err = intra_node_aggregation(ncp, is_incr, num_pairs, offsets, lengths,
                                 bufCount, bufType, buf);
    if (status == NC_NOERR) status = err;

    return status;
}

