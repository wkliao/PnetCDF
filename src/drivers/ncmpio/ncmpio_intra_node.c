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
 * aggregators on that node, i.e. there is one aggregator per group. The rank
 * IDs of processes belong to the same group must be established.
 *
 * Note this subroutine must be called before MPI_File_open().
 *
 * 1. Find the affinity of each MPI process to its compute node. This should
 *    have been done from a call to ncmpii_construct_node_list() during
 *    ncmpio_create() and ncmpio_open().
 * 2. Determine whether self process is an intra-node aggregator.
 * 3. For an aggregator, find the number of non-aggregators assigned to it and
 *    construct a list of rank IDs of non-aggregators assigned to it.
 * 4. For a non-aggregator, find the rank ID of its assigned aggregator.
 *
 * This subroutine should be called only once per file at ncmpio_create() or
 * ncmpio_open(). It sets the following variables.
 *   ncp->my_aggr          rank ID of my aggregator
 *   ncp->num_nonaggrs     number of non-aggregators assigned to this rank
 *   ncp->nonaggr_ranks[]  rank IDs of assigned non-aggregators
 */
int
ncmpio_intra_node_aggr_init(NC *ncp)
{
    int i, naggrs_my_node;
    int my_rank_index, *ranks_my_node, my_node_id, nprocs_my_node;

    /* initialize parameters of local-node aggregation */
    ncp->my_aggr = -1;         /* rank ID of my aggregator */
    ncp->num_nonaggrs = 0;     /* number of non-aggregators assigned */
    ncp->nonaggr_ranks = NULL; /* ranks of assigned non-aggregators */

    /* Here, ncp->num_aggrs_per_node must be > 0, because checking if
     * ncp->num_aggrs_per_node == 0  has been done before entering this
     * subroutine.
     */
    assert(ncp->num_aggrs_per_node > 0);

    /* check ill value of num_aggrs_per_node */
    if (ncp->num_aggrs_per_node * ncp->num_nodes >= ncp->nprocs)
        /* disable intra-node aggregation */
        return NC_NOERR;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
    ncp->aggr_time[0] = ncp->aggr_time[1] = 0.0;
    ncp->aggr_time[2] = ncp->aggr_time[3] = 0.0;
    ncp->aggr_time[4] = ncp->aggr_time[5] = 0.0;
#endif

    /* ncp->node_ids[] has been established in ncmpii_construct_node_list()
     * called in ncmpio_create() or ncmpio_open() before entering this
     * subroutine.
     * my_node_id is this rank's node ID.
     */
    my_node_id = ncp->node_ids[ncp->rank];

    /* Set the following variables:
     *   nprocs_my_node: the number of processes in my nodes
     *   ranks_my_node[]: rank IDs of all processes in my node.
     *   my_rank_index: points to ranks_my_node[] where
     *       ranks_my_node[my_rank_index] == ncp->rank
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

    /* Make sure number of aggregators in my node <= nprocs_my_node. In some
     * cases, the number of processes allocated to the last few nodes can be
     * less than the ones before them.
     */
    naggrs_my_node = MIN(ncp->num_aggrs_per_node, nprocs_my_node);

    /* For the aggregation group this rank belongs to, calculate the number of
     * non-aggregators assigned to this group's aggregator, ncp->num_nonaggrs.
     * Note ncp->num_nonaggrs includes self rank.
     */
    ncp->num_nonaggrs = nprocs_my_node / naggrs_my_node;
    if (nprocs_my_node % naggrs_my_node) ncp->num_nonaggrs++;

    /* Set the number of non-aggregators assigned to this rank. For the last
     * group, make sure it does not go beyond nprocs_my_node.
     */
    int first_rank = my_rank_index - my_rank_index % ncp->num_nonaggrs;
    ncp->num_nonaggrs = MIN(ncp->num_nonaggrs, nprocs_my_node - first_rank);

    if (ncp->num_nonaggrs <= 1) {
        /* Disable aggregation if the number of non-aggregators assigned to
         * this aggregator is 1. Note num_nonaggrs includes self rank. It is
         * possible for aggregation to be enabled or disabled on different
         * groups (thus also among different nodes).
         *
         * Use whether ncp->my_aggr < 0 to tell if aggregation is disabled or
         * enabled.
         */
        ncp->my_aggr = -1;
    }
    else { /* num_nonaggrs > 1 */
        /* Construct ncp->nonaggr_ranks[], the rank IDs of non-aggregators
         * assigned in this group.  Note ncp->nonaggr_ranks[] will be freed
         * when closing the file, if allocated.
         */
        ncp->nonaggr_ranks = (int*)NCI_Malloc(sizeof(int) * ncp->num_nonaggrs);

        memcpy(ncp->nonaggr_ranks, ranks_my_node + first_rank,
               sizeof(int) * ncp->num_nonaggrs);

        /* Assign the first rank as the intra-node aggregator of this group */
        ncp->my_aggr = ncp->nonaggr_ranks[0];

#if 0
        /* It is no longer necessary to make intra-node aggregators also as
         * MPI-IO aggregators, because we later construct a new MPI
         * communicator consisting only the intra-node aggregators and pass it
         * to MPI_File_open(). It means only the intra-node aggregators have
         * non-zero sized data to do I/O.
         */
        if (ncp->fstype != ADIO_FSTYPE_MPIIO && ncp->adio_fh != NULL) {
            /* Check if the aggregator is also an MPI-IO aggregator. If not,
             * change it to the first MPI-IO aggregator in the group. Note that
             * binary search cannot be used, because the MPI rank IDs of all
             * MPI-IO aggregators, ncp->adio_fh->hints->ranklist[], may not be
             * sorted in an increasing order. On the other hand, rank IDs in
             * ranks_my_node[], rank IDs of this rank's intra-node group, are
             * in an increasing order.
             *
             * Note MPI ranks allocated to compute nodes are set at the run
             * time. Two common allocation methods are block (consecutive IDs
             * on a node) and cyclic (round-robin over all nodes). Both result
             * in MPI rank IDs allocated to each compute node in an increasing
             * order.
             */
            int j, cb_nodes = ncp->adio_fh->hints->cb_nodes;
            int *ranklist = ncp->adio_fh->hints->ranklist;

            for (i=0; i<ncp->num_nonaggrs; i++) {
                for (j=0; j<cb_nodes; j++) {
                    if (ncp->nonaggr_ranks[i] == ranklist[j])
                        break; /* found */
                }
                if (j < cb_nodes) {
                    if (i > 0) {
                        /* change intra-node aggregator */
                        ncp->my_aggr = ncp->nonaggr_ranks[i];
                        /* swap aggregator's rank with nonaggr_ranks[0] */
                        int tmp = ncp->nonaggr_ranks[0];
                        ncp->nonaggr_ranks[0] = ncp->my_aggr;
                        ncp->nonaggr_ranks[i] = tmp;
                    }
                    break;
                }
            }
        }
#endif
        if (ncp->my_aggr != ncp->rank) {
            /* free ncp->nonaggr_ranks as it is used by aggregators only */
            NCI_Free(ncp->nonaggr_ranks);
            ncp->nonaggr_ranks = NULL;
        }
    }
    NCI_Free(ranks_my_node);

    /* Note when ncp->my_aggr < 0, it does not mean this rank will not write to
     * the file. It just means intra-node aggregation is disabled for this
     * rank. This rank may still have data to write to the file.
     */

    /* Now only the intra-node aggregators have non-zero sized data to do I/O.
     * We can create a new MPI communicator consisting of all intra-node
     * aggregators and use it to call MPI_File_open().
     *
     * When using the internal ADIO driver, we also must create a new node_ids
     * based on the new communicator and pass it to PnetCDF's internal ADIO
     * file handler, ncp->adio_fh, so to avoid repeating the work of
     * constructing the node_ids again. If using MPI-IO driver, then ROMIO will
     * do this internally again anyway.
     */

    /* create a new MPI communicator containing intra-node aggregators */
    int j, do_io, ina_nprocs, *ina_ranks;
    MPI_Group origin_group, ina_group;

    do_io = 0;
    if (ncp->my_aggr == -1 || ncp->my_aggr == ncp->rank)
        do_io = 1;

    /* construct an array of ranks that access to file */
    ina_ranks = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
    MPI_Allgather(&do_io, 1, MPI_INT, ina_ranks, 1, MPI_INT, ncp->comm);

    /* node_ids[] can be two kinds (nodes=3, nprocs=10, num_aggrs_per_node=2)
     *     block  process allocation: 0,0,0,0,1,1,1,2,2,2
     *     cyclic process allocation: 0,1,2,0,1,2,0,1,2,0
     * Accordingly, ina_ranks[] can be two kinds
     *     block  process allocation: 1,0,1,0,1,0,1,1,0,1
     *     cyclic process allocation: 1,1,1,0,0,0,1,1,1,0
     */

    /* calculate number of intra-node aggregators, the ones that will access to
     * the file
     */
    for (ina_nprocs=0, i=0; i<ncp->nprocs; i++)
        if (ina_ranks[i]) ina_nprocs++;

// if (ncp->nprocs==4) printf("%s at %d: ina_nprocs=%d ina_ranks=%d %d %d %d\n", __func__,__LINE__,ina_nprocs, ina_ranks[0],ina_ranks[1],ina_ranks[2],ina_ranks[3]);

    /* construct the rank IDs of the new MPI communicator, ina_ranks[] */
    for (j=0,i=0; i<ncp->nprocs; i++) {
        if (ina_ranks[i]) {
            ina_ranks[j] = i;
            /* Modify ncp->node_ids[] to store the node IDs of the processes in
             * the new communicator. Note ncp->node_ids[] from now on is used
             * by PnetCDF's ADIO driver only.
             */
            ncp->node_ids[j] = ncp->node_ids[i];
            j++;
        }
    }

    /* create a new communicator containing only ranks that access to file */
    MPI_Comm_group(ncp->comm, &origin_group);
    MPI_Group_incl(origin_group, ina_nprocs, ina_ranks, &ina_group);
    MPI_Comm_create(ncp->comm, ina_group, &ncp->ina_comm);
    MPI_Group_free(&ina_group);
    MPI_Group_free(&origin_group);
    NCI_Free(ina_ranks);

#ifdef WKL_DEBUG
int io_rank;
printf("%s at %d: world rank %d ncp->ina_comm = %s\n", __func__,__LINE__,ncp->rank, (ncp->ina_comm == MPI_COMM_NULL) ? "MPI_COMM_NULL" : "NOT NULL");

if (ncp->ina_comm != MPI_COMM_NULL)
        MPI_Comm_rank(ncp->ina_comm, &io_rank);
else
        io_rank = -1;
printf("%s at %d: World rank = %d I/O comm rank = %d\n", __func__,__LINE__,ncp->rank, io_rank);
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
    MPI_Count prev_end_off;
    *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num);
    *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num);
#else
    MPI_Offset prev_end_off;
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
    prev_end_off = -1;
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
            if (prev_end_off > (*offsets)[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = (*offsets)[idx+j];
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
    MPI_Count prev_end_off;
    *offsets = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs));
    *lengths = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * (*num_pairs));
#else
    MPI_Offset prev_end_off;
    *offsets = (MPI_Offset*)NCI_Malloc(sizeof(MPI_Offset) * (*num_pairs));
    *lengths = (int*)       NCI_Malloc(sizeof(int)        * (*num_pairs));
#endif

    ones = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * max_ndims);
    for (i=0; i<max_ndims; i++) ones[i] = 1;

    idx = 0;
    prev_end_off = -1;
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
            if (prev_end_off > (*offsets)[idx])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = (*offsets)[idx];
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
            if (prev_end_off > (*offsets)[idx+j])
                *is_incr = 0;  /* offsets are not incrementing */
            else
                prev_end_off = (*offsets)[idx+j];
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
                   const NC_req *reqs,      /* IN: [num_reqs] requests */
                   MPI_Aint     *bufLen,    /* OUT: buffer size in bytes */
                   MPI_Datatype *bufType)   /* OUT: buffer datatype */
{
    int i, err, mpireturn, status=NC_NOERR;
    NC_lead_req *lead;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *disps     = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
    MPI_Count *blocklens = (MPI_Count*)NCI_Malloc(sizeof(MPI_Count) * num_reqs);
#else
    MPI_Aint  *disps     = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint)  * num_reqs);
    int       *blocklens = (int*)      NCI_Malloc(sizeof(int)       * num_reqs);
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

/*----< ina_collect_md() >---------------------------------------------------*/
/* Within each intra-node aggregation group, the aggregator collects request
 * metadata from the non-aggregators.
 */
static
int ina_collect_md(NC          *ncp,
                   MPI_Aint    *meta,
#ifdef HAVE_MPI_LARGE_COUNT
                   MPI_Count  **offsets, /* OUT: may be realloc-ed */
                   MPI_Count  **lengths, /* OUT: may be realloc-ed */
#else
                   MPI_Offset **offsets, /* OUT: may be realloc-ed */
                   int        **lengths, /* OUT: may be realloc-ed */
#endif
                   MPI_Aint    *npairs)  /* OUT: total no. off-len pairs */
{
    int i, err, mpireturn, status=NC_NOERR, nreqs;
    MPI_Request *req=NULL;
    MPI_Aint num_pairs=meta[0];

    /* Aggregator collects each non-aggregator's num_pairs and bufLen */
    if (ncp->rank == ncp->my_aggr) {

        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
        nreqs = 0;
        for (i=1; i<ncp->num_nonaggrs; i++)
            MPI_Irecv(meta + i*3, 3, MPI_AINT, ncp->nonaggr_ranks[i], 0,
                      ncp->comm, &req[nreqs++]);

#ifdef HAVE_MPI_STATUSES_IGNORE
        mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               ADIOI_Malloc(nreqs * sizeof(MPI_Status));
        mpireturn = MPI_Waitall(nreqs, req, statuses);
        ADIOI_Free(statuses);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    else /* non-aggregator */
        MPI_Send(meta, 3, MPI_AINT, ncp->my_aggr, 0, ncp->comm);

    /* Secondly, aggregators collect offset-length pairs from all its
     * non-aggregators
     */
    if (ncp->rank == ncp->my_aggr) {
        MPI_Datatype recvType;

        /* calculate the total number of offset-length pairs to receive */
        for (*npairs=0, i=0; i<ncp->num_nonaggrs; i++) *npairs += meta[i*3];

        /* offsets and lengths have been allocated for storing this rank's
         * offsets and lengths, realloc them to receive offsets and lengths
         * from non-aggregators so they can be in a contiguous buffer.
         */
#ifdef HAVE_MPI_LARGE_COUNT
        if (*npairs > num_pairs) {
            *offsets = (MPI_Count*) NCI_Realloc(*offsets, *npairs * sizeof(MPI_Count));
            *lengths = (MPI_Count*) NCI_Realloc(*lengths, *npairs * sizeof(MPI_Count));
        }
#else
        if (*npairs > num_pairs) {
            /* realloc to store all pairs in a contiguous buffer */
            *offsets = (MPI_Offset*) NCI_Realloc(*offsets, *npairs * sizeof(MPI_Offset));
            *lengths = (int*) NCI_Realloc(*lengths, *npairs * sizeof(int));
        }
#endif

        /* To minimize number of MPI recv calls per non-aggregator, below
         * creates a derived datatype, recvType, to combine offsets and lengths
         * into one MPI_Irecv call.
         */
        nreqs = 0;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        MPI_Get_address(*offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Count) * meta[0]);
        MPI_Get_address(*lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(MPI_Count) * meta[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (meta[i*3] == 0) continue;
            bklens[0] = meta[i*3] * sizeof(MPI_Count);
            bklens[1] = meta[i*3] * sizeof(MPI_Count);
            mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                                   &recvType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            MPI_Irecv_c(MPI_BOTTOM, 1, recvType, ncp->nonaggr_ranks[i],
                        0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvType);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#else
        int bklens[2];
        MPI_Aint aint, disps[2];

        MPI_Get_address(*offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Aint) * meta[0]);
        MPI_Get_address(*lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(int) * meta[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (meta[i*3] == 0) continue;
            bklens[0] = meta[i*3] * sizeof(MPI_Aint);
            bklens[1] = meta[i*3] * sizeof(int);
            mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                                 &recvType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else {
                mpireturn = MPI_Type_commit(&recvType);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                    /* return the first encountered error if there is any */
                    if (status == NC_NOERR) status = err;
                }
            }
            /* post to receive offset-length pairs from non-aggregators */
            MPI_Irecv(MPI_BOTTOM, 1, recvType, ncp->nonaggr_ranks[i],
                      0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvType);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#endif
#ifdef HAVE_MPI_STATUSES_IGNORE
        mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               ADIOI_Malloc(nreqs * sizeof(MPI_Status));
        mpireturn = MPI_Waitall(nreqs, req, statuses);
        ADIOI_Free(statuses);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    else if (num_pairs > 0) { /* non-aggregator */
        /* To minimize number of MPI send calls to the aggregator, below
         * creates a derived datatype, sendType, to combine offsets and lengths
         * into one MPI_Send call.
         */
        MPI_Datatype sendType;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Aint aint;
        MPI_Count bklens[2];
        MPI_Count disps[2];

        bklens[0] = meta[0] * sizeof(MPI_Count);
        bklens[1] = bklens[0];
        MPI_Get_address(*offsets, &aint);
        disps[0] = aint;
        MPI_Get_address(*lengths, &aint);
        disps[1] = aint;
        mpireturn = MPI_Type_create_hindexed_c(2, bklens, disps, MPI_BYTE,
                                               &sendType);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed_c");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&sendType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        MPI_Send_c(MPI_BOTTOM, 1, sendType, ncp->my_aggr, 0, ncp->comm);
        MPI_Type_free(&sendType);
#else
        int bklens[2];
        MPI_Aint disps[2];

        bklens[0] = meta[0] * sizeof(MPI_Aint);
        bklens[1] = meta[0] * sizeof(int);
        MPI_Get_address(*offsets, &disps[0]);
        MPI_Get_address(*lengths, &disps[1]);
        mpireturn = MPI_Type_create_hindexed(2, bklens, disps, MPI_BYTE,
                                             &sendType);
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&sendType);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
        }
        MPI_Send(MPI_BOTTOM, 1, sendType, ncp->my_aggr, 0, ncp->comm);
        MPI_Type_free(&sendType);
#endif
    }

    if (req != NULL) NCI_Free(req);

    return status;
}

/*----< ina_put() >----------------------------------------------------------*/
/* This collective subroutine implements the intra-node aggregation for write
 * operations.
 * offsets and lengths will be freed when return.
 */
static
int ina_put(NC           *ncp,
            int           is_incr,   /* if offsets are incremental */
            MPI_Aint      num_pairs, /* number of offset-length pairs */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count    *offsets,
            MPI_Count    *lengths,
#else
            MPI_Offset   *offsets,
            int          *lengths,
#endif
            MPI_Offset    bufCount,  /* number of user buffer data types */
            MPI_Datatype  bufType,   /* user buffer data type */
            void         *buf)       /* user buffer */
{
    int i, j, err, mpireturn, status=NC_NOERR, nreqs;
    char *recv_buf=NULL, *wr_buf = NULL;
    MPI_Aint npairs=0, *meta=NULL, *count=NULL;
    MPI_Offset disp=0, buf_count=0;
    MPI_Datatype saved_fileType, fileType=MPI_BYTE;
    MPI_File fh;
    MPI_Request *req=NULL;
#ifdef PNETCDF_PROFILING
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

#ifdef WKL_DEBUG
MPI_Offset maxm;
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

    /* Firstly, aggregators collect metadata from non-aggregators ------------*
     * This rank tells its aggregator how much metadata to receive from this
     * rank, by sending: the number of offset-length pairs (num_pairs) and user
     * buffer size in bytes (bufLen). This message size to be sent by this rank
     * is 3 MPI_Offset.
     */
    if (ncp->rank == ncp->my_aggr) {
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs * 3);
        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
    }
    else
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 3);

    meta[0] = num_pairs;
    meta[1] = bufLen;
    meta[2] = is_incr;

    /* Each aggregator first collects metadata about the number of
     * offset-length pairs, size of write request, and whether the offsets are
     * in an incremental order.  from non-aggregators assigned to it.
     */
    err = ina_collect_md(ncp, meta, &offsets, &lengths, &npairs);
    if (err != NC_NOERR) goto fn_exit;

    /* For write operation, the non-aggregators now can start sending their
     * write data to the aggregator.
     */
    if (ncp->rank != ncp->my_aggr && num_pairs > 0) { /* non-aggregator */
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

#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

#ifdef PNETCDF_PROFILING
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[2] += endT - startT;
    startT = endT;
#endif

    /* MPI-IO requires the flattened file offsets of a fileview to be in a
     * monotonic non-decreasing order.
     */
    if (ncp->rank == ncp->my_aggr && npairs > 0) {
        /* Now the aggregator has received all offset-length pairs from
         * non-aggregators. At first, check if qsort is necessary. Find the
         * first non-aggregator with non-zero pairs whose offsets are not
         * already sorted.
         */
        int do_sort=0, indv_sorted=1;

        for (i=-1,j=0; j<ncp->num_nonaggrs; j++) {
            if (i == -1 && meta[j*3] > 0) /* find 1st one whose num_pairs > 0 */
                i = j;
            if (meta[j*3+2] == 0) { /* non-aggregator j's offsets are not sorted */
                indv_sorted = 0;
                do_sort = 1;
                break;
            }
        }
        /* i is the first non-aggregator whose num_pairs > 0
         * j is the first non-aggregator whose is_incr is false
         */

        if (i >= 0 && indv_sorted == 1) {
            /* This is when all non-aggregators' offsets are individually
             * sorted. We still need to check if offsets are interleaved among
             * all non-aggregators to determine whether a sorting is needed.
             */
            if (meta[i*3+2] == 0) /* first non-zero offsets are not sorted */
                do_sort = 1;
            else {
                MPI_Aint sum = meta[i*3];
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count prev_end_off = offsets[sum-1];
#else
                MPI_Offset prev_end_off = offsets[sum-1];
#endif
                /* check if the offsets are interleaved */
                for (++i; i<ncp->num_nonaggrs; i++) {
                    if (meta[i*3] == 0) continue;
                    if (meta[i*3+2] == 0 || prev_end_off > offsets[sum]) {
                        do_sort = 1; /* offsets are not incrementing */
                        break;
                    }
                    sum += meta[i*3];
                    prev_end_off = offsets[sum-1];
                }
            }
        }

        if (do_sort && indv_sorted) {
            /* Need to do sorting. But individual offsets are already sorted.
             * For this case, heap_merge() is called to merge all offsets into
             * one single sorted offset list. Note count will be used in
             * heap_merge()
             */
            count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) *
            ncp->num_nonaggrs);
            for (i=0; i<ncp->num_nonaggrs; i++) count[i] = meta[i*3];
        }

        /* construct an array of buffer addresses containing a mapping of the
         * buffer used to receive write data from non-aggregators and the
         * buffer used to write to file.
         */
        MPI_Aint *bufAddr = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * npairs);
        bufAddr[0] = 0;
        for (i=1; i<npairs; i++)
            bufAddr[i] = bufAddr[i-1] + lengths[i-1];

        if (do_sort) {
            /* sort offsets, lengths, bufAddr altogether, based on offsets into
             * an increasing order
             */
            if (indv_sorted) {
                /* heap-merge of already sorted individual lists is much faster
                 * than qsort.  However, it has a much bigger memory footprint.
                 */
                heap_merge(ncp->num_nonaggrs, count, npairs, offsets, lengths,
                           bufAddr);
                NCI_Free(count);
            }
            else
                /* As some individual offsets are not sorted, we cannot use
                 * heap_merge().  Note qsort is an in-place sorting.
                 */
                qsort_off_len_buf(npairs, offsets, lengths, bufAddr);
        }
#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB sort=%s\n",__func__,__LINE__,(float)maxm/1048576.0, (do_sort)?((indv_sorted)?"HEAP":"QSORT"):"NO");
#endif

#ifdef PNETCDF_PROFILING
        endT = MPI_Wtime();
        if (ncp->rank == ncp->my_aggr) {
            ncp->aggr_time[3] += endT - startT;
            ncp->aggr_time[5] = MAX(ncp->aggr_time[5], npairs);
        }
        startT = endT;
#endif

        /* calculate the total amount to be written by this aggregator */
        for (buf_count=0,i=0; i<ncp->num_nonaggrs; i++)
            buf_count += meta[i*3 + 1];

        /* Allocate receive buffer and write buffer. Once write data from
         * non-aggregators have received into recv_buf, recv_buf is packed into
         * wr_buf. wr_buf is used when calling MPI-IO to write to file.
         */
        if (buf_count > 0) {
            recv_buf = (char*) NCI_Malloc(buf_count);
            wr_buf = (char*) NCI_Malloc(buf_count);
        }
#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

        /* First, pack self write data into front of the recv_buf */
        if (bufLen > 0) {
            int is_predef;
            if (bufType == MPI_DATATYPE_NULL)
                is_predef = 0;
            else
                ADIOI_Type_ispredef(bufType, &is_predef);

            if (is_predef)
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

        /* receive write data sent from non-aggregators */
        if (buf_count > 0) {
            char *ptr = recv_buf + bufLen;
            nreqs = 0;
            for (i=1; i<ncp->num_nonaggrs; i++) {
                if (meta[i*3 + 1] == 0) continue;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(ptr, meta[i*3 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                            0, ncp->comm, &req[nreqs++]);
#else
                MPI_Irecv(ptr, meta[i*3 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                          0, ncp->comm, &req[nreqs++]);
#endif
                ptr += meta[i*3 + 1];
            }
#ifdef HAVE_MPI_STATUSES_IGNORE
            mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
#else
            MPI_Status *statuses = (MPI_Status *)
                                   ADIOI_Malloc(nreqs * sizeof(MPI_Status));
            mpireturn = MPI_Waitall(nreqs, req, statuses);
            ADIOI_Free(statuses);
#endif
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
            if (gap >= 0) { /* segments i and j overlap */
                if (bufAddr[i] + lengths[i] == bufAddr[j] + gap) {
                    /* buffers i and j are contiguous, merge j into i */
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

    NCI_Free(meta);
    if (ncp->rank == ncp->my_aggr)
        NCI_Free(req);

ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("==== %s line %d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[4] += endT - startT;
#endif

    /* Only aggregators call MPI-IO functions to write data to the file.
     * Non-aggregators do not participate MPI-IO calls.
     */
    if (ncp->rank != ncp->my_aggr)
        goto fn_exit;

    /* intra-node aggregation only takes effect in collective data mode */
    fh = ncp->collective_fh;

    /* Preserve fd->filetype previously used */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        disp = 0;
        saved_fileType = ncp->adio_fh->filetype;

        /* When using internal Lustre driver, a call to ncmpio_file_set_view()
         * will pass offsets[] and lengths[] to ADIO_File_set_view() to be
         * reused as flattened fileType. This can avoid re-create a datatype
         * and re-flattening it. See ADIOI_Calc_my_off_len()
         *
         * Use MPI_DATATYPE_NULL to indicate this setting of fileview comes
         * from intra-node aggregation, so the ADIO subroutines can reuse
         * offsets[] and lengths[] and avoid freeing them.
         */
        fileType = MPI_DATATYPE_NULL;
    }

    /* Set the MPI-IO fileview (this is a collective call). This call returns
     * disp which points to the very first file offset to be written.
     */
    err = ncmpio_file_set_view(ncp, fh, &disp, fileType, npairs,
                               offsets, lengths);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        buf_count = 0;
    }
    if (ncp->fstype == ADIO_FSTYPE_MPIIO && fileType != MPI_BYTE)
        MPI_Type_free(&fileType);

#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

    /* call MPI_File_write_at_all */
    err = ncmpio_read_write(ncp, NC_REQ_WR, NC_REQ_COLL, disp, buf_count,
                            MPI_BYTE, wr_buf, 1);
    if (status == NC_NOERR) status = err;

    if (wr_buf  != NULL) NCI_Free(wr_buf);
    if (offsets != NULL) NCI_Free(offsets);
    if (lengths != NULL) NCI_Free(lengths);

#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif
    /* ncp->adio_fh->flat_file is allocated in ncmpio_file_set_view() */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        NCI_Free(ncp->adio_fh->flat_file);
        ncp->adio_fh->flat_file = NULL;

        /* restore the original filetype */
        ncp->adio_fh->filetype = saved_fileType;
    }

fn_exit:
    return status;
}

/*----< ina_get() >----------------------------------------------------------*/
/* This collective subroutine implements the intra-node aggregation for read
 * operations.
 * offsets and lengths will be freed when return.
 */
static
int ina_get(NC           *ncp,
            int           is_incr,   /* if offsets are incremental */
            MPI_Aint      num_pairs, /* number of offset-length pairs */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count    *offsets,
            MPI_Count    *lengths,
#else
            MPI_Offset   *offsets,
            int          *lengths,
#endif
            MPI_Offset    bufCount, /* number of user buffer data types */
            MPI_Datatype  bufType,  /* user buffer data type */
            void         *buf)      /* user buffer */
{
    int i, j, err, mpireturn, status=NC_NOERR, nreqs;
    char *send_buf=NULL, *rd_buf = NULL;
    MPI_Aint npairs=0, *meta=NULL, *count=NULL, *bufAddr=NULL;
    MPI_Offset disp=0, buf_count=0;
    MPI_Datatype saved_fileType, fileType=MPI_BYTE;
    MPI_File fh;
    MPI_Request *req=NULL;
#ifdef PNETCDF_PROFILING
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

#ifdef WKL_DEBUG
MPI_Offset maxm;
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

    /* Firstly, aggregators collect metadata from non-aggregators ------------*
     * This rank tells its aggregator how much metadata to receive from this
     * rank, by sending: the number of offset-length pairs (num_pairs) and user
     * buffer size in bytes (bufLen). This message size to be sent by this rank
     * is 3 MPI_Offset.
     */
    if (ncp->rank == ncp->my_aggr) {
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs * 3);
        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
    }
    else
        meta = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * 3);

    meta[0] = num_pairs;
    meta[1] = bufLen;
    meta[2] = is_incr;

    /* Each aggregator first collects metadata about the number of
     * offset-length pairs, size of write request, and whether the offsets are
     * in an incremental order.  from non-aggregators assigned to it.
     */
    err = ina_collect_md(ncp, meta, &offsets, &lengths, &npairs);
    if (err != NC_NOERR) goto fn_exit;

    if (ncp->rank != ncp->my_aggr && num_pairs > 0) {
        /* For read operation, the non-aggregators now can start receiving
         * their read data from the aggregator.
         */
        MPI_Status st;
        void *buf_ptr = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count num = (buf == NULL) ? 1 : bufCount;
        MPI_Recv_c(buf_ptr, num, bufType, ncp->my_aggr, 0, ncp->comm, &st);
#else
        int num = (buf == NULL) ? 1 : bufCount;
        MPI_Recv(buf_ptr, num, bufType, ncp->my_aggr, 0, ncp->comm), &st;
#endif
printf("%s at %d: bufCount=%lld num=%lld\n",__func__,__LINE__,bufCount,num);
        NCI_Free(offsets);
        NCI_Free(lengths);
        offsets = NULL;
        lengths = NULL;
    }

    /* Non-aggregators are now done. Below is for aggregators which must call
     * MPI-IO functions to read data from the file. Non-aggregators do not
     * participate MPI-IO calls.
     */
    if (ncp->rank != ncp->my_aggr)
        goto fn_exit;

#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

#ifdef PNETCDF_PROFILING
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[2] += endT - startT;
    startT = endT;
#endif

    /* MPI-IO requires the flattened file offsets of a fileview to be in a
     * monotonic non-decreasing order.
     */
    if (npairs > 0) {
        /* Now the aggregator has received all offset-length pairs from
         * non-aggregators. At first, check if qsort is necessary. Find the
         * first non-aggregator with non-zero pairs whose offsets are not
         * already sorted.
         */
        int do_sort=0, indv_sorted=1;

        for (i=-1,j=0; j<ncp->num_nonaggrs; j++) {
            if (i == -1 && meta[j*3] > 0) /* find 1st one whose num_pairs > 0 */
                i = j;
            if (meta[j*3+2] == 0) { /* non-aggregator j's offsets are not sorted */
                indv_sorted = 0;
                do_sort = 1;
                break;
            }
        }
        /* i is the first non-aggregator whose num_pairs > 0
         * j is the first non-aggregator whose is_incr is false
         */

        if (i >= 0 && indv_sorted == 1) {
            /* This is when all non-aggregators' offsets are individually
             * sorted. We still need to check if offsets are interleaved among
             * all non-aggregators to determine whether a sorting is needed.
             */
            if (meta[i*3+2] == 0) /* first non-zero offsets are not sorted */
                do_sort = 1;
            else {
                MPI_Aint sum = meta[i*3];
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count prev_end_off = offsets[sum-1];
#else
                MPI_Offset prev_end_off = offsets[sum-1];
#endif
                /* check if the offsets are interleaved */
                for (++i; i<ncp->num_nonaggrs; i++) {
                    if (meta[i*3] == 0) continue;
                    if (meta[i*3+2] == 0 || prev_end_off > offsets[sum]) {
                        do_sort = 1; /* offsets are not incrementing */
                        break;
                    }
                    sum += meta[i*3];
                    prev_end_off = offsets[sum-1];
                }
            }
        }

        if (do_sort && indv_sorted) {
            /* Need to do sorting. But individual offsets are already sorted.
             * For this case, heap_merge() is called to merge all offsets into
             * one single sorted offset list. Note count will be used in
             * heap_merge()
             */
            count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * ncp->num_nonaggrs);
            for (i=0; i<ncp->num_nonaggrs; i++) count[i] = meta[i*3];
        }

        /* construct an array of buffer addresses containing a mapping of the
         * buffer used to read from file and the buffer used to send read data
         * to non-aggregator.
         */
        bufAddr = (MPI_Aint*)NCI_Malloc(sizeof(MPI_Aint) * npairs);
        bufAddr[0] = 0;
        for (i=1; i<npairs; i++)
            bufAddr[i] = bufAddr[i-1] + lengths[i-1];

        if (do_sort) {
            /* sort offsets, lengths, bufAddr altogether, based on offsets into
             * an increasing order
             */
            if (indv_sorted) {
                /* heap-merge of already sorted individual lists is much faster
                 * than qsort.  However, it has a much bigger memory footprint.
                 */
                heap_merge(ncp->num_nonaggrs, count, npairs, offsets, lengths,
                           bufAddr);
                NCI_Free(count);
            }
            else
                /* As some individual offsets are not sorted, we cannot use
                 * heap_merge().  Note qsort is an in-place sorting.
                 */
                qsort_off_len_buf(npairs, offsets, lengths, bufAddr);
        }
#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB sort=%s\n",__func__,__LINE__,(float)maxm/1048576.0, (do_sort)?((indv_sorted)?"HEAP":"QSORT"):"NO");
#endif

#ifdef PNETCDF_PROFILING
        endT = MPI_Wtime();
        if (ncp->rank == ncp->my_aggr) {
            ncp->aggr_time[3] += endT - startT;
            ncp->aggr_time[5] = MAX(ncp->aggr_time[5], npairs);
        }
        startT = endT;
#endif

        /* calculate the total amount to be read by this aggregator */
        for (buf_count=0,i=0; i<ncp->num_nonaggrs; i++)
            buf_count += meta[i*3 + 1];

        /* Allocate read buffer and send buffer. Once data are read from file
         * into rd_buf, rd_buf is unpacked into send_buf for each
         * non-aggregator. send_buf will be used to send the read request data
         * to non-aggregators.
         */
        if (buf_count > 0) {
            send_buf = (char*) NCI_Malloc(buf_count);
            rd_buf = (char*) NCI_Malloc(buf_count);
        }
#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

        if (ncp->fstype == ADIO_FSTYPE_MPIIO) {
            /* When using MPI-IO, we need to create a fileType to set the file
             * view. The fileType is constructed based on the offset-length
             * pairs.
             */
            if (npairs == 1) {
                /* Need not create fileType if reading from a contiguous space,
                 * as when fileType is MPI_BYTE, ncmpio_file_set_view() will
                 * make the entire file visible.
                 */
                disp = offsets[0];
            }
            else {
                /* Must construct a fileType */
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
        }
    }
    /* else case: This aggregation group may not have data to read, but must
     * participate the collective MPI-IO calls.
     */

    /* intra-node aggregation only takes effect in collective data mode */
    fh = ncp->collective_fh;

    /* Preserve fd->filetype previously used */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        disp = 0;
        saved_fileType = ncp->adio_fh->filetype;

        /* When using internal Lustre driver, a call to ncmpio_file_set_view()
         * will pass offsets[] and lengths[] to ADIO_File_set_view() to be
         * reused as flattened fileType. This can avoid re-create a datatype
         * and re-flattening it. See ADIOI_Calc_my_off_len()
         *
         * Use MPI_DATATYPE_NULL to indicate this setting of fileview comes
         * from intra-node aggregation, so the ADIO subroutines can reuse
         * offsets[] and lengths[] and avoid freeing them.
         */
        fileType = MPI_DATATYPE_NULL;
    }

    /* Set the MPI-IO fileview (this is a collective call). This call returns
     * disp which points to the very first file offset to be read.
     */
    err = ncmpio_file_set_view(ncp, fh, &disp, fileType, npairs, offsets,
                               lengths);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        buf_count = 0;
    }
    if (ncp->fstype == ADIO_FSTYPE_MPIIO && fileType != MPI_BYTE)
        MPI_Type_free(&fileType);

#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

    /* call MPI_File_read_at_all */
    err = ncmpio_read_write(ncp, NC_REQ_RD, NC_REQ_COLL, disp, buf_count,
                            MPI_BYTE, rd_buf, 1);
    if (status == NC_NOERR) status = err;

#ifdef WKL_DEBUG
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif
    /* ncp->adio_fh->flat_file is allocated in ncmpio_file_set_view() */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        NCI_Free(ncp->adio_fh->flat_file);
        ncp->adio_fh->flat_file = NULL;

        /* restore the original filetype */
        ncp->adio_fh->filetype = saved_fileType;
    }

    /* offsets[] and lengths[] are used in ADIO read subroutines as flattened
     * filetype. They cannot be freed before the I/O is done.
     */
    if (offsets != NULL) NCI_Free(offsets);
    if (lengths != NULL) NCI_Free(lengths);

    if (npairs > 0) {
        /* Unpack rd_buf, data read from file, into send_buf, to be sent to
         * non-aggregators.
         */
        char *ptr = rd_buf;
        if (npairs > 0) {
            memcpy(send_buf + bufAddr[0], ptr, lengths[0]);
            ptr += lengths[0];
        }
        for (i=0, j=1; j<npairs; j++) {
            memcpy(send_buf + bufAddr[j], ptr, lengths[j]);
            ptr += lengths[j];
        }
        NCI_Free(bufAddr);
        NCI_Free(rd_buf);

        /* First, for self rank, unpack read data into self's buf */
        if (bufLen > 0) {
            int is_predef;
            if (bufType == MPI_DATATYPE_NULL)
                is_predef = 0;
            else
                ADIOI_Type_ispredef(bufType, &is_predef);

            if (is_predef)
                memcpy(buf, send_buf, bufLen);
            else {
                void *outbuf = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count position=0;
                MPI_Count outcount = (buf == NULL) ? 1 : bufCount;
                MPI_Unpack_c(send_buf, bufLen, &position, outbuf, outcount,
                             bufType, MPI_COMM_SELF);
#else
                int position=0;
                int outcount = (buf == NULL) ? 1 : bufCount;
                MPI_Unpack(send_buf, bufLen, &position, outbuf, outcount,
                           bufType, MPI_COMM_SELF);
#endif
            }
        }

        /* send read data to non-aggregators */
        ptr = send_buf + bufLen;
        nreqs = 0;
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (meta[i*3 + 1] == 0) continue;
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Isend_c(ptr, meta[i*3 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                        0, ncp->comm, &req[nreqs++]);
#else
            MPI_Isend(ptr, meta[i*3 + 1], MPI_BYTE, ncp->nonaggr_ranks[i],
                      0, ncp->comm, &req[nreqs++]);
#endif
            ptr += meta[i*3 + 1];
        }
#ifdef HAVE_MPI_STATUSES_IGNORE
        mpireturn = MPI_Waitall(nreqs, req, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               ADIOI_Malloc(nreqs * sizeof(MPI_Status));
        mpireturn = MPI_Waitall(nreqs, req, statuses);
        ADIOI_Free(statuses);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        NCI_Free(send_buf);
    }

#ifdef PNETCDF_PROFILING
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[4] += endT - startT;
#endif

fn_exit:
    if (meta != NULL) NCI_Free(meta);
    if (req != NULL) NCI_Free(req);

    return status;
}

/*----< ncmpio_intra_node_aggregation_nreqs() >------------------------------*/
/* This subroutine is a collective call, and is to handle non-blocking APIs,
 * which contain requests to one or more variable.
 * reqMode: NC_REQ_RD for read request and NC_REQ_WR for write.
 * num_reqs: number of elements in array req_list.
 * req_list[]: stores pending requests from non-blocking API calls, which is
 *     used to construct file offset-length pairs and user buffer datatype.
 * newnumrecs: number of new records
 */
int
ncmpio_intra_node_aggregation_nreqs(NC         *ncp,
                                    int         reqMode,
                                    int         num_reqs,
                                    NC_req     *req_list,
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

#ifdef WKL_DEBUG
MPI_Offset maxm;
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

    /* ensure the intra-node aggregation is enabled in this rank's group */
    assert(ncp->my_aggr >= 0);

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     *     is_incr:   whether offsets are incremental
     */
    if (num_reqs > 0)
        flatten_reqs(ncp, num_reqs, req_list, &is_incr, &num_pairs, &offsets,
                     &lengths);
    else
        num_pairs = 0;

    /* construct user buffer datatype, bufType.
     * bufLen is the buffer size in bytes
     */
    if (num_reqs > 0) {
        construct_buf_type(ncp, num_reqs, req_list, &bufLen, &bufType);
        bufLen = 1;
    }
    else
        bufLen = 0;

    if (req_list != NULL)
        /* All metadata in req_list have been used to construct bufType and
         * bufLen. It is now safe to release the space occupied by req_list.
         */
        NCI_Free(req_list);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (ncp->rank == ncp->my_aggr) ncp->aggr_time[1] += MPI_Wtime() - timing;
#endif

    /* perform intra-node aggregation */
    if (fIsSet(reqMode, NC_REQ_WR)) {
        err = ina_put(ncp, is_incr, num_pairs, offsets, lengths,
                      bufLen, bufType, NULL);
        if (status == NC_NOERR) status = err;
    }
    else {
        err = ina_get(ncp, is_incr, num_pairs, offsets, lengths,
                      bufLen, bufType, NULL);
        if (status == NC_NOERR) status = err;
    }

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
/* This subroutine is a collective call, and is to handle blocking APIs, which
 * contain requests to only one variable.
 * reqMode: NC_REQ_RD for read request and NC_REQ_WR for write.
 * varp: pointer to the variable struct.
 * start[]: starting offsets
 * count[]: counts along each dimension
 * stride[]: stride along each dimension
 * bufCount: number of bufType in the user buffer, pointed by buf
 * bufType: MPI derived datatype describing user buffer layout
 * buf: pointer to the user buffer
 */
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

#ifdef WKL_DEBUG
MPI_Offset maxm;
ncmpi_inq_malloc_max_size(&maxm); if (ncp->rank == 0)  printf("xxxx %s line %4d: maxm=%.2f MB\n",__func__,__LINE__,(float)maxm/1048576.0);
#endif

    /* ensure the intra-node aggregation is enabled in this rank's group */
    assert(ncp->my_aggr >= 0);

    if (buf == NULL) { /* zero-length request */
        /* This rank must still participate intra_node_aggregation() to tell
         * its aggregator that it has no I/O data.
         */
        if (fIsSet(reqMode, NC_REQ_WR))
            return ina_put(ncp, 1, 0, NULL, NULL, 0, MPI_BYTE, NULL);
        else
            return ina_get(ncp, 1, 0, NULL, NULL, 0, MPI_BYTE, NULL);
    }

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     *     is_incr:   whether offsets are incremental
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

    /* perform intra-node aggregation */
    if (fIsSet(reqMode, NC_REQ_WR)) {
        err = ina_put(ncp, is_incr, num_pairs, offsets, lengths,
                      bufCount, bufType, buf);
        if (status == NC_NOERR) status = err;
    }
    else {
        err = ina_get(ncp, is_incr, num_pairs, offsets, lengths,
                      bufCount, bufType, buf);
        if (status == NC_NOERR) status = err;
    }

    return status;
}

