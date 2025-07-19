/*
 *  Copyright (C) 2024, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 * This file contains the implementation of intra-node aggregation feature,
 * which is designed to improve performance for I/O patterns that contain many
 * noncontiguous requests interleaved among processes, with a wide aggregate
 * access region on each process that involves file stripes responsible by
 * almost all the file servers. By reducing the number of processes per node
 * to participate MPI-IO operations, this feature can effectively reduce the
 * communication contention, particularly often happened to jobs that run a
 * large the number of MPI processes per compute node.
 *
 * Users can enable this feature by setting the PnetCDF I/O hint named
 * 'nc_num_aggrs_per_node' to a positive integral value, indicating the desired
 * number of processes per compute node to be selected as the intra-node I/O
 * aggregators. Processes running on the same node are divided into groups.
 * The process with the lowest rank ID is selected as the I/O aggregator of
 * that group. Non-aggregators send their requests to their aggregators, and
 * then the aggregators make I/O requests to the file, i.e. only aggregators
 * make MPI-IO calls.
 *
 * Because communication within a node can be achieved by memory copy operation
 * and thus its cost is expected to be much lower than the inter-node
 * communication, this feature can effectively reduce the communication
 * congestion or exhaustion of message queues, due to many pending asynchronous
 * messages produced in the two-phase I/O, the strategy used to implement
 * MPI collective I/O.
 *
 * The concept of intra-node request aggregation and its performance results
 * are presented in the following paper.
 * Q. Kang, S. Lee, K. Hou, R. Ross, A. Agrawal, A. Choudhary, and W. Liao.
 * Improving MPI Collective I/O for High Volume Non-Contiguous Requests With
 * Intra-Node Aggregation. IEEE Transactions on Parallel and Distributed
 * Systems, 31(11):2682-2695, November 2020.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
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

/* swap values of x and y */
#define SWAP1(x, y, tmp) { tmp = x ; x = y; y = tmp ; }

#ifdef HAVE_MPI_LARGE_COUNT
/* swap elements of arrays x, y, and corresponding lengths and bufAddr */
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    MPI_Count aint; \
    MPI_Count cint; \
    MPI_Count d0 = (x) - offsets; \
    MPI_Count d1 = (y) - offsets; \
    if (d0 != d1) { \
        SWAP1(*(x), *(y), cint); \
        SWAP1(lengths[d0], lengths[d1], cint); \
        if (bufAddr != NULL) \
            SWAP1(bufAddr[d0], bufAddr[d1], aint); \
    } \
}
#else
#define SWAP(offsets, lengths, bufAddr, x, y) { \
    int int4; \
    MPI_Offset aint; \
    MPI_Offset d0 = (x) - offsets; \
    MPI_Offset d1 = (y) - offsets; \
    if (d0 != d1) { \
        SWAP1(*(x), *(y), aint); \
        SWAP1(lengths[d0], lengths[d1], int4); \
        if (bufAddr != NULL) \
            SWAP1(bufAddr[d0], bufAddr[d1], aint); \
    } \
}
#endif

#define MEDIAN(a,b,c) ((*(a) < *(b)) ? \
                      ((*(b) < *(c)) ? (b) : ((*(a) < *(c)) ? (c) : (a))) : \
                      ((*(b) > *(c)) ? (b) : ((*(a) < *(c)) ? (a) : (c))))

/*----< qsort_off_len_buf() >------------------------------------------------*/
/* Sort three arrays of offsets, lengths, and buffer addresses based on array
 * offsets into an increasing order. This code is based on the qsort routine
 * from Bentley & McIlroy's "Engineering a Sort Function".
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
            if (bufAddr != NULL)
                bufAddr = bufAddr + (num - r);
            offsets = pn - r;
            num = r;
        }
        else
            break;
    }
}

/*----< heap_merge() >-------------------------------------------------------*/
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

    /* This heap_merge is not in-place, taking too much memory footprint */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *srt_off = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
    MPI_Count *srt_len = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * nelems);
#else
    MPI_Aint *srt_off = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nelems);
    int      *srt_len = (int*)      NCI_Malloc(sizeof(int)      * nelems);
#endif
    MPI_Aint *srt_addr = NULL;

    if (bufAddr != NULL)
        srt_addr = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nelems);

    a = (heap_struct *) NCI_Calloc(nprocs, sizeof(heap_struct));

    /* there are nprocs number of lists to be merged */
    j = 0;
    sum = 0;
    for (i = 0; i < nprocs; i++) {
        if (count[i]) {
            /* each of a[j].off_list is already sorted */
            a[j].off_list = offsets + sum;
            a[j].len_list = blklens + sum;
            if (bufAddr != NULL)
                a[j].addr_list = bufAddr + sum;
            sum += count[i];
            a[j].count = count[i];
            j++;
        }
    }
    nprocs = j; /* some count[i] may be zero */

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
        if (bufAddr != NULL)
            srt_addr[i] = a[0].addr_list[0];
        a[0].count--;

        if (!a[0].count) {
            a[0] = a[heapsize - 1];
            heapsize--;
        } else {
            a[0].off_list++;
            a[0].len_list++;
            if (bufAddr != NULL)
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
    if (bufAddr != NULL)
        memcpy(bufAddr, srt_addr, sizeof(MPI_Aint) * nelems);

    NCI_Free(a);
    if (bufAddr != NULL) NCI_Free(srt_addr);
    NCI_Free(srt_len);
    NCI_Free(srt_off);
}


/*----< ncmpio_intra_node_aggr_init() >--------------------------------------*/
/* When intra-node write aggregation is enabled, this subroutine initializes
 * the metadata to be used for intra-node communication and I/O requests.
 *
 * Processes on the same node will first be divided into groups. A process with
 * the lowest rank ID in a group is selected as the aggregator. Only the
 * aggregators call the MPI-IO functions to perform I/O to the file. Thus, this
 * subroutine must be called before MPI_File_open() and should be called only
 * once at ncmpio_create() or ncmpio_open().
 *
 * The subroutine performs the following tasks.
 * 1. Make use of the affinity of each MPI process to its compute node,
 *    represented by ncp->num_nodes and ncp->node_ids[]. These two member of
 *    ncp should have been set from a call to ncmpii_construct_node_list()
 *    earlier during ncmpio_create() and ncmpio_open().
 *    + ncp->num_nodes is the number of unique compute nodes.
 *    + ncp->node_ids[ncp->nprocs] contains node IDs for all processes.
 * 2. Divide processes into groups, select aggregators, and determine whether
 *    self process is an intra-node aggregator.
 *    + ncp->my_aggr is rank ID of my aggregator.
 *    + if (ncp->my_aggr == ncp->rank) then this rank is an aggregator.
 * 3. For an aggregator, find the number of non-aggregators assigned to it and
 *    construct a list of rank IDs of non-aggregators of its group.
 *    + ncp->num_nonaggrs is the number of non-aggregators in its group.
 * 4. For a non-aggregator, find the rank ID of its assigned aggregator.
 *    + ncp->my_aggr is rank ID of my aggregator.
 *    + ncp->nonaggr_ranks[] contains the rank IDs of assigned non-aggregators.
 * 5. Create a new MPI communicator consisting of only the aggregators only.
 *    Obtain the rank ID and total process number of the new communicator.
 *    + ncp->ina_comm contains the aggregators across all nodes.
 *    + ncp->ina_nprocs is the number of processes in intra-node communicator.
 *    + ncp->ina_rank is this process's rank ID in intra-node communicator.
 */
int
ncmpio_intra_node_aggr_init(NC *ncp)
{
    int i, j, mpireturn, do_io, ina_nprocs, naggrs_my_node, first_rank;
    int my_rank_index, *ranks_my_node, my_node_id, nprocs_my_node;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
    int nelems = sizeof(ncp->ina_time_put) / sizeof(ncp->ina_time_put[0]);
    ncp->ina_time_init = ncp->ina_time_flatten = 0.0;
    for (i=0; i<nelems; i++) {
        ncp->ina_time_put[i] = ncp->ina_time_get[i] = 0;
        ncp->maxmem_put[i] = ncp->maxmem_get[i] = 0;
    }
    ncp->ina_npairs_put = ncp->ina_npairs_get = 0;
#endif

    /* initialize parameters of intra-node aggregation */
    ncp->my_aggr = -1;         /* rank ID of my aggregator */
    ncp->num_nonaggrs = 0;     /* number of non-aggregators assigned */
    ncp->nonaggr_ranks = NULL; /* ranks of assigned non-aggregators */

    /* Note that ill value of ncp->num_aggrs_per_node has been checked before
     * entering this subroutine. Thus ncp->num_aggrs_per_node must be > 0.
     */

    /* ncp->node_ids[] has been established in ncmpii_construct_node_list()
     * called in ncmpio_create() or ncmpio_open() before entering this
     * subroutine. my_node_id is this rank's node ID.
     */
    my_node_id = ncp->node_ids[ncp->rank];

    /* nprocs_my_node:  the number of processes in my nodes
     * ranks_my_node[]: rank IDs of all processes in my node.
     * my_rank_index:   points to ranks_my_node[] where
     *                  ranks_my_node[my_rank_index] == ncp->rank
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
     * less than others.
     */
    naggrs_my_node = MIN(ncp->num_aggrs_per_node, nprocs_my_node);

    /* For each aggregation group, calculate the number of non-aggregators,
     * ncp->num_nonaggrs. Note ncp->num_nonaggrs includes self rank.
     */
    ncp->num_nonaggrs = nprocs_my_node / naggrs_my_node;
    if (nprocs_my_node % naggrs_my_node) ncp->num_nonaggrs++;

    /* Adjust the number of non-aggregators for the last group of each node,
     * to make sure it does not go beyond nprocs_my_node.
     */
    first_rank = my_rank_index - my_rank_index % ncp->num_nonaggrs;
    ncp->num_nonaggrs = MIN(ncp->num_nonaggrs, nprocs_my_node - first_rank);

    /* Assign the first rank as the intra-node aggregator of this group and
     * set the rank ID of my aggregator for each process.
     */
    ncp->my_aggr = ranks_my_node[first_rank];

    if (ncp->num_nonaggrs == 1) {
        /* When the number of processes in this group is 1, the aggregation
         * is not performed. Note num_nonaggrs includes self rank.
         *
         * Note this does not mean intra-node aggregation is disabled. The
         * indicator of whether intra-node aggregation is enabled or disabled
         * is ncp->num_aggrs_per_node, whose value should be consistent across
         * all processes. It is possible for some groups containing only one
         * process, in which the aggregation is not necessarily performed
         * within that group.
         */
        assert(ncp->my_aggr == ncp->rank);
    }
    else if (ncp->my_aggr == ncp->rank) { /* ncp->num_nonaggrs > 1 */
        /* Construct ncp->nonaggr_ranks[], the rank IDs of non-aggregators of
         * this group. Note ncp->nonaggr_ranks[], if malloc-ed, will only be
         * freed when closing the file.
         */
        ncp->nonaggr_ranks = (int*)NCI_Malloc(sizeof(int) * ncp->num_nonaggrs);

        memcpy(ncp->nonaggr_ranks, ranks_my_node + first_rank,
               sizeof(int) * ncp->num_nonaggrs);
    }
    NCI_Free(ranks_my_node);

    /* Next step is to construct a new MPI communicator consisting of all
     * intra-node aggregators. It will later be used to call MPI_File_open(),
     * so that only aggregators call MPI-IO functions to access the file.
     *
     * When using the PnetCDF's internal ADIO driver, we can pass a list of
     * node_ids of the new communicator to the ADIO file handler, ncp->adio_fh,
     * so to prevent the driver from the repeated work of constructing the list
     * of node IDs, node_ids. If using MPI-IO driver, then ROMIO will do this
     * internally again anyway.
     */

    do_io = (ncp->my_aggr == ncp->rank) ? 1 : 0;

    /* construct an array containing ranks of aggregators */
    ncp->ina_node_list = (int*) NCI_Malloc(sizeof(int) * ncp->nprocs);
    TRACE_COMM(MPI_Allgather)(&do_io, 1, MPI_INT, ncp->ina_node_list, 1,
                              MPI_INT,ncp->comm);

    /* Calculate the total number of intra-node aggregators */
    for (ina_nprocs=0, i=0; i<ncp->nprocs; i++)
        if (ncp->ina_node_list[i]) ina_nprocs++;

    /* Construct ncp->node_ids[] and ncp->ina_node_list[]. Their contents
     * depend on the layout of MPI process allocation to the compute nodes.
     * The common layouts can be two kinds:
     *   + cyclic - MPI ranks are assigned to nodes round-robin-ly,
     *   + block - MPI ranks are assigned to a node and then move on to next.
     *
     * Below uses an example of nodes=3, nprocs=10, * num_aggrs_per_node=2.
     * ncp->node_ids[] should be
     *     block  process allocation: 0,0,0,0,1,1,1,2,2,2
     *     cyclic process allocation: 0,1,2,0,1,2,0,1,2,0
     * Accordingly, ncp->ina_node_list[] can be two kinds
     *     block  process allocation: 1,0,1,0,1,0,1,1,0,1
     *     cyclic process allocation: 1,1,1,0,0,0,1,1,1,0
     */

    /* ncp->node_ids[]: node IDs of processes in the new MPI communicator.
     * ncp->ina_node_list[]: the rank IDs of the new MPI communicator.
     */
    for (j=0,i=0; i<ncp->nprocs; i++) {
        if (ncp->ina_node_list[i]) {
            ncp->ina_node_list[j] = i;
            /* Modify ncp->node_ids[] to store the node IDs of the processes in
             * the new communicator. Note ncp->node_ids[] from now on is used
             * by PnetCDF's ADIO driver only.
             */
            ncp->node_ids[j] = ncp->node_ids[i];
            j++;
        }
    }

    /* Make MPI calls to create a new communicator. */
    MPI_Group origin_group, ina_group;
    TRACE_COMM(MPI_Comm_group)(ncp->comm, &origin_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_group");
    TRACE_COMM(MPI_Group_incl)(origin_group, ina_nprocs, ncp->ina_node_list, &ina_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_incl");
    TRACE_COMM(MPI_Comm_create)(ncp->comm, ina_group, &ncp->ina_comm);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Comm_create");
    TRACE_COMM(MPI_Group_free)(&ina_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_free");
    TRACE_COMM(MPI_Group_free)(&origin_group);
    if (mpireturn != MPI_SUCCESS)
        return ncmpii_error_mpi2nc(mpireturn, "MPI_Group_free");

    /* Non-aggregators will have ncp->ina_comm set to MPI_COMM_NULL */
    if (ncp->ina_comm == MPI_COMM_NULL) {
        ncp->ina_nprocs = 0;
        ncp->ina_rank = -1;
    }
    else {
        MPI_Comm_size(ncp->ina_comm, &ncp->ina_nprocs);
        MPI_Comm_rank(ncp->ina_comm, &ncp->ina_rank);
    }

    /* TODO: automatically determine whether or not to enable intra-node
     * aggregation.
     *
     * The ideal case is it can be determined right before each collective
     * write call, because only at that time, the communication pattern is
     * known. If the pattern can cause contention, then enable it. Otherwise,
     * disable it.
     *
     * Such mechanism may depends on the followings.
     *   1. MPI-IO hint cb_noddes, and striping_unit
     *   2. calculate aggregate access region
     *   3. If the number of senders to each cb_nodes is very large, then
     *      intra-node aggregation should be enabled.
     *   4. Average of nprocs_per_node across all processes may be a factor for
     *      determining whether to enable intra-node aggregation. It indicates
     *      whether the high number of processes are allocated on the same
     *      node.
     */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncp->ina_time_init = MPI_Wtime() - timing;
#endif

    return NC_NOERR;
}

/*----< flatten_subarray() >-------------------------------------------------*/
/* Flatten a subarray request, specified by start[], count[], and stride[] into
 * a list of file offset-length pairs, offsets[] and lengths[].
 */
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

    /* check if the list can be coalesced */
    for (i=0, j=1; j<*npairs; j++) {
        if (offsets[i] + lengths[i] == offsets[j])
            lengths[i] += lengths[j];
        else {
            i++;
            if (i < j) {
                offsets[i] = offsets[j];
                lengths[i] = lengths[j];
            }
        }
    }
    *npairs = i + 1;

    return NC_NOERR;
}

/*----< flatten_req() >------------------------------------------------------*/
/* Flatten one subarray request into offset-length pairs. Arrays offsets and
 * lengths are allocated in this subroutine and need to be freed by the caller.
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

    /* num_pairs may be less than originally calculated, because offset-length
     * pairs are coalesced in the call to flatten_subarray().
     */
    *num_pairs = idx;

    return err;
}

/*----< flatten_reqs() >-----------------------------------------------------*/
/* Flatten multiple subarray requests into file offset-length pairs. Arrays
 * offsets and lengths are allocated here and need to be freed by the caller.
 */
static int
flatten_reqs(NC            *ncp,
             int            reqMode,   /* IN: NC_REQ_RD or NC_REQ_WR */
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
    for (i=0; i<num_reqs; i++) {
        /* reqs[i].npairs is the number of offset-length pairs of this request,
         * calculated in ncmpio_igetput_varm() and igetput_varn()
         */
        *num_pairs += reqs[i].npairs;
        if (fIsSet(reqMode, NC_REQ_WR))
            ndims = ncp->put_lead_list[reqs[i].lead_off].varp->ndims;
        else
            ndims = ncp->get_lead_list[reqs[i].lead_off].varp->ndims;
        max_ndims = MAX(max_ndims, ndims);
    }

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
        NC_lead_req *lead;
        if (fIsSet(reqMode, NC_REQ_WR))
            lead = ncp->put_lead_list + reqs[i].lead_off;
        else
            lead = ncp->get_lead_list + reqs[i].lead_off;

        if (reqs[i].npairs == 1) {
            /* When reqs[i] contains only one offset-length pair, re-use
             * reqs[i].offset_start, which has been generated earlier at a call
             * to ncmpio_intra_node_aggregation_nreqs().
             */
            (*offsets)[idx] = reqs[i].offset_start;
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

    /* num_pairs may be less than originally calculated, because offset-length
     * pairs are coalesced in the call to flatten_subarray().
     */
    *num_pairs = idx;

    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead;
        if (fIsSet(reqMode, NC_REQ_WR))
            lead = ncp->put_lead_list + reqs[i].lead_off;
        else
            lead = ncp->get_lead_list + reqs[i].lead_off;
        if (fIsSet(lead->flag, NC_REQ_TO_FREE)) {
            NCI_Free(lead->start);
            lead->start = NULL;
        }
    }

    return status;
}

/*----< construct_buf_type() >-----------------------------------------------*/
/* Construct an MPI derived datatype for I/O buffers from multiple requests
 * described in the request list, reqs, by concatenate memory addresses of all
 * buffers.
 */
static int
construct_buf_type(const NC     *ncp,
                   int           reqMode,   /* IN: NC_REQ_RD or NC_REQ_WR */
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
        if (fIsSet(reqMode, NC_REQ_WR))
            lead = ncp->put_lead_list + reqs[i].lead_off;
        else
            lead = ncp->get_lead_list + reqs[i].lead_off;
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
 * metadata from the non-aggregators into meta, including:
 *   1. the number of offset-length pairs on each non-aggregator
 *   2. offsets array of each non-aggregator
 *   3. lengths array of each non-aggregator
 *   4. npairs is the total number of offset-length pairs of this group.
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
    if (ncp->my_aggr == ncp->rank) {

        req = (MPI_Request*)NCI_Malloc(sizeof(MPI_Request) * ncp->num_nonaggrs);
        nreqs = 0;
        for (i=1; i<ncp->num_nonaggrs; i++)
            TRACE_COMM(MPI_Irecv)(meta + i*3, 3, MPI_AINT,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);

#ifdef HAVE_MPI_STATUSES_IGNORE
        TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               ADIOI_Malloc(nreqs * sizeof(MPI_Status));
        TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
        ADIOI_Free(statuses);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
    }
    else /* non-aggregator */
        TRACE_COMM(MPI_Send)(meta, 3, MPI_AINT, ncp->my_aggr, 0, ncp->comm);

    /* Secondly, aggregators collect offset-length pairs from all its
     * non-aggregators
     */
    if (ncp->my_aggr == ncp->rank) {
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
            *lengths = (int*)        NCI_Realloc(*lengths, *npairs * sizeof(int));
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
            TRACE_COMM(MPI_Irecv_c)(MPI_BOTTOM, 1, recvType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvType);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#else
        int bklens[2];
        MPI_Aint aint, disps[2];

        MPI_Get_address(*offsets, &aint);
        disps[0] = MPI_Aint_add(aint, sizeof(MPI_Offset) * meta[0]);
        MPI_Get_address(*lengths, &aint);
        disps[1] = MPI_Aint_add(aint, sizeof(int) * meta[0]);
        for (i=1; i<ncp->num_nonaggrs; i++) {
            if (meta[i*3] == 0) continue;
            bklens[0] = meta[i*3] * sizeof(MPI_Offset);
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
            TRACE_COMM(MPI_Irecv)(MPI_BOTTOM, 1, recvType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs]);
            MPI_Type_free(&recvType);

            disps[0] = MPI_Aint_add(disps[0], bklens[0]);
            disps[1] = MPI_Aint_add(disps[1], bklens[1]);
            nreqs++;
        }
#endif
#ifdef HAVE_MPI_STATUSES_IGNORE
        TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               ADIOI_Malloc(nreqs * sizeof(MPI_Status));
        TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
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
        TRACE_COMM(MPI_Send_c)(MPI_BOTTOM, 1, sendType, ncp->my_aggr, 0,
                               ncp->comm);
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
        TRACE_COMM(MPI_Send)(MPI_BOTTOM, 1, sendType, ncp->my_aggr, 0,
                             ncp->comm);
        MPI_Type_free(&sendType);
#endif
    }

    if (req != NULL) NCI_Free(req);

    return status;
}

/*----< ina_put() >----------------------------------------------------------*/
/* This collective subroutine implements the intra-node aggregation for write
 * operations. Arrays offsets and lengths will be freed when return.
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

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double endT, startT = MPI_Wtime();
    MPI_Offset mem_max;
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_put[0] = MAX(ncp->maxmem_put[0], mem_max);
#endif

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufLen;
    MPI_Type_size_c(bufType, &bufLen);
#else
    int bufLen;
    MPI_Type_size(bufType, &bufLen);
#endif
    bufLen *= bufCount;

    /* Firstly, aggregators collect metadata from non-aggregators.
     *
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

    /* Each aggregator first collects metadata about its offset-length pairs,
     * size of write request, and whether the offsets are in an incremental
     * order. The aggregator will gather these metadata from non-aggregators
     * assigned to it.
     * For write operation, keeping the original offset-length pairs is not
     * necessary, as they will later be sorted and coalesced before calling
     * MPI-IO.
     */
    err = ina_collect_md(ncp, meta, &offsets, &lengths, &npairs);
    if (err != NC_NOERR) {
        if (req != NULL) NCI_Free(req);
        NCI_Free(meta);
        return err;
    }

    /* For write operation, the non-aggregators now can start sending their
     * write data to the aggregator.
     */
    if (ncp->rank != ncp->my_aggr && num_pairs > 0) { /* non-aggregator */
        /* Non-aggregators send write data to the aggregator */
        void *buf_ptr = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count num = (buf == NULL) ? 1 : bufCount;
        TRACE_COMM(MPI_Send_c)(buf_ptr, num, bufType, ncp->my_aggr, 0,
                               ncp->comm);
#else
        int num = (buf == NULL) ? 1 : bufCount;
        TRACE_COMM(MPI_Send)(buf_ptr, num, bufType, ncp->my_aggr, 0,
                             ncp->comm);
#endif
        NCI_Free(offsets);
        NCI_Free(lengths);
        offsets = NULL;
        lengths = NULL;
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_put[1] = MAX(ncp->maxmem_put[1], mem_max);
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_put[0] += endT - startT;
    startT = endT;
#endif

    /* MPI-IO has the following requirements about filetype.
     * 1. The (flattened) displacements (of a filetype) are not required to be
     *    distinct, but they cannot be negative, and they must be monotonically
     *    non-decreasing.
     * 2. If the file is opened for writing, neither the etype nor the filetype
     *    is permitted to contain overlapping regions.
     */
    if (ncp->rank == ncp->my_aggr && npairs > 0) {
        /* Now the aggregator has received all offset-length pairs from its
         * non-aggregators. At first, check if qsort is necessary. Find the
         * first non-aggregator with non-zero pairs whose offsets are not
         * already sorted.
         */
        int do_sort=0, indv_sorted=1;

        for (i=-1,j=0; j<ncp->num_nonaggrs; j++) {
            if (i == -1 && meta[j*3] > 0) /* find 1st whose num_pairs > 0 */
                i = j;
            if (meta[j*3+2] == 0) { /* j's offsets are not sorted */
                indv_sorted = 0;
                do_sort = 1;
                break;
            }
        }
        /* i is the first non-aggregator whose num_pairs > 0, and
         * j is the first non-aggregator whose is_incr is false
         */

        if (i >= 0 && indv_sorted == 1) {
            /* This is when all non-aggregators' offsets are individually
             * sorted. We still need to check if offsets are interleaved among
             * all non-aggregators to determine whether a sorting of all
             * offset-length pairs is necessary.
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
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        ncmpi_inq_malloc_size(&mem_max);
        ncp->maxmem_put[2] = MAX(ncp->maxmem_put[2], mem_max);

        endT = MPI_Wtime();
        if (ncp->rank == ncp->my_aggr) {
            ncp->ina_time_put[1] += endT - startT;
            ncp->ina_npairs_put = MAX(ncp->ina_npairs_put, npairs);
        }
        startT = endT;
#endif

        /* calculate the total amount to be written by this aggregator */
        for (buf_count=0,i=0; i<ncp->num_nonaggrs; i++)
            buf_count += meta[i*3 + 1];

        /* Allocate receive buffer and write buffer. Once write data from
         * non-aggregators have received into recv_buf, recv_buf is packed into
         * wr_buf. Then, wr_buf is used when calling MPI-IO to write to file.
         */
        if (buf_count > 0) {
            recv_buf = (char*) NCI_Malloc(buf_count);
            wr_buf = (char*) NCI_Malloc(buf_count);
        }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        ncmpi_inq_malloc_size(&mem_max);
        ncp->maxmem_put[3] = MAX(ncp->maxmem_put[3], mem_max);
#endif

        /* First, pack self write data into front of the recv_buf */
        if (bufLen > 0) {
            int is_predef;
            if (bufType == MPI_DATATYPE_NULL)
                is_predef = 0;
            else
                PNCIO_Type_ispredef(bufType, &is_predef);

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

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        endT = MPI_Wtime();
        if (ncp->rank == ncp->my_aggr)
            ncp->ina_time_put[2] += endT - startT;
        startT = endT;
#endif

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
                TRACE_COMM(MPI_Irecv_c)(ptr, meta[i*3 + 1], MPI_BYTE,
                           ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#else
                TRACE_COMM(MPI_Irecv)(ptr, meta[i*3 + 1], MPI_BYTE,
                           ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#endif
                ptr += meta[i*3 + 1];
            }
#ifdef HAVE_MPI_STATUSES_IGNORE
            TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
            MPI_Status *statuses = (MPI_Status *)
                                   ADIOI_Malloc(nreqs * sizeof(MPI_Status));
            TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
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
    } /* if (ncp->rank == ncp->my_aggr && npairs > 0) */

    NCI_Free(meta);
    if (ncp->rank == ncp->my_aggr)
        NCI_Free(req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_put[3] += endT - startT;
#endif

    /* Only aggregators call MPI-IO functions to write data to the file.
     * Non-aggregators do not participate MPI-IO calls.
     */
    if (ncp->rank != ncp->my_aggr)
        return status;

    /* intra-node aggregation only takes effect in collective data mode */
    fh = ncp->collective_fh;

    /* Preserve fd->filetype previously used */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        disp = 0;
        saved_fileType = ncp->adio_fh->filetype;

        /* When using PnetCDF's internal ADIO driver, a call to
         * ncmpio_file_set_view() will pass offsets[] and lengths[] to
         * PNCIO_File_set_view() to be reused as flattened fileType. This can
         * avoid re-create an MPI datatype and re-flattening it. See such
         * adaptation in ADIOI_Calc_my_off_len()
         *
         * Here we use MPI_DATATYPE_NULL to tell PNCIO_File_set_view() that it
         * is called from the intra-node aggregation subroutine, so the ADIO
         * subroutines can reuse offsets[] and lengths[].
         */
        fileType = MPI_DATATYPE_NULL;
    }

    /* Set the MPI-IO fileview (this is a collective call). This call returns
     * disp which points to the very first file offset to be written.
     */
    err = ncmpio_file_set_view(ncp, fh, &disp, fileType, npairs, offsets,
                               lengths);
    if (err != NC_NOERR) {
        if (status == NC_NOERR) status = err;
        buf_count = 0;
    }
    if (ncp->fstype == ADIO_FSTYPE_MPIIO && fileType != MPI_BYTE)
        MPI_Type_free(&fileType);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_put[4] = MAX(ncp->maxmem_put[4], mem_max);
#endif

    /* call MPI_File_write_at_all or PNCIO_File_write_at_all */
    err = ncmpio_read_write(ncp, NC_REQ_WR, NC_REQ_COLL, disp, buf_count,
                            MPI_BYTE, wr_buf, 1);
    if (status == NC_NOERR) status = err;

    if (wr_buf  != NULL) NCI_Free(wr_buf);
    if (offsets != NULL) NCI_Free(offsets);
    if (lengths != NULL) NCI_Free(lengths);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_put[5] = MAX(ncp->maxmem_put[5], mem_max);
#endif
    /* ncp->adio_fh->flat_file is allocated in ncmpio_file_set_view() */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        NCI_Free(ncp->adio_fh->flat_file);
        ncp->adio_fh->flat_file = NULL;

        /* restore the original filetype */
        ncp->adio_fh->filetype = saved_fileType;
    }

    return status;
}

static
size_t bin_search(
#ifdef HAVE_MPI_LARGE_COUNT
                  MPI_Count key, MPI_Count *base,
#else
                  MPI_Offset key, MPI_Offset *base,
#endif
                  size_t nmemb)
{
    size_t low, high;

    /* only one element */
    if (nmemb == 1)
        return (base[0] <= key) ? 0 : -1;

    /* check the 1st emelemt */
    if (base[0] <= key && key < base[1])
        return 0;

    low = 1;
    high = nmemb - 1;

    while (low <= high) {
        size_t mid = low + (high - low) / 2;
        if (base[mid] == key)
            return mid;
        if (base[mid] < key)
            low = mid + 1;
        else
            high = mid - 1;
    }
    return (low - 1);
}

/*----< ina_get() >----------------------------------------------------------*/
/* This collective subroutine implements the intra-node aggregation for read
 * operations. Arrays offsets and lengths will be freed when return.
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
    int do_sort=0, indv_sorted=1;
    char *rd_buf = NULL;
    MPI_Aint npairs=0, max_npairs, *meta=NULL, *count=NULL;
    MPI_Offset disp=0, buf_count=0;
    MPI_Datatype saved_fileType, fileType=MPI_BYTE;
    MPI_File fh;
    MPI_Request *req=NULL;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double endT, startT = MPI_Wtime();
    MPI_Offset mem_max;
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_get[0] = MAX(ncp->maxmem_get[0], mem_max);
#endif

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufLen, *orig_offsets=NULL, *orig_lengths=NULL;
    MPI_Type_size_c(bufType, &bufLen);
#else
    MPI_Offset *orig_offsets=NULL;
    int bufLen, *orig_lengths=NULL;
    MPI_Type_size(bufType, &bufLen);
#endif
    bufLen *= bufCount;

    /* Firstly, aggregators collect metadata from non-aggregators.
     *
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

    /* Each aggregator first collects metadata about its offset-length pairs,
     * size of read request, and whether the offsets are in an incremental
     * order. The aggregator will gather these metadata from non-aggregators
     * assigned to it.
     * For read operation, the original offsets and lengths must be kept,
     * because the sorting and coalescing later will mess up the original order
     * of offsets and lengths. The original ones are needed to construct a
     * datatype when an aggregator sends read data to its non-aggregators.
     */
    err = ina_collect_md(ncp, meta, &offsets, &lengths, &npairs);
    if (err != NC_NOERR) {
        if (req != NULL) NCI_Free(req);
        NCI_Free(meta);
        return err;
    }

    if (ncp->rank == ncp->my_aggr) {
        /* duplicate to keep a copy of the original offset-length pairs */
#ifdef HAVE_MPI_LARGE_COUNT
        orig_offsets = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
        orig_lengths = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * npairs);
        memcpy(orig_offsets, offsets, sizeof(MPI_Count) * npairs);
        memcpy(orig_lengths, lengths, sizeof(MPI_Count) * npairs);
#else
        orig_offsets = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * npairs);
        orig_lengths = (int*)        NCI_Malloc(sizeof(int) * npairs);
        memcpy(orig_offsets, offsets, sizeof(MPI_Offset) * npairs);
        memcpy(orig_lengths, lengths, sizeof(int) * npairs);
#endif
    }

    if (ncp->rank != ncp->my_aggr && num_pairs > 0) {
        /* For read operation, the non-aggregators now can start receiving
         * their read data from the aggregator.
         */
        MPI_Status st;
        void *buf_ptr = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count num = (buf == NULL) ? 1 : bufCount;
        TRACE_COMM(MPI_Recv_c)(buf_ptr, num, bufType, ncp->my_aggr, 0,
                               ncp->comm, &st);
#else
        int num = (buf == NULL) ? 1 : bufCount;
        TRACE_COMM(MPI_Recv)(buf_ptr, num, bufType, ncp->my_aggr, 0,
                             ncp->comm, &st);
#endif
        NCI_Free(offsets);
        NCI_Free(lengths);
        offsets = NULL;
        lengths = NULL;
    }

    /* Non-aggregators are now done. */
    if (ncp->rank != ncp->my_aggr) {
        NCI_Free(meta);
        return status;
    }

    /* Below are tasks for aggregators only, which must call MPI-IO functions
     * to read data from the file. Non-aggregators do not participate MPI-IO
     * calls.
     */

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_get[1] = MAX(ncp->maxmem_get[1], mem_max);
#endif

    /* MPI-IO has the following requirements about filetype.
     * 1. The (flattened) displacements (of a filetype) are not required to be
     *    distinct, but they cannot be negative, and they must be monotonically
     *    non-decreasing.
     * 2. If the file is opened for writing, neither the etype nor the filetype
     *    is permitted to contain overlapping regions.
     */
    if (npairs > 0) {
        /* Now the aggregator has received all offset-length pairs from
         * non-aggregators.
         *
         * At first, check if qsort is necessary. Find the first non-aggregator
         * with non-zero pairs whose offsets are not already sorted.
         */
        for (i=-1,j=0; j<ncp->num_nonaggrs; j++) {
            if (i == -1 && meta[j*3] > 0) /* find 1st whose num_pairs > 0 */
                i = j;
            if (meta[j*3+2] == 0) { /* j's offsets are not sorted */
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
            count = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint)* ncp->num_nonaggrs);
            for (i=0; i<ncp->num_nonaggrs; i++) count[i] = meta[i*3];
        }

        /* construct an array of buffer addresses containing a mapping of the
         * buffer used to read from file and the buffer used to send read data
         * to non-aggregator.
         */
        if (do_sort) {
            /* sort offsets and lengths together, based on offsets into an
             * increasing order
             */
            if (indv_sorted) {
                /* heap-merge of already sorted individual lists is much faster
                 * than qsort.  However, it has a much bigger memory footprint.
                 */
                heap_merge(ncp->num_nonaggrs, count, npairs, offsets, lengths,
                           NULL);
                NCI_Free(count);
            }
            else
                /* As some individual offsets are not sorted, we cannot use
                 * heap_merge().  Note qsort is an in-place sorting.
                 */
                qsort_off_len_buf(npairs, offsets, lengths, NULL);
        }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        ncmpi_inq_malloc_size(&mem_max);
        ncp->maxmem_get[2] = MAX(ncp->maxmem_get[2], mem_max);
        if (ncp->rank == ncp->my_aggr)
            ncp->ina_npairs_get = MAX(ncp->ina_npairs_get, npairs);
#endif

        /* Coalesce the offset-length pairs and calculate the total amount to
         * be read by this aggregator.
         */
        buf_count = (npairs > 0) ? lengths[0] : 0;
        for (i=0, j=1; j<npairs; j++) {
            if (offsets[i] + lengths[i] >= offsets[j]) {
                /* coalesce j into i */
                MPI_Offset end_off = offsets[j] + lengths[j];
                if (offsets[i] + lengths[i] < end_off) {
                    lengths[i] = end_off - offsets[i];
                    buf_count += offsets[i] + lengths[i] - offsets[j];
                }
                /* else: j is entirely covered by i */
            }
            else { /* j and i are not overlapped */
                buf_count += lengths[j];
                i++;
                if (i < j) {
                    offsets[i] = offsets[j];
                    lengths[i] = lengths[j];
                }
            }
        }
        /* update npairs after coalesce */
        npairs = i+1;

        /* Allocate read buffer and send buffer. Once data are read from file
         * into rd_buf, it is unpacked into send_buf for each non-aggregator.
         * send_buf will be directly used to send the read request data to
         * non-aggregators, without extra malloc.
         */
        if (buf_count > 0)
            rd_buf = (char*) NCI_Malloc(buf_count);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        ncmpi_inq_malloc_size(&mem_max);
        ncp->maxmem_get[3] = MAX(ncp->maxmem_get[3], mem_max);
#endif

        if (ncp->fstype == ADIO_FSTYPE_MPIIO) {
            /* When using MPI-IO, we need to create a fileType to set the file
             * view before read. The fileType is constructed based on the
             * offset-length pairs.
             */
            if (npairs == 1) {
                /* Need not create fileType if reading from a contiguous space
                 * (when fileType is MPI_BYTE), ncmpio_file_set_view() will
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

    /* Save fd->filetype previously used and will restore it later. */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        disp = 0;
        saved_fileType = ncp->adio_fh->filetype;

        /* When using PnetCDF's internal ADIO driver, a call to
         * ncmpio_file_set_view() will pass offsets[] and lengths[] to
         * PNCIO_File_set_view() to be reused as flattened fileType. This can
         * avoid re-create an MPI datatype and re-flattening it. See such
         * adaptation in ADIOI_Calc_my_off_len()
         *
         * Here we use MPI_DATATYPE_NULL to tell PNCIO_File_set_view() that it
         * is called from the intra-node aggregation subroutine, so the ADIO
         * subroutines can reuse offsets[] and lengths[].
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

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_get[4] = MAX(ncp->maxmem_get[4], mem_max);
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_get[0] += endT - startT;
#endif

    /* call MPI_File_read_at_all */
    err = ncmpio_read_write(ncp, NC_REQ_RD, NC_REQ_COLL, disp, buf_count,
                            MPI_BYTE, rd_buf, 1);
    if (status == NC_NOERR) status = err;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    ncmpi_inq_malloc_size(&mem_max);
    ncp->maxmem_get[5] = MAX(ncp->maxmem_get[5], mem_max);
    startT = MPI_Wtime();
#endif

    /* ncp->adio_fh->flat_file is allocated in ncmpio_file_set_view() */
    if (ncp->fstype != ADIO_FSTYPE_MPIIO) {
        NCI_Free(ncp->adio_fh->flat_file);
        ncp->adio_fh->flat_file = NULL;

        /* restore the original filetype */
        ncp->adio_fh->filetype = saved_fileType;
    }

    /* If sorting has been performed, the orders of offsets[] and lengths[] may
     * no longer be the same as the original ones. We use binary search to find
     * the aggregated offset-length pair containing each non-aggregator's
     * offset-length pair to construct a send buffer datatype, a view layout to
     * the read buffer, rd_buf, so the data can be directly sent from
     * rd_buf, without copying it to a separate buffer.
     */

    /* First, copy data for self rank, also the aggregator. Note offsets[] is
     * in an incremental order.
     */
    int buftype_is_contig;
    char *ptr=NULL, *tmp_buf=NULL;

    /* Check if this rank's user buftype is contiguous. If not, a temporary
     * buffer is allocated and the read data is copied over, before unpacking
     * it to user's buffer.
     */
    PNCIO_Datatype_iscontig(bufType, &buftype_is_contig);

    if (buf != NULL && buftype_is_contig)
        ptr = buf;
    else if (bufLen > 0)
        ptr = tmp_buf = (char*) NCI_Malloc(bufLen);

    size_t m=0, k, scan_off=0;
    for (j=0; j<num_pairs; j++) {
        /* Find the offset-length pair in rd_buf containing this pair.
         * Note that if the offset-length pairs are not already sorted, i.e.
         * is_incr == 1, this bin_search() below can be very expensive!
         */
        if (!is_incr) m = 0;
        k = bin_search(orig_offsets[j], &offsets[m], npairs-m);
        assert(k >= 0); /* k returned from bin_search is relative to m */
        k += m;

        /* When is_incr is 1, the orig_offsets[] are in an incremental order,
         * we can continue binary search using the index from the previous
         * search.  When is_incr is 0, the orig_offsets[] are NOT in an
         * incremental order, we must do binary search on the entire offsets[].
         */
        if (!is_incr) scan_off = 0;
        for (; m<k; m++)
            scan_off += lengths[m];

        /* Note orig_offsets[j] and lengths[j] must entirely covered by
         * offsets[k] and lengths[k], because offsets[] and lengths[]
         * have been coalesced.
         */

        memcpy(ptr, rd_buf + (scan_off + orig_offsets[j] - offsets[k]),
                    orig_lengths[j]);
        ptr += orig_lengths[j];
    }

    /* unpack read data to user read buffer */
    if (bufLen > 0 && (buf == NULL || !buftype_is_contig)) {

        void *buf_ptr = (buf == NULL) ? MPI_BOTTOM : buf;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count pos=0;
        MPI_Count num = (buf == NULL) ? 1 : bufCount;
        MPI_Unpack_c(tmp_buf, bufLen, &pos, buf_ptr, num, bufType, MPI_COMM_SELF);
#else
        int pos=0;
        int num = (buf == NULL) ? 1 : bufCount;
        MPI_Unpack(tmp_buf, bufLen, &pos, buf_ptr, num, bufType, MPI_COMM_SELF);
#endif
        NCI_Free(tmp_buf);
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_get[1] += endT - startT;
    startT = endT;
#endif

    /* allocate array_of_blocklengths[] and array_of_displacements[] */
    for (max_npairs=0, i=1; i<ncp->num_nonaggrs; i++)
        max_npairs = MAX(meta[3*i], max_npairs);

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *blks = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * max_npairs);
    MPI_Count *disps = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * max_npairs);
#else
    int *blks = (int*) NCI_Malloc(sizeof(int) * max_npairs);
    MPI_Aint *disps = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * max_npairs);
#endif

    /* Now, send data to each non-aggregator */
    nreqs = 0;
    MPI_Offset off_start = num_pairs;
    for (i=1; i<ncp->num_nonaggrs; i++) {
        /* populate disps[] and blks[] */
        MPI_Aint remote_num_pairs = meta[3*i];
        MPI_Aint remote_is_incr = meta[3*i+2];

        if (remote_num_pairs == 0) continue; /* zero sized request */

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count *off = orig_offsets + off_start;
        MPI_Count *len = orig_lengths + off_start;
#else
        MPI_Offset *off = orig_offsets + off_start;
        int        *len = orig_lengths + off_start;
#endif
        m = 0;
        scan_off = 0;
        for (j=0; j<remote_num_pairs; j++) {
            MPI_Aint addr;

            /* Find the offset-length pair in rd_buf containing this pair.
             * Note that if the offset-length pairs are not already sorted,
             * i.e. remote_is_incr == 1, this bin_search() below can be very
             * expensive!
             */
            if (!remote_is_incr) m = 0;

            k = bin_search(off[j], &offsets[m], npairs-m);
            assert(k >= 0); /* k returned from bin_search is relative to m */
            k += m;
            assert(offsets[k] <= off[j] && off[j] < offsets[k] + lengths[k]);

            /* When is_incr is 1, the orig_offsets[] are in an incremental
             * order, we can continue binary search using the index from the
             * previous search.  When is_incr is 0, the orig_offsets[] are NOT
             * in an incremental order, we must do binary search on the entire
             * offsets[].
             */
            if (!remote_is_incr) scan_off = 0;
            for (; m<k; m++)
                scan_off += lengths[m];
            /* Note orig_offsets[j] and lengths[j] must entirely covered by
             * offsets[k] and lengths[k], because offsets[] and lengths[] have
             * been coalesced.
             */
            ptr = rd_buf + (scan_off + off[j] - offsets[k]);
            MPI_Get_address(ptr, &addr);
            disps[j] = addr;
            blks[j] = len[j];
        }
        off_start += remote_num_pairs;

        /* Construct a send buffer MPI datatype */
        MPI_Datatype sendType;
#ifdef HAVE_MPI_LARGE_COUNT
        mpireturn = MPI_Type_create_hindexed_c(remote_num_pairs, blks, disps,
                                               MPI_BYTE, &sendType);
#else
        mpireturn = MPI_Type_create_hindexed(remote_num_pairs, blks, disps,
                                             MPI_BYTE, &sendType);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            MPI_Type_commit(&sendType);

#ifdef HAVE_MPI_LARGE_COUNT
            TRACE_COMM(MPI_Isend_c)(MPI_BOTTOM, 1, sendType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#else
            TRACE_COMM(MPI_Isend)(MPI_BOTTOM, 1, sendType,
                       ncp->nonaggr_ranks[i], 0, ncp->comm, &req[nreqs++]);
#endif
            MPI_Type_free(&sendType);
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_get[2] += endT - startT;
    startT = endT;
#endif

#ifdef HAVE_MPI_STATUSES_IGNORE
    TRACE_COMM(MPI_Waitall)(nreqs, req, MPI_STATUSES_IGNORE);
#else
    MPI_Status *statuses = (MPI_Status *)
                           ADIOI_Malloc(nreqs * sizeof(MPI_Status));
    TRACE_COMM(MPI_Waitall)(nreqs, req, statuses);
    ADIOI_Free(statuses);
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn,"MPI_Waitall");
        /* return the first encountered error if there is any */
        if (status == NC_NOERR) status = err;
    }
    NCI_Free(blks);
    NCI_Free(disps);

    /* offsets[] and lengths[] are used in ADIO read subroutines as flattened
     * filetype. They cannot be freed before the I/O is done.
     */
    if (offsets != NULL) NCI_Free(offsets);
    if (lengths != NULL) NCI_Free(lengths);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    endT = MPI_Wtime();
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_get[3] += endT - startT;
#endif

    if (rd_buf != NULL) NCI_Free(rd_buf);
    if (orig_lengths != NULL) NCI_Free(orig_lengths);
    if (orig_offsets != NULL) NCI_Free(orig_offsets);
    if (req != NULL) NCI_Free(req);
    if (meta != NULL) NCI_Free(meta);

    return status;
}

/*----< req_compare() >------------------------------------------------------*/
/* used to sort the the string file offsets of reqs[] */
static int
req_compare(const void *a, const void *b)
{
    if (((NC_req*)a)->offset_start > ((NC_req*)b)->offset_start) return (1);
    if (((NC_req*)a)->offset_start < ((NC_req*)b)->offset_start) return (-1);
    return (0);
}

/*----< ncmpio_intra_node_aggregation_nreqs() >------------------------------*/
/* This subroutine is a collective call, which handles PnetCDF's requests made
 * from non-blocking APIs, which contain multiple requests to one or more
 * variable. The input arguments are described below.
 *    reqMode: NC_REQ_RD for read request and NC_REQ_WR for write.
 *    num_reqs: number of elements in array req_list.
 *    req_list[]: stores pending requests from non-blocking API calls, which is
 *                used to construct file offset-length pairs and user buffer
 *                datatype.
 *    newnumrecs: number of new records
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

    /* ensure the intra-node aggregation is enabled in this rank's group */
    assert(ncp->num_aggrs_per_node > 0);

    /* populate reqs[].offset_start, starting offset of each request */
    NC_req *reqs = req_list;
    int i, descreasing=0;
    for (i=0; i<num_reqs; i++) {
        NC_lead_req *lead;
        NC_var *varp;

        lead = (reqMode == NC_REQ_RD) ? ncp->get_lead_list
                                      : ncp->put_lead_list;
        lead += reqs[i].lead_off;
        varp = lead->varp;

        if (varp->ndims == 0) { /* scalar variable */
            reqs[i].offset_start += varp->begin;
        }
        else if (reqs[i].npairs == 1) { /* only one offset-length pair */
            MPI_Offset off = varp->begin;

            if (IS_RECVAR(varp)) off += reqs[i].start[0] * ncp->recsize;

            reqs[i].offset_start += off;
        }
        else {
            /* start/count/stride have been allocated in a contiguous array */
            MPI_Offset *count, *stride, offset_end;
            count  = reqs[i].start + varp->ndims;
            stride = (fIsSet(lead->flag, NC_REQ_STRIDE_NULL)) ? NULL :
                     count + varp->ndims;

            /* calculate access range of this request */
            ncmpio_calc_start_end(ncp, varp, reqs[i].start, count, stride,
                                  &reqs[i].offset_start, &offset_end);
        }
        /* check if offset_start are in a monotonic nondecreasing order */
        if (i > 0 && reqs[i].offset_start < reqs[i-1].offset_start)
            descreasing = 1;
    }

    /* If a decreasing order is found, sort reqs[] based on reqs[].offset_start
     * into an increasing order.
     */
    if (descreasing)
        qsort(reqs, (size_t)num_reqs, sizeof(NC_req), req_compare);

    /* construct file offset-length pairs
     *     num_pairs: total number of off-len pairs
     *     offsets:   array of flattened offsets
     *     lengths:   array of flattened lengths
     *     is_incr:   whether offsets are incremental
     */
    if (num_reqs > 0)
        flatten_reqs(ncp, reqMode, num_reqs, reqs, &is_incr, &num_pairs,
                     &offsets, &lengths);
    else
        num_pairs = 0;

    /* construct user buffer datatype, bufType.
     * bufLen is the buffer size in bytes
     */
    if (num_reqs > 0) {
        construct_buf_type(ncp, reqMode, num_reqs, reqs, &bufLen,
                           &bufType);
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
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_flatten += MPI_Wtime() - timing;
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
/* This subroutine is a collective call, which to handles a single request made
 * by blocking APIs, which involves only one variable. Below describe the
 * subroutine arguments.
 *    reqMode: NC_REQ_RD for read request and NC_REQ_WR for write.
 *    varp: pointer to the variable struct.
 *    start[]: starting offsets
 *    count[]: counts along each dimension
 *    stride[]: stride along each dimension
 *    bufCount: number of bufType in the user buffer, pointed by buf
 *    bufType: MPI derived datatype describing user buffer layout
 *    buf: pointer to the user buffer
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

    /* ensure the intra-node aggregation is enabled in this rank's group */
    assert(ncp->num_aggrs_per_node > 0);

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
        if (offsets != NULL) NCI_Free(offsets);
        offsets = NULL;
        if (lengths != NULL) NCI_Free(lengths);
        lengths = NULL;
    }
    status = err;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (ncp->rank == ncp->my_aggr) ncp->ina_time_flatten += MPI_Wtime() - timing;
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

