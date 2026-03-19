/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <limits.h> /* INT_MAX */

#include <pncio.h>

/* This file contains four functions:
 *
 * PNCIO_Calc_aggregator()
 * PNCIO_Calc_file_domains()
 * PNCIO_Calc_my_req()
 * PNCIO_Free_my_req()
 * PNCIO_Calc_others_req()
 * PNCIO_Free_others_req()
 */

/*----< PNCIO_Calc_aggregator() >--------------------------------------------*/
/* This subroutine returns the rank ID of aggregator who is responsible for the
 * request represented by (off, *len).  The "len" parameter may be modified to
 * indicate the amount of data actually available in this file domain.
 */
int PNCIO_Calc_aggregator(int               striping_unit,
                          int               cb_nodes,
                          const int        *cb_node_list, /* IN: [cb_nodes] */
                          MPI_Offset        min_st_off,
                          MPI_Offset        fd_size,
                          const MPI_Offset *fd_end,       /* IN: [cb_nodes] */
                          MPI_Offset        off,
                          MPI_Offset       *len)          /* IN/OUT: */
{
    int rank_index, rank;
    MPI_Offset avail_bytes;

    /* get an index into array of aggregators */
    rank_index = (int) ((off - min_st_off + fd_size) / fd_size - 1);

    if (striping_unit > 0) {
        /* Implementation for file domain alignment. Note fd_end[] have been
         * aligned with file system lock boundaries when it was produced by
         * PNCIO_Calc_file_domains().
         */
        rank_index = 0;
        while (off > fd_end[rank_index])
            rank_index++;
    }

    /* we index into fd_end with rank_index, and fd_end was allocated to be no
     * bigger than cb_nodes.   If we ever violate that, we're
     * overrunning arrays.  Obviously, we should never ever hit this abort */
    if (rank_index >= cb_nodes || rank_index < 0) {
        fprintf(stderr,
                "Error %s(): rank_index(%d) >= cb_nodes(%d) fd_size="OFFFMT" off="OFFFMT"\n",
                __func__,rank_index, cb_nodes, fd_size, off);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* fd_end[] is to make sure that we know how much data this aggregator is
     * working with. The +1 is to take into account the end vs. length issue.
     */
    avail_bytes = fd_end[rank_index] + 1 - off;
    if (avail_bytes < *len) {
        /* this file domain only has part of the requested contiguous region */
        *len = avail_bytes;
    }

    /* map our index to a rank */
    rank = cb_node_list[rank_index];

    return rank;
}

/*----<PNCIO_Calc_file_domains() >-------------------------------------------*/
/* Divide the aggregate access region of a collective read/write call into a
 * set of contiguous but disjoined file regions, called file domain (denoted as
 * 'fd'). Each of them is assigned to an I/O aggregator. An aggregator is
 * responsible for carrying out the file I/O for other processes whose requests
 * fall into its file domain.
 *
 * fd_end[cb_nodes]   - end location of file domains, inclusive offsets.
 * fd_size - average size (ceiling) of file domain among cb_nodes.
 */
void PNCIO_Calc_file_domains(int          cb_nodes,
                             int          striping_unit,
                             MPI_Offset   min_st_off,
                             MPI_Offset   max_end_off,
                             MPI_Offset **fd_end,    /* OUT: [cb_nodes] */
                             MPI_Offset  *fd_size)   /* OUT: */
{
    int i, rem_front, rem_back;
    MPI_Offset end_off;

    *fd_end = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * cb_nodes);

    /* partition the aggregate access region equally among I/O aggregators */
    *fd_size  = (max_end_off - min_st_off + 1 + cb_nodes - 1) / cb_nodes;

    /* Align file domain to the nearest file lock boundary (as specified by
     * striping_unit hint).
     */
#ifdef PNETCDF_DEBUG
    assert(striping_unit > 0);
#endif

    /* align fd_end[0] to the nearest file lock boundary */
    end_off = min_st_off + *fd_size;
    rem_front = end_off % striping_unit;
    rem_back = striping_unit - rem_front;
    if (rem_front < rem_back)
        end_off -= rem_front;
    else
        end_off += rem_back;
    (*fd_end)[0] = end_off - 1;

    /* align (*fd_end)[i] to the nearest file lock boundary */
    for (i=1; i<cb_nodes; i++) {
        end_off = min_st_off + *fd_size * (i + 1);
        rem_front = end_off % striping_unit;
        rem_back = striping_unit - rem_front;
        if (rem_front < rem_back)
            end_off -= rem_front;
        else
            end_off += rem_back;
        (*fd_end)[i] = end_off - 1;
    }
    (*fd_end)[cb_nodes - 1] = max_end_off;

    /* Take care of cases in which the aggregate access region is not divisible
     * by the number of aggregators. In such cases, the last process, or the
     * last few processes, may have less load (even 0). For example, a region
     * of 97 divided among 16 processes.  Note that the division is ceiling
     * division.
     */
    for (i=0; i<cb_nodes; i++) {
        if ((*fd_end)[i] > max_end_off)
            (*fd_end)[i] = max_end_off;
    }
}

/*----< PNCIO_Calc_my_req() >------------------------------------------------*/
/* PNCIO_Calc_my_req() calculates every portions of this rank's requests that
 * fall into each aggregator's file domain. When returned, it set the following
 * variables:
 * my_req_naggr - number of aggregators for which this rank has requests that
 *      fall into their file domains
 * count_per_aggr - count of requests for each aggregator, indexed by rank of
 *      the process
 * my_req[nprocs] - array of data structures describing the requests to be
 *      performed by each aggregator.
 * buf_idx[nprocs] - array of locations into which data in the user buffer that
 *      can be directly used to perform read/write; this is only valid when
 *      the user buffer is contiguous.
 */
void PNCIO_Calc_my_req(PNCIO_File         *fh,
                       MPI_Offset          min_st_off,
                       const MPI_Offset   *fd_end,
                       MPI_Offset          fd_size,
                       MPI_Count          *my_req_naggr,   /* OUT: */
                       MPI_Count         **count_per_aggr, /* OUT: */
                       PNCIO_Access      **my_req,         /* OUT: */
                       MPI_Aint          **buf_idx)        /* OUT: */
{
    size_t memLen, alloc_sz;
    int i, nprocs, aggr;
    MPI_Count j, l;
    MPI_Offset fd_len, rem_len, curr_idx, off, *off_ptr;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset *len_ptr;
#else
    int *len_ptr;
#endif

    MPI_Comm_size(fh->comm, &nprocs);

    *my_req_naggr = 0;

    *count_per_aggr = NCI_Calloc(nprocs, sizeof(MPI_Count));

    *my_req = (PNCIO_Access*) NCI_Calloc(nprocs, sizeof(PNCIO_Access));

    /* buf_idx is relevant only if the user buffer is contiguous. buf_idx[i]
     * stores the index into user_buf where data received from rank i should be
     * placed. This allows receives to be done without extra buffer. This can't
     * be done if buftype is not contiguous.
     */
    *buf_idx = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * nprocs);
    for (i=0; i<nprocs; i++) /* initialize buf_idx to -1 */
        (*buf_idx)[i] = -1;

    if (fh->file_view.size == 0) /* zero-sized request */
        return;

#ifdef PNETCDF_DEBUG
    /* For non-zero sized requests, fh->file_view.count has been checked and
     * adjusted to a positive number at the beginning of PNCIO_UFS_read_coll()
     * and PNCIO_UFS_write_coll().
     */
    assert(fh->file_view.count > 0);

    /* fh->file_view's offset-length pairs has been coalesced in ncmpio's
     * ina_put() and ina_get().
     */
    for (j=0; j<fh->file_view.count; j++)
        assert(fh->file_view.len[j] > 0);
#endif

    /* one pass just to calculate how much space to allocate for my_req */
    memLen = 0;
    for (i=0; i<fh->file_view.count; i++) {
        off = fh->file_view.off[i];
        fd_len = fh->file_view.len[i];
        /* Note: we set fd_len to be the total size of the access, then
         * PNCIO_Calc_aggregator() will modify the value to return the
         * amount that was available from the file domain that holds the
         * first part of the access.
         */
        aggr = PNCIO_Calc_aggregator(fh->hints->striping_unit,
                                     fh->hints->cb_nodes, fh->hints->aggr_ranks,
                                     min_st_off, fd_size, fd_end, off,
                                     &fd_len);
        (*count_per_aggr)[aggr]++;
        memLen++;

        /* figure out how much data is remaining in the access (i.e. wasn't
         * part of the file domain that had the starting byte); we'll take
         * care of this data (if there is any) in the while loop below.
         */
        rem_len = fh->file_view.len[i] - fd_len;

        while (rem_len != 0) {
            off += fd_len;      /* point to first remaining byte */
            fd_len = rem_len;   /* save remaining size, pass to calc */
            aggr = PNCIO_Calc_aggregator(fh->hints->striping_unit,
                                         fh->hints->cb_nodes,
                                         fh->hints->aggr_ranks, min_st_off,
                                         fd_size, fd_end, off, &fd_len);

            (*count_per_aggr)[aggr]++;
            memLen++;
            rem_len -= fd_len;  /* reduce remaining length by amount from fd */
        }
    }

#ifdef HAVE_MPI_LARGE_COUNT
    alloc_sz = sizeof(MPI_Offset) * 2;
    (*my_req)[0].offsets = (MPI_Offset *) NCI_Malloc(alloc_sz * memLen);
    (*my_req)[0].lens = (*my_req)[0].offsets + memLen;
#else
    alloc_sz = sizeof(MPI_Offset) + sizeof(int);
    (*my_req)[0].offsets = (MPI_Offset *) NCI_Malloc(alloc_sz * memLen);
    (*my_req)[0].lens = (int*) ((*my_req)[0].offsets + memLen);
#endif

    off_ptr = (*my_req)[0].offsets;
    len_ptr = (*my_req)[0].lens;
    for (i=0; i<nprocs; i++) {
        if ((*count_per_aggr)[i]) {
            (*my_req)[i].offsets = off_ptr;
            off_ptr += (*count_per_aggr)[i];
            (*my_req)[i].lens = len_ptr;
            len_ptr += (*count_per_aggr)[i];
            (*my_req_naggr)++;
        }
        (*my_req)[i].count = 0; /* will be incremented in loop j below */
    }

    /* now fill in my_req */
    curr_idx = 0;
    for (j=0; j<fh->file_view.count; j++) {
        off = fh->file_view.off[j];
        fd_len = fh->file_view.len[j];

        aggr = PNCIO_Calc_aggregator(fh->hints->striping_unit,
                                     fh->hints->cb_nodes, fh->hints->aggr_ranks,
                                     min_st_off, fd_size, fd_end, off,
                                     &fd_len);

        /* for each separate contiguous access from this process */
        if ((*buf_idx)[aggr] == -1)
            (*buf_idx)[aggr] = (MPI_Aint) curr_idx;

        l = (*my_req)[aggr].count;
        curr_idx += fd_len;

        rem_len = fh->file_view.len[j] - fd_len;

        /* store aggr, offset, and len in an array of structures, my_req. Each
         * structure contains the offsets and lengths located in that process's
         * FD, and the associated count.
         */
        (*my_req)[aggr].offsets[l] = off;
        (*my_req)[aggr].lens[l] = fd_len;
        (*my_req)[aggr].count++;

        while (rem_len != 0) {
            off += fd_len;
            fd_len = rem_len;
            aggr = PNCIO_Calc_aggregator(fh->hints->striping_unit,
                                         fh->hints->cb_nodes,
                                         fh->hints->aggr_ranks, min_st_off,
                                         fd_size, fd_end, off, &fd_len);

            if ((*buf_idx)[aggr] == -1)
                (*buf_idx)[aggr] = (MPI_Aint) curr_idx;

            l = (*my_req)[aggr].count;
            curr_idx += fd_len;
            rem_len -= fd_len;

            (*my_req)[aggr].offsets[l] = off;
            (*my_req)[aggr].lens[l] = fd_len;
            (*my_req)[aggr].count++;
        }
    }
}

void PNCIO_Free_my_req(MPI_Count    *count_per_aggr,
                       PNCIO_Access *my_req,
                       MPI_Aint     *buf_idx)
{
    NCI_Free(count_per_aggr);
    NCI_Free(my_req[0].offsets);
    NCI_Free(my_req);
    NCI_Free(buf_idx);
}

/*----< PNCIO_Calc_others_req() >--------------------------------------------*/
/* PNCIO_Calc_others_req() is only relevant to the I/O aggregators. Based
 * on everyone's my_req, PNCIO_Calc_others_req() calculates what requests
 * of all processes fall into this aggregator's file domain. This
 * subroutine sets the following variables:
 * count_others_req_procs - number of processes whose requests fall into
 *      this aggregator's file domain (including this rank itself)
 * count_others_req_per_proc[i] - how many non-contiguous requests of rank
 *      i fall into this aggregator's file domain.
 */
void PNCIO_Calc_others_req(PNCIO_File    *fh,
                           MPI_Count      my_req_naggr,
                           MPI_Count     *count_per_aggr,
                           PNCIO_Access  *my_req,
                           MPI_Count     *count_others_req_procs,    /* OUT: */
                           MPI_Count    **count_others_req_per_proc, /* OUT: */
                           PNCIO_Access **others_req)                /* OUT: */
{
    size_t alloc_sz, memLen;
    int i, j, nprocs, myrank;
    MPI_Request *reqs;
    MPI_Offset *off_ptr;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset *len_ptr;
    MPI_Count *mem_ptr;
#else
    int *len_ptr;
    MPI_Aint *mem_ptr;
#endif

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

    /* first find out how much to send/recv and from/to whom */
    *count_others_req_per_proc = NCI_Malloc(sizeof(MPI_Count) * nprocs);

    MPI_Alltoall(count_per_aggr, 1, MPI_COUNT, *count_others_req_per_proc,
                                 1, MPI_COUNT, fh->comm);

    *others_req = (PNCIO_Access*) NCI_Malloc(sizeof(PNCIO_Access) * nprocs);

    memLen = 0;
    for (i=0; i<nprocs; i++)
        memLen += (*count_others_req_per_proc)[i];

#ifdef HAVE_MPI_LARGE_COUNT
    alloc_sz = sizeof(MPI_Offset) * 2 + sizeof(MPI_Count);
    (*others_req)[0].offsets = (MPI_Offset *) NCI_Malloc(alloc_sz * memLen);
    (*others_req)[0].lens = (*others_req)[0].offsets + memLen;
    (*others_req)[0].mem_ptrs = (MPI_Count*) ((*others_req)[0].lens + memLen);
#else
    alloc_sz = sizeof(MPI_Offset) + sizeof(int) + sizeof(MPI_Aint);
    (*others_req)[0].offsets = (MPI_Offset *) NCI_Malloc(alloc_sz * memLen);
    (*others_req)[0].lens = (int *) ((*others_req)[0].offsets + memLen);
    (*others_req)[0].mem_ptrs = (MPI_Aint*) ((*others_req)[0].lens + memLen);
#endif
    off_ptr = (*others_req)[0].offsets;
    len_ptr = (*others_req)[0].lens;
    mem_ptr = (*others_req)[0].mem_ptrs;

    *count_others_req_procs = 0;
    for (i=0; i<nprocs; i++) {
        if ((*count_others_req_per_proc)[i]) {
            (*others_req)[i].count = (*count_others_req_per_proc)[i];
            (*others_req)[i].offsets = off_ptr;
            off_ptr += (*count_others_req_per_proc)[i];
            (*others_req)[i].lens = len_ptr;
            len_ptr += (*count_others_req_per_proc)[i];
            (*others_req)[i].mem_ptrs = mem_ptr;
            mem_ptr += (*count_others_req_per_proc)[i];
            (*count_others_req_procs)++;
        } else
            (*others_req)[i].count = 0;
    }

    /* now send the calculated offsets and lengths to respective processes */
    reqs = (MPI_Request*) NCI_Malloc(sizeof(MPI_Request) *
           (my_req_naggr + *count_others_req_procs) * 2);

    j = 0;
    for (i=0; i<nprocs; i++) {
        if ((*others_req)[i].count == 0)
            continue;
        if (i == myrank) {
            /* send to self by using memcpy().
             * Note (*others_req)[i].count == my_req[i].count
             */
            memcpy((*others_req)[i].offsets, my_req[i].offsets,
                   my_req[i].count * sizeof(MPI_Offset));
#ifdef HAVE_MPI_LARGE_COUNT
            memcpy((*others_req)[i].lens, my_req[i].lens,
                   my_req[i].count * sizeof(MPI_Offset));
#else
            memcpy((*others_req)[i].lens, my_req[i].lens,
                   my_req[i].count * sizeof(int));
#endif
        }
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Irecv_c((*others_req)[i].offsets, (*others_req)[i].count,
                        MPI_OFFSET, i, i + myrank, fh->comm, &reqs[j++]);
            MPI_Irecv_c((*others_req)[i].lens, (*others_req)[i].count,
                        MPI_OFFSET, i, i + myrank, fh->comm, &reqs[j++]);
#else
            /* check overflow 4-byte int */
            assert((*others_req)[i].count <= 2147483647);

            MPI_Irecv((*others_req)[i].offsets, (int)(*others_req)[i].count,
                      MPI_OFFSET, i, i + myrank, fh->comm, &reqs[j++]);
            MPI_Irecv((*others_req)[i].lens, (int)(*others_req)[i].count,
                      MPI_INT, i, i + myrank, fh->comm, &reqs[j++]);
#endif
        }
    }

    for (i=0; i<nprocs; i++) {
        if (my_req[i].count && i != myrank) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Isend_c(my_req[i].offsets, my_req[i].count,
                        MPI_OFFSET, i, i + myrank, fh->comm, &reqs[j++]);
            MPI_Isend_c(my_req[i].lens, my_req[i].count,
                        MPI_OFFSET, i, i + myrank, fh->comm, &reqs[j++]);
#else
            assert(my_req[i].count <= INT_MAX); /* overflow 4-byte int */
            MPI_Isend(my_req[i].offsets, (int)my_req[i].count,
                      MPI_OFFSET, i, i + myrank, fh->comm, &reqs[j++]);
            MPI_Isend(my_req[i].lens, (int)my_req[i].count,
                      MPI_INT, i, i + myrank, fh->comm, &reqs[j++]);
#endif
        }
    }

    if (j) {
#ifdef HAVE_MPI_STATUSES_IGNORE
        MPI_Waitall(j, reqs, MPI_STATUSES_IGNORE);
#else
        MPI_Status *sts = (MPI_Status*) NCI_Malloc(sizeof(MPI_Status) * j);
        MPI_Waitall(j, reqs, sts);
        NCI_Free(sts);
#endif
    }

    NCI_Free(reqs);
}

void PNCIO_Free_others_req(MPI_Count    *count_others_req_per_proc,
                           PNCIO_Access *others_req)
{
    NCI_Free(count_others_req_per_proc);
    NCI_Free(others_req[0].offsets);
    NCI_Free(others_req);
}

