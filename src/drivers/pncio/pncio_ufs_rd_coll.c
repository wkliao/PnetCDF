/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdbool.h> /* type bool */

#include <pncio.h>

#define BUF_INCR {                                                  \
    while (buf_incr) {                                              \
        size_in_buf = MIN(buf_incr, buf_rem);                       \
        user_buf_idx += size_in_buf;                                \
        buf_rem -= size_in_buf;                                     \
        buf_incr -= size_in_buf;                                    \
        if (buf_incr > 0 && buf_rem == 0) {                         \
            buf_indx++;                                             \
            user_buf_idx = buf_view.off[buf_indx];                  \
            buf_rem = buf_view.len[buf_indx];                       \
        }                                                           \
    }                                                               \
}


#define BUF_COPY {                                                  \
    while (size) {                                                  \
        size_in_buf = MIN(size, buf_rem);                           \
        memcpy((char*)buf + user_buf_idx,                           \
               &recv_buf[aggr][recv_buf_idx[aggr]], size_in_buf);   \
        recv_buf_idx[aggr] += size_in_buf;                          \
        user_buf_idx += size_in_buf;                                \
        buf_rem -= size_in_buf;                                     \
        size -= size_in_buf;                                        \
        buf_incr -= size_in_buf;                                    \
        if (size > 0 && buf_rem == 0) {                             \
            buf_indx++;                                             \
            user_buf_idx = buf_view.off[buf_indx];                  \
            buf_rem = buf_view.len[buf_indx];                       \
        }                                                           \
    }                                                               \
    BUF_INCR                                                        \
}

/*----< fill_user_buffer() >-------------------------------------------------*/
/* This subroutine is only called when buf_view is not contiguous. */
static
void fill_user_buffer(PNCIO_File       *fh,
                      void             *buf,
                      PNCIO_View        buf_view,
                      MPI_Offset        min_st_off,
                      MPI_Offset        fd_size,
                      const MPI_Offset *fd_end,         /* IN: [cb_nodes] */
                      const MPI_Count  *recv_size,      /* IN: [nprocs] */
                      MPI_Count        *recd_from_proc, /* IN/OUT: [nprocs] */
                      char *const      *recv_buf)       /* IN: [nprocs] */
{
    int i, nprocs, aggr, buf_indx;
    MPI_Offset buf_rem, size_in_buf, buf_incr, size;
    MPI_Offset off, user_buf_idx;
    MPI_Offset len, rem_len;
    MPI_Count j, *curr_from, *done_from, *recv_buf_idx;

    MPI_Comm_size(fh->comm, &nprocs);

    /* curr_from[nprocs] - amount of data received from each rank that has
     *      already been accounted for so far.
     * done_from[nprocs] - amount of data already received from each rank and
     *      filled into user buffer in previous round.
     * user_buf_idx - current location in user buffer
     * recv_buf_idx[nprocs] = current location in recv_buf of each rank
     */
    curr_from = NCI_Malloc(sizeof(MPI_Count) * nprocs * 3);
    done_from = curr_from + nprocs;
    recv_buf_idx = done_from + nprocs;

    for (i=0; i<nprocs; i++) {
        recv_buf_idx[i] = curr_from[i] = 0;
        done_from[i] = recd_from_proc[i];
    }

    /* buf_indx - index buf_view's offset-length pairs being processed
     * buf_rem - remaining length of the current offset-length pair
     */
    user_buf_idx = buf_view.off[0];
    buf_indx = 0;
    buf_rem = buf_view.len[0];

    for (j=0; j<fh->file_view.count; j++) {
        off = fh->file_view.off[j];
        rem_len = fh->file_view.len[j];

        /* this request may span file domains of more than one aggregator */
        while (rem_len != 0) {
            len = rem_len;

            /* NOTE: len value will be modified by PNCIO_Calc_aggregator() to
             * be no more than the single file domain that aggregator 'aggr'
             * is responsible for.
             */
            aggr = PNCIO_Calc_aggregator(fh->hints->striping_unit,
                                         fh->hints->cb_nodes,
                                         fh->hints->aggr_ranks, min_st_off,
                                         fd_size, fd_end, off, &len);

            if (recv_buf_idx[aggr] < recv_size[aggr]) {
                if (curr_from[aggr] + len > done_from[aggr]) {
                    if (done_from[aggr] > curr_from[aggr]) {
                        size = MIN(curr_from[aggr] + len - done_from[aggr],
                                   recv_size[aggr] - recv_buf_idx[aggr]);
                        buf_incr = done_from[aggr] - curr_from[aggr];
                        BUF_INCR
                        buf_incr = curr_from[aggr] + len - done_from[aggr];
                        curr_from[aggr] = done_from[aggr] + size;
                        BUF_COPY
                    } else {
                        size = MIN(len, recv_size[aggr] - recv_buf_idx[aggr]);
                        buf_incr = len;
                        curr_from[aggr] += size;
                        BUF_COPY
                    }
                } else {
                    curr_from[aggr] += len;
                    buf_incr = len;
                    BUF_INCR
                }
            } else {
                buf_incr = len;
                BUF_INCR
            }
            off += len;
            rem_len -= len;
        }
    }
    for (i=0; i<nprocs; i++)
        if (recv_size[i])
            recd_from_proc[i] = curr_from[i];

    NCI_Free(curr_from);
}

static
void R_Exchange_data(PNCIO_File   *fh,
                     void         *buf,
                     PNCIO_View    buf_view,
                     MPI_Count    *send_size,      /* IN: [nprocs] */
                     MPI_Count    *recv_size,      /* OUT: [nprocs] */
                     MPI_Count    *count,          /* IN: [nprocs] */
                     MPI_Count    *start_pos,      /* IN: [nprocs] */
                     MPI_Count    *partial_send,   /* IN: [nprocs] */
                     MPI_Count    *recd_from_proc, /* IN/OUT: [nprocs] */
                     MPI_Offset    min_st_off,
                     MPI_Offset    fd_size,
                     MPI_Offset   *fd_end,         /* IN: [cb_nodes] */
                     PNCIO_Access *others_req,     /* IN: [nprocs] */
                     MPI_Aint     *buf_idx,        /* IN: [nprocs] */
                     MPI_Aint     *recved_bytes)   /* OUT: */
{
    char **recv_buf = NULL;
    int i, nprocs, myrank, nrecvs, nsends;
    MPI_Request *reqs;
    MPI_Datatype send_type;
    MPI_Status *sts;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double curT = MPI_Wtime();
#endif

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

    /* Exchange send_size info so that each aggregator knows how much to
     * receive from whom and how much memory to allocate.
     */
    MPI_Alltoall(send_size, 1, MPI_COUNT, recv_size, 1, MPI_COUNT, fh->comm);

    reqs = (MPI_Request*) NCI_Malloc(sizeof(MPI_Request) * 2 * nprocs);

    /* Post nonblocking receive calls. If buf_view is contiguous, data can be
     * directly received into user buf at location given by buf_idx. Otherwise,
     * allocate recv_buf and use it to receive.
     */
    nrecvs = 0;
    if (buf_view.count <= 1) {
        for (i=0; i<nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c((char*)buf + buf_idx[i], recv_size[i], MPI_BYTE, i,
                            0, fh->comm, &reqs[nrecvs++]);
#else
                MPI_Irecv((char*)buf + buf_idx[i], recv_size[i], MPI_BYTE, i,
                           0, fh->comm, &reqs[nrecvs++]);
#endif
                buf_idx[i] += recv_size[i];
            }
        }
    } else {
        size_t memLen = 0;
        for (i=0; i<nprocs; i++)
            memLen += recv_size[i];

        /* allocate memory for recv_buf */
        recv_buf = (char **) NCI_Malloc(sizeof(char*) * nprocs);
        recv_buf[0] = (char *) NCI_Malloc(memLen);
        for (i=1; i<nprocs; i++)
            recv_buf[i] = recv_buf[i - 1] + recv_size[i - 1];

        /* post receives */
        for (i=0; i<nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(recv_buf[i], recv_size[i], MPI_BYTE, i,
                            0, fh->comm, &reqs[nrecvs++]);
#else
                MPI_Irecv(recv_buf[i], recv_size[i], MPI_BYTE, i,
                            0, fh->comm, &reqs[nrecvs++]);
#endif
            }
        }
    }

    /* Construct derived datatypes and use them to send data */
    nsends = 0;
    for (i=0; i<nprocs; i++) {
        if (send_size[i]) {
            /* take care the last off-len pair if is a partial send */
            MPI_Offset tmp = 0;
            MPI_Count k = 0;
            if (partial_send[i]) {
                k = start_pos[i] + count[i] - 1;
                tmp = others_req[i].lens[k];
                others_req[i].lens[k] = partial_send[i];
            }
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(count[i],
                                       &others_req[i].lens[start_pos[i]],
                                       &others_req[i].mem_ptrs[start_pos[i]],
                                       MPI_BYTE, &send_type);
#else
            MPI_Type_create_hindexed(count[i],
                                     &others_req[i].lens[start_pos[i]],
                                     &others_req[i].mem_ptrs[start_pos[i]],
                                     MPI_BYTE, &send_type);
#endif
            /* absolute displacement; use MPI_BOTTOM in send */
            MPI_Type_commit(&send_type);
            MPI_Isend(MPI_BOTTOM, 1, send_type, i, 0, fh->comm,
                      reqs + nrecvs + nsends);
            MPI_Type_free(&send_type);
            if (partial_send[i])
                others_req[i].lens[k] = tmp;
            nsends++;
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->read_timing[4] += MPI_Wtime() - curT;
#endif

    sts = (MPI_Status*) NCI_Malloc(sizeof(MPI_Status) * (nsends + nrecvs));

    /* wait on the receives */
    if (nrecvs) {
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        curT = MPI_Wtime();
#endif
        MPI_Waitall(nrecvs, reqs, sts);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        if (fh->is_agg) fh->read_timing[3] += MPI_Wtime() - curT;
#endif

        *recved_bytes = 0;
        for (i=0; i<nrecvs; i++) {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count count_recved;
            MPI_Get_count_c(&sts[i], MPI_BYTE, &count_recved);
#else
            int count_recved;
            MPI_Get_count(&sts[i], MPI_BYTE, &count_recved);
#endif
            *recved_bytes += count_recved;
        }

        /* When buf is noncontiguous, copy data from recv_buf to buf */
        if (buf_view.count > 1)
            fill_user_buffer(fh, buf, buf_view, min_st_off, fd_size, fd_end,
                             recv_size, recd_from_proc, recv_buf);
    }

    /* wait on the sends */
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    curT = MPI_Wtime();
#endif
#ifdef HAVE_MPI_STATUSES_IGNORE
    MPI_Waitall(nsends, reqs + nrecvs, MPI_STATUSES_IGNORE);
#else
    MPI_Waitall(nsends, reqs + nrecvs, sts + nrecvs);
#endif
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->read_timing[3] += MPI_Wtime() - curT;
#endif

    NCI_Free(sts);
    NCI_Free(reqs);

    if (buf_view.count > 1) { /* buf_view is noncontiguous */
        NCI_Free(recv_buf[0]);
        NCI_Free(recv_buf);
    }
}

/*----< Read_and_exch() >----------------------------------------------------*/
/* Each aggregator reads in sizes of no more than coll_bufsize, an I/O info,
 * sends read data to the requesting processes. All processes place nonblocking
 * receive calls to receive read data into user buffer. The idea is to reduce
 * the amount of extra memory required for collective I/O. If all data were
 * read all at once, which is much easier, it would require allocating a
 * temporary space more than the size of user_buf, which is often unacceptable.
 * For example, to read a distributed array from a file, where each local array
 * is 8 MiB, requiring at least another 8 MiB of space is unacceptable.
 */
static
MPI_Offset Read_and_exch(PNCIO_File   *fh,
                         void         *buf,
                         PNCIO_View    buf_view,
                         PNCIO_Access *others_req, /* IN: [nprocs] */
                         MPI_Offset    min_st_off,
                         MPI_Offset    fd_size,
                         MPI_Offset   *fd_end,   /* IN: [cb_nodes] */
                         MPI_Aint     *buf_idx)  /* IN: [nprocs] */
{
    char *read_buf = NULL;
    int i, m, ntimes, max_ntimes, nprocs, myrank;
    MPI_Offset st_loc = -1, end_loc = -1, rem_off, done, real_off;
    MPI_Count j, *curr_offlen_ptr, *count, *send_size, *recv_size;
    MPI_Count *partial_send, *recd_from_proc, *start_pos;
    /* Not convinced end_loc-st_loc couldn't be > int, so make these offsets */
    MPI_Offset real_size, rem_size, for_curr_round, for_next_round;
    MPI_Aint coll_bufsize, actual_recved_bytes = 0;
    MPI_Offset r_len;

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

    /* Calculate the first and last file offsets (st_loc and end_loc,
     * respectively) this aggregator will read from the file.
     */
    for (i=0; i<nprocs; i++) {
        /* Some processes may not have data for this aggregator */
        if (others_req[i].count) {
            st_loc = others_req[i].offsets[0];
            end_loc = others_req[i].offsets[0];
            break;
        }
    }
    for (; i<nprocs; i++) {
        for (j=0; j<others_req[i].count; j++) {
            st_loc = MIN(st_loc, others_req[i].offsets[j]);
            end_loc = MAX(end_loc, (others_req[i].offsets[j]
                                  + others_req[i].lens[j] - 1));
        }
    }

    /* Calculate the number of rounds of two-phase read, ntimes, each round
     * reads an amount of no more than coll_bufsize. Then, a call to
     * MPI_Allreduce() to obtain the max number of rounds, max_ntimes, among
     * all processes.
     */
    coll_bufsize = fh->hints->cb_buffer_size;
    if ((st_loc == -1) && (end_loc == -1)) {
        /* this process does no I/O. */
        ntimes = 0;
    } else {
        /* ntimes is a ceiling */
        ntimes = (int) ((end_loc - st_loc + coll_bufsize) / coll_bufsize);
    }

    MPI_Allreduce(&ntimes, &max_ntimes, 1, MPI_INT, MPI_MAX, fh->comm);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fh->read_counter[0] = MAX(fh->read_counter[0], max_ntimes);
#endif

    curr_offlen_ptr = NCI_Calloc(nprocs * 7, sizeof(*curr_offlen_ptr));
    /* curr_offlen_ptr is explained below, it must be initialized to 0s. */

    count = curr_offlen_ptr + nprocs;
    /* count is the number of off-len pairs per processes that are satisfied in
     * a round. Must be initialized to 0s.
     */

    partial_send = count + nprocs;
    /* If only a portion of the last off-len pair is sent to an aggregator in a
     * round, then partial_send is the amount sent. Must be initialized to 0s.
     */

    send_size = partial_send + nprocs;
    /* Total size of data to be sent to each aggregator in a round */

    recv_size = send_size + nprocs;
    /* Total size of data to be received from each process in a round. Make it
     * an array of size nprocs, so that I can use MPI_Alltoall later.
     */

    recd_from_proc = recv_size + nprocs;
    /* Amount of data received so far from each process. It is only used in
     * fill_user_buffer(). Must be initialized to 0s.
     */

    start_pos = recd_from_proc + nprocs;
    /* start_pos stores the starting value of curr_offlen_ptr[i] in a round */

    done = 0;
    rem_off = st_loc;
    for_curr_round = for_next_round = 0;

    for (m=0; m<ntimes; m++) {
        /* Each of ntimes rounds, an aggregator reads into buf of amount no
         * more than coll_bufsize bytes. Each aggregator goes through all
         * others_req[] and check if any are satisfied by the current round of
         * read.
         */

        /* Since MPI standard requires that displacements in filetypes are
         * sorted in a monotonically non-decreasing order, I can maintain a
         * pointer, curr_offlen_ptr, to the current off-len pair being
         * processed for each process in others_req[] and scan further only
         * from there. There is still a problem of filetypes such as described
         * below: (1, 2, 3 are not process ranks. They are just three
         * offset-length pairs in a filetype.)
         *
         * pair 1  -------!--
         * pair 2    -----!----
         * pair 3       --!-----
         *
         * where '!' indicates where the current read_size limitation cuts
         * through the filetype. I resolve this by reading up to '!', but
         * filling the communication buffer only for pair 1. Move the portion
         * left over for 2 to the front of read_buf for use in the next round.
         * i.e., pair 2 and pair 3 will be satisfied in the next round. This
         * simplifies filling in the user's buf at the other end, as only one
         * offset-length pair with incomplete data will be sent. I also don't
         * need to send the individual offsets and lengths along with the data,
         * as the data is being sent in a particular order.
         *
         * rem_off   = start file offset for data actually read in this round
         * rem_size  = size of data to be read, corresponding to rem_off
         * real_off  = rem_off minus whatever data was retained in read_buf
         *             from the previous round for cases like pair 2, or pair
         *             3 illustrated above
         * real_size = size plus the extra corresponding to real_off
         * req_off   = file offset for a particular offset-length pair minus
         *             what has been satisfied in previous round
         * req_size  = size corresponding to req_off
         */

        /* fh->io_buf has been allocated at file open time and may be
         * re-allocated at the end of each round.
         */
        read_buf = fh->io_buf;

        rem_size = MIN(coll_bufsize, end_loc - st_loc + 1 - done);

        MPI_Offset round_end = rem_off + rem_size;

        for (i=0; i<nprocs; i++) {
            if (others_req[i].count == 0)
                continue;

            /* This should be only reachable by I/O aggregators only */
            for (j=curr_offlen_ptr[i]; j<others_req[i].count; j++) {
                if (others_req[i].offsets[j] + partial_send[i] < round_end) {
#ifdef PNETCDF_DEBUG
                    assert(for_curr_round + rem_size <= coll_bufsize);
#endif
                    r_len = PNCIO_UFS_read_contig(fh, read_buf + for_curr_round,
                                                  rem_size, rem_off);
                    if (r_len < 0) {
                        actual_recved_bytes = r_len;
                        goto err_out;
                    }
                    rem_size = r_len;
                    break;
                }
            }
        }
        real_off = rem_off - for_curr_round;
        real_size = rem_size + for_curr_round;
        for_next_round = 0;

        for (i=0; i<nprocs; i++) {
            count[i] = send_size[i] = 0;
            if (others_req[i].count == 0)
                continue;

            start_pos[i] = curr_offlen_ptr[i];
            for (j=curr_offlen_ptr[i]; j<others_req[i].count; j++) {
                MPI_Aint addr;
                MPI_Offset req_off, rem_len;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Offset req_len;
#else
                int req_len;
#endif
                if (partial_send[i]) {
                    /* This request may have been partially satisfied in the
                     * previous round.
                     */
                    req_off = others_req[i].offsets[j] + partial_send[i];
                    req_len = others_req[i].lens[j]    - partial_send[i];
                    partial_send[i] = 0;
                    /* modify the off-len pair to reflect this change */
                    others_req[i].offsets[j] = req_off;
                    others_req[i].lens[j]    = req_len;
                } else {
                    req_off = others_req[i].offsets[j];
                    req_len = others_req[i].lens[j];
                }

                rem_len = real_off + real_size - req_off;
                if (rem_len <= 0)
                    break;

                /* now req_off < real_off + real_size */
                count[i]++;

#ifdef PNETCDF_DEBUG
                assert(req_off - real_off <= coll_bufsize);
#endif
                MPI_Get_address(read_buf + req_off - real_off, &addr);
                others_req[i].mem_ptrs[j] = addr;
                send_size[i] += MIN(rem_len, req_len);

                if (rem_len < req_len) {
                    partial_send[i] = rem_len;
                    if (j + 1 < others_req[i].count &&
                        others_req[i].offsets[j + 1] < real_off + real_size) {
                        /* This is the case illustrated in the figure above. */
                        for_next_round = MAX(for_next_round,
                                             real_off + real_size -
                                             others_req[i].offsets[j + 1]);
                        /* max because it must cover requests from different
                         * processes */
                    }
                    break;
                }
            }
            curr_offlen_ptr[i] = j;
        }

        for_curr_round = for_next_round;

        /* carry out the communication phase */
        MPI_Aint recved_bytes = 0;
        R_Exchange_data(fh, buf, buf_view, send_size, recv_size, count,
                        start_pos, partial_send, recd_from_proc, min_st_off,
                        fd_size, fd_end, others_req, buf_idx, &recved_bytes);
        actual_recved_bytes += recved_bytes;

        if (for_next_round) {
            /* move remaining data to the front of fh->io_buf for next round */
#ifdef PNETCDF_DEBUG
            assert(real_size - for_next_round <= coll_bufsize);
#endif
            fh->io_buf = (char*) NCI_Malloc(coll_bufsize);
            memcpy(fh->io_buf, read_buf + real_size - for_next_round,
                   for_next_round);
            NCI_Free(read_buf);
        }

        rem_off += rem_size;
        done += rem_size;
    }

    /* This process is done with its I/O, and must run the remaining rounds to
     * participate the collective communication in R_Exchange_data().
     */
    for (i=0; i<nprocs; i++)
        count[i] = send_size[i] = 0;

    for (m=ntimes; m<max_ntimes; m++) {
        /* nothing to send, but check for recv. */
        MPI_Aint recved_bytes = 0;
        R_Exchange_data(fh, buf, buf_view, send_size, recv_size, count,
                        start_pos, partial_send, recd_from_proc, min_st_off,
                        fd_size, fd_end, others_req, buf_idx, &recved_bytes);
        actual_recved_bytes += recved_bytes;
    }

err_out:
    NCI_Free(curr_offlen_ptr);

    return actual_recved_bytes;
}

/*----< PNCIO_UFS_read_coll() >----------------------------------------------*/
MPI_Offset PNCIO_UFS_read_coll(PNCIO_File *fh,
                               void       *buf,
                               PNCIO_View  buf_view)
{
    /* Uses a generalized version of the extended two-phase method described in
     * "An Extended Two-Phase Method for Accessing Sections of Out-of-Core
     * Arrays", Rajeev Thakur and Alok Choudhary, Scientific Programming,
     * (5)4:301--317, Winter 1996.
     * http://www.mcs.anl.gov/home/thakur/ext2ph.ps
     */

    /* my_req contains access structures of this rank, describing the request
     * offset-length pairs that fall into each aggregator's file domain.
     */
    PNCIO_Access *my_req;

    /* others_req contains access structures of all processes whose requests
     * fall into this aggregator's file domain. It is only relevant of this
     * rank is an I/O aggregator.
     */
    PNCIO_Access *others_req;

    int i, nprocs, rank, interleave_count = 0;
    MPI_Aint *buf_idx = NULL;
    MPI_Count *count_my_req_per_proc, count_my_req_procs;
    MPI_Count *count_others_req_per_proc, count_others_req_procs;
    MPI_Offset fd_size, min_st_off, max_end_off;
    MPI_Offset *fd_end=NULL, r_len, total_r_len=0;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
double curT = MPI_Wtime();
#endif

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &rank);

    /* PnetCDF never reuses a fileview across two or more PNCIO calls. As
     * fh->file_view will be reset right after this subroutine returns, it
     * can be modify within this subroutine.
     */

#ifdef PNETCDF_DEBUG
    if (fh->file_view.size > 0) {
        assert(fh->file_view.count > 0);
        assert(fh->file_view.off != NULL);
        assert(fh->file_view.len != NULL);
    }
    assert(fh->file_view.size == buf_view.size);
#endif

    /* only check for interleaving if romio_cb_read isn't disabled */
    if (fh->hints->romio_cb_read != PNCIO_HINT_DISABLE) {
        MPI_Offset *st_end_all;

        /* For this process's request, calculate its aggregate access region,
         * representing a range from starting offset till end_offset. Note:
         * end offset points to the last byte-offset that will be accessed,
         * e.g., if starting offset=0 and 100 bytes to be read, end_offset=99.
         *
         * Note file_view.off[] is always relative to beginning of file.
         *
         * All processes gather the aggregate access regions of all other
         * processes in order to tell whether there is an interleaving access
         * among all.
         */
        st_end_all = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * 2 * nprocs);
        if (fh->file_view.size == 0)
            /* set to -1 to indicate zero-sized request */
            st_end_all[2*rank] = st_end_all[2*rank+1] = -1;
        else {
            st_end_all[2*rank]   = fh->file_view.off[0];
            st_end_all[2*rank+1] = fh->file_view.off[fh->file_view.count-1]
                                 + fh->file_view.len[fh->file_view.count-1] - 1;
        }

        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, st_end_all, 2,
                      MPI_OFFSET, fh->comm);

        /* Are the accesses of different processes interleaved? Below is a
         * rudimentary check for interleaving, but should suffice for the
         * moment.
         */
        for (i=0; i<2*nprocs; i+=2) { /* find the 1st non-zero sized */
            if (st_end_all[i] >= 0) {
                min_st_off  = st_end_all[i];
                max_end_off = st_end_all[i+1];
                break;
            }
        }
        for (i+=2; i<2*nprocs; i+=2) {
            if (st_end_all[i] == -1) /* skip zero-sized request */
                continue;
            if (st_end_all[i] <  st_end_all[i-1] &&
                st_end_all[i] <= st_end_all[i+1])
                interleave_count++;
            min_st_off  = MIN(min_st_off,  st_end_all[i]);
            max_end_off = MAX(max_end_off, st_end_all[i+1]);
        }
        NCI_Free(st_end_all);
    }

    if (fh->hints->romio_cb_read == PNCIO_HINT_DISABLE ||
        (!interleave_count && fh->hints->romio_cb_read == PNCIO_HINT_AUTO)) {
        /* switch to perform independent read */

        if (buf_view.size == 0) /* zero-size request */
            return 0;

        if (buf_view.count <= 1 && fh->file_view.count <= 1)
            /* Both buf_view and file_view are contiguous. */
            return PNCIO_UFS_read_contig(fh, buf, buf_view.size,
                                         fh->file_view.off[0]);
        else
            return PNCIO_UFS_read_indep(fh, buf, buf_view);
    }

    /* We now proceed to perform two-phase I/O. At first, a call to
     * PNCIO_Calc_file_domains() to calculate the file domains assigned to each
     * I/O aggregator. fh->hints->cb_nodes is the number of aggregators.
     * The aggregate access region of this collective read call is divided
     * among all aggregators into a set of disjoined file domains. A file
     * domain (denoted as 'fd') is the set of file regions an aggregator is
     * responsible for their file access. Thus, a file domain is only relevant
     * to I/O aggregators.
     *
     * PNCIO_Calc_file_domains() set the following 3 variables:
     * fd_end[cb_nodes]   - holds the ending byte of file domains.
     * min_st_off - starting file offset of the aggregate access region
     * max_end_off - last byte offset of the aggregate access region
     *
     * fd_end[] are indexed by an aggregator number; this needs to be mapped
     * to an actual rank in the communicator later.
     */
    PNCIO_Calc_file_domains(fh->hints->cb_nodes, fh->hints->striping_unit,
                            min_st_off, max_end_off, &fd_end, &fd_size);

    /* PNCIO_Calc_my_req() calculates where the portions of this rank's
     * requests fall into aggregator's file domains. When returned, it set the
     * following variables:
     * count_my_req_procs - number of aggregators for which this rank has
     *      requests fall into their file domains
     * count_my_req_per_proc - count of requests for each aggregator, indexed
     *      by rank of the process
     * my_req[nprocs] - array of data structures describing the requests to be
     *      performed by each aggregator.
     * buf_idx[nprocs] - array of locations into which data in the user buffer
     *      that can be directly used to perform read; this is only valid when
     *      the user buffer is contiguous.
     */
    PNCIO_Calc_my_req(fh, min_st_off, fd_end, fd_size, &count_my_req_procs,
                      &count_my_req_per_proc, &my_req, &buf_idx);

    /* PNCIO_Calc_others_req() is only relevant to the I/O aggregators. Based
     * on everyone's my_req, PNCIO_Calc_others_req() calculates what requests
     * of all processes fall into this aggregator's file domain. This
     * subroutine sets the following variables:
     * count_others_req_procs - number of processes whose requests fall into
     *      this aggregator's file domain (including this rank itself)
     * count_others_req_per_proc[i] - how many non-contiguous requests of rank
     *      i fall into this aggregator's file domain.
     */
    PNCIO_Calc_others_req(fh, count_my_req_procs, count_my_req_per_proc,
                          my_req, &count_others_req_procs,
                          &count_others_req_per_proc, &others_req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->read_timing[1] += MPI_Wtime() - curT;
#endif

    /* read data in sizes of no more than collective buffer size, communicate
     * to exchange read data, and fill user buf.
     */
    r_len = Read_and_exch(fh, buf, buf_view, others_req, min_st_off, fd_size,
                          fd_end, buf_idx);
    if (r_len > 0) total_r_len += r_len;

    /* free all memory allocated for collective I/O */
    PNCIO_Free_my_req(count_my_req_per_proc, my_req, buf_idx);
    PNCIO_Free_others_req(count_others_req_per_proc, others_req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->read_timing[0] += MPI_Wtime() - curT;
#endif

    return (r_len < 0) ? r_len : total_r_len;
}

