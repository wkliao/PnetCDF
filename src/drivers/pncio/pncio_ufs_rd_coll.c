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
        memcpy(((char *) buf) + user_buf_idx,                       \
               &(recv_buf[aggr][recv_buf_idx[aggr]]), size_in_buf); \
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

static
void Fill_user_buffer(PNCIO_File  *fh,
                      void        *buf,
                      PNCIO_View   buf_view,
                      char       **recv_buf,        /* IN: [nprocs] */
                      MPI_Count   *recv_size,       /* IN: [nprocs] */
                      MPI_Count   *recd_from_proc,  /* IN/OUT: [nprocs] */
                      MPI_Offset   min_st_off,
                      MPI_Offset   fd_size,
                      MPI_Offset  *fd_start,        /* IN: [cb_nodes] */
                      MPI_Offset  *fd_end)          /* IN: [cb_nodes] */
{
    /* This subroutine is only called when buf_view is not contiguous. */

    int i, nprocs, aggr, buf_indx;
    MPI_Offset buf_rem, size_in_buf, buf_incr, size;
    MPI_Offset off, user_buf_idx;
    MPI_Offset len, rem_len;
    MPI_Count j, *curr_from, *done_from, *recv_buf_idx;

    MPI_Comm_size(fh->comm, &nprocs);

    /* curr_from[nprocs] - amount of data received from each rank that has
     *      already been accounted for so far.
     * done_from[nprocs] - amount of data already received from each rank
     *      and filled into user buffer in previous iterations.
     * user_buf_idx - current location in user buffer
     * recv_buf_idx[nprocs] = current location in recv_buf of each rank */
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

        /* this request may span the file domains of more than one process */
        while (rem_len != 0) {
            len = rem_len;

            /* NOTE: len value will be modified by PNCIO_Calc_aggregator() to
             * be no more than the single file domain region that aggregator
             * 'aggr' is responsible for.
             */
            aggr = PNCIO_Calc_aggregator(fh, off, min_st_off, &len, fd_size,
                                         fd_end);

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
                     MPI_Count    *recd_from_proc, /* [nprocs] */
                     MPI_Offset    min_st_off,
                     MPI_Offset    fd_size,
                     MPI_Offset   *fd_start,       /* IN: [cb_nodes] */
                     MPI_Offset   *fd_end,         /* IN: [cb_nodes] */
                     PNCIO_Access *others_req,     /* IN: [nprocs] */
                     MPI_Aint     *buf_idx,        /* IN: [nprocs] */
                     MPI_Aint     *recved_bytes)   /* OUT: */
{
    int i, nprocs, myrank, nprocs_recv, nprocs_send;
    char **recv_buf = NULL;
    size_t memLen;
    MPI_Count j;
    MPI_Request *requests;
    MPI_Datatype send_type;
    MPI_Status *statuses;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double curT = MPI_Wtime();
#endif

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

/* exchange send_size info so that each process knows how much to
   receive from whom and how much memory to allocate. */

    MPI_Alltoall(send_size, 1, MPI_COUNT, recv_size, 1, MPI_COUNT, fh->comm);

    nprocs_recv = 0;
    nprocs_send = 0;
    memLen = 0;
    for (i=0; i<nprocs; i++) {
        memLen += recv_size[i];
        if (recv_size[i])
            nprocs_recv++;
        if (send_size[i])
            nprocs_send++;
    }

    requests = (MPI_Request *)
        NCI_Malloc((nprocs_send + nprocs_recv + 1) * sizeof(MPI_Request));
/* +1 to avoid a 0-size malloc */

    /* Post recvs. If buf_view is contiguous, data can be directly received
     * into user buf at location given by buf_idx. Otherwise, use recv_buf to
     * receive.
     */
    j = 0; // think of this as a counter of non-zero sends/recs
    if (buf_view.count <= 1) {
        for (i=0; i<nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(((char *) buf) + buf_idx[i], recv_size[i],
                            MPI_BYTE, i, 0, fh->comm, requests + j);
#else
                MPI_Irecv(((char *) buf) + buf_idx[i], recv_size[i],
                            MPI_BYTE, i, 0, fh->comm, requests + j);
#endif
                j++;
                buf_idx[i] += recv_size[i];
            }
        }
    } else {
        /* allocate memory for recv_buf and post receives */
        recv_buf = (char **) NCI_Malloc(nprocs * sizeof(char *));
        recv_buf[0] = (char *) NCI_Malloc(memLen);
        for (i=1; i<nprocs; i++)
            recv_buf[i] = recv_buf[i - 1] + recv_size[i - 1];

        j = 0;
        for (i=0; i<nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Irecv_c(recv_buf[i], recv_size[i], MPI_BYTE, i,
                            0, fh->comm, requests + j);
#else
                MPI_Irecv(recv_buf[i], recv_size[i], MPI_BYTE, i,
                            0, fh->comm, requests + j);
#endif
                j++;
            }
        }
    }

/* create derived datatypes and send data */

    j = 0;
    for (i=0; i<nprocs; i++) {
        if (send_size[i]) {
            /* take care if the last off-len pair is a partial send */
            MPI_Offset tmp = 0;
            MPI_Count k = 0;
            if (partial_send[i]) {
                k = start_pos[i] + count[i] - 1;
                tmp = others_req[i].lens[k];
                others_req[i].lens[k] = partial_send[i];
            }
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(count[i],
                                       &(others_req[i].lens[start_pos[i]]),
                                       &(others_req[i].mem_ptrs[start_pos[i]]),
                                       MPI_BYTE, &send_type);
#else
            MPI_Type_create_hindexed(count[i],
                                     &(others_req[i].lens[start_pos[i]]),
                                     &(others_req[i].mem_ptrs[start_pos[i]]),
                                     MPI_BYTE, &send_type);
#endif
            /* absolute displacement; use MPI_BOTTOM in send */
            MPI_Type_commit(&send_type);
            MPI_Isend(MPI_BOTTOM, 1, send_type, i, 0,
                      fh->comm, requests + nprocs_recv + j);
            MPI_Type_free(&send_type);
            if (partial_send[i])
                others_req[i].lens[k] = tmp;
            j++;
        }
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->read_timing[4] += MPI_Wtime() - curT;
#endif


    /* +1 to avoid a 0-size malloc */
    statuses = (MPI_Status *) NCI_Malloc((nprocs_send + nprocs_recv + 1) * sizeof(MPI_Status));

    /* wait on the receives */
    if (nprocs_recv) {
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        curT = MPI_Wtime();
#endif
        MPI_Waitall(nprocs_recv, requests, statuses);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        if (fh->is_agg) fh->read_timing[3] += MPI_Wtime() - curT;
#endif

        *recved_bytes = 0;
        j = 0;
        for (i=0; i<nprocs; i++) {
            if (recv_size[i]) {
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count count_recved;
                MPI_Get_count_c(&statuses[j], MPI_BYTE, &count_recved);
#else
                int count_recved;
                MPI_Get_count(&statuses[j], MPI_BYTE, &count_recved);
#endif
                *recved_bytes += count_recved;
                j++;
            }
        }

        /* if noncontiguous, to the copies from the recv buffers */
        if (buf_view.count > 1)
            Fill_user_buffer(fh, buf, buf_view, recv_buf, recv_size,
                             recd_from_proc, min_st_off, fd_size, fd_start,
                             fd_end);
    }

    /* wait on the sends */
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    curT = MPI_Wtime();
#endif
#ifdef HAVE_MPI_STATUSES_IGNORE
    MPI_Waitall(nprocs_send, requests + nprocs_recv, MPI_STATUSES_IGNORE);
#else
    MPI_Waitall(nprocs_send, requests + nprocs_recv, statuses + nprocs_recv);
#endif
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->read_timing[3] += MPI_Wtime() - curT;
#endif

    NCI_Free(statuses);
    NCI_Free(requests);

    if (buf_view.count > 1) {
        NCI_Free(recv_buf[0]);
        NCI_Free(recv_buf);
    }
}

static
MPI_Offset Read_and_exch(PNCIO_File   *fh,
                         void         *buf,
                         PNCIO_View    buf_view,
                         PNCIO_Access *others_req, /* IN: [nprocs] */
                         MPI_Offset    min_st_off,
                         MPI_Offset    fd_size,
                         MPI_Offset   *fd_start, /* IN: [cb_nodes] */
                         MPI_Offset   *fd_end,   /* IN: [cb_nodes] */
                         MPI_Aint     *buf_idx)  /* IN: [nprocs] */
{
/* Read in sizes of no more than coll_bufsize, an info parameter.
   Send data to appropriate processes.
   Place recd. data in user buf.
   The idea is to reduce the amount of extra memory required for
   collective I/O. If all data were read all at once, which is much
   easier, it would require temp space more than the size of user_buf,
   which is often unacceptable. For example, to read a distributed
   array from a file, where each local array is 8Mbytes, requiring
   at least another 8Mbytes of temp space is unacceptable. */

    char *read_buf = NULL, *tmp_buf;
    int i, m, ntimes, max_ntimes, nprocs, myrank;
    MPI_Offset st_loc = -1, end_loc = -1, off, done, real_off;
    MPI_Count j, *curr_offlen_ptr, *count, *send_size, *recv_size;
    MPI_Count *partial_send, *recd_from_proc, *start_pos;
    /* Not convinced end_loc-st_loc couldn't be > int, so make these offsets */
    MPI_Offset real_size, size, for_curr_iter, for_next_iter;
    MPI_Aint coll_bufsize, actual_recved_bytes = 0;
    MPI_Offset r_len;

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

/* calculate the number of reads of size coll_bufsize
   to be done by each process and the max among all processes.
   That gives the no. of communication phases as well.
   coll_bufsize is obtained from the hints object. */

    coll_bufsize = fh->hints->cb_buffer_size;

    /* grab some initial values for st_loc and end_loc */
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

    /* calculate ntimes, the number of times this process must perform I/O
     * operations in order to complete all the requests it has received.
     * the need for multiple I/O operations comes from the restriction that
     * we only use coll_bufsize bytes of memory for internal buffering.
     */
    if ((st_loc == -1) && (end_loc == -1)) {
        /* this process does no I/O. */
        ntimes = 0;
    } else {
        /* ntimes=ceiling_div(end_loc - st_loc + 1, coll_bufsize) */
        ntimes = (int) ((end_loc - st_loc + coll_bufsize) / coll_bufsize);
    }

    MPI_Allreduce(&ntimes, &max_ntimes, 1, MPI_INT, MPI_MAX, fh->comm);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fh->read_counter[0] = MAX(fh->read_counter[0], max_ntimes);
#endif

    read_buf = fh->io_buf; /* has been allocated at open time */

    curr_offlen_ptr = NCI_Calloc(nprocs * 7, sizeof(*curr_offlen_ptr));
    /* its use is explained below. calloc initializes to 0. */

    count = curr_offlen_ptr + nprocs;
    /* to store count of how many off-len pairs per proc are satisfied
     * in an iteration. */

    partial_send = count + nprocs;
    /* if only a portion of the last off-len pair is sent to a process
     * in a particular iteration, the length sent is stored here.
     * calloc initializes to 0. */

    send_size = partial_send + nprocs;
    /* total size of data to be sent to each proc. in an iteration */

    recv_size = send_size + nprocs;
    /* total size of data to be recd. from each proc. in an iteration.
     * Of size nprocs so that I can use MPI_Alltoall later. */

    recd_from_proc = recv_size + nprocs;
    /* amount of data recd. so far from each proc. Used in Fill_user_buffer.
     * initialized to 0 here. */

    start_pos = recd_from_proc + nprocs;
    /* used to store the starting value of curr_offlen_ptr[i] in
     * this iteration */

    done = 0;
    off = st_loc;
    for_curr_iter = for_next_iter = 0;

    for (m = 0; m < ntimes; m++) {
        /* read buf of size coll_bufsize (or less) */
        /* go through all others_req and check if any are satisfied
         * by the current read */

        /* since MPI guarantees that displacements in filetypes are in
         * monotonically nondecreasing order, I can maintain a pointer
         * (curr_offlen_ptr) to current off-len pair for each process in
         * others_req and scan further only from there. There is still a
         * problem of filetypes such as:  (1, 2, 3 are not process nos. They
         * are just numbers for three chunks of data, specified by a filetype.)
         *
         * 1  -------!--
         * 2    -----!----
         * 3       --!-----
         *
         * where ! indicates where the current read_size limitation cuts
         * through the filetype.  I resolve this by reading up to !, but
         * filling the communication buffer only for 1. I copy the portion left
         * over for 2 into a tmp_buf for use in the next iteration. i.e., 2 and
         * 3 will be satisfied in the next iteration. This simplifies filling
         * in the user's buf at the other end, as only one off-len pair with
         * incomplete data will be sent. I also don't need to send the
         * individual offsets and lens along with the data, as the data is
         * being sent in a particular order.
         */

        /* off = start offset in the file for the data actually read in this
         *       iteration
         * size = size of data read corresponding to off
         * real_off = off minus whatever data was retained in memory from
         *            previous iteration for cases like 2, 3 illustrated above
         * real_size = size plus the extra corresponding to real_off
         * req_off = off in file for a particular contiguous request minus
         *           what was satisfied in previous iteration
         * req_size = size corresponding to req_off
         */

        size = MIN(coll_bufsize, end_loc - st_loc + 1 - done);
        bool flag = false;
        for (i = 0; i < nprocs; i++) {
            if (others_req[i].count) {
                for (j=curr_offlen_ptr[i]; j<others_req[i].count; j++) {
                    MPI_Offset req_off;
                    if (partial_send[i]) {
                        req_off = others_req[i].offsets[j] + partial_send[i];
                    } else {
                        req_off = others_req[i].offsets[j];
                    }
                    if (req_off < off + size) {
                        flag = true;
                    }
                }
            }
        }
        if (flag) {
            /* This should be only reached by I/O aggregators only */
            r_len = PNCIO_UFS_read_contig(fh, read_buf + for_curr_iter, size,
                                          off);
            if (r_len < 0) return r_len;
            size = r_len;
        }

        real_off = off - for_curr_iter;
        real_size = size + for_curr_iter;

        for (i = 0; i < nprocs; i++)
            count[i] = send_size[i] = 0;
        for_next_iter = 0;

        for (i = 0; i < nprocs; i++) {
            if (others_req[i].count) {
                start_pos[i] = curr_offlen_ptr[i];
                for (j = curr_offlen_ptr[i]; j < others_req[i].count; j++) {
                    MPI_Offset req_off;
#ifdef HAVE_MPI_LARGE_COUNT
                    MPI_Offset req_len;
#else
                    int req_len;
#endif
                    if (partial_send[i]) {
                        /* this request may have been partially
                         * satisfied in the previous iteration. */
                        req_off = others_req[i].offsets[j] + partial_send[i];
                        req_len = others_req[i].lens[j] - partial_send[i];
                        partial_send[i] = 0;
                        /* modify the off-len pair to reflect this change */
                        others_req[i].offsets[j] = req_off;
                        others_req[i].lens[j] = req_len;
                    } else {
                        req_off = others_req[i].offsets[j];
                        req_len = others_req[i].lens[j];
                    }
                    if (req_off < real_off + real_size) {
                        count[i]++;
                        MPI_Aint addr;
                        MPI_Get_address(read_buf + req_off - real_off, &addr);
                        others_req[i].mem_ptrs[j] = addr;
                        send_size[i] += (MIN(real_off + real_size - req_off, req_len));

                        if (real_off + real_size - req_off < req_len) {
                            partial_send[i] = (real_off + real_size - req_off);
                            if ((j + 1 < others_req[i].count) &&
                                (others_req[i].offsets[j + 1] < real_off + real_size)) {
                                /* this is the case illustrated in the
                                 * figure above. */
                                for_next_iter = MAX(for_next_iter,
                                                        real_off + real_size -
                                                        others_req[i].offsets[j + 1]);
                                /* max because it must cover requests
                                 * from different processes */
                            }
                            break;
                        }
                    } else
                        break;
                }
                curr_offlen_ptr[i] = j;
            }
        }

        for_curr_iter = for_next_iter;

        MPI_Aint recved_bytes = 0;
        R_Exchange_data(fh, buf, buf_view, send_size, recv_size, count,
                        start_pos, partial_send, recd_from_proc,
                        min_st_off, fd_size, fd_start, fd_end,
                        others_req, buf_idx, &recved_bytes);
        actual_recved_bytes += recved_bytes;


        if (for_next_iter) {
            tmp_buf = (char *) NCI_Malloc(for_next_iter);
            memcpy(tmp_buf, read_buf + real_size - for_next_iter, for_next_iter);
            NCI_Free(fh->io_buf);
            fh->io_buf = (char *) NCI_Malloc(for_next_iter + coll_bufsize);
            memcpy(fh->io_buf, tmp_buf, for_next_iter);
            read_buf = fh->io_buf;
            NCI_Free(tmp_buf);
        }

        off += size;
        done += size;
    }

    for (i = 0; i < nprocs; i++)
        count[i] = send_size[i] = 0;
    for (m = ntimes; m < max_ntimes; m++) {
        /* nothing to send, but check for recv. */
        MPI_Aint recved_bytes = 0;
        R_Exchange_data(fh, buf, buf_view, send_size, recv_size, count,
                        start_pos, partial_send, recd_from_proc,
                        min_st_off, fd_size, fd_start, fd_end,
                        others_req, buf_idx, &recved_bytes);
        actual_recved_bytes += recved_bytes;
    }

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
    MPI_Offset *fd_start=NULL, *fd_end=NULL, r_len, total_r_len=0;

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
     * fd_start[cb_nodes] - holds the starting byte location of file domains.
     * fd_end[cb_nodes]   - holds the ending byte of file domains.
     * min_st_off - holds the minimum byte location that will be accessed.
     *
     * Both fd_start[] and fd_end[] are indexed by an aggregator number; this
     * needs to be mapped to an actual rank in the communicator later.
     */
    PNCIO_Calc_file_domains(min_st_off, max_end_off, fh->hints->cb_nodes,
                            &fd_start, &fd_end, &fd_size,
                            fh->hints->striping_unit);

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
                          fd_start, fd_end, buf_idx);
    if (r_len > 0) total_r_len += r_len;

    /* free all memory allocated for collective I/O */
    PNCIO_Free_my_req(count_my_req_per_proc, my_req, buf_idx);
    PNCIO_Free_others_req(count_others_req_per_proc, others_req);

    NCI_Free(fd_start);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->read_timing[0] += MPI_Wtime() - curT;
#endif

    return (r_len < 0) ? r_len : total_r_len;
}

