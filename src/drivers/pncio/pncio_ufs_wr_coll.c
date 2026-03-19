/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "pncio.h"

#define BUF_INCR {                                      \
    while (buf_incr) {                                  \
        size_in_buf = MIN(buf_incr, buf_rem);           \
        user_buf_idx += size_in_buf;                    \
        buf_rem -= size_in_buf;                         \
        buf_incr -= size_in_buf;                        \
        if (buf_incr > 0 && buf_rem == 0) {             \
            buf_indx++;                                 \
            user_buf_idx = buf_view.off[buf_indx];      \
            buf_rem = buf_view.len[buf_indx];           \
        }                                               \
    }                                                   \
}

#define BUF_COPY {                                      \
    while (size) {                                      \
        size_in_buf = MIN(size, buf_rem);               \
        memcpy(&send_buf[aggr][send_buf_idx[aggr]],     \
               (char*)buf + user_buf_idx, size_in_buf); \
        send_buf_idx[aggr] += size_in_buf;              \
        user_buf_idx += size_in_buf;                    \
        buf_rem -= size_in_buf;                         \
        size -= size_in_buf;                            \
        buf_incr -= size_in_buf;                        \
        if (size > 0 && buf_rem == 0) {                 \
            buf_indx++;                                 \
            user_buf_idx = buf_view.off[buf_indx];      \
            buf_rem = buf_view.len[buf_indx];           \
        }                                               \
    }                                                   \
    BUF_INCR                                            \
}

/*----< fill_send_buffer() >-------------------------------------------------*/
/* This subroutine is only called when buf_view is not contiguous. */
static
void fill_send_buffer(PNCIO_File       *fh,
                      const void       *buf,
                      PNCIO_View        buf_view,
                      MPI_Offset        min_st_off,
                      MPI_Offset        fd_size,
                      const MPI_Offset *fd_end,       /* IN: [cb_nodes] */
                      const MPI_Count  *send_size,    /* IN: [nprocs] */
                      MPI_Count        *sent_to_proc, /* IN/OUT: [nprocs] */
                      char *const      *send_buf,     /* OUT: [nprocs] */
                      MPI_Request      *reqs)         /* OUT: [nprocs] */
{
    int i, k, nprocs, myrank, aggr;
    MPI_Offset buf_indx, buf_rem, size_in_buf, buf_incr, size;
    MPI_Offset off, len, rem_len, user_buf_idx;
    MPI_Count j, *curr_to, *done_to, *send_buf_idx;

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

    /* curr_to[nprocs] - amount of data sent to each rank that has already
     *      been accounted for so far.
     * done_to[nprocs] - amount of data already sent to each rank in previous
     *      round.
     * user_buf_idx - current location in user buffer
     * send_buf_idx[nprocs] = current location in send_buf of each rank
     */
    curr_to = NCI_Malloc(sizeof(MPI_Count) * nprocs * 3);
    done_to = curr_to + nprocs;
    send_buf_idx = done_to + nprocs;

    for (i=0; i<nprocs; i++) {
        send_buf_idx[i] = curr_to[i] = 0;
        done_to[i] = sent_to_proc[i];
    }

    /* buf_indx - index buf_view's offset-length pairs being processed
     * buf_rem - remaining length of the current offset-length pair
     */
    user_buf_idx = buf_view.off[0];
    buf_indx = 0;
    buf_rem = buf_view.len[0];

    k = 0;
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

            if (send_buf_idx[aggr] < send_size[aggr]) {
                if (curr_to[aggr] + len > done_to[aggr]) {
                    if (done_to[aggr] > curr_to[aggr]) {
                        size = MIN(curr_to[aggr] + len - done_to[aggr],
                                   send_size[aggr] - send_buf_idx[aggr]);
                        buf_incr = done_to[aggr] - curr_to[aggr];
                        BUF_INCR
                        buf_incr = curr_to[aggr] + len - done_to[aggr];
                        curr_to[aggr] = done_to[aggr] + size;
                        BUF_COPY
                    } else {
                        size = MIN(len, send_size[aggr] - send_buf_idx[aggr]);
                        buf_incr = len;
                        curr_to[aggr] += size;
                        BUF_COPY
                    }
                    if (send_buf_idx[aggr] == send_size[aggr] &&
                        aggr != myrank) {
#if MPI_VERSION >= 4
                        MPI_Isend_c(send_buf[aggr], send_size[aggr], MPI_BYTE,
                                    aggr, 0, fh->comm, &reqs[k++]);
#else
                        MPI_Isend(send_buf[aggr], send_size[aggr], MPI_BYTE,
                                  aggr, 0, fh->comm, &reqs[k++]);
#endif
                    }
                } else {
                    curr_to[aggr] += len;
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
        if (send_size[i])
            sent_to_proc[i] = curr_to[i];

    NCI_Free(curr_to);
}

static MPI_Offset
W_Exchange_data(PNCIO_File       *fh,
                const void       *buf,
                char             *write_buf,
                PNCIO_View        buf_view,
                MPI_Count        *send_size,    /* OUT: [nprocs] */
                const MPI_Count  *recv_size,    /* IN: [nprocs] */
                MPI_Offset        rem_off,
                MPI_Count         rem_size,
                const MPI_Count  *count,        /* IN: [nprocs] */
                const MPI_Count  *start_pos,    /* IN: [nprocs] */
                const MPI_Count  *partial_recv, /* IN: [nprocs] */
                MPI_Count        *sent_to_proc, /* IN/OUT: [nprocs] */
                MPI_Offset        min_st_off,
                MPI_Offset        fd_size,
                const MPI_Offset *fd_end,       /* IN: [cb_nodes] */
                PNCIO_Access     *others_req,   /* IN/OUT: [nprocs] */
                MPI_Aint         *buf_idx)      /* IN/OUT: [nprocs] */

{
    char **send_buf = NULL;
    int i, j, nprocs, myrank, err=NC_NOERR;
    int nrecvs, nsends, num_rtypes, nreqs, hole;
    MPI_Request *reqs, *send_req;
    MPI_Datatype *recv_types, self_recv_type=MPI_DATATYPE_NULL;
    MPI_Count sum, *srt_len=NULL, *tmp_len;
    MPI_Offset *srt_off=NULL;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double curT = MPI_Wtime();
#endif

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

    /* Exchange recv_size info so that each aggregator knows how much to send
     * to whom and how much memory to allocate.
     */
    MPI_Alltoall(recv_size, 1, MPI_COUNT, send_size, 1, MPI_COUNT, fh->comm);

    /* construct derived datatypes for recv */
    recv_types = (MPI_Datatype*) NCI_Malloc(sizeof(MPI_Datatype) * nprocs);

    tmp_len = NCI_Malloc(sizeof(MPI_Count) * nprocs);
    j = 0;
    nsends = 0;
    nrecvs = 0;
    sum = 0;
    for (i=0; i<nprocs; i++) {
        sum += count[i];
        if (send_size[i])
            nsends++;

        if (recv_size[i]) {
            MPI_Datatype *dtype;

            nrecvs++;
            dtype = (i != myrank) ? (recv_types + j) : (&self_recv_type);

            if (partial_recv[i]) {
                /* take care if the last off-len pair is a partial recv */
                MPI_Count k = start_pos[i] + count[i] - 1;
                tmp_len[i] = others_req[i].lens[k];
                others_req[i].lens[k] = partial_recv[i];
            }
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(count[i],
                                       &(others_req[i].lens[start_pos[i]]),
                                       &(others_req[i].mem_ptrs[start_pos[i]]),
                                       MPI_BYTE, dtype);
#else
            MPI_Type_create_hindexed(count[i],
                                     &(others_req[i].lens[start_pos[i]]),
                                     &(others_req[i].mem_ptrs[start_pos[i]]),
                                     MPI_BYTE, dtype);
#endif
            /* absolute displacements; use MPI_BOTTOM in recv */
            MPI_Type_commit(dtype);
            if (i != myrank)
                j++;
        }
    }
    num_rtypes = j;     /* number of non-self receive datatypes created */

    /* To avoid a read-modify-write, check if there are holes in the data to
     * be written. For this, merge the (sorted) offset lists others_req using
     * a heap-merge sort.
     */

    if (sum) {
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        double timing = MPI_Wtime();
#endif
        srt_off = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * sum);
        srt_len = (MPI_Count*) NCI_Malloc(sizeof(MPI_Count) * sum);

/* TODO: PNCIO_Heap_merge is expensive, borrow codes from ad_lustre_wrcoll.c to skip it when possible */

    /* Skip hole checking if there is no write data by this aggregator */
        PNCIO_Heap_merge(others_req, count, srt_off, srt_len, start_pos,
                         nprocs, nrecvs, sum);
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
        if (fh->is_agg) fh->write_timing[5] += MPI_Wtime() - timing;
#endif
    }

    /* for partial recvs, restore original lengths */
    for (i=0; i<nprocs; i++)
        if (partial_recv[i])
            others_req[i].lens[start_pos[i] + count[i] - 1] = tmp_len[i];

    NCI_Free(tmp_len);

    /* Check if there are any holes. If yes, must do read-modify-write.
     * Holes can occur in three places. 'middle' is what you'd expect: the
     * processes are operating on noncontiguous data. However, holes can also
     * show up at the beginning or end of the file domain. Missing these holes
     * would result in us writing more data than received by everyone else.
     */
    hole = 0;
    if (sum) {
        if (rem_off != srt_off[0])  /* hole at the front */
            hole = 1;
        else {  /* coalesce the sorted offset-length pairs */
            MPI_Count k, new_len;
            for (k=1; k<sum; k++) {
                if (srt_off[k] <= srt_off[0] + srt_len[0]) {
                    new_len = srt_off[k] + srt_len[k] - srt_off[0];
                    if (new_len > srt_len[0])
                        srt_len[0] = new_len;
                } else
                    break;
            }
            if (i < sum || rem_size != srt_len[0])  /* hole in middle or end */
                hole = 1;
        }

        NCI_Free(srt_off);
        NCI_Free(srt_len);
    }

    if (nrecvs && hole) {
        MPI_Offset r_len;
        r_len = PNCIO_UFS_read_contig(fh, write_buf, rem_size, rem_off);
        if (r_len < 0) return r_len;
    }

    if (fh->atomicity) {
        /* nreqs is the number of Isend and Irecv to be posted */
        nreqs = (send_size[myrank]) ? (nsends - 1) : nsends;
        reqs = (MPI_Request*) NCI_Malloc(sizeof(MPI_Request) * nreqs);
        send_req = reqs;
    } else {
        nreqs = nsends + nrecvs;
        if (send_size[myrank])  /* NO send to and recv from self */
            nreqs -= 2;
        reqs = (MPI_Request*) NCI_Malloc(sizeof(MPI_Request) * nreqs);

        /* post receives */
        j = 0;
        for (i=0; i<nprocs; i++) {
            if (recv_size[i] == 0)
                continue;
            if (i != myrank) {
                MPI_Irecv(MPI_BOTTOM, 1, recv_types[j], i, 0, fh->comm,
                          &reqs[j]);
                j++;
            } else if (buf_view.count <= 1) {
                /* sen/recv to/from self uses MPI_Unpack() */
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count pos=0;
                MPI_Unpack_c((char*)buf + buf_idx[i], recv_size[i], &pos,
                             write_buf, 1, self_recv_type, MPI_COMM_SELF);
#else
                int pos = 0;
                assert(recv_size[i] < INT_MAX);
                MPI_Unpack((char*)buf + buf_idx[i], (int)recv_size[i], &pos,
                           write_buf, 1, self_recv_type, MPI_COMM_SELF);
#endif
                buf_idx[i] += recv_size[i];
            }
        }
        send_req = reqs + j;
    }

    /* Post nonblocking send calls. If buf_view is contiguous, data can be
     * directly sent from user buf at location given by buf_idx. Otherwise,
     * allocate send_buf and use it to send.
     */
    if (buf_view.count <= 1) {
        j = 0;
        for (i=0; i<nprocs; i++)
            if (send_size[i] && i != myrank) {
                assert(buf_idx[i] != -1);
#if MPI_VERSION >= 4
                MPI_Isend_c((char*)buf + buf_idx[i], send_size[i], MPI_BYTE,
                            i, 0, fh->comm, &send_req[j++]);
#else
                MPI_Isend((char*)buf + buf_idx[i], send_size[i], MPI_BYTE,
                          i, 0, fh->comm, &send_req[j++]);
#endif
                buf_idx[i] += send_size[i];
            }
    } else if (nsends) {
        /* buf_view is not contiguous */
        size_t msgLen = 0;
        for (i=0; i<nprocs; i++)
            msgLen += send_size[i];
        send_buf = (char**) NCI_Malloc(sizeof(char*) * nprocs);
        send_buf[0] = (char*) NCI_Malloc(msgLen);
        for (i=1; i<nprocs; i++)
            send_buf[i] = send_buf[i - 1] + send_size[i - 1];

        fill_send_buffer(fh, buf, buf_view, min_st_off, fd_size, fd_end,
                         send_size, sent_to_proc, send_buf, send_req);

        /* the send is done in fill_send_buffer() */
    }

    if (fh->atomicity) {
        /* In atomic mode, we must use blocking receives to receive data in the
         * same increasing order of MPI process rank IDs.
         */
        j = 0;
        for (i=0; i<nprocs; i++) {
            if (recv_size[i] == 0)
                continue;
            if (i != myrank) {
                MPI_Status st;
                MPI_Recv(MPI_BOTTOM, 1, recv_types[j++], i, 0, fh->comm, &st);
            } else {
                /* sen/recv to/from self uses MPI_Unpack() */
                char *ptr = (buf_view.count <= 1) ? (char*)buf + buf_idx[i]
                                                  : send_buf[i];
assert(self_recv_type != MPI_DATATYPE_NULL);

#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count pos=0;
                MPI_Unpack_c(ptr, recv_size[i], &pos, write_buf, 1,
                             self_recv_type, MPI_COMM_SELF);
#else
                int pos = 0;
                assert(recv_size[i] < INT_MAX);
                MPI_Unpack(ptr, (int)recv_size[i], &pos, write_buf, 1,
                           self_recv_type, MPI_COMM_SELF);
#endif
                buf_idx[i] += recv_size[i];
            }
        }
    } else if (buf_view.count > 1 && recv_size[myrank]) {

assert(self_recv_type != MPI_DATATYPE_NULL);

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count pos=0;
        MPI_Unpack_c(send_buf[myrank], recv_size[myrank], &pos, write_buf,
                     1, self_recv_type, MPI_COMM_SELF);
#else
        int pos = 0;
        assert(recv_size[myrank] < INT_MAX);
        MPI_Unpack(send_buf[myrank], (int)recv_size[myrank], &pos, write_buf,
                   1, self_recv_type, MPI_COMM_SELF);
#endif
    }

    for (i=0; i<num_rtypes; i++)
        MPI_Type_free(recv_types + i);
    NCI_Free(recv_types);

    if (self_recv_type != MPI_DATATYPE_NULL)
        MPI_Type_free(&self_recv_type);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->write_timing[4] += MPI_Wtime() - curT;
    curT = MPI_Wtime();
#endif

#ifdef HAVE_MPI_STATUSES_IGNORE
    MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
#else
    MPI_Status *sts = (MPI_Status*) NCI_Malloc(sizeof(MPI_Status) * nreqs);
    MPI_Waitall(nreqs, reqs, sts);
    NCI_Free(sts);
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->write_timing[3] += MPI_Wtime() - curT;
#endif

    NCI_Free(reqs);
    if (buf_view.count > 1 && nsends) {
        NCI_Free(send_buf[0]);
        NCI_Free(send_buf);
    }

    return err;
}

/* If successful, it returns the amount written. Otherwise a NetCDF error code
 * (negative value) is returned.
 */
static MPI_Offset
Exch_and_write(PNCIO_File       *fh,
               const void       *buf,
               PNCIO_View        buf_view,
               MPI_Offset        min_st_off,
               MPI_Offset        fd_size,
               const MPI_Offset *fd_end,     /* IN: [cb_nodes] */
               PNCIO_Access     *others_req, /* IN/OUT: [nprocs] */
               MPI_Aint         *buf_idx)    /* IN/OUT: [nprocs] */
{
/* Send data to appropriate processes and write in sizes of no more
   than coll_bufsize.
   The idea is to reduce the amount of extra memory required for
   collective I/O. If all data were written all at once, which is much
   easier, it would require temp space more than the size of user_buf,
   which is often unacceptable. For example, to write a distributed
   array to a file, where each local array is 8 MiB, requiring
   at least another 8 MiB of temp space is unacceptable. */

    /* Not convinced end_loc-st_loc couldn't be > int, so make these offsets */

    char *value, *write_buf = NULL;
    int i, m, ntimes, max_ntimes, nprocs, myrank, do_write, info_flag;
    MPI_Offset rem_size=0, w_len, total_w_len=0;
    MPI_Offset st_loc = -1, end_loc = -1, rem_off, done, req_off;
    MPI_Count j, *curr_offlen_ptr, *send_size, *count, req_len, *recv_size;
    MPI_Count *partial_recv, *sent_to_proc, *start_pos;
    MPI_Aint coll_bufsize;

    /* only I/O errors are currently reported */

    MPI_Comm_size(fh->comm, &nprocs);
    MPI_Comm_rank(fh->comm, &myrank);

/* calculate the number of writes of size coll_bufsize
   to be done by each process and the max among all processes.
   That gives the no. of communication phases as well. */

    value = (char *) NCI_Malloc(MPI_MAX_INFO_VAL + 1);
    MPI_Info_get(fh->info, "cb_buffer_size", MPI_MAX_INFO_VAL, value, &info_flag);
    coll_bufsize = atoi(value);
    NCI_Free(value);

    for (i=0; i<nprocs; i++) {
        if (others_req[i].count) {
            st_loc = others_req[i].offsets[0];
            end_loc = others_req[i].offsets[0];
            break;
        }
    }

    for (i=0; i<nprocs; i++)
        for (j=0; j<others_req[i].count; j++) {
            st_loc = MIN(st_loc, others_req[i].offsets[j]);
            end_loc = MAX(end_loc, (others_req[i].offsets[j]
                                  + others_req[i].lens[j] - 1));
        }

/* ntimes=ceiling_div(end_loc - st_loc + 1, coll_bufsize)*/

    ntimes = (int) ((end_loc - st_loc + coll_bufsize) / coll_bufsize);

    if ((st_loc == -1) && (end_loc == -1)) {
        ntimes = 0;     /* this process does no writing. */
    }

    MPI_Allreduce(&ntimes, &max_ntimes, 1, MPI_INT, MPI_MAX, fh->comm);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fh->write_counter[0] = MAX(fh->write_counter[0], max_ntimes);
#endif

    write_buf = fh->io_buf;

    curr_offlen_ptr = NCI_Calloc(nprocs * 7, sizeof(MPI_Count));
    /* its use is explained below. calloc initializes to 0. */

    count = curr_offlen_ptr + nprocs;
    /* to store count of how many off-len pairs per proc are satisfied
     * in an iteration. */

    partial_recv = count + nprocs;
    /* if only a portion of the last off-len pair is recd. from a process
     * in a particular iteration, the length recd. is stored here.
     * calloc initializes to 0. */

    send_size = partial_recv + nprocs;
    /* total size of data to be sent to each proc. in an iteration.
     * Of size nprocs so that I can use MPI_Alltoall later. */

    recv_size = send_size + nprocs;
    /* total size of data to be recd. from each proc. in an iteration. */

    sent_to_proc = recv_size + nprocs;
    /* amount of data sent to each proc so far. Used in
     * fill_send_buffer(). initialized to 0 here. */

    start_pos = sent_to_proc + nprocs;
    /* used to store the starting value of curr_offlen_ptr[i] in
     * this iteration */

    done = 0;
    rem_off = st_loc;

    for (m=0; m<ntimes; m++) {
        /* go through all others_req and check which will be satisfied
         * by the current write */

        /* Note that MPI guarantees that displacements in filetypes are in
         * monotonically non-decreasing order and that, for writes, the
         * filetypes cannot specify overlapping regions in the file. This
         * simplifies implementation a bit compared to reads. */

        /* rem_off  = start file offset for data actually written in this round
         * rem_size = size of data to be written, corresponding to rem_off
         * req_off  = file offset for a particular offset-length pair minus
         *            what has been satisfied in previous round
         * req_size = size corresponding to req_off
         */

        /* first calculate what should be communicated */

        for (i=0; i<nprocs; i++)
            count[i] = recv_size[i] = 0;

        do_write = 0; /* will change to 1 if any of count[i] becomes > 0 */

        rem_size = MIN(coll_bufsize, end_loc - st_loc + 1 - done);

        MPI_Offset round_end = rem_off + rem_size;

        for (i=0; i<nprocs; i++) {
            if (others_req[i].count == 0)
                continue;

            start_pos[i] = curr_offlen_ptr[i];
            for (j=curr_offlen_ptr[i]; j<others_req[i].count; j++) {
                if (partial_recv[i]) {
                    /* this request may have been partially satisfied in the
                     * previous iteration.
                     */
                    req_off = others_req[i].offsets[j] + partial_recv[i];
                    req_len = others_req[i].lens[j]    - partial_recv[i];
                    partial_recv[i] = 0;
                    /* modify the off-len pair to reflect this change */
                    others_req[i].offsets[j] = req_off;
                    others_req[i].lens[j]    = req_len;
                } else {
                    req_off = others_req[i].offsets[j];
                    req_len = others_req[i].lens[j];
                }

                if (req_off >= round_end)
                    break;

                /* Now req_off < round_end */
                count[i]++;
                do_write = 1;

                if (myrank != i) {
                    MPI_Aint addr;
                    MPI_Get_address(write_buf + req_off - rem_off, &addr);
                    others_req[i].mem_ptrs[j] = addr;
                }
                else
                    others_req[i].mem_ptrs[j] = req_off - rem_off;
                recv_size[i] += MIN(round_end - req_off, req_len);

                if (round_end - req_off < req_len) {
                    partial_recv[i] = (round_end - req_off);

                    if (j + 1 < others_req[i].count &&
                        others_req[i].offsets[j + 1] < round_end) {
                        /* This error should not happen to PnetCDF, as
                         * fileview is checked before entering this
                         * subroutine.
                         */
                        fprintf(stderr, "Filetype specifies overlapping write regions (which is illegal according to the MPI-2 specification\n");
                        /* allow to continue since additional communication
                         * might have to occur
                         */
                        total_w_len = NC_EFILE;
                        goto err_out;
                    }
                    break;
                }
            }
            curr_offlen_ptr[i] = j;
        }

        w_len = W_Exchange_data(fh, buf, write_buf, buf_view, send_size,
                                recv_size, rem_off, rem_size, count, start_pos,
                                partial_recv, sent_to_proc, min_st_off,
                                fd_size, fd_end, others_req, buf_idx);

        if (w_len < 0) {
            total_w_len = w_len;
            goto err_out;
        }
        else
            total_w_len += w_len;

        if (do_write) {
            w_len = PNCIO_UFS_write_contig(fh, write_buf, rem_size, rem_off, 1);
            if (w_len < 0) {
                total_w_len = w_len;
                goto err_out;
            }
            else
                total_w_len += w_len;
        }

        rem_off += rem_size;
        done += rem_size;
    }

    for (i=0; i<nprocs; i++)
        count[i] = recv_size[i] = 0;
    for (m=ntimes; m<max_ntimes; m++) {
        /* nothing to recv, but check for send. */
        w_len = W_Exchange_data(fh, buf, write_buf, buf_view, send_size,
                                recv_size, rem_off, rem_size, count, start_pos,
                                partial_recv, sent_to_proc, min_st_off,
                                fd_size, fd_end, others_req, buf_idx);
        if (w_len < 0) {
            total_w_len = w_len;
            goto err_out;
        }
        else
            total_w_len += w_len;
    }

err_out:
    NCI_Free(curr_offlen_ptr);

    return total_w_len;
}

/*----< PNCIO_UFS_write_coll() >---------------------------------------------*/
MPI_Offset PNCIO_UFS_write_coll(PNCIO_File *fh,
                                const void *buf,
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

    int i, nprocs, rank, interleave_count=0;
    MPI_Aint *buf_idx = NULL;
    MPI_Count *count_my_req_per_proc, count_my_req_procs;
    MPI_Count *count_others_req_per_proc, count_others_req_procs;
    MPI_Offset fd_size, min_st_off, max_end_off;
    MPI_Offset *fd_end=NULL, w_len=0;

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

    /* only check for interleaving if romio_cb_write isn't disabled */
    if (fh->hints->romio_cb_write != PNCIO_HINT_DISABLE) {
        MPI_Offset *st_end_all;

        /* For this process's request, calculate its aggregate access region,
         * representing a range from starting offset till end offset. Note:
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
        if (fh->file_view.size == 0) /* -1 to indicate zero-sized request */
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
            if (st_end_all[i] == -1) /* zero-sized request */
                continue;
            if (st_end_all[i] <  st_end_all[i-1] &&
                st_end_all[i] <= st_end_all[i+1])
                interleave_count++;
            min_st_off  = MIN(min_st_off,  st_end_all[i]);
            max_end_off = MAX(max_end_off, st_end_all[i+1]);
        }
        NCI_Free(st_end_all);
    }

    if (fh->hints->romio_cb_write == PNCIO_HINT_DISABLE ||
        (!interleave_count && fh->hints->romio_cb_write == PNCIO_HINT_AUTO)) {
        /* switch to perform independent write */

        if (buf_view.size == 0) /* zero_sized request */
            return 0;

        if (buf_view.count <= 1 && fh->file_view.count <= 1)
            /* Both buf_view and file_view are contiguous */
            w_len = PNCIO_UFS_write_contig(fh, buf, buf_view.size,
                                           fh->file_view.off[0], 0);
        else
            w_len = PNCIO_UFS_write_indep(fh, buf, buf_view);

        return w_len;
    }

    /* We now proceed to perform two-phase I/O. At first, a call to
     * PNCIO_Calc_file_domains() to calculate the file domains assigned to each
     * I/O aggregator. fh->hints->cb_nodes is the number of aggregators.
     * The aggregate access region of this collective write call is divided
     * among all aggregators into a set of disjoined file domains. A file
     * domain (denoted as 'fd') is the set of file regions an aggregator is
     * responsible for their file access. Thus, a file domain is only relevant
     * to I/O aggregators.
     *
     * PNCIO_Calc_file_domains() set the following 3 variables:
     * fd_end     - holds the ending byte location.
     * min_st_off - holds the minimum byte location that will be accessed.
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
     * my_req[] - array of data structures describing the requests to be
     *      performed by each aggregator.
     * buf_idx[] - array of locations into which data in the user buffer that
     *      can be directly used to perform write; this is only valid when the
     *      user buffer is contiguous.
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
    if (fh->is_agg) fh->write_timing[1] += MPI_Wtime() - curT;
#endif

    /* exchange data and write in sizes of no more than coll_bufsize. */
    w_len = Exch_and_write(fh, buf, buf_view, min_st_off, fd_size, fd_end,
                           others_req, buf_idx);

    /* If this collective write is followed by an independent write, it is
     * possible to have those subsequent writes on other processes race ahead
     * and sneak in before the read-modify-write completes. Below, we carry out
     * a collective communication at the end of this subroutine, so no one can
     * start independent I/O before collective I/O completes.
     *
     * need to do some gymnastics with the error codes so that if something
     * went wrong, all processes report error, but if a process has a more
     * specific error code, we can still have that process report the
     * additional information */

    if (fh->hints->cb_nodes == 1)
        /* If there is only one aggregator, we can perform a less-expensive
         * Bcast().
         */
        MPI_Bcast(&w_len, 1, MPI_OFFSET, fh->hints->aggr_ranks[0], fh->comm);
    else
        MPI_Allreduce(MPI_IN_PLACE, &w_len, 1, MPI_OFFSET, MPI_MIN, fh->comm);

    /* free all memory allocated for collective I/O */
    PNCIO_Free_my_req(count_my_req_per_proc, my_req, buf_idx);
    PNCIO_Free_others_req(count_others_req_per_proc, others_req);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    if (fh->is_agg) fh->write_timing[0] += MPI_Wtime() - curT;
#endif

    /* w_len may not be the same as buf_view.size, because data sieving may
     * write more than requested.
     */
    return buf_view.size;
}

