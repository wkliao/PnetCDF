/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "ad_lustre.h"
#include "adio_extern.h"

#define MEMCPY_UNPACK(x, inbuf, start, count, outbuf) {          \
    int _k;                                                      \
    char *_ptr = (inbuf);                                        \
    MPI_Count   *mem_ptrs = others_req[x].mem_ptrs + (start)[x]; \
    ADIO_Offset *mem_lens = others_req[x].lens     + (start)[x]; \
    for (_k=0; _k<(count)[x]; _k++) {                            \
        memcpy((outbuf) + mem_ptrs[_k], _ptr, mem_lens[_k]);     \
        _ptr += mem_lens[_k];                                    \
    }                                                            \
}

typedef struct {
    MPI_Count    num; /* number of elements in the above off-len list */
    MPI_Count   *len; /* list of write lengths by this rank in round m */
    ADIO_Offset *off; /* list of write offsets by this rank in round m */
} off_len_list;

typedef struct {
    MPI_Count  count; /* number displacement-length pairs */
    MPI_Count *len;   /* [count]: size in bytes */
    MPI_Count *disp;  /* [count]: displacement */
} disp_len_list;

typedef struct {
    MPI_Datatype type;      /* MPI derived datatype */
    MPI_Count    count;     /* number of off-len pairs (blocks) */
    MPI_Offset  *off;       /* array of byte offsets of each block */
    MPI_Offset  *len;       /* array of contiguous block lengths (bytes) */
    MPI_Count    rnd;       /* number of whole type already consumed */
    MPI_Count    idx;       /* index of off-len pairs consumed so far */
    MPI_Aint     rem;       /* remaining amount in the pair to be consumed */
    MPI_Aint     extent;    /* data type extent */
    int          is_contig; /* for fileview, whether file access is contiguous
                             * for buffer, whether user buffer is contiguous
                             * Note this is not whether filetype or buftype
                             * is contiguous or not.
                             */
} Flat_list;

/* prototypes of functions used for collective writes only. */
static void ADIOI_LUSTRE_Exch_and_write(ADIO_File fd,
                                        const void *buf,
                                        Flat_list *flat_bview,
                                        ADIOI_Access *others_req,
                                        ADIOI_Access *my_req,
                                        Flat_list *flat_fview,
                                        ADIO_Offset min_st_loc,
                                        ADIO_Offset max_end_loc,
                                        ADIO_Offset **buf_idx,
                                        int *error_code);
static void ADIOI_LUSTRE_Fill_send_buffer(ADIO_File fd, const void *buf,
                                          Flat_list *flat_fview,
                                          Flat_list *flat_bview,
                                          char **send_buf,
                                          size_t send_total_size,
                                          const MPI_Count *send_size,
                                          char **self_buf,
                                          disp_len_list *send_list);
static void ADIOI_LUSTRE_W_Exchange_data(ADIO_File fd, const void *buf,
                                         char *write_buf,
                                         char **recve_buf,
                                         char **send_buf_ptr,
                                         Flat_list *flat_fview,
                                         Flat_list *flat_bview,
                                         const MPI_Count *send_size,
                                         const MPI_Count *recv_size,
                                         ADIO_Offset range_off,
                                         MPI_Count range_size,
                                         const MPI_Count *count,
                                         const MPI_Count *start_pos,
                                         const ADIOI_Access *others_req,
                                         const ADIO_Offset *buf_idx,
                                         off_len_list *srt_off_len,
                                         disp_len_list *send_list,
                                         disp_len_list *recv_list,
                                         int *error_code);
static
int ADIOI_LUSTRE_Calc_aggregator(ADIO_File fd,
                                 ADIO_Offset off,
                                 ADIO_Offset *len)
{
    ADIO_Offset avail_bytes, stripe_id;

    stripe_id = off / fd->hints->striping_unit;

    avail_bytes = (stripe_id + 1) * fd->hints->striping_unit - off;
    if (avail_bytes < *len) {
        /* The request [off, off+len) has only [off, off+avail_bytes) part
         * falling into aggregator's file domain */
        *len = avail_bytes;
    }
    /* return the index to ranklist[] */
    return (stripe_id % fd->hints->cb_nodes);
}

/*----< ADIOI_LUSTRE_Calc_my_req() >-----------------------------------------*/
/* calculates what portions of the read/write requests of this process fall
 * into the file domains of all I/O aggregators.
 *   IN: flat_fview: this rank's flattened write requests
 *       flat_fview.count: number of noncontiguous offset-length file requests
 *       flat_fview.off[flat_fview.count] file offsets of individual
 *       noncontiguous requests.
 *       flat_fview.len[flat_fview.count] lengths of individual noncontiguous
 *       requests.
 *   IN: buf_is_contig: whether the write buffer is contiguous or not
 *   OUT: my_req_ptr[cb_nodes] offset-length pairs of this process's requests
 *        fall into the file domain of each aggregator
 *   OUT: buf_idx_ptr[cb_nodes] index pointing to the starting location in
 *        user_buf for data to be sent to each aggregator.
 */
static
void ADIOI_LUSTRE_Calc_my_req(ADIO_File fd,
                              Flat_list flat_fview,
                              int buf_is_contig,
                              ADIOI_Access **my_req_ptr,
                              ADIO_Offset **buf_idx)
{
    int aggr, *aggr_ranks, cb_nodes;
    MPI_Count i, l;
    size_t nelems;
    ADIO_Offset avail_len, rem_len, curr_idx, off, *ptr;
    ADIO_Offset *avail_lens;
    ADIOI_Access *my_req;

    cb_nodes = fd->hints->cb_nodes;

    /* my_req[i].count gives the number of contiguous requests of this process
     * that fall in aggregator i's file domain (not process MPI rank i).
     */
    my_req = (ADIOI_Access *) ADIOI_Calloc(cb_nodes, sizeof(ADIOI_Access));
    *my_req_ptr = my_req;

    /* First pass is just to calculate how much space is needed to allocate
     * my_req. Note that flat_fview.count has been calculated way back in
     * ADIOI_Calc_my_off_len()
     */
    aggr_ranks = (int *) ADIOI_Malloc(flat_fview.count * sizeof(int));
    avail_lens = (ADIO_Offset *) ADIOI_Malloc(flat_fview.count * sizeof(ADIO_Offset));

    /* nelems will be the number of offset-length pairs for my_req[] */
    nelems = 0;
    for (i = 0; i < flat_fview.count; i++) {
        /* short circuit offset/len processing if zero-byte read/write. */
        if (flat_fview.len[i] == 0)
            continue;

        off = flat_fview.off[i];
        avail_len = flat_fview.len[i];
        /* ADIOI_LUSTRE_Calc_aggregator() modifies the value of 'avail_len' to
         * the amount that is only covered by the aggr's file domain. The
         * remaining (tail) will continue to be processed to determine to whose
         * file domain it belongs. As ADIOI_LUSTRE_Calc_aggregator() can be
         * expensive for large value of flat_fview.count, we keep a copy of
         * the returned values of 'aggr' and 'avail_len' in aggr_ranks[] and
         * avail_lens[] to be used in the next for loop (not next iteration).
         *
         * Note the returned value in 'aggr' is the index to ranklist[], i.e.
         * the 'aggr'th element of array ranklist[], rather than the
         * aggregator's MPI rank ID in fd->comm.
         */
        aggr = ADIOI_LUSTRE_Calc_aggregator(fd, off, &avail_len);
        aggr_ranks[i] = aggr;          /* first aggregator ID of this request */
        avail_lens[i] = avail_len;     /* length covered, may be < flat_fview.len[i] */
#ifdef WKL_DEBUG
assert(aggr >= 0 && aggr <= cb_nodes);
#endif
        my_req[aggr].count++; /* increment for aggregator aggr */
        nelems++;             /* true number of noncontiguous requests
                               * in terms of file domains */

        /* rem_len is the amount of ith offset-length pair that is not covered
         * by aggregator aggr's file domain.
         */
        rem_len = flat_fview.len[i] - avail_len;

#ifdef WKL_DEBUG
assert(rem_len >= 0);
#endif
        while (rem_len > 0) {
            off += avail_len;    /* move forward to first remaining byte */
            avail_len = rem_len; /* save remaining size, pass to calc */
            aggr = ADIOI_LUSTRE_Calc_aggregator(fd, off, &avail_len);
            my_req[aggr].count++;
            nelems++;
            rem_len -= avail_len;/* reduce remaining length by amount from fd */
        }
    }

    /* allocate space for buf_idx.
     * buf_idx is relevant only if buftype is contiguous. buf_idx[i] gives the
     * starting index in user_buf where data will be sent to aggregator 'i'.
     * This allows sends to be done without extra buffer.
     */
    if (buf_idx != NULL && buf_is_contig) {
        buf_idx[0] = (ADIO_Offset *) ADIOI_Malloc(nelems * sizeof(ADIO_Offset));
        for (i = 1; i < cb_nodes; i++)
            buf_idx[i] = buf_idx[i - 1] + my_req[i - 1].count;

#ifdef WKL_DEBUG
int wkl=0; for (i=0; i<cb_nodes; i++) wkl+=my_req[i].count;
assert(wkl == nelems);
for (i=0; i<nelems; i++) buf_idx[0][i] = -1;
#endif
    }

    /* allocate space for my_req and its members offsets and lens */
    my_req[0].offsets = (ADIO_Offset *) ADIOI_Malloc(nelems * 2 * sizeof(ADIO_Offset));

    ptr = my_req[0].offsets;
    for (i = 0; i < cb_nodes; i++) {
        if (my_req[i].count) {
            my_req[i].offsets = ptr;
            ptr += my_req[i].count;
            my_req[i].lens = ptr;
            ptr += my_req[i].count;
        }
        my_req[i].count = 0; /* reset, will be incremented where needed later */
    }

    /* now fill in my_req */
    curr_idx = 0;
    for (i = 0; i < flat_fview.count; i++) {
        /* short circuit offset/len processing if zero-byte read/write. */
        if (flat_fview.len[i] == 0)
            continue;

        off = flat_fview.off[i];
        aggr = aggr_ranks[i];
#ifdef WKL_DEBUG
assert(aggr >= 0 && aggr <= cb_nodes);
#endif
        avail_len = avail_lens[i];

        l = my_req[aggr].count;
        if (buf_idx != NULL && buf_is_contig) {
            buf_idx[aggr][l] = curr_idx;
            curr_idx += avail_len;
        }
        rem_len = flat_fview.len[i] - avail_len;

        /* Each my_req[i] contains the number of this process's noncontiguous
         * requests that fall into aggregator aggr's file domain.
         * my_req[aggr].offsets[] and my_req[aggr].lens store the offsets and
         * lengths of the requests.
         */
        my_req[aggr].offsets[l] = off;
        my_req[aggr].lens[l] = avail_len;
        my_req[aggr].count++;

        while (rem_len != 0) {
            off += avail_len;
            avail_len = rem_len;
            aggr = ADIOI_LUSTRE_Calc_aggregator(fd, off, &avail_len);

#ifdef WKL_DEBUG
assert(aggr >= 0 && aggr <= cb_nodes);
#endif
            l = my_req[aggr].count;
            if (buf_idx != NULL && buf_is_contig) {
                buf_idx[aggr][l] = curr_idx;
                curr_idx += avail_len;
            }
            rem_len -= avail_len;

            my_req[aggr].offsets[l] = off;
            my_req[aggr].lens[l] = avail_len;
            my_req[aggr].count++;
        }
    }
    ADIOI_Free(aggr_ranks);
    ADIOI_Free(avail_lens);

#ifdef AGG_DEBUG
    for (i = 0; i < cb_nodes; i++) {
        if (my_req[i].count > 0) {
            FPRINTF(stdout, "data needed from %d (count = %d):\n", i, my_req[i].count);
            for (l = 0; l < my_req[i].count; l++) {
                FPRINTF(stdout, "   off[%d] = %lld, len[%d] = %d\n",
                        l, (long long) my_req[i].offsets[l], l, (long long) my_req[i].lens[l]);
            }
        }
    }
#endif
}

/* ADIOI_LUSTRE_Calc_others_req() calculates what requests of other processes
 * lie in this aggregator's file domain.
 *   IN: my_req[cb_nodes]: offset-length pairs of this rank's requests fall into
 *       each aggregator
 *   OUT: count_others_req_per_proc[i]: number of noncontiguous requests of
 *        rank i that falls in this rank's file domain.
 *   OUT: others_req_ptr[nprocs]: requests of all other ranks fall into this
 *        aggregator's file domain.
 */
static
void ADIOI_LUSTRE_Calc_others_req(ADIO_File fd,
                                  const ADIOI_Access *my_req,
                                  ADIOI_Access **others_req_ptr)
{
    int i, myrank, nprocs;
    MPI_Count *count_my_req_per_proc, *count_others_req_per_proc;
    ADIOI_Access *others_req;
    size_t memLen;
    ADIO_Offset *ptr;
    MPI_Count *mem_ptrs;

/* first find out how much to send/recv and from/to whom */
#ifdef AGGREGATION_PROFILE
    MPE_Log_event(5026, 0, NULL);
#endif

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    /* Use my_req[].count to set count_others_req_per_proc[], a contiguous
     * buffer, so it can be used in th call to MPI_Alltoall().
     * my_req[i].count is the number of contiguous requests of this rank that
     * fall in aggregator i's file domain.
     */
    count_others_req_per_proc = (MPI_Count *) ADIOI_Malloc(nprocs * sizeof(MPI_Count));
    count_my_req_per_proc = (MPI_Count*) ADIOI_Calloc(nprocs, sizeof(MPI_Count));
    for (i=0; i<fd->hints->cb_nodes; i++)
        count_my_req_per_proc[fd->hints->ranklist[i]] = my_req[i].count;

    MPI_Alltoall(count_my_req_per_proc, 1, MPI_COUNT,
                 count_others_req_per_proc, 1, MPI_COUNT, fd->comm);
    ADIOI_Free(count_my_req_per_proc);

#ifdef WKL_DEBUG
if (!fd->is_agg) {
    /* non-aggregator should not receive any request */
    for (i=0; i<nprocs; i++) assert(count_others_req_per_proc[i] == 0);
}
#endif

    others_req = (ADIOI_Access *) ADIOI_Malloc(nprocs * sizeof(ADIOI_Access));
    *others_req_ptr = others_req;

    memLen = 0;
    for (i = 0; i < nprocs; i++)
        memLen += count_others_req_per_proc[i];
    ptr = (ADIO_Offset *) ADIOI_Malloc(memLen * 2 * sizeof(ADIO_Offset));
    mem_ptrs = (MPI_Count *) ADIOI_Malloc(memLen * sizeof(MPI_Count));
    others_req[0].offsets = ptr;
    others_req[0].mem_ptrs = mem_ptrs;

    for (i = 0; i < nprocs; i++) {
        if (count_others_req_per_proc[i]) {
            others_req[i].count = count_others_req_per_proc[i];
            others_req[i].offsets = ptr;
            ptr += count_others_req_per_proc[i];
            others_req[i].lens = ptr;
            ptr += count_others_req_per_proc[i];
            others_req[i].mem_ptrs = mem_ptrs;
            mem_ptrs += count_others_req_per_proc[i];
        } else
            others_req[i].count = 0;
    }
    ADIOI_Free(count_others_req_per_proc);

#ifdef WKL_DEBUG
if (!fd->is_agg) {
    /* non-aggregator should not receive any request */
    for (i=0; i<nprocs; i++) assert(others_req[i].count == 0);
}
MPI_Offset recv_amnt=0;
double timing, max_timing;
MPI_Barrier(fd->comm);
timing = MPI_Wtime();
#endif

    /* now send the calculated offsets and lengths to respective processes */

#define USE_ALLTOALLV_
#ifdef USE_ALLTOALLV_
    MPI_Offset *r_buf=NULL, *s_buf=NULL;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *sendCounts, *recvCounts;
    MPI_Aint *sdispls, *rdispls;
    sendCounts = (MPI_Count*) ADIOI_Calloc(nprocs * 2, sizeof(MPI_Count));
    sdispls = (MPI_Aint*) ADIOI_Calloc(nprocs * 2, sizeof(MPI_Aint));
#else
    int *sendCounts, *recvCounts, *sdispls, *rdispls;
    sendCounts = (int*) ADIOI_Calloc(nprocs * 2, sizeof(int));
    sdispls = (int*) ADIOI_Calloc(nprocs * 2, sizeof(int));
#endif
    recvCounts = sendCounts + nprocs;
    rdispls = sdispls + nprocs;

    /* prepare receive side */
    for (i=0; i<nprocs; i++) {
        if (others_req[i].count == 0) continue;
        if (i == myrank) {
            /* send to self uses memcpy(), here
             * others_req[i].count == my_req[fd->my_cb_nodes_index].count
             */
            memcpy(others_req[i].offsets, my_req[fd->my_cb_nodes_index].offsets,
                   2 * my_req[fd->my_cb_nodes_index].count * sizeof(ADIO_Offset));
        }
        else {
            recvCounts[i] = 2 * others_req[i].count;
            if (r_buf == NULL) r_buf = others_req[i].offsets;
            /* Note all others_req[*].offsets are allocated in a single malloc(). */
            rdispls[i] = others_req[i].offsets - r_buf;
        }
    }

    /* prepare send side */
    for (i = 0; i < fd->hints->cb_nodes; i++) {
        if (my_req[i].count && i != fd->my_cb_nodes_index) {
            int dest = fd->hints->ranklist[i];
            sendCounts[dest] = 2 * my_req[i].count;
            if (s_buf == NULL) s_buf = my_req[i].offsets;
            /* Note all my_req[*].offsets are allocated in a single malloc(). */
            sdispls[dest] = my_req[i].offsets - s_buf;
        }
    }

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Alltoallv_c(s_buf, sendCounts, sdispls, MPI_OFFSET,
                    r_buf, recvCounts, rdispls, MPI_OFFSET, fd->comm);
#else
    MPI_Alltoallv(s_buf, sendCounts, sdispls, MPI_OFFSET,
                  r_buf, recvCounts, rdispls, MPI_OFFSET, fd->comm);
#endif

    ADIOI_Free(sendCounts);
    ADIOI_Free(sdispls);

#else /* #ifdef USE_ALLTOALLV_ */
    int nreqs;
    MPI_Request *requests = (MPI_Request *)
        ADIOI_Malloc((nprocs + fd->hints->cb_nodes) * sizeof(MPI_Request));

    nreqs = 0;
    for (i = 0; i < nprocs; i++) {
        if (others_req[i].count == 0)
            continue;
#ifdef WKL_DEBUG
if (i == myrank) assert(fd->my_cb_nodes_index >= 0 && fd->my_cb_nodes_index <= fd->hints->cb_nodes);
else recv_amnt += 2 * others_req[i].count;
#endif
        if (i == myrank)
            /* send to self uses memcpy(), here
             * others_req[i].count == my_req[fd->my_cb_nodes_index].count
             */
            memcpy(others_req[i].offsets, my_req[fd->my_cb_nodes_index].offsets,
                   2 * my_req[fd->my_cb_nodes_index].count * sizeof(ADIO_Offset));
        else
            MPI_Irecv(others_req[i].offsets, 2 * others_req[i].count,
                      ADIO_OFFSET, i, 0, fd->comm, &requests[nreqs++]);
    }

#ifdef WKL_DEBUG
/* WRF hangs below when calling MPI_Waitall(), at running 16 nodes, 128 ranks
 * per node on Perlmutter, when these 3 env variables are set:
 *    FI_UNIVERSE_SIZE        = 2048
 *    FI_CXI_DEFAULT_CQ_SIZE  = 524288
 *    FI_CXI_RX_MATCH_MODE    = software
 *
 * Using MPI_Alltoallv seems to be able to avoid such hanging problem. (above)
 */
// MPI_Barrier(fd->comm); /* This barrier prevents the MPI_Waitall below from hanging !!! */
#endif

    for (i = 0; i < fd->hints->cb_nodes; i++) {
        if (my_req[i].count && i != fd->my_cb_nodes_index)
            MPI_Issend(my_req[i].offsets, 2 * my_req[i].count,
                      ADIO_OFFSET, fd->hints->ranklist[i], 0, fd->comm,
                      &requests[nreqs++]);
    }

    if (nreqs) {
#ifdef MPI_STATUSES_IGNORE
        MPI_Waitall(nreqs, requests, MPI_STATUSES_IGNORE);
#else
        MPI_Status *statuses = (MPI_Status *)
                               ADIOI_Malloc(nreqs * sizeof(MPI_Status));
        MPI_Waitall(nreqs, requests, statuses);
        ADIOI_Free(statuses);
#endif
    }
    ADIOI_Free(requests);
#endif

#ifdef AGGREGATION_PROFILE
    MPE_Log_event(5027, 0, NULL);
#endif

#ifdef WKL_DEBUG
timing = MPI_Wtime() - timing;
MPI_Offset max_recv_amnt=0;
MPI_Reduce(&recv_amnt, &max_recv_amnt, 1, MPI_OFFSET, MPI_MAX, 0, fd->comm);
MPI_Reduce(&timing, &max_timing, 1, MPI_DOUBLE, MPI_MAX, 0, fd->comm);
if (myrank == 0) printf("%s --- max_recv_amnt=%lld time=%.4f\n",__func__,max_recv_amnt,max_timing);
#endif
}

void ADIOI_LUSTRE_WriteStridedColl(ADIO_File fd, const void *buf, MPI_Aint count,
                                   MPI_Datatype buftype,
                                   int file_ptr_type, ADIO_Offset offset,
                                   ADIO_Status *status, int *error_code)
{
    /* Uses a generalized version of the extended two-phase method described in
     * "An Extended Two-Phase Method for Accessing Sections of Out-of-Core
     * Arrays", Rajeev Thakur and Alok Choudhary, Scientific Programming,
     * (5)4:301--317, Winter 1996.
     * http://www.mcs.anl.gov/home/thakur/ext2ph.ps
     */

    int i, j, nprocs, nonzero_nprocs, myrank, old_error, tmp_error;
    int do_collect = 0;
    ADIO_Offset orig_fp, start_offset, end_offset;
    ADIO_Offset min_st_loc = -1, max_end_loc = -1;
    MPI_Aint lb;
    ADIOI_Flatlist_node *flat_view = NULL;
    Flat_list flat_fview;
    Flat_list flat_bview;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    orig_fp = fd->fp_ind;

fd->lustre_write_metrics[1] = fd->lustre_write_metrics[2] = 0;
fd->lustre_write_metrics[0] = MPI_Wtime();

    /* Using user buffer datatype, count, and fileview to construct a list of
     * starting file offsets and write lengths of this rank and store them in
     * flat_fview.off[] and flat_fview.len[], respectively. Note a rank's
     * fileview is changed only when MPI_File_set_view is called and thus The
     * same fileview can be used by multiple collective writes. On the other
     * hand, user buffer datatype and count can be different at each collective
     * write call. Thus constructing a new flat_fview is necessary at each
     * collective write.
     *
     * Note flat_fview here is NOT about the fileview datatype itself. It
     * actually contains this rank's file access information of this collective
     * write call. flat_fview.count is the number of noncontiguous file
     * offset-length pairs, thus the size of both flat_fview.off[] and
     * flat_fview.len[]. flat_fview.count has taken into account of argument
     * 'count', i.e. the number of user buffer datatypes in this request.
     *
     * TODO: In the current implementation, even for a small fileview type, the
     *       flat_fview.count can still be large, when the write amount is
     *       larger than the file type size. In order to reduce the memory
     *       footprint, flat_fview should be modified to describe only one file
     *       type and use flat_fview.rnd, flat_fview.idx, flat_fview.rem to
     *       keep track the latest processed offset-length pairs, just like the
     *       way flat_bview is used.
     *
     * This rank's aggregate file access region is from start_offset to
     * end_offset. Note: end_offset points to the last byte-offset to be
     * accessed. E.g., if start_offset=0 and end_offset=99, then the aggregate
     * file access region is of size 100 bytes. If this rank has no data to
     * write, end_offset == (start_offset - 1)
     *
     * ADIOI_Calc_my_off_len() requires no inter-process communication.
     */
    ADIOI_Calc_my_off_len(fd, count, buftype, file_ptr_type, offset,
                          &flat_fview.off, &flat_fview.len,
                          &start_offset, &end_offset, &flat_fview.count);
    flat_fview.idx = 0;
    flat_fview.rnd = 0; /* currently for flat_fview, rnd is not used at all */
    flat_fview.rem = (flat_fview.count > 0) ? flat_fview.len[0] : 0;

    /* When this rank's fileview datatype is not contiguous, it is still
     * possible that flat_fview.count == 1. It happens when this write is small
     * enough to fall into one of the contiguous segment of the fileview.
     * Thus, flat_fview.count being 1 or not is the true indicator of whether
     * or not this write is contiguous in the file. This is because
     * ADIOI_Calc_my_off_len has taken into account of both fileview and user
     * buffer type.
     *
     * Note in ADIOI_Calc_my_off_len(), checking whether a fileview is
     * contiguous by calling
     *     ADIOI_Datatype_iscontig(fd->filetype, &filetype_is_contig);
     * But whether a filetype is contiguous is not necessary requal to whether
     * this collective write is contiguous in file.
     *
     * Note flat_fview.is_contig is whether or not the file access region of
     * this collective write call is contiguous in file. It is not whether the
     * filetype is contiguous or not. See test program noncontig_filetype.c
     */
    flat_fview.is_contig = (flat_fview.count > 1) ? 0 : 1;

    /* flatten user buffer datatype, buftype */
    flat_view = ADIOI_Flatten_and_find(buftype);

    flat_bview.type  = buftype;
    flat_bview.count = flat_view->count;
    flat_bview.off   = flat_view->indices;
    flat_bview.len   = flat_view->blocklens;
    flat_bview.rnd   = 0;
    flat_bview.idx   = 0;
    flat_bview.rem   = (flat_bview.count > 0) ? flat_bview.len[0] : 0;
    MPI_Type_get_extent(buftype, &lb, &flat_bview.extent);

    /* Check if the user buffer is truly contiguous or not.
     * flat_bview.count being 1 is not the true indicator of whether or not
     * buftype is contiguous. For example, a buftype may contain only one
     * offset-length pair, but its ub and lb may have been resized to make
     * buftype noncontiguous. In this case, if the value of argument 'count' in
     * this collective write call is larger than 1, then the user buffer is not
     * contiguous. On the other hand, when 'count' == 1, the user buffer is
     * contiguous.
     *
     * Note flat_bview.is_contig is whether or not the user buffer of this
     * collective write call is contiguous in memory or not. It is not whether
     * or not the buftype is contiguous. See test program resize_buftype.c
     */
    flat_bview.is_contig = 0;
    if (flat_bview.count == 1) {
        if (count == 1) /* write amount is less than one buftype */
            flat_bview.is_contig = 1;
        else if (flat_bview.extent == flat_bview.len[0])
            /* buftype extent is the same as the only length */
            flat_bview.is_contig = 1;
    }

    /* Check if collective write is actually necessary, if cb_write hint isn't
     * disabled by users.
     */
    if (fd->hints->cb_write != ADIOI_HINT_DISABLE) {
        int is_interleaved;
        ADIO_Offset st_end[2], *st_end_all = NULL;

        /* Gather starting and ending file offsets of requests from all ranks
         * into st_end_all[]. Even indices of st_end_all[] are start offsets,
         * odd indices are end offsets. st_end_all[] is used below to tell
         * whether access across all ranks is interleaved.
         */
        st_end[0] = start_offset;
        st_end[1] = end_offset;
        st_end_all = (ADIO_Offset *) ADIOI_Malloc(nprocs * 2 * sizeof(ADIO_Offset));
        MPI_Allgather(st_end, 2, ADIO_OFFSET, st_end_all, 2, ADIO_OFFSET, fd->comm);

        /* check if the request pattern is non-interleaved among all processes
         * and each process writes a large amount. Here, "large" means a
         * process's write range is > striping_factor * striping_unit. In this
         * case, independent write will perform faster than collective.
         */
        int large_indv_req = 1;
        MPI_Offset striping_range = fd->hints->striping_unit * fd->hints->striping_factor;

        /* Find the starting and ending file offsets of aggregate access region
         * and the number of ranks that have non-zero length write requests.
         * Also, check whether accesses are interleaved across ranks. Below is
         * a rudimentary check for interleaving, but should suffice for the
         * moment.
         */
        is_interleaved = 0;
        nonzero_nprocs = 0;
        for (i = 0; i < nprocs * 2; i += 2) {
            if (st_end_all[i] > st_end_all[i + 1]) {
                /* process rank (i/2) has no data to write */
                continue;
            }
            min_st_loc = st_end_all[i];
            max_end_loc = st_end_all[i + 1];
            if (st_end_all[i+1] - st_end_all[i] < striping_range)
                large_indv_req = 0;
            nonzero_nprocs = 1;
            j = i; /* j is the rank of making first non-zero request */
            i += 2;
            break;
        }
        for (; i < nprocs * 2; i += 2) {
            if (st_end_all[i] > st_end_all[i + 1]) {
                /* process rank (i/2) has no data to write */
                continue;
            }
            if (st_end_all[i] < st_end_all[j+1]) {
                /* start offset of process rank (i/2) is less than the end
                 * offset of process rank (i/2-1)
                 */
                is_interleaved = 1;
            }
            min_st_loc = MPL_MIN(st_end_all[i], min_st_loc);
            max_end_loc = MPL_MAX(st_end_all[i + 1], max_end_loc);
            nonzero_nprocs++;
            if (st_end_all[i+1] - st_end_all[i] < striping_range)
                large_indv_req = 0;
            j = i;
        }
        ADIOI_Free(st_end_all);

        /* Two typical access patterns can benefit from collective write.
         *   1) access file regions among all processes are interleaved, and
         *   2) the individual request sizes are not too big, i.e. no bigger
         *      than hint coll_threshold.  Large individual requests may cause
         *      a high communication cost for redistributing requests to the
         *      I/O aggregators.
         */
        if (is_interleaved > 0)
            do_collect = 1;
        else if (nprocs == 1)
            do_collect = 0;
        else if (large_indv_req &&
                 fd->hints->cb_nodes <= fd->hints->striping_factor)
            /* do independent write, if every rank's write range >
             * striping_range and writes are not interleaved in file space */
            do_collect = 0;
    }

    /* If collective I/O is not necessary, use independent I/O */
    if ((!do_collect && fd->hints->cb_write == ADIOI_HINT_AUTO) ||
        fd->hints->cb_write == ADIOI_HINT_DISABLE) {

#ifdef WKL_DEBUG
if (do_collect == 0) printf("%s --- SWITCH to independent write !!!\n",__func__);
#endif

        fd->fp_ind = orig_fp;

        if (flat_fview.is_contig && flat_bview.is_contig) {
            /* both buffer and fileview are contiguous */
            ADIO_Offset off = 0;
            if (file_ptr_type == ADIO_EXPLICIT_OFFSET)
                off = flat_fview.off[0];
                /* In ADIOI_Calc_my_off_len(), (offset * fd->etype_size) has
                 * been counted into flat_fview.off[]. Similarly, fd->disp has
                 * been counted into flat_fview.off[].
                 */

            if (flat_fview.off != NULL) ADIOI_Free(flat_fview.off);

static int wkl=0; if (wkl==0 && fd->disp>0) { printf("%s %d: --- SWITCH to independent write ADIO_WriteContig fd->disp=%lld off=%lld !!!\n",__func__,__LINE__,fd->disp,off); wkl++; }

            ADIO_WriteContig(fd, buf, count, buftype, file_ptr_type, off, status, error_code);
        } else {
            if (flat_fview.off != NULL) ADIOI_Free(flat_fview.off);

            ADIO_WriteStrided(fd, buf, count, buftype, file_ptr_type, offset, status, error_code);
        }

        return;
    }

    /* Now we are using collective I/O (two-phase I/O strategy) */

    /* adjust striping_unit when striping_factor is twice or more than the
     * number of compute nodes. Note cb_node is set to at least
     * striping_factor, if nprocs >= striping_factor. Adjustment below is to
     * let each aggregator to write to two or more consecutive OSTs, which can
     * most likely improve the performance. This will still yield an effect of
     * any one OST receiving write requests from aggregators running on only
     * one compute node.
     */
    int orig_striping_unit = fd->hints->striping_unit;
    if (fd->hints->striping_factor >= fd->num_nodes * 2) {
        fd->hints->striping_unit *= (fd->hints->striping_factor / fd->num_nodes);
#ifdef WKL_DEBUG
        if (myrank == 0)
            printf("Warning: %s line %d: Change striping_unit from %d to %d\n",
                   __func__, __LINE__, orig_striping_unit, fd->hints->striping_unit);
#endif
    }

    /* my_req[cb_nodes] is an array of access info, one for each I/O
     * aggregator whose file domain has this rank's request.
     */
    ADIOI_Access *my_req;

    /* others_req[nprocs] is an array of access info, one for each process
     * whose write requests fall into this process's file domain.
     */
    ADIOI_Access *others_req;
    ADIO_Offset **buf_idx = NULL;

    /* Calculate the portions of this process's write requests that fall
     * into the file domains of each I/O aggregator. No inter-process
     * communication is needed.
     */
    if (flat_bview.is_contig)
        buf_idx = (ADIO_Offset **) ADIOI_Malloc(fd->hints->cb_nodes * sizeof(ADIO_Offset*));

    ADIOI_LUSTRE_Calc_my_req(fd, flat_fview, flat_bview.is_contig,
                             &my_req, buf_idx);

    /* Calculate the portions of all other ranks' requests fall into this
     * process's file domain (note only I/O aggregators are assigned file
     * domains). Inter-process communication is required to construct
     * others_req[], including MPI_Alltoall, MPI_Issend, MPI_Irecv, and
     * MPI_Waitall.
     */
    ADIOI_LUSTRE_Calc_others_req(fd, my_req, &others_req);

    /* Two-phase I/O: first communication phase to exchange write data from all
     * processes to the I/O aggregators, followed by the write phase where only
     * I/O aggregators write to the file. There is no collective MPI
     * communication in ADIOI_LUSTRE_Exch_and_write(), only MPI_Issend,
     * MPI_Irecv, and MPI_Waitall.
     */
    ADIOI_LUSTRE_Exch_and_write(fd, buf, &flat_bview, others_req,
                                my_req, &flat_fview, min_st_loc,
                                max_end_loc, buf_idx, error_code);

    /* free all memory allocated */
    ADIOI_Free(others_req[0].offsets);
    ADIOI_Free(others_req[0].mem_ptrs);
    ADIOI_Free(others_req);

    if (buf_idx != NULL) {
        ADIOI_Free(buf_idx[0]);
        ADIOI_Free(buf_idx);
    }
    ADIOI_Free(my_req[0].offsets);
    ADIOI_Free(my_req);

    /* restore the original striping_unit */
    fd->hints->striping_unit = orig_striping_unit;

    if (flat_fview.off != NULL)
        ADIOI_Free(flat_fview.off);

    /* If this collective write is followed by an independent write, it's
     * possible to have those subsequent writes on other processes race ahead
     * and sneak in before the read-modify-write completes.  We carry out a
     * collective communication at the end here so no one can start independent
     * I/O before collective I/O completes.
     *
     * need to do some gymnastics with the error codes so that if something
     * went wrong, all processes report error, but if a process has a more
     * specific error code, we can still have that process report the
     * additional information
     */
    old_error = *error_code;
    if (*error_code != MPI_SUCCESS)
        *error_code = MPI_ERR_IO;

    /* optimization: if only one process performing I/O, we can perform
     * a less-expensive Bcast. */
#ifdef ADIOI_MPE_LOGGING
    MPE_Log_event(ADIOI_MPE_postwrite_a, 0, NULL);
#endif
    if (fd->hints->cb_nodes == 1)
        MPI_Bcast(error_code, 1, MPI_INT, fd->hints->ranklist[0], fd->comm);
    else {
        tmp_error = *error_code;
        MPI_Allreduce(&tmp_error, error_code, 1, MPI_INT, MPI_MAX, fd->comm);
    }
#ifdef ADIOI_MPE_LOGGING
    MPE_Log_event(ADIOI_MPE_postwrite_b, 0, NULL);
#endif

    if ((old_error != MPI_SUCCESS) && (old_error != MPI_ERR_IO))
        *error_code = old_error;

    if (status) {
        /* This is a temporary way of filling in status. The right way is to
         * keep track of how much data was actually written during collective
         * I/O.
         */
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
        MPI_Status_set_elements_x(status, buftype, count);
#else
        MPI_Status_set_elements(status, buftype, count);
#endif
    }

    fd->fp_sys_posn = -1;       /* set it to null. */

#ifdef WKL_DEBUG
fd->lustre_write_metrics[0] = MPI_Wtime() - fd->lustre_write_metrics[0];
double max_t[3]; MPI_Reduce(fd->lustre_write_metrics, max_t, 3, MPI_DOUBLE, MPI_MAX, 0, fd->comm);
if (myrank == 0) printf("%s line %d: MAX lustre time write=%.4f pwrite=%.4f all-to-many senders=%ld (nprocs=%d)\n", __func__, __LINE__, max_t[0], max_t[1], (long)max_t[2], nprocs);
#endif
}

#define CAST_INT32(count, bklen, disp, dType, newType) {                     \
    int kk, iCount, *iBklen;                                                 \
    MPI_Aint *iDisp;                                                         \
    ADIOI_Assert(count <= 2147483647); /* overflow 4-byte int */             \
    iCount = (int)count;                                                     \
    iBklen = (int*) ADIOI_Malloc(sizeof(int) * iCount);                      \
    for (kk=0; kk<iCount; kk++) {                                            \
        ADIOI_Assert(bklen[kk] <= 2147483647); /* overflow 4-byte int */     \
        iBklen[kk] = (int)bklen[kk];                                         \
    }                                                                        \
    if (sizeof(MPI_Aint) < sizeof(MPI_Count)) {                              \
        iDisp = (MPI_Aint*) ADIOI_Malloc(sizeof(MPI_Aint) * iCount);         \
        for (kk=0; kk<iCount; kk++) {                                        \
            ADIOI_Assert(disp[kk] <= 2147483647); /* overflow 4-byte int */  \
            iDisp[kk] = (MPI_Aint)disp[kk];                                  \
        }                                                                    \
    }                                                                        \
    else                                                                     \
        iDisp = (MPI_Aint*)disp;                                             \
    MPI_Type_create_hindexed(iCount, iBklen, iDisp, dType, newType);         \
    ADIOI_Free(iBklen);                                                      \
    if (sizeof(MPI_Aint) != sizeof(MPI_Count))                               \
        ADIOI_Free(iDisp);                                                   \
}

static
void commit_comm_phase(ADIO_File      fd,
                       disp_len_list *send_list,  /* [cb_nodes] */
                       disp_len_list *recv_list)  /* [nprocs] */
{
    /* This subroutine creates a datatype combining all displacement-length
     * pairs in each element of send_list[]. The datatype is used when calling
     * MPI_Issend to send write data to the I/O aggregators. Similarly, it
     * creates a datatype combining all displacement-length pairs in each
     * element of recv_list[] and uses it when calling MPI_Irecv or MPI_Recv
     * to receive write data from all processes.
     */
    int i, nprocs;

    MPI_Comm_size(fd->comm, &nprocs);

#ifdef USE_ALLTOALLW_
    MPI_Datatype *sendTypes, *recvTypes;
    sendTypes = (MPI_Datatype*)ADIOI_Malloc(sizeof(MPI_Datatype) * nprocs * 2);
    recvTypes = sendTypes + nprocs;

    for (i=0; i<nprocs; i++)
        sendTypes[i] = recvTypes[i] = MPI_BYTE;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count *sendCounts, *recvCounts;
    MPI_Aint *sdispls, *rdispls;
    sendCounts = (MPI_Count*) ADIOI_Calloc(nprocs * 2, sizeof(MPI_Count));
    sdispls = (MPI_Aint*) ADIOI_Calloc(nprocs * 2, sizeof(MPI_Aint));
#else
    int *sendCounts, *recvCounts, *sdispls, *rdispls;
    sendCounts = (int*) ADIOI_Calloc(nprocs * 2, sizeof(int));
    sdispls = (int*) ADIOI_Calloc(nprocs * 2, sizeof(int));
#endif
    recvCounts = sendCounts + nprocs;
    rdispls = sdispls + nprocs;

    /* prepare receive side */
    if (fd->is_agg) {
        for (i=0; i<nprocs; i++) {
            if (recv_list[i].count == 0) continue;

            recvCounts[i] = 1;

            /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(recv_list[i].count, recv_list[i].len,
                                       recv_list[i].disp, MPI_BYTE,
                                       &recvTypes[i]);
#else
            CAST_INT32(recv_list[i].count, recv_list[i].len,
                                       recv_list[i].disp, MPI_BYTE,
                                       &recvTypes[i]);
#endif
            MPI_Type_commit(&recvTypes[i]);
        }
    }

    /* prepare send side */
    for (i=0; i<fd->hints->cb_nodes; i++) {
        if (send_list[i].count == 0) continue;

        int dest = fd->hints->ranklist[i];
        sendCounts[dest] = 1;

        /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Type_create_hindexed_c(send_list[i].count, send_list[i].len,
                                   send_list[i].disp, MPI_BYTE,
                                   &sendTypes[dest]);
#else
        CAST_INT32(send_list[i].count, send_list[i].len,
                                   send_list[i].disp, MPI_BYTE,
                                   &sendTypes[dest]);
#endif
        MPI_Type_commit(&sendTypes[dest]);
    }

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Alltoallw_c(MPI_BOTTOM, sendCounts, sdispls, sendTypes,
                    MPI_BOTTOM, recvCounts, rdispls, recvTypes, fd->comm);
#else
    MPI_Alltoallw(MPI_BOTTOM, sendCounts, sdispls, sendTypes,
                  MPI_BOTTOM, recvCounts, rdispls, recvTypes, fd->comm);
#endif

    for (i=0; i<nprocs; i++) {
        if (sendTypes[i] != MPI_BYTE)
            MPI_Type_free(&sendTypes[i]);
        if (recvTypes[i] != MPI_BYTE)
            MPI_Type_free(&recvTypes[i]);
    }
    ADIOI_Free(sendCounts);
    ADIOI_Free(sdispls);
    ADIOI_Free(sendTypes);

#else /* #ifdef USE_ALLTOALLW_ */
    int nreqs;
    MPI_Request *reqs;
    MPI_Status status;
    MPI_Datatype sendType, recvType;

    nreqs = fd->hints->cb_nodes;
    nreqs += (fd->is_agg) ? nprocs : 0;
    reqs = (MPI_Request *)ADIOI_Malloc(sizeof(MPI_Request) * nreqs);
    nreqs = 0;

    /* receiving part */
int nrecvs = 0;
    if (fd->is_agg) {
        for (i = 0; i < nprocs; i++) {
            if (recv_list[i].count > 0) {
                /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Type_create_hindexed_c(recv_list[i].count, recv_list[i].len,
                                           recv_list[i].disp, MPI_BYTE,
                                           &recvType);
#else
                CAST_INT32(recv_list[i].count, recv_list[i].len,
                                           recv_list[i].disp, MPI_BYTE,
                                           &recvType);
#endif
                MPI_Type_commit(&recvType);

                if (fd->atomicity) /* Blocking Recv */
                    MPI_Recv(MPI_BOTTOM, 1, recvType, i, 0, fd->comm, &status);
                else
                    MPI_Irecv(MPI_BOTTOM, 1, recvType, i, 0, fd->comm,
                              &reqs[nreqs++]);
                MPI_Type_free(&recvType);
nrecvs++;
            }
        }
    }

    /* send reqs */
    for (i = 0; i < fd->hints->cb_nodes; i++) {
        if (send_list[i].count > 0) {
            /* combine reqs using new datatype */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Type_create_hindexed_c(send_list[i].count, send_list[i].len,
                                       send_list[i].disp, MPI_BYTE, &sendType);
#else
            CAST_INT32(send_list[i].count, send_list[i].len,
                                       send_list[i].disp, MPI_BYTE, &sendType);
#endif
            MPI_Type_commit(&sendType);

            MPI_Issend(MPI_BOTTOM, 1, sendType, fd->hints->ranklist[i], 0,
                       fd->comm, &reqs[nreqs++]);
            MPI_Type_free(&sendType);
        }
    }
fd->lustre_write_metrics[2] = (fd->lustre_write_metrics[2] > nrecvs) ? fd->lustre_write_metrics[2] : nrecvs;

    if (nreqs > 0)
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);

    ADIOI_Free(reqs);
#endif

    /* clear send_list and recv_list for future reuse */
    for (i = 0; i < fd->hints->cb_nodes; i++)
        send_list[i].count = 0;

    if (fd->is_agg)
        for (i = 0; i < nprocs; i++)
            recv_list[i].count = 0;
}

/* If successful, error_code is set to MPI_SUCCESS.  Otherwise an error code is
 * created and returned in error_code.
 */
static void ADIOI_LUSTRE_Exch_and_write(ADIO_File fd,
                                        const void *buf,
                                        Flat_list *flat_bview,
                                        ADIOI_Access *others_req,
                                        ADIOI_Access *my_req,
                                        Flat_list *flat_fview,
                                        ADIO_Offset min_st_loc,
                                        ADIO_Offset max_end_loc,
                                        ADIO_Offset **buf_idx,
                                        int *error_code)
{
    /* Each process sends all its write requests to I/O aggregators based on
     * the file domain assignment to the aggregators. In this implementation,
     * a file is first divided into stripes and all its stripes are assigned to
     * I/O aggregators in a round-robin fashion. The collective write is
     * carried out in 'ntimes' rounds of two-phase I/O. Each round covers an
     * aggregate file region of size equal to the file stripe size times the
     * number of I/O aggregators. In other words, the 'collective buffer size'
     * used in each aggregator is always set equally to the file stripe size,
     * ignoring the MPI-IO hint 'cb_buffer_size'. There are other algorithms
     * allowing an aggregator to write more than a file stripe size in each
     * round, up to the cb_buffer_size hint. For those, refer to the paper:
     * Wei-keng Liao, and Alok Choudhary. "Dynamically Adapting File Domain
     * Partitioning Methods for Collective I/O Based on Underlying Parallel
     * File System Locking Protocols", in The Supercomputing Conference, 2008.
     */

    char **write_buf = NULL, **recv_buf = NULL, **send_buf = NULL;
    int i, nprocs, myrank, nbufs;
    MPI_Count j, m, ntimes;
    MPI_Count *recv_curr_offlen_ptr, **recv_size=NULL, **recv_count=NULL;
    MPI_Count **recv_start_pos=NULL;
    MPI_Count *send_curr_offlen_ptr, *send_size;
    int batch_idx = 0, cb_nodes, striping_unit;
    ADIO_Offset end_loc, req_off, iter_end_off, *off_list, step_size;
    ADIO_Offset *this_buf_idx;
    off_len_list *srt_off_len = NULL;
    disp_len_list *send_list = NULL, *recv_list = NULL;

#ifdef WKL_DEBUG
double timing[6]={0,0,0,0,0,0}, s_time, e_time;
MPI_Barrier(fd->comm);
timing[0] = s_time = MPI_Wtime();
#endif

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    /* The aggregate access region (across all processes) of this collective
     * write starts from min_st_loc and ends at max_end_loc. The collective
     * write is carried out in 'ntimes' rounds of two-phase I/O. Each round
     * covers an aggregate file region of size 'step_size' written only by
     * cb_nodes number of I/O aggregators. Note
     * non-aggregators must also participate all ntimes rounds to send their
     * requests to I/O aggregators.
     */

    cb_nodes = fd->hints->cb_nodes;
    striping_unit = fd->hints->striping_unit;

    /* step_size is the size of aggregate access region covered by each round
     * of two-phase I/O
     */
    step_size = (ADIO_Offset)cb_nodes * striping_unit;

    /* align min_st_loc downward to the nearest file stripe boundary */
    min_st_loc -= min_st_loc % (ADIO_Offset) striping_unit;

    /* ntimes is the number of rounds of two-phase I/O */
    ntimes = (max_end_loc - min_st_loc + 1) / step_size;
    if ((max_end_loc - min_st_loc + 1) % step_size)
        ntimes++;

    /* off_list[m] is the starting file offset of this process's write region
     * in iteration m (file domain of iteration m). This offset may not be
     * aligned with file stripe boundaries. end_loc is the ending file offset
     * of this process's file domain.
     */
    off_list = (ADIO_Offset *) ADIOI_Malloc(ntimes * sizeof(ADIO_Offset));
    end_loc = -1;
    for (m = 0; m < ntimes; m++)
        off_list[m] = max_end_loc;
    for (i = 0; i < nprocs; i++) {
        for (j = 0; j < others_req[i].count; j++) {
            req_off = others_req[i].offsets[j];
            m = (int) ((req_off - min_st_loc) / step_size);
            off_list[m] = MPL_MIN(off_list[m], req_off);
            end_loc = MPL_MAX(end_loc, (others_req[i].offsets[j] + others_req[i].lens[j] - 1));
        }
    }

    /* collective buffer size is divided by 2, one for write buffer and the
     * other for receive messages from other processes. nbufs is the number of
     * write buffers, each of size equal to Lustre stripe size. Write requests
     * are sent to aggregators for nbufs amount before writing to the file.
     */
    nbufs = fd->hints->cb_buffer_size / striping_unit / 2;
    if (fd->hints->cb_buffer_size % striping_unit) nbufs++;
    nbufs = (nbufs > ntimes) ? ntimes : nbufs;
    if (nbufs == 0) nbufs = 1; /* must at least 1 */

    /* Allocate displacement-length pair arrays, describing the send buffer.
     * send_list[i].count: number displacement-length pairs.
     * send_list[i].len: length in bytes.
     * send_list[i].disp: displacement (send buffer address).
     */
    send_list = (disp_len_list*) ADIOI_Malloc(sizeof(disp_len_list) * cb_nodes);
    for (i = 0; i < cb_nodes; i++) {
        send_list[i].count = 0;
        send_list[i].len   = (MPI_Count*) ADIOI_Malloc(sizeof(MPI_Count) * nbufs);
        send_list[i].disp  = (MPI_Count*) ADIOI_Malloc(sizeof(MPI_Count) * nbufs);
    }

    /* Allocate displacement-length pair arrays, describing the recv buffer.
     * recv_list[i].count: number displacement-length pairs.
     * recv_list[i].len: length in bytes.
     * recv_list[i].disp: displacement (recv buffer address).
     */
    if (fd->is_agg) {
        recv_list = (disp_len_list*) ADIOI_Malloc(sizeof(disp_len_list) * nprocs);
        for (i = 0; i < nprocs; i++) {
            recv_list[i].count = 0;
            recv_list[i].len   = (MPI_Count*) ADIOI_Malloc(sizeof(MPI_Count) * nbufs);
            recv_list[i].disp  = (MPI_Count*) ADIOI_Malloc(sizeof(MPI_Count) * nbufs);
        }
    }

    /* end_loc >= 0 indicates this process has something to write. Only I/O
     * aggregators can have end_loc > 0. write_buf is the collective buffer and
     * only matter for I/O aggregators. recv_buf is the buffer used only in
     * aggregators to receive requests from other processes. Its size may be
     * larger then the file stripe size. In this case, it will be realloc-ed in
     * ADIOI_LUSTRE_W_Exchange_data(). The received data is later copied over
     * to write_buf, whose contents will be written to file.
     */
    if (end_loc >= 0) {
        write_buf = (char **) ADIOI_Malloc(nbufs * sizeof(char*));

        /* collective buffer was allocated in ADIO_Open(). For this Lustre ADIO
         * driver, its size must be at least (nbufs * striping_unit)
         */
        if (fd->hints->cb_buffer_size < nbufs * striping_unit)
            fd->io_buf = (char*) ADIOI_Realloc(fd->io_buf, nbufs * striping_unit);

        write_buf[0] = fd->io_buf;
        for (j = 1; j < nbufs; j++)
            write_buf[j] = write_buf[j-1] + striping_unit;

        recv_buf = (char **) ADIOI_Malloc(nbufs * sizeof(char*));
        for (j = 0; j < nbufs; j++)
            /* recv_buf[j] may be realloc in ADIOI_LUSTRE_W_Exchange_data() */
            recv_buf[j] = (char *) ADIOI_Malloc(striping_unit);
    }

    /* send_buf[] will be allocated in ADIOI_LUSTRE_W_Exchange_data(), when the
     * use buffer is not contiguous.
     */
    send_buf = (char **) ADIOI_Malloc(nbufs * sizeof(char*));

    /* this_buf_idx contains indices to user buffer for sending this rank's
     * write data to remote aggregators. It is used only when user buffer is
     * contiguous.
     */
    if (flat_bview->is_contig)
        this_buf_idx = (ADIO_Offset *) ADIOI_Malloc(cb_nodes * sizeof(ADIO_Offset));

    /* Allocate multiple buffers of type int altogether at once in a single
     * calloc call. Their use is explained below. calloc initializes to 0.
     */
    recv_curr_offlen_ptr = (MPI_Count *) ADIOI_Calloc((nprocs + 2 * cb_nodes), sizeof(MPI_Count));
    send_curr_offlen_ptr = recv_curr_offlen_ptr + nprocs;

    /* array of data sizes to be sent to each aggregator in a 2-phase round */
    send_size = send_curr_offlen_ptr + cb_nodes;

    /* I need to check if there are any outstanding nonblocking writes to
     * the file, which could potentially interfere with the writes taking
     * place in this collective write call. Since this is not likely to be
     * common, let me do the simplest thing possible here: Each process
     * completes all pending nonblocking operations before completing.
     */
    /*ADIOI_Complete_async(error_code);
     * if (*error_code != MPI_SUCCESS) return;
     * MPI_Barrier(fd->comm);
     */

    /* min_st_loc has been downward aligned to the nearest file stripe
     * boundary, iter_end_off is the ending file offset of aggregate write
     * region of iteration m, upward aligned to the file stripe boundary.
     */
    iter_end_off = min_st_loc + step_size;

    /* the number of off-len pairs to be received from each proc in a round. */
    recv_count = (MPI_Count**) ADIOI_Malloc(3 * nbufs * sizeof(MPI_Count*));
    recv_count[0] = (MPI_Count*) ADIOI_Malloc(3 * nbufs * nprocs * sizeof(MPI_Count));
    for (i = 1; i < nbufs; i++)
        recv_count[i] = recv_count[i-1] + nprocs;

    /* recv_size is array of data sizes to be received from each proc in a
     * round. */
    recv_size = recv_count + nbufs;
    recv_size[0] = recv_count[0] + nbufs * nprocs;
    for (i = 1; i < nbufs; i++)
        recv_size[i] = recv_size[i-1] + nprocs;

    /* recv_start_pos[j][i] stores the starting value of
     * recv_curr_offlen_ptr[i] for remote rank i in round j
     */
    recv_start_pos = recv_size + nbufs;
    recv_start_pos[0] = recv_size[0] + nbufs * nprocs;
    for (i = 1; i < nbufs; i++)
        recv_start_pos[i] = recv_start_pos[i-1] + nprocs;

    srt_off_len = (off_len_list*) ADIOI_Malloc(nbufs * sizeof(off_len_list));

    int ibuf = 0;

#ifdef WKL_DEBUG
e_time = MPI_Wtime();
timing[3] += e_time - s_time;
#endif
    for (m = 0; m < ntimes; m++) {
        int range_size;
        ADIO_Offset range_off;

#ifdef WKL_DEBUG
s_time = MPI_Wtime();
#endif
        /* Note that MPI standard (MPI3.1 Chapter 13.1.1 and MPI 4.0 Chapter
         * 14.1.1) requires that the typemap displacements of etype and
         * filetype are non-negative and monotonically non-decreasing. This
         * simplifies implementation a bit compared to reads.
         *
         * iter_end_off  = ending file offset of aggregate write region of this
         *                 round, and upward aligned to the file stripe
         *                 boundary. Note the aggregate write region of this
         *                 round starts from (iter_end_off-step_size) to
         *                 iter_end_off, aligned with file stripe boundaries.
         * send_size[i]  = total size in bytes of this process's write data
         *                 fall into aggregator i's FD in this round.
         * recv_size[j][i] = total size in bytes of write data to be received
         *                 by this process (aggregator) in round j.
         * recv_count[j][i] = number of noncontiguous offset-length pairs from
         *                 process i fall into this aggregator's write region
         *                 in round j.
         */

        /* Calculate what should be communicated. First, calculate the amount
         * to be send to each aggregator for this round of m, by going through
         * all my_req.
         */
        for (i = 0; i < cb_nodes; i++) {
            /* reset communication metadata to all 0s for this round */
            send_size[i] = 0;

            if (my_req[i].count == 0) continue;
            /* my_req[i].count is the number of this rank's offset-length pairs
             * to be sent to aggregator i
             */

            if (send_curr_offlen_ptr[i] == my_req[i].count)
                continue; /* done with aggregator i */

            if (flat_bview->is_contig)
                /* buf_idx is used only when user buffer is contiguous.
                 * this_buf_idx[i] points to the starting offset of user
                 * buffer, buf, for amount of send_size[i] to be sent to
                 * aggregator i in this round of m.
                 */
                this_buf_idx[i] = buf_idx[i][send_curr_offlen_ptr[i]];

            /* calculate the send amount from this rank to aggregator i */
            for (j = send_curr_offlen_ptr[i]; j < my_req[i].count; j++) {
                if (my_req[i].offsets[j] < iter_end_off)
                    send_size[i] += my_req[i].lens[j];
                else
                    break;
            }
            send_curr_offlen_ptr[i] = j;
            /* point to the jth offset-length pair of my_req[i], which will
             * be used as the first pair in the next round of iteration.
             */
        }

        /* range_off is the starting file offset of this aggregator's write
         * region for this round (may not be aligned to stripe boundary).
         * range_size is the size (in bytes) of this aggregator's write region
         * for this round (whose size is always <= striping_unit).
         */
        range_off = off_list[m];
        range_size = (int) MPL_MIN(striping_unit - range_off % striping_unit,
                                   end_loc - range_off + 1);

        /* Calculate the amount to be received from each process for this round
         * of m, by going through others_req.
         */
        for (i = 0; i < nprocs; i++) {
            /* reset communication metadata to all 0s for this round */
            recv_count[ibuf][i] = recv_size[ibuf][i] = recv_start_pos[ibuf][i] = 0;

            if (others_req[i].count == 0) continue;

            recv_start_pos[ibuf][i] = recv_curr_offlen_ptr[i];
            for (j = recv_curr_offlen_ptr[i]; j < others_req[i].count; j++) {
                if (others_req[i].offsets[j] < iter_end_off) {
                    recv_count[ibuf][i]++;
                    others_req[i].mem_ptrs[j] =
                        (MPI_Count) (others_req[i].offsets[j] - range_off);
                    recv_size[ibuf][i] += others_req[i].lens[j];
                } else {
                    break;
                }
            }
            recv_curr_offlen_ptr[i] = j;
        }
        iter_end_off += step_size;

#ifdef WKL_DEBUG
e_time = MPI_Wtime();
timing[4] += e_time - s_time;
s_time = e_time;
#endif
        /* redistribute (exchange) this process's write requests to I/O
         * aggregators. Communication are Issend and Irecv only. No collective
         * communication. Only aggregators have non-NULL write_buf and
         * recv_buf. All processes have non-NULL send_buf.
         */
        char *wbuf = (write_buf == NULL) ? NULL : write_buf[ibuf];
        char *rbuf = (recv_buf  == NULL) ? NULL :  recv_buf[ibuf];
        send_buf[ibuf] = NULL;

        /* ADIOI_LUSTRE_W_Exchange_data performs one round of communication
         * phase and there are ntimes rounds
         */
        ADIOI_LUSTRE_W_Exchange_data(fd,
                                     buf,                /* IN: user buffer */
                                     wbuf,               /* OUT: write buffer */
                                     &rbuf,              /* OUT: receive buffer */
                                     &send_buf[ibuf],    /* OUT: send buffer */
                                     flat_fview,
                                     flat_bview,
                                     send_size,            /* IN: changed each round */
                                     recv_size[ibuf],      /* IN: changed each round */
                                     range_off,            /* IN: changed each round */
                                     range_size,           /* IN: changed each round */
                                     recv_count[ibuf],     /* IN: changed each round */
                                     recv_start_pos[ibuf], /* IN: changed each round */
                                     others_req,           /* IN: changed each round */
                                     this_buf_idx,         /* IN: changed each round */
                                     &srt_off_len[ibuf],   /* OUT: list of write request off-len pairs */
                                     send_list,            /* OUT: send displacement-length pairs */
                                     recv_list,            /* OUT: send displacement-length pairs */
                                     error_code);

        if (*error_code != MPI_SUCCESS)
            goto over;

        /* rbuf might be realloc-ed */
        if (recv_buf != NULL) recv_buf[ibuf] = rbuf;

#ifdef WKL_DEBUG
e_time = MPI_Wtime();
timing[5] += e_time - s_time;
#endif
        /* commit communication and write for this batch of numBufs */
        if (m % nbufs == nbufs - 1 || m == ntimes - 1) {
            int numBufs = ibuf + 1;

            /* reset ibuf to the first element of nbufs */
            ibuf = 0;

#ifdef WKL_DEBUG
MPI_Barrier(fd->comm);
s_time = MPI_Wtime();
#endif
            /* communication phase */
            commit_comm_phase(fd, send_list, recv_list);
#ifdef WKL_DEBUG
e_time = MPI_Wtime();
timing[1] += e_time - s_time;
s_time = e_time;
#endif

            /* free send_buf allocated in ADIOI_LUSTRE_W_Exchange_data() */
            for (j = 0; j < numBufs; j++) {
                if (send_buf[j] != NULL) {
                    ADIOI_Free(send_buf[j]);
                    send_buf[j] = NULL;
                }
            }
            if (!fd->is_agg) /* non-aggregators are done for this batch */
                continue;    /* next run of loop m */

            /* unpack the data in recv_buf[] into write_buf */
            if (end_loc >= 0) {
                for (j = 0; j < numBufs; j++) {
                    char *buf_ptr = recv_buf[j];
                    for (i = 0; i < nprocs; i++) {
                        if (recv_count[j][i] > 1 && i != myrank) {
                            /* When recv_count[j][i] == 1, this case has
                             * been taken care of earlier by receiving the
                             * message directly into write_buf.
                             */
                            MEMCPY_UNPACK(i, buf_ptr, recv_start_pos[j], recv_count[j],
                                          write_buf[j]);
                            buf_ptr += recv_size[j][i];
                        }
                    }
                }
            }

            /* write to numBufs number of stripes */
            for (j=0; j<numBufs; j++) {

                /* if there is no data to write in round (batch_idx + j) */
                if (srt_off_len[j].num == 0)
                    continue;

                /* range_off  starting file offset of this aggregator's write
                 *            region for this round (may not be aligned to
                 *            stripe boundary)
                 * range_size size (in bytes) of this rank's write region for
                 *            this round, <= striping_unit
                 */
                range_off = off_list[batch_idx + j];
                range_size = MPL_MIN(striping_unit - range_off % striping_unit,
                                     end_loc - range_off + 1);

                /* When srt_off_len[j].num == 1, either there is no hole in the
                 * write buffer or the file domain has been read into write
                 * buffer and updated with the received write data. When
                 * srt_off_len[j].num > 1, holes have been found and the list
                 * of sorted offset-length pairs describing noncontiguous
                 * writes have been constructed. Call writes for each
                 * offset-length pair. Note the offset-length pairs
                 * (represented by srt_off_len[j].off, srt_off_len[j].len, and
                 * srt_off_len[j].num) have been coalesced in
                 * ADIOI_LUSTRE_W_Exchange_data().
                 */

                for (i = 0; i < srt_off_len[j].num; i++) {
                    MPI_Status status;

                    /* all write requests in this round should fall into this
                     * range of [range_off, range_off+range_size). This
                     * assertion should never fail.
                     */
                    ADIOI_Assert(srt_off_len[j].off[i] < range_off + range_size &&
                                 srt_off_len[j].off[i] >= range_off);

                    ADIO_WriteContig(fd,
                                     write_buf[j] + (srt_off_len[j].off[i] - range_off),
                                     srt_off_len[j].len[i],
                                     MPI_BYTE,
                                     ADIO_EXPLICIT_OFFSET,
                                     srt_off_len[j].off[i],
                                     &status, error_code);
                    if (*error_code != MPI_SUCCESS)
                        goto over;
                }
                if (srt_off_len[j].num > 0) {
                    ADIOI_Free(srt_off_len[j].off);
                    ADIOI_Free(srt_off_len[j].len);
                    srt_off_len[j].num = 0;
                }
            }
            batch_idx += numBufs; /* only matters for aggregators */
#ifdef WKL_DEBUG
e_time = MPI_Wtime();
timing[2] += e_time - s_time;
#endif
        }
        else
            ibuf++;
    }

  over:
    if (srt_off_len)
        ADIOI_Free(srt_off_len);
    if (write_buf != NULL)
        ADIOI_Free(write_buf);
    if (recv_buf != NULL) {
        for (j = 0; j < nbufs; j++)
            ADIOI_Free(recv_buf[j]);
        ADIOI_Free(recv_buf);
    }
    if (recv_count != NULL) {
        ADIOI_Free(recv_count[0]);
        ADIOI_Free(recv_count);
    }
    ADIOI_Free(recv_curr_offlen_ptr);
    ADIOI_Free(off_list);
    if (flat_bview->is_contig)
        ADIOI_Free(this_buf_idx);
    if (send_buf != NULL)
        ADIOI_Free(send_buf);
    if (send_list != NULL) {
        for (i = 0; i < cb_nodes; i++) {
            ADIOI_Free(send_list[i].len);
            ADIOI_Free(send_list[i].disp);
        }
        ADIOI_Free(send_list);
    }
    if (recv_list != NULL) {
        for (i = 0; i < nprocs; i++) {
            ADIOI_Free(recv_list[i].len);
            ADIOI_Free(recv_list[i].disp);
        }
        ADIOI_Free(recv_list);
    }
#ifdef WKL_DEBUG
timing[0] = MPI_Wtime() - timing[0];

/* check any pending messages to be received */
MPI_Status probe_st;
int probe_flag;
MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, fd->comm, &probe_flag, &probe_st);
if (probe_flag) {
printf("ERROR ++++ MPI_Iprobe rank=%4d is_agg=%d: ---- cb_nodes=%d ntimes=%lld nbufs=%d\n",myrank,fd->is_agg,cb_nodes,ntimes,nbufs);
fflush(stdout);
}

double max_t[6];
MPI_Reduce(timing, max_t, 6, MPI_DOUBLE, MPI_MAX, 0, fd->comm);
if (myrank == 0) printf("%s ---- time collective=%.4f comm=%.4f pwrite=%.4f others=%.4f %.4f %.4f\n",__func__,max_t[0],max_t[1],max_t[2],max_t[3],max_t[4],max_t[5]);
#endif
}

/* This subroutine is copied from ADIOI_Heap_merge(), but modified to coalesce
 * sorted offset-length pairs whenever possible.
 *
 * Heapify(a, i, heapsize); Algorithm from Cormen et al. pg. 143 modified for a
 * heap with smallest element at root. The recursion has been removed so that
 * there are no function calls. Function calls are too expensive.
 */
static
void heap_merge(const ADIOI_Access *others_req, const MPI_Count *count,
                ADIO_Offset *srt_off, MPI_Count *srt_len,
                const MPI_Count *start_pos, int nprocs, int nprocs_recv,
                MPI_Count *total_elements)
{
    typedef struct {
        ADIO_Offset *off_list;
        ADIO_Offset *len_list;
        MPI_Count nelem;
    } heap_struct;

    heap_struct *a, tmp;
    int i, j, heapsize, l, r, k, smallest;

    a = (heap_struct *) ADIOI_Malloc((nprocs_recv + 1) * sizeof(heap_struct));

    j = 0;
    for (i = 0; i < nprocs; i++) {
        if (count[i]) {
            a[j].off_list = others_req[i].offsets + start_pos[i];
            a[j].len_list = others_req[i].lens + start_pos[i];
            a[j].nelem = count[i];
            j++;
        }
    }

#define SWAP(x, y, tmp) { tmp = x ; x = y ; y = tmp ; }

    heapsize = nprocs_recv;

    /* Build a heap out of the first element from each list, with the smallest
     * element of the heap at the root. The first for loop is to find and move
     * the smallest a[*].off_list[0] to a[0].
     */
    for (i = heapsize / 2 - 1; i >= 0; i--) {
        k = i;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                SWAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }

    /* The heap keeps the smallest element in its first element, i.e.
     * a[0].off_list[0].
     */
    j = 0;
    for (i = 0; i < *total_elements; i++) {
        /* extract smallest element from heap, i.e. the root */
        if (j == 0 || srt_off[j - 1] + srt_len[j - 1] < *(a[0].off_list)) {
            srt_off[j] = *(a[0].off_list);
            srt_len[j] = *(a[0].len_list);
            j++;
        } else {
            /* this offset-length pair can be coalesced into the previous one */
            srt_len[j - 1] = *(a[0].off_list) + *(a[0].len_list) - srt_off[j - 1];
        }
        (a[0].nelem)--;

        if (a[0].nelem) {
            (a[0].off_list)++;
            (a[0].len_list)++;
        } else {
            a[0] = a[heapsize - 1];
            heapsize--;
        }

        /* Heapify(a, 0, heapsize); */
        k = 0;
        for (;;) {
            r = 2 * (k + 1);
            l = r - 1;
            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                SWAP(a[k], a[smallest], tmp);
                k = smallest;
            } else
                break;
        }
    }
    ADIOI_Free(a);
    *total_elements = j;
}

#define CACHE_REQ(list, nelems, buf) {             \
    MPI_Aint buf_addr;                             \
    list.len[list.count] = nelems;                 \
    MPI_Get_address(buf, &buf_addr);               \
    list.disp[list.count] = (MPI_Aint)buf_addr;    \
    list.count++;                                  \
}

static void ADIOI_LUSTRE_W_Exchange_data(
            ADIO_File       fd,
      const void           *buf,          /* user buffer */
            char           *write_buf,    /* OUT: internal buffer used to write to file */
            char          **recv_buf,     /* OUT: [nbufs] internal buffer used to receive */
            char          **send_buf_ptr, /* OUT: [cb_nodes] internal buffer used to send */
            Flat_list      *flat_fview,   /* IN/OUT: flattened file offset-length pairs */
            Flat_list      *flat_bview,   /* IN/OUT: flattened buffer offset-length pairs */
      const MPI_Count      *send_size,    /* [cb_nodes] send_size[i] is amount of this rank sent to aggregator i */
      const MPI_Count      *recv_size,    /* [nprocs] recv_size[i] is amount of this aggregator recv from rank i */
            ADIO_Offset     range_off,    /* starting file offset of this aggregator's write region */
            MPI_Count       range_size,   /* amount of this aggregator's write region */
      const MPI_Count      *recv_count,   /* [nprocs] recv_count[i] is No. offset-length pairs received from rank i */
      const MPI_Count      *start_pos,    /* [nprocs] start_pos[i] starting value of recv_curr_offlen_ptr[i] */
      const ADIOI_Access   *others_req,   /* [nprocs] others_req[i] is rank i's write requests fall into this aggregator's file domain */
      const ADIO_Offset    *buf_idx,      /* [cb_nodes] indices to user buffer offsets for sending this rank's write data to aggregator i */
            off_len_list   *srt_off_len,  /* OUT: list of write offset-length pairs of this aggregator */
            disp_len_list  *send_list,    /* OUT: displacement-length pairs of send buffer */
            disp_len_list  *recv_list,    /* OUT: displacement-length pairs of recv buffer */
            int            *error_code)   /* OUT: */
{
    char *buf_ptr, *contig_buf;
    int i, nprocs, myrank, nprocs_recv, err;
    int hole, check_hole, cb_nodes, striping_unit;
    MPI_Count sum_recv;
    MPI_Status status;
    static char myname[] = "ADIOI_LUSTRE_W_EXCHANGE_DATA";

    *send_buf_ptr = NULL;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);

    cb_nodes = fd->hints->cb_nodes;
    striping_unit = fd->hints->striping_unit;

    /* srt_off_len->num   OUT: number of elements in the above off-len list */
    /* srt_off_len->off[] OUT: list of write offsets by this rank in this round */
    /* srt_off_len->len[] OUT: list of write lengths by this rank in this round */

    /* calculate send receive metadata */
    srt_off_len->num = 0;
    srt_off_len->off = NULL;
    sum_recv = 0;
    nprocs_recv = 0;
    for (i = 0; i < nprocs; i++) {
        srt_off_len->num += recv_count[i];
        sum_recv += recv_size[i];
        if (recv_size[i])
            nprocs_recv++;
    }

    /* determine whether checking holes is necessary */
    check_hole = 1;
    if (srt_off_len->num == 0) {
        /* this process has nothing to receive and hence no hole */
        check_hole = 0;
        hole = 0;
    } else if (fd->hints->ds_write == ADIOI_HINT_AUTO) {
        if (srt_off_len->num > fd->hints->ds_wr_lb) {
            /* Number of offset-length pairs is too large, making merge sort
             * expensive. Skip the sorting in hole checking and proceed with
             * read-modify-write.
             */
            check_hole = 0;
            hole = 1;
        }
        /* else: merge sort is less expensive, proceed to check_hole */
    }
    else if (fd->hints->ds_write == ADIOI_HINT_ENABLE) {
        check_hole = 0;
        hole = 1;
    }
    /* else: fd->hints->ds_write == ADIOI_HINT_DISABLE,
     * proceed to check_hole, as we must construct srt_off_len->off and
     * srt_off_len->len.
     */

    if (check_hole) {
        /* merge the offset-length pairs of all others_req[] (already sorted
         * individually) into a single list of offset-length pairs
         * (srt_off_len->off and srt_off_len->len) in an increasing order of
         * file offsets using a heap-merge sorting algorithm.
         */
        srt_off_len->off = (ADIO_Offset *) ADIOI_Malloc(srt_off_len->num * sizeof(ADIO_Offset));
        srt_off_len->len = (MPI_Count *) ADIOI_Malloc(srt_off_len->num * sizeof(MPI_Count));

        heap_merge(others_req, recv_count, srt_off_len->off, srt_off_len->len,
                   start_pos, nprocs, nprocs_recv, &srt_off_len->num);

        /* srt_off_len->num has been updated in heap_merge() such that
         * srt_off_len->off and srt_off_len->len were coalesced
         */
        hole = (srt_off_len->num > 1);
    }

    /* data sieving */
    if (fd->hints->ds_write != ADIOI_HINT_DISABLE && hole) {
        ADIO_ReadContig(fd, write_buf, range_size, MPI_BYTE, ADIO_EXPLICIT_OFFSET,
                        range_off, &status, &err);
        if (err != MPI_SUCCESS) {
            *error_code = MPIO_Err_create_code(err,
                                               MPIR_ERR_RECOVERABLE,
                                               myname, __LINE__, MPI_ERR_IO, "**ioRMWrdwr", 0);
            return;
        }

        /* Once read, holes have been filled and thus the number of
         * offset-length pairs, srt_off_len->num, becomes one.
         */
        srt_off_len->num = 1;
        if (srt_off_len->off == NULL) {
            srt_off_len->off = (ADIO_Offset *) ADIOI_Malloc(sizeof(ADIO_Offset));
            srt_off_len->len = (MPI_Count *) ADIOI_Malloc(sizeof(MPI_Count));
        }
        srt_off_len->off[0] = range_off;
        srt_off_len->len[0] = range_size;
    }

    /* It is possible sum_recv (sum of message sizes to be received) is larger
     * than the size of collective buffer, write_buf, if writes from multiple
     * remote processes overlap. Receiving messages into overlapped regions of
     * the same write_buffer may cause a problem. To avoid it, we allocate a
     * temporary buffer big enough to receive all messages into disjointed
     * regions. Earlier in ADIOI_LUSTRE_Exch_and_write(), write_buf is already
     * allocated with twice amount of the file stripe size, with the second half
     * to be used to receive messages. If sum_recv is smaller than file stripe
     * size, we can reuse that space. But if sum_recv is bigger (an overlap
     * case, which is rare), we allocate a separate buffer of size sum_recv.
     */
    sum_recv -= recv_size[myrank];
    if (sum_recv > striping_unit)
        *recv_buf = (char *) ADIOI_Realloc(*recv_buf, sum_recv);
    contig_buf = *recv_buf;

    /* cache displacement-length pairs of receive buffer */
    buf_ptr = contig_buf;
    for (i = 0; i < nprocs; i++) {
        if (recv_size[i] == 0)
            continue;
        if (i != myrank) {
            if (recv_count[i] > 1) {
                CACHE_REQ(recv_list[i], recv_size[i], buf_ptr)
                buf_ptr += recv_size[i];
            } else {
                /* recv_count[i] is the number of noncontiguous offset-length
                 * pairs describing the write requests of rank i that fall
                 * into this aggregator's file domain. When recv_count[i] is 1,
                 * there is only one such pair, meaning the receive message is
                 * to be stored contiguously. Such message can be received
                 * directly into write_buf.
                 */
                CACHE_REQ(recv_list[i], recv_size[i],
                          write_buf + others_req[i].mem_ptrs[start_pos[i]])
            }
        } else if (flat_bview->is_contig && recv_count[i] > 0) {
            /* send/recv to/from self uses memcpy()
             * buftype is not contiguous is handled at the send time below.
             */
            char *fromBuf = (char *) buf + buf_idx[fd->my_cb_nodes_index];
            MEMCPY_UNPACK(i, fromBuf, start_pos, recv_count, write_buf);
        }
    }

    if (flat_bview->is_contig) {
        /* If buftype is contiguous, data can be directly sent from user buf
         * at location given by buf_idx.
         */
        for (i = 0; i < cb_nodes; i++) {
            if (send_size[i] && i != fd->my_cb_nodes_index)
                CACHE_REQ(send_list[i], send_size[i], (char*)buf + buf_idx[i]);
        }
    } else {
        char **send_buf, *self_buf;

        /* total send size of this round */
        size_t send_total_size = 0;
        for (i = 0; i < cb_nodes; i++)
            send_total_size += send_size[i];

        if (send_total_size == 0) return;

        /* The user buffer to be used to send in this round is not contiguous,
         * allocate send_buf[], a contiguous space, copy data to send_buf,
         * including ones to be sent to self, and then use send_buf to send.
         */
        send_buf = (char **) ADIOI_Malloc(cb_nodes * sizeof(char *));
        send_buf[0] = (char *) ADIOI_Malloc(send_total_size);
        for (i = 1; i < cb_nodes; i++)
            send_buf[i] = send_buf[i - 1] + send_size[i - 1];

        ADIOI_LUSTRE_Fill_send_buffer(fd, buf, flat_fview, flat_bview,
                                      send_buf, send_total_size, send_size,
                                      &self_buf, send_list);
        /* Send buffers must not be touched before MPI_Waitall() is completed,
         * and thus send_buf will be freed in ADIOI_LUSTRE_Exch_and_write()
         */

        if (fd->my_cb_nodes_index >= 0 && send_size[fd->my_cb_nodes_index] > 0) {
            /* contents of user buf that must be sent to self has been copied
             * into send_buf[fd->my_cb_nodes_index]. Now unpack it into
             * write_buf.
             */
            if (self_buf == NULL) self_buf = send_buf[fd->my_cb_nodes_index];
            MEMCPY_UNPACK(myrank, self_buf, start_pos, recv_count, write_buf);
        }

        *send_buf_ptr = send_buf[0];
        ADIOI_Free(send_buf);
    }
}

static void ADIOI_LUSTRE_Fill_send_buffer(ADIO_File fd,
                                          const void *buf,
                                          Flat_list *flat_fview,
                                          Flat_list *flat_bview,
                                          char **send_buf,
                                          size_t send_total_size,
                                          const MPI_Count *send_size,
                                          char **self_buf,
                                          disp_len_list *send_list)
{
    /* this function is only called if buftype is not contiguous */
    int q, first_q=-1, isUserBuf;
    MPI_Count send_size_rem=0, size, copy_size;
    char *user_buf_ptr, *send_buf_ptr, *same_buf_ptr;
    ADIO_Offset off, len, rem_len, user_buf_idx;

#ifdef WKL_DEBUG
int num_memcpy=0;
#endif

    *self_buf = NULL;

    /* user_buf_idx is to the index offset to buf, indicating the starting
     * location to be copied.
     *
     * flat_bview stores the offset-length pairs of the flattened user buffer
     *     data type. Note this stores offset-length pairs of the data type,
     *     and write amount can be a multiple of the data type.
     * flat_bview->count: the number of pairs
     * flat_bview->off[i]: the ith pair's byte offset to buf. Note the
     *     flattened offsets of user buffer type may not be sorted in an
     *     increasing order, unlike fileview which is required by MPI to be
     *     sorted in a monotonically non-decreasing order.
     * flat_bview->len[i]: length of the ith pair
     * flat_bview->rnd stores the current number of data types being processed.
     * flat_bview->idx: index to the offset-length pair currently being
     *     processed, incremented each round.
     * flat_bview->rem: amount of data in the pair that has not been copied
     *     over, changed each round.
     * flat_bview->extent: extent size of user buffer data type.
     */
    user_buf_idx = flat_bview->extent * flat_bview->rnd
                 + flat_bview->off[flat_bview->idx]
                 + flat_bview->len[flat_bview->idx]
                 - flat_bview->rem;
                 /* in case data left to be copied from previous round */

    /* flat_fview->count: the number of noncontiguous file segments this
     *     rank writes to. Each segment i is described by flat_fview->offs[i]
     *     and flat_fview->len[i].
     * flat_fview->idx: the index to the flat_fview->offs[], flat_fview->len[]
     *     that have been processed in the previous round.
     * The while loop below packs write data into send buffers, send_buf[],
     * based on this rank's off-len pairs in its file view,
     */
    off     = flat_fview->off[flat_fview->idx]
            + flat_fview->len[flat_fview->idx]
            - flat_fview->rem;
    rem_len = flat_fview->rem;

    while (send_total_size > 0) {
        /* this off-len request may span to more than one I/O aggregator */
        while (rem_len != 0) {
            len = rem_len;
            q = ADIOI_LUSTRE_Calc_aggregator(fd, off, &len);
            /* NOTE: len will be modified by ADIOI_Calc_aggregator() to be no
             * more than a file stripe unit size that aggregator "q" is
             * responsible for. Note q is not the MPI rank ID, It is the array
             * index to fd->hints->ranklist[].
             *
             * Now len is the amount of data in ith off-len pair that should be
             * sent to aggregator q. Note q can also be self. In this case,
             * data is also packed into send_buf[q] or pointed to a segment of
             * buf when the data to be packed is contiguous. send_buf[q] will
             * later be copied to write buffer in MEMCPY_UNPACK, instead of
             * calling MPI_Issend to send.
             *
             * send_size[q]: data amount of this rank needs to send to
             * aggregator q in this round.
             *
             * len and send_size[q] are all always <= striping_unit
             */

            if (first_q != q) {
                ADIOI_Assert(send_size_rem == 0);
                first_q = q;
                isUserBuf = 1;
                send_size_rem = send_size[q];
                copy_size = 0;
                same_buf_ptr = (char*)buf + user_buf_idx; /* no increment */
                user_buf_ptr = same_buf_ptr; /* increment after each memcpy */
                if (send_buf != NULL)
                    send_buf_ptr = send_buf[q]; /* increment after each memcpy */
            }

            /* copy len amount of data from buf to send_buf[q] */
            size = len;

            while (size) {
                MPI_Count size_in_buf = MPL_MIN(size, flat_bview->rem);
                copy_size += size_in_buf;
                user_buf_idx += size_in_buf;
                send_size_rem -= size_in_buf;
                flat_bview->rem -= size_in_buf;
                if (flat_bview->rem == 0) { /* move on to next off-len pair */
                    if (! flat_bview->is_contig) {
                        /* user buffer type is not contiguous */
                        if (send_size_rem) {
                            /* after this copy send_buf[q] is still not full */
                            isUserBuf = 0;
                            memcpy(send_buf_ptr, user_buf_ptr, copy_size);
#ifdef WKL_DEBUG
num_memcpy++;
#endif
                            send_buf_ptr += copy_size;
                            copy_size = 0;
                        } else if (isUserBuf == 0) {
                            /* send_buf[q] is full and not using user buf,
                             * copy the remaining delayed data */
                            memcpy(send_buf_ptr, user_buf_ptr, copy_size);
#ifdef WKL_DEBUG
num_memcpy++;
#endif
                        }
                    }
                    /* update flat_bview->idx, flat_bview->rem,
                     * flat_bview->rnd, and user_buf_idx
                     */
                    if (flat_bview->idx < (flat_bview->count - 1))
                        flat_bview->idx++;
                    else {
                        flat_bview->idx = 0;
                        flat_bview->rnd++;
                    }
                    user_buf_idx = flat_bview->off[flat_bview->idx] +
                                   flat_bview->rnd * flat_bview->extent;
                    flat_bview->rem = flat_bview->len[flat_bview->idx];
                    user_buf_ptr = (char*) buf + user_buf_idx;
                }
                else if (send_size_rem == 0 && isUserBuf == 0) {
                    /* flat_bview->rem > 0, send_buf[q] is full, and not using
                     * user buf to send, copy the remaining delayed data
                     */
                    memcpy(send_buf_ptr, user_buf_ptr, copy_size);
#ifdef WKL_DEBUG
num_memcpy++;
#endif
                    user_buf_ptr += copy_size;
                }
                size -= size_in_buf;
            }

            if (send_size_rem == 0) { /* data to q is fully packed */
                first_q = -1;

                if (q != fd->my_cb_nodes_index) { /* send only if not self rank */
                    if (isUserBuf)
                        CACHE_REQ(send_list[q], send_size[q], same_buf_ptr)
                    else
                        CACHE_REQ(send_list[q], send_size[q], send_buf[q])
                }
                else if (isUserBuf) {
                    /* send buffer is also (part of) user's buf. Return the
                     * buffer pointer, so the self send data can be directly
                     * unpack from user buf to write buffer.
                     */
                    *self_buf = same_buf_ptr;
                }
            }
            /* len is the amount of data copied */
            off += len;
            rem_len -= len;
            flat_fview->rem -= len;
            send_total_size -= len;
            if (send_total_size == 0) break;
        }

        /* done with this off-len pair, move on to the next */
        if (flat_fview->rem == 0) {
            if (flat_fview->idx == flat_fview->count-1) {
                /* Note flat_fview->rnd is never used ! This is because
                 * flat_fview consists of all offset-length pairs of this
                 * collective write call. Unlike flat_bview whose offset-length
                 * pairs can be used multiple rounds when 'count' is > 1.
                 */
                ADIOI_Assert(flat_fview->rnd <= 1);
                flat_fview->idx = 0;
                flat_fview->rnd++;
            }
            else
                flat_fview->idx++;
            flat_fview->rem = flat_fview->len[flat_fview->idx];
        }
        off = flat_fview->off[flat_fview->idx];
        rem_len = flat_fview->rem;
    }

#ifdef WKL_DEBUG
if (num_memcpy> 0) printf("---- flat_fview->count=%lld flat_bview->rnd=%lld flat_fview->idx=%lld flat_bview->count=%lld num_memcpy=%d\n",flat_fview->count,flat_bview->rnd,flat_fview->idx,flat_bview->count,num_memcpy);
#endif
}

