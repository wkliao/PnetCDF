/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/errno.h>
#include <unistd.h>   /* pwrite() */

#include <mpi.h>

#include "pncio.h"


/*----< PNCIO_UFS_write_contig() >-------------------------------------------*/
MPI_Offset PNCIO_UFS_write_contig(PNCIO_File *fh,
                                  const void *buf,
                                  MPI_Offset  w_size,
                                  MPI_Offset  offset,
                                  int         is_coll)
{
    char *p;
    ssize_t err = 0;
    size_t w_count;
    MPI_Offset bytes_xfered = 0;

    if (w_size == 0) return NC_NOERR;

#if defined(PNETCDF_DEBUG) && (defined(HAVE_LUSTRE) || defined(MIMIC_LUSTRE))
    if (is_coll) { /* This call arrived from a collective write call */
        int rank, ost_blk;
        static int striping_factor=-1, cb_nodes=-1;
        static MPI_Offset first_ost_id = -1;
        MPI_Offset ost_id;

        MPI_Comm_rank(fh->comm, &rank);

        /* reset first_ost_id when striping_factor or cb_nodes changes */
        if (striping_factor != fh->hints->striping_factor ||
            cb_nodes != fh->hints->cb_nodes)
            first_ost_id = -1;

        if (fh->hints->striping_factor > fh->hints->cb_nodes)
            ost_blk = fh->hints->striping_factor / fh->hints->cb_nodes;
        else
            ost_blk = 1;

        ost_id = (offset / fh->hints->striping_unit) % fh->hints->cb_nodes;
        ost_id /= ost_blk;

        if (first_ost_id == -1)
            first_ost_id = ost_id;
        else if (ost_id != first_ost_id) {
            printf("Warning in %s rank %d: pwrite offset=%lld w_size=%lld ost_id=%lld != 1st ost %lld\n",__func__,rank,offset,w_size,ost_id,first_ost_id);
            first_ost_id = ost_id;
        }
    }
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
    p = (char *) buf;
    while (bytes_xfered < w_size) {
        w_count = w_size - bytes_xfered;
        err = pwrite(fh->fd_sys, p, w_count, offset + bytes_xfered);
        if (err == -1)
            goto err_out;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fh->write_timing[2] += MPI_Wtime() - timing;
#endif

err_out:
    if (err == -1)
        bytes_xfered = ncmpii_error_posix2nc("pwrite");

    return bytes_xfered;
}

/*----< PNCIO_UFS_write_indep() >--------------------------------------------*/
/* This subroutine implements independent write. It consists of two major code
 * segments. The first one is for when data sieving is disabled and the second
 * one enabled.
 *
 * Note in PnetCDF, the file_view and buf_view are never used for more than
 * one round, which greatly simplifies the implementation of this subroutine.
 */
MPI_Offset PNCIO_UFS_write_indep(PNCIO_File *fh,
                                 const void *buf,
                                 PNCIO_View  buf_view)
{
    char *ptr, *cpy_ptr, *tmp_buf=NULL;
    MPI_Count i, j, k, ntimes;
    MPI_Offset lock_off, lock_len, len, total_len=0;

    MPI_Offset tmp_buf_size, file_off, zero=0;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset file_rem, buf_rem, buf_size;
#else
    int file_rem, buf_rem, buf_size;
#endif

#ifdef PNETCDF_DEBUG
    /* When both file_view and buf_view are contiguous, file_write() calls
     * PNCIO_UFS_write_contig().
     */
    assert(!(buf_view.count <= 1 && fh->file_view.count <= 1));

    /* fh->file_view.off and fh->file_view.len should not be NULL */
    assert(fh->file_view.count > 0);
    assert(fh->file_view.off != NULL);
    assert(fh->file_view.len != NULL);

    if (fh->file_view.count == 1)
        assert(fh->file_view.size == fh->file_view.len[0]);

    if (buf_view.count == 1)
        assert(buf_view.size == buf_view.len[0]);

    /* In PnetCDF, fh->file_view.size always == buf_view.size. */
    assert(fh->file_view.size == buf_view.size);
#endif

    if (fh->file_view.size == 0) /* zero-sized request */
        return 0; /* independent I/O can return now */

    lock_off = fh->file_view.off[0];
    if (fh->file_view.count > 1)
        lock_len = fh->file_view.off[fh->file_view.count-1]
                 + fh->file_view.len[fh->file_view.count-1]
                 - lock_off;
    else
        lock_len = fh->file_view.size;

    /* if atomicity is true, lock (exclusive) the whole region */
    if (fh->atomicity)
        PNCIO_WRITE_LOCK(fh, lock_off, SEEK_SET, lock_len);

    if (buf_view.count == 0) {
        /* In this case, the buffer view is contiguous, adjust buf_view.count,
         * buf_view.off, and buf_view.len, so they can be used in the while
         * loop below.
         */
        buf_size = buf_view.size;
        buf_view.count = 1;
        buf_view.off = &zero;
        buf_view.len = &buf_size;
    }

    /* When data sieving write is disabled or file_view is contiguous, write is
     * carried out in multiple rounds. In each round, write data is first
     * aggregated into a contiguous buffer, tmp_buf, and then write it to file
     * at the end of each round. This can improve performance when the write
     * buffer is non-contiguous and file_view is contiguous, i.e. by reducing
     * the number of file writes.
     */
    if (fh->hints->romio_ds_write == PNCIO_HINT_DISABLE ||
        fh->file_view.count <= 1) {

        if (buf_view.count <= 1) { /* use buf directly to write */
            tmp_buf = (char*)buf;
            buf_rem = buf_view.size;
            ntimes = 1;
            tmp_buf_size = buf_view.size;
        }
        else { /* buf is noncontiguous */
            tmp_buf_size = MIN(buf_view.size, fh->hints->ind_wr_buffer_size);
            tmp_buf = (char*) NCI_Malloc(tmp_buf_size);
            buf_rem = buf_view.len[0];
            ntimes = buf_view.size / tmp_buf_size;
            if (buf_view.size % tmp_buf_size)
                ntimes++;
        }

        file_off = fh->file_view.off[0];
        file_rem = fh->file_view.len[0];

        /* pointer to buf, starting location to copy over to tmp_buf */
        cpy_ptr = (char*)buf;

        k = (buf_view.count <= 1) ? 1 : 0; /* whether to skip while loop k */
        j = 0;
        for (i=0; i<ntimes; i++) { /* perform write in ntimes rounds */
            MPI_Offset req_len, tmp_buf_rem;

            /* copy data from buf to tmp_buf */
            tmp_buf_rem = tmp_buf_size;
            ptr = tmp_buf;
            while (k < buf_view.count) {
                req_len = MIN(tmp_buf_rem, buf_rem);
                memcpy(ptr, cpy_ptr, req_len);

                ptr += req_len;
                tmp_buf_rem -= req_len;
                if (tmp_buf_rem == 0) break;

                if (buf_rem == req_len) { /* done with pair k */
                    k++;
                    cpy_ptr = (char*)buf + buf_view.off[k];
                    buf_rem = buf_view.len[k];
                }
                else { /* there is still data remained in pair k */
                    cpy_ptr += req_len;
                    buf_rem -= req_len;
                }
            }

            /* using tmp_buf to write to the file */
            tmp_buf_rem = tmp_buf_size;
            ptr = tmp_buf;
            while (j < fh->file_view.count) {
                req_len = MIN(tmp_buf_rem, file_rem);
                /* write from offset file_off of length req_len */
                len = PNCIO_UFS_write_contig(fh, ptr, req_len, file_off, 0);
                if (len < 0) return len;
                total_len += len;

                ptr += req_len;
                tmp_buf_rem -= req_len;
                if (tmp_buf_rem == 0) break;

                if (file_rem == req_len) { /* done with pair j */
                    j++;
                    file_off = fh->file_view.off[j];
                    file_rem = fh->file_view.len[j];
                }
                else { /* there is still data remained in pair j */
                    file_off += req_len;
                    file_rem -= req_len;
                }
            }
        }

        /* free tmp_buf if allocated */
        if (tmp_buf != buf) NCI_Free(tmp_buf);
    }
    else {
        /* file_view is noncontiguous and data sieving is not disabled */
        MPI_Offset lock_rem, disp, cpy_len;

        /* allocate read-modify-write buffer */
        tmp_buf_size = MIN(lock_len, fh->hints->ind_wr_buffer_size);
        tmp_buf = (char*) NCI_Malloc(tmp_buf_size);

        /* lock_rem is the amount remained to be locked for the entire
         * read-modify-write region.
         */
        lock_rem = lock_len;

        /* perform ntimes rounds of read-modify-write */
        ntimes = lock_len / tmp_buf_size;
        if (lock_len % tmp_buf_size)
            ntimes++;

#ifdef PNETCDF_DEBUG
        /* file_view's offsets should have already sorted into a monotonically
         * non-decreasing order without overlaps. In addition, earlier checks
         * have ensured all fh->file_view.len[] > 0.
         */
        assert(fh->file_view.len[0] > 0);
        for (i=1; i<fh->file_view.count; i++) {
            assert(fh->file_view.off[i-1] < fh->file_view.off[i]);
            assert(fh->file_view.off[i-1] + fh->file_view.len[i-1] <
                   fh->file_view.off[i]);
            assert(fh->file_view.len[i] > 0);
        }
#endif

        /* initialize loop local variables with the 1st pair of file_view and
         * buf_view
         */
        file_off = fh->file_view.off[0];
        file_rem = fh->file_view.len[0];
        buf_rem  = buf_view.len[0];

        /* pointer to buf, starting location to copy over to tmp_buf */
        cpy_ptr = (char*)buf;

        disp=0;
        k = 0; /* index pointed to buf_view's  offset-length pairs */
        j = 0; /* index pointed to file_view's offset-length pairs */
        for (i=0; i<ntimes; i++) { /* perform write in ntimes rounds */
            MPI_Offset req_len, tmp_buf_rem, gap;

            if (disp >= tmp_buf_size) {
                /* This displacement at the beginning of read-modify-write
                 * region of this round i is too large, containing no data to
                 * be copied from the user write buffer. This allows to skip a
                 * few rounds.
                 */
                MPI_Offset skip = disp / tmp_buf_size;
                i        += skip - 1;
                disp     -= skip * tmp_buf_size;
                lock_rem -= skip * tmp_buf_size;
                file_off += skip * tmp_buf_size;
                continue;
            }

            /* read a chunk from the file into tmp_buf */
            req_len = MIN(tmp_buf_size, lock_rem);

            if (!fh->atomicity) /* lock the read-modify-write region */
                PNCIO_WRITE_LOCK(fh, file_off, SEEK_SET, req_len);

            len = PNCIO_UFS_read_contig(fh, tmp_buf, req_len, file_off);
            if (len < 0) return len;

            /* Copy data from buf to tmp_buf. Skip 'disp' bytes at the front
             * for both buffers.
             */
            tmp_buf_rem = MIN(tmp_buf_size, lock_rem) - disp;
            ptr = tmp_buf + disp;

            /* disp will be reset at the end of each round */
            disp = 0;

            while (tmp_buf_rem > 0) {
                /* cpy_len is the amount to be copied over. It should be no
                 * more than remaining of temporary buffer size of this round
                 * i, remaining of the file_view offset-length pair j, and
                 * remaining of buf_view's offset-length pair k.
                 */
                cpy_len = MIN(file_rem, buf_rem);
                cpy_len = MIN(cpy_len, tmp_buf_rem);
                memcpy(ptr, cpy_ptr, cpy_len);
                total_len += cpy_len;

                /* Deduct remaining of temp buffer. Note even if tmp_buf_rem ==
                 * 0, we still need continue to calculate disp for the next
                 * round. */
                tmp_buf_rem -= cpy_len;

                /* advance buffer pointer */
                ptr += cpy_len;

                if (buf_rem == cpy_len) { /* done with pair k */
                    k++;
                    if (k == buf_view.count) /* all data have been copied */
                        break;
                    cpy_ptr = (char*)buf + buf_view.off[k];
                    buf_rem = buf_view.len[k];
                }
                else { /* there is still data remained in pair k */
                    cpy_ptr += cpy_len;
                    buf_rem -= cpy_len;
                }

                if (file_rem == cpy_len) { /* done with pair j */
                    /* When j is the last pair of this round i, the end offset
                     * of tmp_buf may fall between the end of pair j and the
                     * beginning of pair j+1. In this case, the beginning of
                     * file offset region for the next round should advance
                     * disp bytes.
                     */

                    j++;
                    assert(j < fh->file_view.count);
                    /* Note j should never become fh->file_view.count, as the
                     * above check of if (k == buf_view.count) should
                     * short-circuit the while loop. This is because PnetCDF
                     * ensures file_view.size == buf_view.size.
                     */

                    file_rem = fh->file_view.len[j];

                    /* calculate the empty size between pairs j-1 and j */
                    gap = fh->file_view.off[j]
                        - (fh->file_view.off[j-1] + fh->file_view.len[j-1]);

                    if (tmp_buf_rem <= gap) {
                        /* j-1 is last pair of this round */
                        disp = gap - tmp_buf_rem;
                        break;
                    }
                    tmp_buf_rem -= gap;
                    ptr += gap;
                }
                else { /* there is still data remained in pair j */
                    file_rem -= cpy_len;
                }
            }

            /* write the modified chunk to the file */
            len = PNCIO_UFS_write_contig(fh, tmp_buf, req_len, file_off, 0);
            if (len < 0) return len;

            if (!fh->atomicity) /* unlock the read-modify-write region */
                PNCIO_UNLOCK(fh, file_off, SEEK_SET, req_len);

            /* reduce remaining size to be locked */
            lock_rem -= req_len;

            /* update file offset for the next round of read-modify-write */
            file_off += req_len;
        }

        /* free tmp_buf if allocated */
        if (tmp_buf != NULL) NCI_Free(tmp_buf);
    }

    /* if atomicity is true, unlock (exclusive) the whole region */
    if (fh->atomicity)
        PNCIO_UNLOCK(fh, lock_off, SEEK_SET, lock_len);

    return total_len;
}

