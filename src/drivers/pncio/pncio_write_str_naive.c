/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

MPI_Offset PNCIO_GEN_WriteStrided_naive(PNCIO_File *fd,
                                        const void *buf,
                                        PNCIO_View  buf_view)
{
    int b_index;
    MPI_Count bufsize;

    /* bwr == buffer write; fwr == file write */
    MPI_Offset bwr_size, fwr_size = 0, size;
    MPI_Offset req_len, userbuf_off;
    MPI_Offset off, req_off, end_offset = 0, start_off;
    MPI_Offset w_len, total_w_len=0;

    /* Contiguous both in buf_view and file_view should have been handled in a
     * call to PNCIO_WriteContig() earlier.
     */
#ifdef PNETCDF_DEBUG
    assert(!(buf_view.count <= 1 && fd->file_view.count <= 1));
#endif

    bufsize = buf_view.size;

    if (buf_view.count > 1 && fd->file_view.count <= 1) {
        /* noncontiguous in memory, contiguous in file. */

        off = fd->file_view.off[0];

        start_off = off;
        end_offset = off + bufsize - 1;

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        /* for each region in the buffer, grab the data and put it in place */
        for (b_index = 0; b_index < buf_view.count; b_index++) {
            userbuf_off = buf_view.off[b_index];
            req_off = off;
            req_len = buf_view.len[b_index];

            w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                      req_len, req_off);
            if (w_len < 0) return w_len;
            total_w_len += w_len;

            /* off is (potentially) used to save the final offset later */
            off += buf_view.len[b_index];
        }

        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);
    }
    else {      /* noncontiguous in file */
        int f_index, st_index = 0;
        MPI_Offset st_fwr_size;

        /* First we're going to calculate a set of values for use in all
         * the noncontiguous in file cases:
         * start_off - starting byte position of data in file
         * end_offset - last byte offset to be accessed in the file
         * st_index - index of block in first file_view that we will be
         *            starting in (?)
         * st_fwr_size - size of the data in the first file_view block
         *               that we will write (accounts for being part-way
         *               into writing this block of the file_view
         */
#if 0
        MPI_Offset size_in_file_view = fd->file_view.off[0];
        MPI_Offset abs_off_in_file_view = 0;

        MPI_Offset sum = 0;
        for (f_index = 0; f_index < fd->file_view.count; f_index++) {
            sum += fd->file_view.len[f_index];
            if (sum > size_in_file_view) {
                st_index = f_index;
                fwr_size = sum - size_in_file_view;
                abs_off_in_file_view = fd->file_view.off[f_index] +
                    size_in_file_view - (sum - fd->file_view.len[f_index]);
                break;
            }
        }
#endif
        st_index = 0;
        fwr_size = fd->file_view.len[0];

        /* abs. offset in bytes in the file */
        start_off = fd->file_view.off[0];

        st_fwr_size = fwr_size;

        /* start_off, st_index, and st_fwr_size are
         * all calculated at this point
         */

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_off=0 and 100 bytes to be written, end_offset=99
         */
        f_index = st_index;
        fwr_size = MIN(st_fwr_size, bufsize);
        userbuf_off = fwr_size;
        end_offset = start_off + fwr_size - 1;
        while (userbuf_off < bufsize) {
            f_index++;
            fwr_size = MIN(fd->file_view.len[f_index],
                               bufsize - userbuf_off);
            userbuf_off += fwr_size;
            end_offset = fd->file_view.off[f_index] + fwr_size - 1;
        }

        /* End of calculations.  At this point the following values have
         * been calculated and are ready for use:
         * - start_off
         * - end_offset
         * - st_index
         * - st_fwr_size
         */

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        if (buf_view.count <= 1 && fd->file_view.count > 1) {
            /* contiguous in memory, noncontiguous in file. should be the
             * most common case.
             */

            userbuf_off = 0;
            f_index = st_index;
            off = start_off;
            fwr_size = MIN(st_fwr_size, bufsize);

            /* while there is still space in the buffer, write more data */
            while (userbuf_off < bufsize) {
                if (fwr_size) {
                    /* TYPE_UB and TYPE_LB can result in
                     * fwr_size = 0. save system call in such cases */
                    req_off = off;
                    req_len = fwr_size;

                    w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                              req_len, req_off);
                    if (w_len < 0) return w_len;
                    total_w_len += w_len;
                }
                userbuf_off += fwr_size;
                if (userbuf_off >= bufsize) break;

                if (off + fwr_size < fd->file_view.off[f_index] +
                    fd->file_view.len[f_index]) {
                    /* important that this value be correct, as it is
                     * used to set the offset in the fd near the end of
                     * this function.
                     */
                    off += fwr_size;
                }
                /* did not reach end of contiguous block in file_view.
                 * no more I/O needed. off is incremented by fwr_size.
                 */
                else {
                    f_index++;
#ifdef PNETCDF_DEBUG
assert(f_index < fd->file_view.count);
#endif
                    off = fd->file_view.off[f_index];
                    fwr_size = MIN(fd->file_view.len[f_index],
                                       bufsize - userbuf_off);
                }
            }
        } else {
            MPI_Offset i_offset, tmp_bufsize = 0;
            /* noncontiguous in memory as well as in file */

            b_index = 0;
            i_offset = buf_view.off[0];
            f_index = st_index;
            off = start_off;
            fwr_size = st_fwr_size;
            bwr_size = buf_view.len[0];

            /* while we haven't read size * count bytes, keep going */
            while (tmp_bufsize < bufsize) {
                MPI_Offset new_bwr_size = bwr_size, new_fwr_size = fwr_size;

                size = MIN(fwr_size, bwr_size);
                /* keep max of a single read amount <= INT_MAX */
                size = MIN(size, INT_MAX);

                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;

                    w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                              req_len, req_off);
                    if (w_len < 0) return w_len;
                    total_w_len += w_len;
                }

                tmp_bufsize += size;
                if (tmp_bufsize >= bufsize) break;

                if (size == fwr_size) {
                    f_index++;
#ifdef PNETCDF_DEBUG
assert(f_index < fd->file_view.count);
#endif
                    off = fd->file_view.off[f_index];
                    new_fwr_size = fd->file_view.len[f_index];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
                    /* reached end of contiguous block in memory */
                    b_index++;
#ifdef PNETCDF_DEBUG
assert(b_index < buf_view.count);
#endif
                    i_offset = buf_view.off[b_index];
                    new_bwr_size = buf_view.len[b_index];
                    if (size != fwr_size) {
                        off += size;
                        new_fwr_size -= size;
                    }
                }
                fwr_size = new_fwr_size;
                bwr_size = new_bwr_size;
            }
        }

        /* unlock the file region if we locked it */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

    }   /* end of (else noncontiguous in file) */

    return total_w_len;
}
