/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <adio.h>

MPI_Offset PNCIO_GEN_WriteStrided_naive(PNCIO_File *fd,
                                        const void *buf,
                                        PNCIO_Flat_list buf_view,
                                        MPI_Offset offset)
{
    /* offset is in units of etype relative to the filetype. */

    PNCIO_Flatlist_node *flat_buf;
    /* bwr == buffer write; fwr == file write */
    MPI_Offset bwr_size, fwr_size = 0, sum, size_in_filetype;
    int b_index;
    MPI_Count bufsize;
    MPI_Offset n_etypes_in_filetype;
    MPI_Offset size, n_filetypes, etype_in_filetype;
    MPI_Offset abs_off_in_filetype = 0, req_len;
    MPI_Count filetype_size;
    MPI_Aint lb, filetype_extent;
    int buf_count, filetype_is_contig;
    MPI_Offset userbuf_off;
    MPI_Offset off, req_off, disp, end_offset = 0, start_off;
    MPI_Offset w_len, total_w_len=0;

assert(fd->filetype == MPI_DATATYPE_NULL || fd->filetype == MPI_BYTE);

    if (fd->filetype == MPI_DATATYPE_NULL) {
        // assert(fd->flat_file != NULL);
        MPI_Count n;
        filetype_is_contig = (fd->flat_file.count <= 1);
        filetype_size = 0;
        for (n=0; n<fd->flat_file.count; n++)
            filetype_size += fd->flat_file.blocklens[n];
        filetype_extent = fd->flat_file.indices[fd->flat_file.count-1]
                        + fd->flat_file.blocklens[fd->flat_file.count-1]
                        - fd->flat_file.indices[0];
    }
    else if (fd->filetype == MPI_BYTE) {
        filetype_is_contig = 1;
        filetype_size = 1;
        filetype_extent = 1;
    }
    else {
        PNCIO_Datatype_iscontig(fd->filetype, &filetype_is_contig);
        MPI_Type_size_x(fd->filetype, &filetype_size);
        if (filetype_size == 0)
            return 0;
        MPI_Type_get_extent(fd->filetype, &lb, &filetype_extent);
    }

    bufsize = buf_view.size;

    /* Contiguous both in buftype and filetype should have been handled in a
     * call to PNCIO_WriteContig() earlier.
     */
    assert(!(buf_view.is_contig && filetype_is_contig));

PNCIO_Flatlist_node tmp_buf;
    flat_buf = &tmp_buf;
    flat_buf->count = buf_view.count;
    flat_buf->indices = buf_view.off;
    flat_buf->blocklens = buf_view.len;

    if (!buf_view.is_contig && filetype_is_contig) {
        /* noncontiguous in memory, contiguous in file. */

        off = fd->disp + offset;

        start_off = off;
        end_offset = off + bufsize - 1;

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS)) {
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
        }

        /* for each region in the buffer, grab the data and put it in
         * place
         */
            for (b_index = 0; b_index < flat_buf->count; b_index++) {
                userbuf_off = flat_buf->indices[b_index];
                req_off = off;
                req_len = flat_buf->blocklens[b_index];

                w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                          req_len, req_off);
                if (w_len < 0) return w_len;
                total_w_len += w_len;

                /* off is (potentially) used to save the final offset later */
                off += flat_buf->blocklens[b_index];
            }

        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS)) {
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
        }
    }

    else {      /* noncontiguous in file */
        int f_index, st_index = 0;
        MPI_Offset st_fwr_size, st_n_filetypes;

        /* First we're going to calculate a set of values for use in all
         * the noncontiguous in file cases:
         * start_off - starting byte position of data in file
         * end_offset - last byte offset to be accessed in the file
         * st_n_filetypes - how far into the file we start in terms of
         *                  whole filetypes
         * st_index - index of block in first filetype that we will be
         *            starting in (?)
         * st_fwr_size - size of the data in the first filetype block
         *               that we will write (accounts for being part-way
         *               into writing this block of the filetype
         *
         */

        disp = fd->disp;

        n_etypes_in_filetype = filetype_size;
        n_filetypes = offset / n_etypes_in_filetype;
        etype_in_filetype = offset % n_etypes_in_filetype;
        size_in_filetype = etype_in_filetype;

        sum = 0;
        for (f_index = 0; f_index < fd->flat_file.count; f_index++) {
            sum += fd->flat_file.blocklens[f_index];
            if (sum > size_in_filetype) {
                st_index = f_index;
                fwr_size = sum - size_in_filetype;
                abs_off_in_filetype = fd->flat_file.indices[f_index] +
                    size_in_filetype - (sum - fd->flat_file.blocklens[f_index]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        start_off = disp + n_filetypes * (MPI_Offset) filetype_extent + abs_off_in_filetype;

        st_fwr_size = fwr_size;
        st_n_filetypes = n_filetypes;

        /* start_off, st_n_filetypes, st_index, and st_fwr_size are
         * all calculated at this point
         */

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_off=0 and 100 bytes to be written, end_offset=99
         */
        userbuf_off = 0;
        f_index = st_index;
        off = start_off;
        fwr_size = MPL_MIN(st_fwr_size, bufsize);
        while (userbuf_off < bufsize) {
            userbuf_off += fwr_size;
            end_offset = off + fwr_size - 1;

            if (f_index < (fd->flat_file.count - 1))
                f_index++;
            else {
                f_index = 0;
                n_filetypes++;
            }

            off = disp + fd->flat_file.indices[f_index] + n_filetypes * (MPI_Offset) filetype_extent;
            fwr_size = MPL_MIN(fd->flat_file.blocklens[f_index], bufsize - userbuf_off);
        }

        /* End of calculations.  At this point the following values have
         * been calculated and are ready for use:
         * - start_off
         * - end_offset
         * - st_n_filetypes
         * - st_index
         * - st_fwr_size
         */

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS)) {
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
        }

        if (buf_view.is_contig && !filetype_is_contig) {
            /* contiguous in memory, noncontiguous in file. should be the
             * most common case.
             */

            userbuf_off = 0;
            f_index = st_index;
            off = start_off;
            n_filetypes = st_n_filetypes;
            fwr_size = MPL_MIN(st_fwr_size, bufsize);

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

                if (off + fwr_size < disp + fd->flat_file.indices[f_index] +
                    fd->flat_file.blocklens[f_index] + n_filetypes * (MPI_Offset) filetype_extent) {
                    /* important that this value be correct, as it is
                     * used to set the offset in the fd near the end of
                     * this function.
                     */
                    off += fwr_size;
                }
                /* did not reach end of contiguous block in filetype.
                 * no more I/O needed. off is incremented by fwr_size.
                 */
                else {
                    if (f_index < (fd->flat_file.count - 1))
                        f_index++;
                    else {
                        f_index = 0;
                        n_filetypes++;
                    }
                    off = disp + fd->flat_file.indices[f_index] +
                        n_filetypes * (MPI_Offset) filetype_extent;
                    fwr_size = MPL_MIN(fd->flat_file.blocklens[f_index], bufsize - userbuf_off);
                }
            }
        } else {
            MPI_Offset i_offset, tmp_bufsize = 0;
            /* noncontiguous in memory as well as in file */

            b_index = buf_count = 0;
            i_offset = flat_buf->indices[0];
            f_index = st_index;
            off = start_off;
            n_filetypes = st_n_filetypes;
            fwr_size = st_fwr_size;
            bwr_size = flat_buf->blocklens[0];

            /* while we haven't read size * count bytes, keep going */
            while (tmp_bufsize < bufsize) {
                MPI_Offset new_bwr_size = bwr_size, new_fwr_size = fwr_size;

                size = MPL_MIN(fwr_size, bwr_size);
                /* keep max of a single read amount <= INT_MAX */
                size = MPL_MIN(size, INT_MAX);

                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;

                    w_len = PNCIO_WriteContig(fd, (char *) buf + userbuf_off,
                                              req_len, req_off);
                    if (w_len < 0) return w_len;
                    total_w_len += w_len;
                }

                if (size == fwr_size) {
                    /* reached end of contiguous block in file */
                    if (f_index < (fd->flat_file.count - 1))
                        f_index++;
                    else {
                        f_index = 0;
                        n_filetypes++;
                    }

                    off = disp + fd->flat_file.indices[f_index] +
                        n_filetypes * (MPI_Offset) filetype_extent;

                    new_fwr_size = fd->flat_file.blocklens[f_index];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
                    /* reached end of contiguous block in memory */

                    b_index = (b_index + 1) % flat_buf->count;
                    buf_count++;
                    i_offset = flat_buf->indices[b_index];
                    new_bwr_size = flat_buf->blocklens[b_index];
                    if (size != fwr_size) {
                        off += size;
                        new_fwr_size -= size;
                    }
                }
                tmp_bufsize += size;
                fwr_size = new_fwr_size;
                bwr_size = new_bwr_size;
            }
        }

        /* unlock the file region if we locked it */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS)) {
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
        }
    }   /* end of (else noncontiguous in file) */

    return total_w_len;
}
