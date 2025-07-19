/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <adio.h>

void ADIOI_GEN_WriteStrided_naive(ADIO_File fd, const void *buf, MPI_Aint count,
                                  MPI_Datatype buftype, MPI_Offset offset,
                                  ADIO_Status *status, int *error_code)
{
    /* offset is in units of etype relative to the filetype. */

    ADIOI_Flatlist_node *flat_buf, *flat_file;
    /* bwr == buffer write; fwr == file write */
    MPI_Offset bwr_size, fwr_size = 0, sum, size_in_filetype;
    int b_index;
    MPI_Count bufsize;
    MPI_Offset n_etypes_in_filetype;
    MPI_Offset size, n_filetypes, etype_in_filetype;
    MPI_Offset abs_off_in_filetype = 0, req_len;
    MPI_Count filetype_size, buftype_size;
    MPI_Aint lb, filetype_extent, buftype_extent;
    int buf_count, buftype_is_contig, filetype_is_contig;
    MPI_Offset userbuf_off;
    MPI_Offset off, req_off, disp, end_offset = 0, start_off;
    ADIO_Status status1;

    *error_code = MPI_SUCCESS;  /* changed below if error */

    if (fd->filetype == MPI_DATATYPE_NULL && fd->flat_file != NULL) {
        MPI_Count n;
        filetype_is_contig = (fd->flat_file->count <= 1);
        filetype_size = 0;
        for (n=0; n<fd->flat_file->count; n++)
            filetype_size += fd->flat_file->blocklens[n];
        filetype_extent = fd->flat_file->indices[fd->flat_file->count-1]
                        + fd->flat_file->blocklens[fd->flat_file->count-1]
                        - fd->flat_file->indices[0];
    }
    else if (fd->filetype == MPI_BYTE) {
        filetype_is_contig = 1;
        filetype_size = 1;
        filetype_extent = 1;
    }
    else {
        ADIOI_Datatype_iscontig(fd->filetype, &filetype_is_contig);
        MPI_Type_size_x(fd->filetype, &filetype_size);
        if (!filetype_size) {
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
            MPI_Status_set_elements_x(status, MPI_BYTE, 0);
#else
            MPI_Status_set_elements(status, MPI_BYTE, 0);
#endif
            *error_code = MPI_SUCCESS;
            return;
        }
        MPI_Type_get_extent(fd->filetype, &lb, &filetype_extent);
    }

    ADIOI_Datatype_iscontig(buftype, &buftype_is_contig);
    MPI_Type_size_x(buftype, &buftype_size);
    MPI_Type_get_extent(buftype, &lb, &buftype_extent);

    bufsize = buftype_size * count;

    /* Contiguous both in buftype and filetype should have been handled in a
     * call to PNCIO_WriteContig() earlier.
     */
    ADIOI_Assert(!(buftype_is_contig && filetype_is_contig));

    if (!buftype_is_contig && filetype_is_contig) {
        int b_count;
        /* noncontiguous in memory, contiguous in file. */

        flat_buf = ADIOI_Flatten_and_find(buftype);

        off = fd->disp + offset;

        start_off = off;
        end_offset = off + bufsize - 1;

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS)) {
            ADIOI_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
        }

        /* for each region in the buffer, grab the data and put it in
         * place
         */
        for (b_count = 0; b_count < count; b_count++) {
            for (b_index = 0; b_index < flat_buf->count; b_index++) {
                userbuf_off = (MPI_Offset) b_count *(MPI_Offset) buftype_extent +
                    flat_buf->indices[b_index];
                req_off = off;
                req_len = flat_buf->blocklens[b_index];

                PNCIO_WriteContig(fd, (char *) buf + userbuf_off, req_len, MPI_BYTE,
                                 req_off, &status1, error_code);
                if (*error_code != MPI_SUCCESS)
                    return;

                /* off is (potentially) used to save the final offset later */
                off += flat_buf->blocklens[b_index];
            }
        }

        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS)) {
            ADIOI_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
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

        flat_file = fd->flat_file;
        disp = fd->disp;

        n_etypes_in_filetype = filetype_size;
        n_filetypes = offset / n_etypes_in_filetype;
        etype_in_filetype = offset % n_etypes_in_filetype;
        size_in_filetype = etype_in_filetype;

        sum = 0;
        for (f_index = 0; f_index < flat_file->count; f_index++) {
            sum += flat_file->blocklens[f_index];
            if (sum > size_in_filetype) {
                st_index = f_index;
                fwr_size = sum - size_in_filetype;
                abs_off_in_filetype = flat_file->indices[f_index] +
                    size_in_filetype - (sum - flat_file->blocklens[f_index]);
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

            if (f_index < (flat_file->count - 1))
                f_index++;
            else {
                f_index = 0;
                n_filetypes++;
            }

            off = disp + flat_file->indices[f_index] + n_filetypes * (MPI_Offset) filetype_extent;
            fwr_size = MPL_MIN(flat_file->blocklens[f_index], bufsize - userbuf_off);
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
        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS)) {
            ADIOI_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
        }

        if (buftype_is_contig && !filetype_is_contig) {
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

                    PNCIO_WriteContig(fd, (char *) buf + userbuf_off, req_len, MPI_BYTE,
                                     req_off, &status1, error_code);
                    if (*error_code != MPI_SUCCESS)
                        return;
                }
                userbuf_off += fwr_size;

                if (off + fwr_size < disp + flat_file->indices[f_index] +
                    flat_file->blocklens[f_index] + n_filetypes * (MPI_Offset) filetype_extent) {
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
                    if (f_index < (flat_file->count - 1))
                        f_index++;
                    else {
                        f_index = 0;
                        n_filetypes++;
                    }
                    off = disp + flat_file->indices[f_index] +
                        n_filetypes * (MPI_Offset) filetype_extent;
                    fwr_size = MPL_MIN(flat_file->blocklens[f_index], bufsize - userbuf_off);
                }
            }
        } else {
            MPI_Offset i_offset, tmp_bufsize = 0;
            /* noncontiguous in memory as well as in file */

            flat_buf = ADIOI_Flatten_and_find(buftype);

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

                    PNCIO_WriteContig(fd, (char *) buf + userbuf_off, req_len, MPI_BYTE,
                                     req_off, &status1, error_code);
                    if (*error_code != MPI_SUCCESS)
                        return;
                }

                if (size == fwr_size) {
                    /* reached end of contiguous block in file */
                    if (f_index < (flat_file->count - 1))
                        f_index++;
                    else {
                        f_index = 0;
                        n_filetypes++;
                    }

                    off = disp + flat_file->indices[f_index] +
                        n_filetypes * (MPI_Offset) filetype_extent;

                    new_fwr_size = flat_file->blocklens[f_index];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
                    /* reached end of contiguous block in memory */

                    b_index = (b_index + 1) % flat_buf->count;
                    buf_count++;
                    i_offset =
                        (MPI_Offset) buftype_extent *(MPI_Offset) (buf_count / flat_buf->count) +
                        flat_buf->indices[b_index];
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
        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS)) {
            ADIOI_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
        }
    }   /* end of (else noncontiguous in file) */

#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
    MPI_Status_set_elements_x(status, buftype, count);
#else
    MPI_Status_set_elements(status, buftype, count);
#endif
    /* This is a temporary way of filling in status. The right way is to
     * keep track of how much data was actually written and placed in buf
     */
}
