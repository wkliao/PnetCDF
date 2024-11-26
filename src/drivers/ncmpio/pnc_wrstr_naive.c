/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy() */
#include <assert.h>
#include <sys/errno.h>

#include <mpi.h>

#include "pnc_lustre.h"

int PNC_GEN_WriteStrided_naive(ADIO_File fd,
                               const void *buf,
                               MPI_Aint count,
                               MPI_Datatype buftype,
                               int file_ptr_type,
                               ADIO_Offset offset,
                               ADIO_Status *status)
{
    /* offset is in units of etype relative to the filetype. */

    ADIOI_Flatlist_node *flat_buf, *flat_file;
    /* bwr == buffer write; fwr == file write */
    ADIO_Offset bwr_size, fwr_size = 0, sum, size_in_filetype;
    int err=NC_NOERR, b_index;
    MPI_Count bufsize;
    ADIO_Offset n_etypes_in_filetype;
    ADIO_Offset size, n_filetypes, etype_in_filetype;
    ADIO_Offset abs_off_in_filetype = 0, req_len;
    int buf_count, buftype_is_contig, filetype_is_contig;
    ADIO_Offset userbuf_off;
    ADIO_Offset off, req_off, disp, start_off;
    ADIO_Status status1;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count buftype_size, buftype_extent;
    MPI_Count lb;
    MPI_Type_size_c(buftype, &buftype_size);
    MPI_Type_get_extent_c(buftype, &lb, &buftype_extent);
#else
    int buftype_size;
    MPI_Aint buftype_extent, lb;
    MPI_Type_size(buftype, &buftype_size);
    MPI_Type_get_extent(buftype, &lb, &buftype_extent);
#endif
    if (!fd->ftype_size) {
        MPI_Status_set_elements(status, buftype, 0);
        return NC_NOERR;
    }

    PNC_Datatype_iscontig(buftype, &buftype_is_contig);
    PNC_Datatype_iscontig(fd->filetype, &filetype_is_contig);

    bufsize = buftype_size * count;

    /* contiguous in buftype and filetype is handled elsewhere */

    if (!buftype_is_contig && filetype_is_contig) {
        int b_count;
        /* noncontiguous in memory, contiguous in file. */

        flat_buf = PNC_Flatten_and_find(buftype);

        off = (file_ptr_type == ADIO_INDIVIDUAL) ? fd->fp_ind :
            fd->disp + (ADIO_Offset) fd->etype_size *offset;

        start_off = off;

        /* for each region in the buffer, grab the data and put it in
         * place
         */
        for (b_count = 0; b_count < count; b_count++) {
            for (b_index = 0; b_index < flat_buf->count; b_index++) {
                userbuf_off = (ADIO_Offset) b_count * buftype_extent +
                    flat_buf->indices[b_index];
                req_off = off;
                req_len = flat_buf->blocklens[b_index];

                PNC_WriteContig(fd, (char *) buf + userbuf_off, req_len, MPI_BYTE,
                                 ADIO_EXPLICIT_OFFSET, req_off, &status1);
                if (err != NC_NOERR)
                    return err;

                /* off is (potentially) used to save the final offset later */
                off += flat_buf->blocklens[b_index];
            }
        }

        if (file_ptr_type == ADIO_INDIVIDUAL)
            fd->fp_ind = off;
    }
    else {      /* noncontiguous in file */
        int f_index, st_index = 0;
        ADIO_Offset st_fwr_size, st_n_filetypes;
        int flag;

        /* First we're going to calculate a set of values for use in all
         * the noncontiguous in file cases:
         * start_off - starting byte position of data in file
         * st_n_filetypes - how far into the file we start in terms of
         *                  whole filetypes
         * st_index - index of block in first filetype that we will be
         *            starting in (?)
         * st_fwr_size - size of the data in the first filetype block
         *               that we will write (accounts for being part-way
         *               into writing this block of the filetype
         *
         */

        flat_file = PNC_Flatten_and_find(fd->filetype);
        disp = fd->disp;

        if (file_ptr_type == ADIO_INDIVIDUAL) {
            start_off = fd->fp_ind;     /* in bytes */
            n_filetypes = -1;
            flag = 0;
            while (!flag) {
                n_filetypes++;
                for (f_index = 0; f_index < flat_file->count; f_index++) {
                    if (disp + flat_file->indices[f_index] +
                        n_filetypes * fd->ftype_extent +
                        flat_file->blocklens[f_index] >= start_off) {
                        /* this block contains our starting position */

                        st_index = f_index;
                        fwr_size = disp + flat_file->indices[f_index] +
                            n_filetypes * fd->ftype_extent +
                            flat_file->blocklens[f_index] - start_off;
                        flag = 1;
                        break;
                    }
                }
            }
        } else {
            n_etypes_in_filetype = fd->ftype_size / fd->etype_size;
            n_filetypes = offset / n_etypes_in_filetype;
            etype_in_filetype = offset % n_etypes_in_filetype;
            size_in_filetype = etype_in_filetype * fd->etype_size;

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
            start_off = disp + n_filetypes * fd->ftype_extent + abs_off_in_filetype;
        }

        st_fwr_size = fwr_size;
        st_n_filetypes = n_filetypes;

        /* start_off, st_n_filetypes, st_index, and st_fwr_size are
         * all calculated at this point
         */

        userbuf_off = 0;
        f_index = st_index;
        off = start_off;
        fwr_size = MIN(st_fwr_size, bufsize);
        while (userbuf_off < bufsize) {
            userbuf_off += fwr_size;

            if (f_index < (flat_file->count - 1))
                f_index++;
            else {
                f_index = 0;
                n_filetypes++;
            }

            off = disp + flat_file->indices[f_index] + n_filetypes * fd->ftype_extent;
            fwr_size = MIN(flat_file->blocklens[f_index], bufsize - userbuf_off);
        }

        /* End of calculations.  At this point the following values have
         * been calculated and are ready for use:
         * - start_off
         * - st_n_filetypes
         * - st_index
         * - st_fwr_size
         */
        if (buftype_is_contig && !filetype_is_contig) {
            /* contiguous in memory, noncontiguous in file. should be the
             * most common case.
             */

            userbuf_off = 0;
            f_index = st_index;
            off = start_off;
            n_filetypes = st_n_filetypes;
            fwr_size = MIN(st_fwr_size, bufsize);

            /* while there is still space in the buffer, write more data */
            while (userbuf_off < bufsize) {
                if (fwr_size) {
                    /* TYPE_UB and TYPE_LB can result in
                     * fwr_size = 0. save system call in such cases */
                    req_off = off;
                    req_len = fwr_size;

                    PNC_WriteContig(fd, (char *) buf + userbuf_off, req_len, MPI_BYTE,
                                     ADIO_EXPLICIT_OFFSET, req_off, &status1);
                    if (err != NC_NOERR)
                        return err;
                }
                userbuf_off += fwr_size;

                if (off + fwr_size < disp + flat_file->indices[f_index] +
                    flat_file->blocklens[f_index] + n_filetypes * fd->ftype_extent) {
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
                        n_filetypes * fd->ftype_extent;
                    fwr_size = MIN(flat_file->blocklens[f_index], bufsize - userbuf_off);
                }
            }
        } else {
            ADIO_Offset i_offset, tmp_bufsize = 0;
            /* noncontiguous in memory as well as in file */

            flat_buf = PNC_Flatten_and_find(buftype);

            b_index = buf_count = 0;
            i_offset = flat_buf->indices[0];
            f_index = st_index;
            off = start_off;
            n_filetypes = st_n_filetypes;
            fwr_size = st_fwr_size;
            bwr_size = flat_buf->blocklens[0];

            /* while we haven't read size * count bytes, keep going */
            while (tmp_bufsize < bufsize) {
                ADIO_Offset new_bwr_size = bwr_size, new_fwr_size = fwr_size;

                size = MIN(fwr_size, bwr_size);
                /* keep max of a single read amount <= INT_MAX */
                size = MIN(size, INT_MAX);

                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;

                    PNC_WriteContig(fd, (char *) buf + userbuf_off, req_len, MPI_BYTE,
                                     ADIO_EXPLICIT_OFFSET, req_off, &status1);
                    if (err != NC_NOERR)
                        return err;
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
                        n_filetypes * fd->ftype_extent;

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
                        (ADIO_Offset) buftype_extent *(ADIO_Offset) (buf_count / flat_buf->count) +
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

        if (file_ptr_type == ADIO_INDIVIDUAL)
            fd->fp_ind = off;
    }   /* end of (else noncontiguous in file) */

    /* This is a temporary way of filling in status. The right way is to
     * keep track of how much data was actually written and placed in buf
     */
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
    MPI_Status_set_elements_x(status, buftype, bufsize);
#else
    MPI_Status_set_elements(status, buftype, bufsize);
#endif

    return NC_NOERR;
}
