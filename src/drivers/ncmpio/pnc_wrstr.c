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

#define ADIOI_BUFFERED_WRITE {                                                \
    if (writebuf_len) {                                                       \
        err = PNC_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,           \
                              ADIO_EXPLICIT_OFFSET, writebuf_off,             \
                              &status1);                                      \
        if (err != NC_NOERR) {                                                \
            ADIOI_Free(writebuf);                                             \
            return ncmpii_error_posix2nc("PNC_WriteContig");                  \
        }                                                                     \
        writebuf_off = req_off;                                               \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) *                 \
                           stripe_size - writebuf_off);                       \
        err = PNC_ReadContig(fd, writebuf, writebuf_len, MPI_BYTE,            \
                             ADIO_EXPLICIT_OFFSET,                            \
                             writebuf_off, &status1);                         \
        if (err != NC_NOERR) {                                                \
            ADIOI_Free(writebuf);                                             \
            return ncmpii_error_posix2nc("PNC_ReadContig");                   \
        }                                                                     \
    }                                                                         \
    write_sz = MIN(req_len, writebuf_off + writebuf_len - req_off);           \
    memcpy(writebuf + req_off - writebuf_off,                                 \
           (char *)buf +userbuf_off, write_sz);                               \
    while (write_sz != req_len) {                                             \
        err = PNC_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,           \
                              ADIO_EXPLICIT_OFFSET, writebuf_off, &status1);  \
        if (err != NC_NOERR) {                                                \
            ADIOI_Free(writebuf);                                             \
            return ncmpii_error_posix2nc("PNC_WriteContig");                  \
        }                                                                     \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) *                 \
                           stripe_size - writebuf_off);                       \
        err = PNC_ReadContig(fd, writebuf, writebuf_len, MPI_BYTE,            \
                        ADIO_EXPLICIT_OFFSET, writebuf_off, &status1 );       \
        if (err != NC_NOERR) {                                                \
            ADIOI_Free(writebuf);                                             \
            return ncmpii_error_posix2nc("PNC_ReadContig");                   \
        }                                                                     \
        write_sz = MIN(req_len, writebuf_len);                            \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}

/* this macro is used when filetype is contig and buftype is not contig.
   it does not do a read-modify-write and does not lock*/
#define ADIOI_BUFFERED_WRITE_WITHOUT_READ {                                   \
    if (req_off >= writebuf_off + writebuf_len) {                             \
        err = PNC_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,           \
                              ADIO_EXPLICIT_OFFSET, writebuf_off, &status1);  \
        if (err != NC_NOERR) {                                                \
            ADIOI_Free(writebuf);                                             \
            return ncmpii_error_posix2nc("PNC_WriteContig");                  \
        }                                                                     \
        writebuf_off = req_off;                                               \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) *                 \
                           stripe_size - writebuf_off);                       \
    }                                                                         \
    write_sz = MIN(req_len, writebuf_off + writebuf_len - req_off);       \
    memcpy(writebuf + req_off - writebuf_off,                                 \
           (char *)buf + userbuf_off, write_sz);                              \
    while (write_sz != req_len) {                                             \
        err = PNC_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,           \
                         ADIO_EXPLICIT_OFFSET, writebuf_off, &status1);       \
        if (err != NC_NOERR) {                                                \
            ADIOI_Free(writebuf);                                             \
            return ncmpii_error_posix2nc("PNC_WriteContig");                  \
        }                                                                     \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) *                 \
                           stripe_size - writebuf_off);                       \
        write_sz = MIN(req_len, writebuf_len);                                \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}

int PNC_WriteStrided(ADIO_File fd,
                     const void *buf,
                     MPI_Aint count,
                     MPI_Datatype datatype,
                     int file_ptr_type,
                     ADIO_Offset offset,
                     ADIO_Status *status)
{
    /* offset is in units of etype relative to the filetype. */
    ADIOI_Flatlist_node *flat_buf, *flat_file;
    ADIO_Offset i_offset, sum, size_in_filetype;
    int i, j, k, err=NC_NOERR, st_index = 0;
    MPI_Count n_etypes_in_filetype;
    ADIO_Offset num, size, n_filetypes, etype_in_filetype, st_n_filetypes;
    ADIO_Offset abs_off_in_filetype = 0;
    int buf_count, buftype_is_contig, filetype_is_contig;
    ADIO_Offset userbuf_off;
    ADIO_Offset off, req_off, disp, end_offset = 0, writebuf_off, start_off;
    char *writebuf;
    MPI_Count bufsize, writebuf_len, write_sz;
    ADIO_Status status1;
    ADIO_Offset new_bwr_size, new_fwr_size, st_fwr_size, fwr_size = 0, bwr_size, req_len;
    int stripe_size;

    if (fd->hints->ds_write == ADIOI_HINT_DISABLE) {
        /* if user has disabled data sieving on writes, use naive
         * approach instead.
         */
        return PNC_GEN_WriteStrided_naive(fd, buf, count, datatype,
                                          file_ptr_type, offset, status);
    }

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count buftype_size, buftype_extent;
    MPI_Count lb;
    MPI_Type_size_c(datatype, &buftype_size);
    MPI_Type_get_extent_c(datatype, &lb, &buftype_extent);
#else
    int buftype_size;
    MPI_Aint buftype_extent, lb;
    MPI_Type_size(datatype, &buftype_size);
    MPI_Type_get_extent(datatype, &lb, &buftype_extent);
#endif
    if (!fd->ftype_size) {
        MPI_Status_set_elements(status, datatype, 0);
        return NC_NOERR;
    }

    PNC_Datatype_iscontig(datatype, &buftype_is_contig);
    PNC_Datatype_iscontig(fd->filetype, &filetype_is_contig);

    bufsize = buftype_size * count;

    /* get striping info */
    stripe_size = fd->hints->striping_unit;

    /* Different buftype to different filetype */
    if (!buftype_is_contig && filetype_is_contig) {
        /* noncontiguous in memory, contiguous in file. */
        flat_buf = PNC_Flatten_and_find(datatype);

        off = (file_ptr_type == ADIO_INDIVIDUAL) ? fd->fp_ind :
              fd->disp + offset * fd->etype_size;

        start_off = off;
        end_offset = start_off + bufsize - 1;
        /* write stripe size buffer each time */
        writebuf = (char *) ADIOI_Malloc(MIN(bufsize, stripe_size));
        writebuf_off = 0;
        writebuf_len = 0;

        for (j = 0; j < count; j++) {
            for (i = 0; i < flat_buf->count; i++) {
                userbuf_off = (ADIO_Offset) j * buftype_extent
                            + flat_buf->indices[i];
                req_off = off;
                req_len = flat_buf->blocklens[i];
                ADIOI_BUFFERED_WRITE_WITHOUT_READ;
                off += flat_buf->blocklens[i];
            }
        }

        /* write the buffer out finally */
        err = PNC_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,
                              ADIO_EXPLICIT_OFFSET, writebuf_off, &status1);
        ADIOI_Free(writebuf);
        if (err != NC_NOERR)
            return err;

        if (file_ptr_type == ADIO_INDIVIDUAL)
            fd->fp_ind = off;
    }
    else {
        /* noncontiguous in file */
        flat_file = PNC_Flatten_and_find(fd->filetype);
        disp = fd->disp;

        if (file_ptr_type == ADIO_INDIVIDUAL) {
            /* Wei-keng reworked type processing to be a bit more efficient */
            offset = fd->fp_ind - disp;
            n_filetypes = (offset - flat_file->indices[0]) / fd->ftype_extent;
            offset -= n_filetypes * fd->ftype_extent;
            /* now offset is local to this extent */

            /* find the block where offset is located, skip blocklens[i]==0 */
            for (i = 0; i < flat_file->count; i++) {
                ADIO_Offset dist;
                if (flat_file->blocklens[i] == 0)
                    continue;
                dist = flat_file->indices[i] + flat_file->blocklens[i] - offset;
                /* fwr_size is from offset to the end of block i */
                if (dist == 0) {
                    i++;
                    offset = flat_file->indices[i];
                    fwr_size = flat_file->blocklens[i];
                    break;
                }
                if (dist > 0) {
                    fwr_size = dist;
                    break;
                }
            }
            st_index = i;       /* starting index in flat_file->indices[] */
            offset += disp + n_filetypes * fd->ftype_extent;
        } else {
            n_etypes_in_filetype = fd->ftype_size / fd->etype_size;
            n_filetypes = offset / n_etypes_in_filetype;
            etype_in_filetype = offset % n_etypes_in_filetype;
            size_in_filetype = etype_in_filetype * fd->etype_size;

            sum = 0;
            for (i = 0; i < flat_file->count; i++) {
                sum += flat_file->blocklens[i];
                if (sum > size_in_filetype) {
                    st_index = i;
                    fwr_size = sum - size_in_filetype;
                    abs_off_in_filetype = flat_file->indices[i] +
                        size_in_filetype - (sum - flat_file->blocklens[i]);
                    break;
                }
            }

            /* abs. offset in bytes in the file */
            offset = disp + n_filetypes * fd->ftype_extent + abs_off_in_filetype;
        }

        start_off = offset;

        /* Wei-keng Liao:write request is within single flat_file
         * contig block*/
        /* this could happen, for example, with subarray types that are
         * actually fairly contiguous */
        if (buftype_is_contig && bufsize <= fwr_size) {
            req_off = start_off;
            req_len = bufsize;
            end_offset = start_off + bufsize - 1;
            writebuf = (char *) ADIOI_Malloc(MIN(bufsize, stripe_size));
            memset(writebuf, -1, MIN(bufsize, stripe_size));
            writebuf_off = 0;
            writebuf_len = 0;
            userbuf_off = 0;
            ADIOI_BUFFERED_WRITE_WITHOUT_READ;
            /* write the buffer out finally */
            err = PNC_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,
                                  ADIO_EXPLICIT_OFFSET, writebuf_off, &status1);
            if (err != NC_NOERR) {
                NCI_Free(writebuf);
                return err;
            }

            if (file_ptr_type == ADIO_INDIVIDUAL) {
                /* update MPI-IO file pointer to point to the first byte
                 * that can be accessed in the fileview. */
                fd->fp_ind = offset + bufsize;
                if (bufsize == fwr_size) {
                    do {
                        st_index++;
                        if (st_index == flat_file->count) {
                            st_index = 0;
                            n_filetypes++;
                        }
                    } while (flat_file->blocklens[st_index] == 0);
                    fd->fp_ind = disp + flat_file->indices[st_index]
                               + n_filetypes * fd->ftype_extent;
                }
            }
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
            MPI_Status_set_elements_x(status, datatype, bufsize);
#else
            MPI_Status_set_elements(status, datatype, bufsize);
#endif
            ADIOI_Free(writebuf);
            return NC_NOERR;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be write, end_offset=99 */

        st_fwr_size = fwr_size;
        st_n_filetypes = n_filetypes;
        i_offset = 0;
        j = st_index;
        off = offset;
        fwr_size = MIN(st_fwr_size, bufsize);
        while (i_offset < bufsize) {
            i_offset += fwr_size;
            end_offset = off + fwr_size - 1;

            j = (j + 1) % flat_file->count;
            n_filetypes += (j == 0) ? 1 : 0;
            while (flat_file->blocklens[j] == 0) {
                j = (j + 1) % flat_file->count;
                n_filetypes += (j == 0) ? 1 : 0;
            }

            off = disp + flat_file->indices[j] + n_filetypes * fd->ftype_extent;
            fwr_size = MIN(flat_file->blocklens[j], bufsize - i_offset);
        }

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) ADIOI_Malloc(stripe_size);
        memset(writebuf, -1, stripe_size);

        if (buftype_is_contig && !filetype_is_contig) {
            /* contiguous in memory, noncontiguous in file. should be the most
             * common case.
             */
            i_offset = 0;
            j = st_index;
            off = offset;
            n_filetypes = st_n_filetypes;
            fwr_size = MIN(st_fwr_size, bufsize);
            while (i_offset < bufsize) {
                if (fwr_size) {
                    /* TYPE_UB and TYPE_LB can result in
                     * fwr_size = 0. save system call in such cases */
                    /* lseek(fd->fd_sys, off, SEEK_SET);
                     * err = write(fd->fd_sys, ((char *) buf) + i_offset, fwr_size); */

                    req_off = off;
                    req_len = fwr_size;
                    userbuf_off = i_offset;
                    if (req_off >= writebuf_off + writebuf_len)
                        ADIOI_BUFFERED_WRITE;
                }
                i_offset += fwr_size;

                if (off + fwr_size < disp + flat_file->indices[j] +
                    flat_file->blocklens[j] + n_filetypes * fd->ftype_extent)
                    off += fwr_size;
                /* did not reach end of contiguous block in filetype.
                 * no more I/O needed. off is incremented by fwr_size. */
                else {
                    j = (j + 1) % flat_file->count;
                    n_filetypes += (j == 0) ? 1 : 0;
                    while (flat_file->blocklens[j] == 0) {
                        j = (j + 1) % flat_file->count;
                        n_filetypes += (j == 0) ? 1 : 0;
                    }
                    off = disp + flat_file->indices[j]
                        + n_filetypes * fd->ftype_extent;
                    fwr_size = MIN(flat_file->blocklens[j], bufsize - i_offset);
                }
            }
        } else {
            /* noncontiguous in memory as well as in file */
            flat_buf = PNC_Flatten_and_find(datatype);

            k = num = buf_count = 0;
            i_offset = flat_buf->indices[0];
            j = st_index;
            off = offset;
            n_filetypes = st_n_filetypes;
            fwr_size = st_fwr_size;
            bwr_size = flat_buf->blocklens[0];

            while (num < bufsize) {
                size = MIN(fwr_size, bwr_size);
                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;
                    if (req_off >= writebuf_off + writebuf_len)
                        ADIOI_BUFFERED_WRITE;
                }

                new_fwr_size = fwr_size;
                new_bwr_size = bwr_size;

                if (size == fwr_size) {
                    /* reached end of contiguous block in file */
                    j = (j + 1) % flat_file->count;
                    n_filetypes += (j == 0) ? 1 : 0;
                    while (flat_file->blocklens[j] == 0) {
                        j = (j + 1) % flat_file->count;
                        n_filetypes += (j == 0) ? 1 : 0;
                    }

                    off = disp + flat_file->indices[j]
                        + n_filetypes * fd->ftype_extent;

                    new_fwr_size = flat_file->blocklens[j];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
                    /* reached end of contiguous block in memory */
                    k = (k + 1) % flat_buf->count;
                    buf_count++;
                    i_offset = (ADIO_Offset) buftype_extent *
                          (buf_count / flat_buf->count) + flat_buf->indices[k];
                    new_bwr_size = flat_buf->blocklens[k];
                    if (size != fwr_size) {
                        off += size;
                        new_fwr_size -= size;
                    }
                }
                num += size;
                fwr_size = new_fwr_size;
                bwr_size = new_bwr_size;
            }
        }

        /* write the buffer out finally */
        if (writebuf_len)
            err = PNC_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,
                                  ADIO_EXPLICIT_OFFSET, writebuf_off, &status1);
        ADIOI_Free(writebuf);
        if (err != NC_NOERR)
            return err;

        if (file_ptr_type == ADIO_INDIVIDUAL)
            fd->fp_ind = off;
    }

    /* This is a temporary way of filling in status. The right way is to
     * keep track of how much data was actually written by
     * ADIOI_BUFFERED_WRITE.
     */
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
    MPI_Status_set_elements_x(status, datatype, bufsize);
#else
    MPI_Status_set_elements(status, datatype, bufsize);
#endif

    return NC_NOERR;
}

