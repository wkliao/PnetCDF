/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <adio.h>

#define BUFFERED_WRITE                                                  \
    {                                                                   \
        if (req_off >= writebuf_off + writebuf_len) {                   \
            if (writebuf_len) {                                         \
                PNCIO_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE, \
                                  writebuf_off, &status1, error_code);  \
                if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                    PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
                if (*error_code != MPI_SUCCESS) {                       \
                    *error_code = PNCIO_Err_create_code(*error_code,    \
                        MPIR_ERR_RECOVERABLE, myname, __LINE__,         \
                        MPI_ERR_IO, "**iowswc", 0);                     \
                    goto fn_exit;                                       \
                }                                                       \
            }                                                           \
            writebuf_off = req_off;                                     \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            PNCIO_ReadContig(fd, writebuf, writebuf_len, MPI_BYTE,      \
                             writebuf_off, &status1, error_code);       \
            if (*error_code != MPI_SUCCESS) {                           \
                *error_code = PNCIO_Err_create_code(*error_code,        \
                    MPIR_ERR_RECOVERABLE, myname, __LINE__, MPI_ERR_IO, \
                    "**iowsrc", 0);                                     \
                goto fn_exit;                                           \
            }                                                           \
        }                                                               \
        write_sz = (MPI_Aint) (MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        assert((MPI_Offset)write_sz == MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        memcpy(writebuf+req_off-writebuf_off, (char *)buf +userbuf_off, write_sz); \
        while (write_sz != req_len) {                                   \
            PNCIO_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,     \
                              writebuf_off, &status1, error_code);      \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            if (*error_code != MPI_SUCCESS) {                           \
                *error_code = PNCIO_Err_create_code(*error_code,        \
                    MPIR_ERR_RECOVERABLE, myname, __LINE__, MPI_ERR_IO, \
                    "**iowswc", 0);                                     \
                goto fn_exit;                                           \
            }                                                           \
            req_len -= write_sz;                                        \
            userbuf_off += write_sz;                                    \
            writebuf_off += writebuf_len;                               \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            PNCIO_ReadContig(fd, writebuf, writebuf_len, MPI_BYTE,      \
                             writebuf_off, &status1, error_code);       \
            if (*error_code != MPI_SUCCESS) {                           \
                *error_code = PNCIO_Err_create_code(*error_code,        \
                    MPIR_ERR_RECOVERABLE, myname, __LINE__, MPI_ERR_IO, \
                    "**iowsrc", 0);                                     \
                goto fn_exit;                                           \
            }                                                           \
            write_sz = MPL_MIN(req_len, writebuf_len);                  \
            memcpy(writebuf, (char *)buf + userbuf_off, write_sz);      \
        }                                                               \
    }


/* this macro is used when filetype is contig and buftype is not contig.
   it does not do a read-modify-write and does not lock*/
#define BUFFERED_WRITE_WITHOUT_READ                                     \
    {                                                                   \
        if (req_off >= writebuf_off + writebuf_len) {                   \
            PNCIO_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,      \
                             writebuf_off, &status1, error_code);       \
            if (*error_code != MPI_SUCCESS) {                           \
                *error_code = PNCIO_Err_create_code(*error_code,        \
                    MPIR_ERR_RECOVERABLE, myname, __LINE__, MPI_ERR_IO, \
                    "**iowswc", 0);                                     \
                goto fn_exit;                                           \
            }                                                           \
            writebuf_off = req_off;                                     \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
        }                                                               \
        write_sz = (MPI_Aint) (MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        assert((MPI_Offset)write_sz == MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        memcpy(writebuf+req_off-writebuf_off, (char *)buf +userbuf_off, write_sz); \
        while (write_sz != req_len) {                                   \
            PNCIO_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,      \
                             writebuf_off, &status1, error_code);       \
            if (*error_code != MPI_SUCCESS) {                           \
                *error_code = PNCIO_Err_create_code(*error_code,        \
                    MPIR_ERR_RECOVERABLE, myname, __LINE__, MPI_ERR_IO, \
                    "**iowswc", 0);                                     \
                goto fn_exit;                                           \
            }                                                           \
            req_len -= write_sz;                                        \
            userbuf_off += write_sz;                                    \
            writebuf_off += writebuf_len;                               \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
            write_sz = MPL_MIN(req_len, writebuf_len);                  \
            memcpy(writebuf, (char *)buf + userbuf_off, write_sz);      \
        }                                                               \
    }
void PNCIO_GEN_WriteStrided(PNCIO_File *fd, const void *buf, MPI_Aint count,
                            MPI_Datatype datatype, MPI_Offset offset,
                            MPI_Status * status, int *error_code)
{

/* offset is in units of etype relative to the filetype. */

    PNCIO_Flatlist_node *flat_buf, *flat_file;
    MPI_Offset i_offset, sum, size_in_filetype;
    int i, j, k, st_index = 0;
    MPI_Offset num, size, n_filetypes, etype_in_filetype, st_n_filetypes;
    MPI_Offset n_etypes_in_filetype, abs_off_in_filetype = 0;
    MPI_Count filetype_size, buftype_size;
    MPI_Aint lb, filetype_extent, buftype_extent;
    int buf_count, buftype_is_contig, filetype_is_contig;
    MPI_Offset userbuf_off;
    MPI_Offset off, req_off, disp, end_offset = 0, writebuf_off, start_off;
    char *writebuf = NULL;
    MPI_Aint writebuf_len, max_bufsize, write_sz;
    MPI_Aint bufsize;
    MPI_Status status1;
    MPI_Offset new_bwr_size, new_fwr_size, st_fwr_size, fwr_size = 0, bwr_size, req_len;
    static char myname[] = "PNCIO_GEN_WriteStrided";

    if (fd->hints->ds_write == PNCIO_HINT_DISABLE) {
        /* if user has disabled data sieving on reads, use naive
         * approach instead.
         */

        PNCIO_GEN_WriteStrided_naive(fd,
                                     buf,
                                     count, datatype, offset, status, error_code);
        return;
    }


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
        PNCIO_Datatype_iscontig(fd->filetype, &filetype_is_contig);
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

    PNCIO_Datatype_iscontig(datatype, &buftype_is_contig);
    MPI_Type_size_x(datatype, &buftype_size);
    MPI_Type_get_extent(datatype, &lb, &buftype_extent);

    bufsize = buftype_size * count;

/* get max_bufsize from the info object. */

    max_bufsize = fd->hints->ind_wr_buffer_size;

    /* Contiguous both in buftype and filetype should have been handled in a
     * call to PNCIO_WriteContig() earlier.
     */
    assert(!(buftype_is_contig && filetype_is_contig));

    if (!buftype_is_contig && filetype_is_contig) {

/* noncontiguous in memory, contiguous in file. */

        flat_buf = PNCIO_Flatten_and_find(datatype);

        off = fd->disp + offset;

        start_off = off;
        end_offset = off + bufsize - 1;
        writebuf_off = off;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        writebuf_len = MPL_MIN(max_bufsize, end_offset - writebuf_off + 1);

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        for (j = 0; j < count; j++) {
            for (i = 0; i < flat_buf->count; i++) {
                userbuf_off = (MPI_Offset) j *(MPI_Offset) buftype_extent + flat_buf->indices[i];
                req_off = off;
                req_len = flat_buf->blocklens[i];
                BUFFERED_WRITE_WITHOUT_READ;
                off += flat_buf->blocklens[i];
            }
        }

        /* write the buffer out finally */
        if (writebuf_len) {
            PNCIO_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,
                             writebuf_off, &status1, error_code);
        }

        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        if (*error_code != MPI_SUCCESS)
            goto fn_exit;
    }

    else {      /* noncontiguous in file */

        flat_file = fd->flat_file;
        disp = fd->disp;

        n_etypes_in_filetype = filetype_size;
        n_filetypes = offset / n_etypes_in_filetype;
        etype_in_filetype = offset % n_etypes_in_filetype;
        size_in_filetype = etype_in_filetype;

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
        offset = disp + (MPI_Offset) n_filetypes *filetype_extent + abs_off_in_filetype;

        start_off = offset;

        /* Wei-keng Liao:write request is within single flat_file contig block */
        /* this could happen, for example, with subarray types that are
         * actually fairly contiguous */
        if (buftype_is_contig && bufsize <= fwr_size) {
            /* though MPI api has an integer 'count' parameter, derived
             * datatypes might describe more bytes than can fit into an integer.
             * if we've made it this far, we can pass a count of original
             * datatypes, instead of a count of bytes (which might overflow)
             * Other WriteContig calls in this path are operating on data
             * sieving buffer */
            PNCIO_WRITE_LOCK(fd, offset, SEEK_SET, bufsize);
            PNCIO_WriteContig(fd, buf, count, datatype, offset, status,
                             error_code);
            PNCIO_UNLOCK(fd, offset, SEEK_SET, bufsize);

            /* bufsize is in bytes */
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
            MPI_Status_set_elements_x(status, MPI_BYTE, bufsize);
#else
            MPI_Status_set_elements(status, MPI_BYTE, bufsize);
#endif
            goto fn_exit;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be write, end_offset=99 */

        st_fwr_size = fwr_size;
        st_n_filetypes = n_filetypes;
        i_offset = 0;
        j = st_index;
        off = offset;
        fwr_size = MPL_MIN(st_fwr_size, bufsize);
        while (i_offset < bufsize) {
            i_offset += fwr_size;
            end_offset = off + fwr_size - 1;

            j = (j + 1) % flat_file->count;
            n_filetypes += (j == 0) ? 1 : 0;
            while (flat_file->blocklens[j] == 0) {
                j = (j + 1) % flat_file->count;
                n_filetypes += (j == 0) ? 1 : 0;
            }

            off = disp + flat_file->indices[j] + n_filetypes * (MPI_Offset) filetype_extent;
            fwr_size = MPL_MIN(flat_file->blocklens[j], bufsize - i_offset);
        }

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        memset(writebuf, -1, max_bufsize);

        if (buftype_is_contig && !filetype_is_contig) {

/* contiguous in memory, noncontiguous in file. should be the most
   common case. */

            i_offset = 0;
            j = st_index;
            off = offset;
            n_filetypes = st_n_filetypes;
            fwr_size = MPL_MIN(st_fwr_size, bufsize);
            while (i_offset < bufsize) {
                if (fwr_size) {
                    /* TYPE_UB and TYPE_LB can result in
                     * fwr_size = 0. save system call in such cases */
                    /* lseek(fd->fd_sys, off, SEEK_SET);
                     * err = write(fd->fd_sys, ((char *) buf) + i_offset, fwr_size); */

                    req_off = off;
                    req_len = fwr_size;
                    userbuf_off = i_offset;
                    BUFFERED_WRITE;
                }
                i_offset += fwr_size;

                if (off + fwr_size < disp + flat_file->indices[j] +
                    flat_file->blocklens[j] + n_filetypes * (MPI_Offset) filetype_extent)
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
                    off = disp + flat_file->indices[j] +
                        n_filetypes * (MPI_Offset) filetype_extent;
                    fwr_size = MPL_MIN(flat_file->blocklens[j], bufsize - i_offset);
                }
            }
        } else {
/* noncontiguous in memory as well as in file */

            flat_buf = PNCIO_Flatten_and_find(datatype);

            k = num = buf_count = 0;
            i_offset = flat_buf->indices[0];
            j = st_index;
            off = offset;
            n_filetypes = st_n_filetypes;
            fwr_size = st_fwr_size;
            bwr_size = flat_buf->blocklens[0];

            while (num < bufsize) {
                size = MPL_MIN(fwr_size, bwr_size);
                if (size) {
                    /* lseek(fd->fd_sys, off, SEEK_SET);
                     * err = write(fd->fd_sys, ((char *) buf) + i_offset, size); */

                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;
                    BUFFERED_WRITE;
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

                    off = disp + flat_file->indices[j] +
                        n_filetypes * (MPI_Offset) filetype_extent;

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
                    i_offset =
                        (MPI_Offset) buftype_extent *(MPI_Offset) (buf_count / flat_buf->count) +
                        flat_buf->indices[k];
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
        if (writebuf_len) {
            PNCIO_WriteContig(fd, writebuf, writebuf_len, MPI_BYTE,
                             writebuf_off, &status1, error_code);
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);
            if (*error_code != MPI_SUCCESS)
                goto fn_exit;
        }
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
    }

    /* bufsize is in bytes */
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
    MPI_Status_set_elements_x(status, MPI_BYTE, bufsize);
#else
    MPI_Status_set_elements(status, MPI_BYTE, bufsize);
#endif
    /* This is a temporary way of filling in status. The right way is to keep
     * track of how much data was actually written by BUFFERED_WRITE.
     */

  fn_exit:
    if (writebuf != NULL)
        NCI_Free(writebuf);
}
