/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <adio.h>

#define ADIOI_BUFFERED_READ                                             \
    {                                                                   \
        if (req_off >= readbuf_off + readbuf_len) {                     \
            readbuf_off = req_off;                                      \
            readbuf_len = (MPI_Aint) (MPL_MIN(max_bufsize, end_offset-readbuf_off+1)); \
            PNCIO_ReadContig(fd, readbuf, readbuf_len, MPI_BYTE,         \
                            readbuf_off, &status1, error_code);         \
            if (*error_code != MPI_SUCCESS) {                           \
                *error_code = PNCIO_Err_create_code(*error_code,        \
                    MPIR_ERR_RECOVERABLE, myname, __LINE__, MPI_ERR_IO, \
                    "**iorsrc", 0);                                     \
                return;                                                 \
            }                                                           \
        }                                                               \
        while (req_len > readbuf_off + readbuf_len - req_off) {         \
            ADIOI_Assert((readbuf_off + readbuf_len - req_off) ==       \
                         (MPI_Aint) (readbuf_off + readbuf_len - req_off)); \
            partial_read = (MPI_Aint) (readbuf_off + readbuf_len - req_off); \
            tmp_buf = (char *) ADIOI_Malloc(partial_read);              \
            memcpy(tmp_buf, readbuf+readbuf_len-partial_read, partial_read); \
            ADIOI_Free(readbuf);                                        \
            readbuf = (char *) ADIOI_Malloc(partial_read + max_bufsize); \
            memcpy(readbuf, tmp_buf, partial_read);                     \
            ADIOI_Free(tmp_buf);                                        \
            readbuf_off += readbuf_len-partial_read;                    \
            readbuf_len = (MPI_Aint) (partial_read +                    \
                          MPL_MIN(max_bufsize, end_offset-readbuf_off+1)); \
            PNCIO_ReadContig(fd, readbuf+partial_read, readbuf_len-partial_read, \
                            MPI_BYTE, readbuf_off+partial_read, \
                            &status1, error_code);                      \
            if (*error_code != MPI_SUCCESS) {                           \
                *error_code = PNCIO_Err_create_code(*error_code,        \
                    MPIR_ERR_RECOVERABLE, myname, __LINE__, MPI_ERR_IO, \
                    "**iorsrc", 0);                                     \
                return;                                                 \
            }                                                           \
        }                                                               \
        ADIOI_Assert(req_len == (size_t)req_len);                       \
        memcpy((char *)buf + userbuf_off, readbuf+req_off-readbuf_off, req_len); \
    }


void PNCIO_GEN_ReadStrided(ADIO_File fd, void *buf, MPI_Aint count,
                           MPI_Datatype datatype, MPI_Offset offset,
                           ADIO_Status * status, int *error_code)
{


/* offset is in units of etype relative to the filetype. */

    ADIOI_Flatlist_node *flat_buf, *flat_file;
    MPI_Offset i_offset, new_brd_size, brd_size, size;
    int i, j, k, st_index = 0;
    MPI_Count num, bufsize;
    MPI_Count n_etypes_in_filetype;
    MPI_Offset n_filetypes, etype_in_filetype, st_n_filetypes, size_in_filetype;
    MPI_Offset abs_off_in_filetype = 0, new_frd_size, frd_size = 0, st_frd_size;
    MPI_Count filetype_size, buftype_size, partial_read;
    MPI_Aint lb, filetype_extent, buftype_extent;
    int buf_count, buftype_is_contig, filetype_is_contig;
    MPI_Offset userbuf_off, req_len, sum;
    MPI_Offset off, req_off, disp, end_offset = 0, readbuf_off, start_off;
    char *readbuf, *tmp_buf, *value;
    int info_flag;
    MPI_Aint max_bufsize, readbuf_len;
    ADIO_Status status1;
    static char myname[] = "PNCIO_GEN_ReadStrided";

    if (fd->hints->ds_read == ADIOI_HINT_DISABLE) {
        /* if user has disabled data sieving on reads, use naive
         * approach instead.
         */
        PNCIO_GEN_ReadStrided_naive(fd,
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
        if (fd->flat_file->count > 0) {
            n = fd->flat_file->count - 1;
            filetype_extent = fd->flat_file->indices[n]
                            + fd->flat_file->blocklens[n]
                            - fd->flat_file->indices[0];
        }
        else
            filetype_extent = 0;
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
            MPI_Status_set_elements(status, datatype, 0);
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

    value = (char *) ADIOI_Malloc((MPI_MAX_INFO_VAL + 1) * sizeof(char));
    ADIOI_Info_get(fd->info, "ind_rd_buffer_size", MPI_MAX_INFO_VAL, value, &info_flag);
    max_bufsize = atoi(value);
    ADIOI_Free(value);


    if (!buftype_is_contig && filetype_is_contig) {

/* noncontiguous in memory, contiguous in file. */

        flat_buf = PNCIO_Flatten_and_find(datatype);

        off = fd->disp + offset;

        start_off = off;
        end_offset = off + bufsize - 1;
        readbuf_off = off;
        readbuf = (char *) ADIOI_Malloc(max_bufsize);
        readbuf_len = (MPI_Aint) (MPL_MIN(max_bufsize, end_offset - readbuf_off + 1));

/* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        PNCIO_ReadContig(fd, readbuf, readbuf_len, MPI_BYTE, readbuf_off,
                        &status1, error_code);
        if (*error_code != MPI_SUCCESS)
            return;

        for (j = 0; j < count; j++) {
            for (i = 0; i < flat_buf->count; i++) {
                userbuf_off = (MPI_Offset) j *(MPI_Offset) buftype_extent + flat_buf->indices[i];
                req_off = off;
                req_len = flat_buf->blocklens[i];
                ADIOI_BUFFERED_READ off += flat_buf->blocklens[i];
            }
        }

        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        ADIOI_Free(readbuf);
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
                frd_size = sum - size_in_filetype;
                abs_off_in_filetype = flat_file->indices[i] +
                    size_in_filetype - (sum - flat_file->blocklens[i]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        offset = disp + (MPI_Offset) n_filetypes *filetype_extent + abs_off_in_filetype;

        start_off = offset;

        /* Wei-keng Liao: read request is within a single flat_file contig
         * block e.g. with subarray types that actually describe the whole
         * array */
        if (buftype_is_contig && bufsize <= frd_size) {
            /* a count of bytes can overflow. operate on original type instead */
            PNCIO_ReadContig(fd, buf, count, datatype, offset, status,
                            error_code);

#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
            MPI_Status_set_elements_x(status, datatype, count);
#else
            MPI_Status_set_elements(status, datatype, count);
#endif
            return;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be read, end_offset=99 */

        st_frd_size = frd_size;
        st_n_filetypes = n_filetypes;
        i_offset = 0;
        j = st_index;
        off = offset;
        frd_size = MPL_MIN(st_frd_size, bufsize);
        while (i_offset < bufsize) {
            i_offset += frd_size;
            end_offset = off + frd_size - 1;

            j = (j + 1) % flat_file->count;
            n_filetypes += (j == 0) ? 1 : 0;
            while (flat_file->blocklens[j] == 0) {
                j = (j + 1) % flat_file->count;
                n_filetypes += (j == 0) ? 1 : 0;
            }
            off = disp + flat_file->indices[j] + n_filetypes * (MPI_Offset) filetype_extent;
            frd_size = MPL_MIN(flat_file->blocklens[j], bufsize - i_offset);
        }

/* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        readbuf_off = 0;
        readbuf_len = 0;
        readbuf = (char *) ADIOI_Malloc(max_bufsize);

        if (buftype_is_contig && !filetype_is_contig) {

/* contiguous in memory, noncontiguous in file. should be the most
   common case. */

            i_offset = 0;
            j = st_index;
            off = offset;
            n_filetypes = st_n_filetypes;
            frd_size = MPL_MIN(st_frd_size, bufsize);
            while (i_offset < bufsize) {
                if (frd_size) {
                    /* TYPE_UB and TYPE_LB can result in
                     * frd_size = 0. save system call in such cases */
                    /* lseek(fd->fd_sys, off, SEEK_SET);
                     * err = read(fd->fd_sys, ((char *) buf) + i, frd_size); */

                    req_off = off;
                    req_len = frd_size;
                    userbuf_off = i_offset;
                ADIOI_BUFFERED_READ}
                i_offset += frd_size;

                if (off + frd_size < disp + flat_file->indices[j] +
                    flat_file->blocklens[j] + n_filetypes * (MPI_Offset) filetype_extent)
                    off += frd_size;
                /* did not reach end of contiguous block in filetype.
                 * no more I/O needed. off is incremented by frd_size. */
                else {
                    j = (j + 1) % flat_file->count;
                    n_filetypes += (j == 0) ? 1 : 0;
                    while (flat_file->blocklens[j] == 0) {
                        j = (j + 1) % flat_file->count;
                        n_filetypes += (j == 0) ? 1 : 0;
                    }
                    off = disp + flat_file->indices[j] +
                        n_filetypes * (MPI_Offset) filetype_extent;
                    frd_size = MPL_MIN(flat_file->blocklens[j], bufsize - i_offset);
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
            frd_size = st_frd_size;
            brd_size = flat_buf->blocklens[0];

            while (num < bufsize) {
                size = MPL_MIN(frd_size, brd_size);
                if (size) {
                    /* lseek(fd->fd_sys, off, SEEK_SET);
                     * err = read(fd->fd_sys, ((char *) buf) + i, size); */

                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;
                ADIOI_BUFFERED_READ}

                new_frd_size = frd_size;
                new_brd_size = brd_size;

                if (size == frd_size) {
/* reached end of contiguous block in file */
                    j = (j + 1) % flat_file->count;
                    n_filetypes += (j == 0) ? 1 : 0;
                    while (flat_file->blocklens[j] == 0) {
                        j = (j + 1) % flat_file->count;
                        n_filetypes += (j == 0) ? 1 : 0;
                    }
                    off = disp + flat_file->indices[j] +
                        n_filetypes * (MPI_Offset) filetype_extent;

                    new_frd_size = flat_file->blocklens[j];
                    if (size != brd_size) {
                        i_offset += size;
                        new_brd_size -= size;
                    }
                }

                if (size == brd_size) {
/* reached end of contiguous block in memory */

                    k = (k + 1) % flat_buf->count;
                    buf_count++;
                    i_offset =
                        ((MPI_Offset) buftype_extent *
                         (MPI_Offset) (buf_count / flat_buf->count) + flat_buf->indices[k]);
                    new_brd_size = flat_buf->blocklens[k];
                    if (size != frd_size) {
                        off += size;
                        new_frd_size -= size;
                    }
                }
                ADIOI_Assert(((MPI_Offset) num + size) == (MPI_Aint) (num + size));
                num += size;
                frd_size = new_frd_size;
                brd_size = new_brd_size;
            }
        }

        if ((fd->atomicity) && ADIO_Feature(fd, ADIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        ADIOI_Free(readbuf);    /* malloced in the buffered_read macro */
    }

#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
    MPI_Status_set_elements_x(status, datatype, count);
#else
    MPI_Status_set_elements(status, datatype, count);
#endif
/* This is a temporary way of filling in status. The right way is to
   keep track of how much data was actually read and placed in buf
   by ADIOI_BUFFERED_READ. */
}
