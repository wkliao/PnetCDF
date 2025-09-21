/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <adio.h>

#define BUFFERED_READ                                                        \
    {                                                                        \
        if (req_off >= readbuf_off + readbuf_len) {                          \
            readbuf_off = req_off;                                           \
            readbuf_len = MPL_MIN(max_bufsize, end_offset-readbuf_off+1);    \
            r_len = PNCIO_ReadContig(fd, readbuf, readbuf_len, readbuf_off); \
            if (r_len < 0) return r_len;                                     \
            total_r_len += r_len;                                            \
        }                                                                    \
        while (req_len > readbuf_off + readbuf_len - req_off) {              \
            partial_read = readbuf_off + readbuf_len - req_off;              \
            tmp_buf = (char *) NCI_Malloc(partial_read);                     \
            memcpy(tmp_buf, readbuf+readbuf_len-partial_read, partial_read); \
            NCI_Free(readbuf);                                               \
            readbuf = (char *) NCI_Malloc(partial_read + max_bufsize);       \
            memcpy(readbuf, tmp_buf, partial_read);                          \
            NCI_Free(tmp_buf);                                               \
            readbuf_off += readbuf_len-partial_read;                         \
            readbuf_len = partial_read +                                     \
                          MPL_MIN(max_bufsize, end_offset-readbuf_off+1);    \
            r_len = PNCIO_ReadContig(fd, readbuf+partial_read,               \
                             readbuf_len-partial_read,                       \
                             readbuf_off+partial_read);                      \
            if (r_len < 0) return r_len;                                     \
            total_r_len += r_len;                                            \
        }                                                                    \
        memcpy((char*)buf+userbuf_off, readbuf+req_off-readbuf_off, req_len); \
    }


MPI_Offset PNCIO_GEN_ReadStrided(PNCIO_File *fd,
                                 void *buf,
                                 PNCIO_Flat_list buf_view,
                                 MPI_Offset offset)
{


/* offset is in units of etype relative to the filetype. */

    PNCIO_Flatlist_node *flat_buf;
    MPI_Offset i_offset, new_brd_size, brd_size, size;
    int i, j, k, st_index = 0;
    MPI_Count num, bufsize;
    MPI_Offset abs_off_in_filetype = 0, new_frd_size, frd_size = 0, st_frd_size;
    MPI_Count partial_read;
    MPI_Aint buftype_extent;
    int buf_count, filetype_is_contig;
    MPI_Offset userbuf_off, req_len, sum;
    MPI_Offset off, req_off, disp, end_offset = 0, readbuf_off, start_off;
    char *readbuf, *tmp_buf, *value;
    int info_flag;
    MPI_Aint max_bufsize, readbuf_len;
    MPI_Offset r_len, total_r_len=0;

// printf("%s at %d:\n",__func__,__LINE__);

    if (fd->hints->ds_read == PNCIO_HINT_DISABLE) {
        /* if user has disabled data sieving on reads, use naive
         * approach instead.
         */
        return PNCIO_GEN_ReadStrided_naive(fd, buf, buf_view, offset);
    }

/* This subroutine is entered with filetype being non-contiguous only */
assert(fd->filetype == MPI_BYTE);

    filetype_is_contig = (fd->flat_file.count <= 1);
if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

assert(filetype_is_contig == 0);

    if (buf_view.is_contig) buftype_extent = buf_view.size;
    else buftype_extent = buf_view.off[buf_view.count-1] + buf_view.len[buf_view.count-1] - buf_view.off[0];

    bufsize = buf_view.size;

/* get max_bufsize from the info object. */

    value = (char *) NCI_Malloc((MPI_MAX_INFO_VAL + 1) * sizeof(char));
    MPI_Info_get(fd->info, "ind_rd_buffer_size", MPI_MAX_INFO_VAL, value, &info_flag);
    max_bufsize = atoi(value);
    NCI_Free(value);


    PNCIO_Flatlist_node flat_tmp;
    flat_buf = &flat_tmp;
        flat_buf->count = buf_view.count;
        flat_buf->indices = buf_view.off;
        flat_buf->blocklens = buf_view.len;
        // flat_buf->lb_idx = 0;
        // flat_buf->ub_idx = 0;
        // flat_buf->refct = 0;
        // flat_buf->flag = 0;

    if (!buf_view.is_contig && filetype_is_contig) {
assert(0);

/* noncontiguous in memory, contiguous in file. */

        off = fd->disp + offset;

        start_off = off;
        end_offset = off + bufsize - 1;
        readbuf_off = off;
        readbuf = (char *) NCI_Malloc(max_bufsize);
        readbuf_len = MPL_MIN(max_bufsize, end_offset - readbuf_off + 1);

/* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        r_len = PNCIO_ReadContig(fd, readbuf, readbuf_len, readbuf_off);
        if (r_len < 0) return r_len;

        for (i = 0; i < flat_buf->count; i++) {
            userbuf_off = (MPI_Offset)j * buftype_extent + flat_buf->indices[i];
            req_off = off;
            req_len = flat_buf->blocklens[i];
            BUFFERED_READ
            off += flat_buf->blocklens[i];
        }

        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        NCI_Free(readbuf);
    }

    else {      /* noncontiguous in file */
        MPI_Offset size_in_filetype = offset;

        disp = fd->disp;

        sum = 0;
        for (i = 0; i < fd->flat_file.count; i++) {
            sum += fd->flat_file.blocklens[i];
            if (sum > size_in_filetype) {
                st_index = i;
                frd_size = sum - size_in_filetype;
                abs_off_in_filetype = fd->flat_file.indices[i] +
                    size_in_filetype - (sum - fd->flat_file.blocklens[i]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        offset = disp + abs_off_in_filetype;

        start_off = offset;

        /* Wei-keng Liao: read request is within a single flat_file contig
         * block e.g. with subarray types that actually describe the whole
         * array */
        if (buf_view.is_contig && bufsize <= frd_size) {
            /* a count of bytes can overflow. operate on original type instead */
            r_len = PNCIO_ReadContig(fd, buf, buf_view.size, offset);

assert(buf_view.size == r_len);
            return r_len;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be read, end_offset=99 */

        st_frd_size = frd_size;
        i_offset = 0;
        j = st_index;
        off = offset;
        frd_size = MPL_MIN(st_frd_size, bufsize);
        while (i_offset < bufsize) {
            i_offset += frd_size;
            end_offset = off + frd_size - 1;

            j = (j + 1) % fd->flat_file.count;
            while (fd->flat_file.blocklens[j] == 0) {
                j = (j + 1) % fd->flat_file.count;
            }
            off = disp + fd->flat_file.indices[j];
            frd_size = MPL_MIN(fd->flat_file.blocklens[j], bufsize - i_offset);
        }

/* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        readbuf_off = 0;
        readbuf_len = 0;
        readbuf = (char *) NCI_Malloc(max_bufsize);

        if (buf_view.is_contig && !filetype_is_contig) {

/* contiguous in memory, noncontiguous in file. should be the most
   common case. */

            i_offset = 0;
            j = st_index;
            off = offset;
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
                    BUFFERED_READ
                }
                i_offset += frd_size;

                if (off + frd_size < disp + fd->flat_file.indices[j] +
                    fd->flat_file.blocklens[j])
                    off += frd_size;
                /* did not reach end of contiguous block in filetype.
                 * no more I/O needed. off is incremented by frd_size. */
                else {
                    j = (j + 1) % fd->flat_file.count;
                    while (fd->flat_file.blocklens[j] == 0) {
                        j = (j + 1) % fd->flat_file.count;
                    }
                    off = disp + fd->flat_file.indices[j];
                    frd_size = MPL_MIN(fd->flat_file.blocklens[j], bufsize - i_offset);
                }
            }
        } else {
/* noncontiguous in memory as well as in file */

            // flat_buf = PNCIO_Flatten_and_find(datatype);

            k = num = buf_count = 0;
            i_offset = flat_buf->indices[0];
            j = st_index;
            off = offset;
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
                    BUFFERED_READ
                }

                new_frd_size = frd_size;
                new_brd_size = brd_size;

                if (size == frd_size) {
/* reached end of contiguous block in file */
                    j = (j + 1) % fd->flat_file.count;
                    while (fd->flat_file.blocklens[j] == 0) {
                        j = (j + 1) % fd->flat_file.count;
                    }
                    off = disp + fd->flat_file.indices[j];

                    new_frd_size = fd->flat_file.blocklens[j];
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
                         (buf_count / flat_buf->count) + flat_buf->indices[k]);
                    new_brd_size = flat_buf->blocklens[k];
                    if (size != frd_size) {
                        off += size;
                        new_frd_size -= size;
                    }
                }
                num += size;
                frd_size = new_frd_size;
                brd_size = new_brd_size;
            }
        }

        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        NCI_Free(readbuf);    /* malloced in the buffered_read macro */
    }

    return buf_view.size;
}
