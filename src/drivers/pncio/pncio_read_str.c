/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

#define BUFFERED_READ {                                                       \
    if (req_off >= readbuf_off + readbuf_len) {                               \
        readbuf_off = req_off;                                                \
        readbuf_len = MIN(max_bufsize, end_offset-readbuf_off+1);             \
        r_len = PNCIO_ReadContig(fd, readbuf, readbuf_len, readbuf_off);      \
        if (r_len < 0) return r_len;                                          \
        total_r_len += r_len;                                                 \
    }                                                                         \
    while (req_len > readbuf_off + readbuf_len - req_off) {                   \
        partial_read = readbuf_off + readbuf_len - req_off;                   \
        tmp_buf = (char *) NCI_Malloc(partial_read);                          \
        memcpy(tmp_buf, readbuf+readbuf_len-partial_read, partial_read);      \
        NCI_Free(readbuf);                                                    \
        readbuf = (char *) NCI_Malloc(partial_read + max_bufsize);            \
        memcpy(readbuf, tmp_buf, partial_read);                               \
        NCI_Free(tmp_buf);                                                    \
        readbuf_off += readbuf_len-partial_read;                              \
        readbuf_len = partial_read +                                          \
                      MIN(max_bufsize, end_offset-readbuf_off+1);             \
        r_len = PNCIO_ReadContig(fd, readbuf+partial_read,                    \
                                 readbuf_len-partial_read,                    \
                                 readbuf_off+partial_read);                   \
        if (r_len < 0) return r_len;                                          \
        total_r_len += r_len;                                                 \
    }                                                                         \
    memcpy((char*)buf+userbuf_off, readbuf+req_off-readbuf_off, req_len);     \
}


MPI_Offset PNCIO_GEN_ReadStrided(PNCIO_File *fd,
                                 void       *buf,
                                 PNCIO_View  buf_view)
{
    char *readbuf, *tmp_buf, *value;
    int i, j, k, st_index=0, info_flag;

    MPI_Aint max_bufsize, readbuf_len;
    MPI_Offset i_offset, new_brd_size, brd_size, size;
    MPI_Offset new_frd_size, frd_size=0, st_frd_size, userbuf_off, req_len;
    MPI_Offset off, req_off, end_offset=0, readbuf_off, start_off;
    MPI_Offset r_len, total_r_len=0;
    MPI_Count num, bufsize, partial_read;

    /* This subroutine is entered with file_view being non-contiguous only */

    if (fd->hints->romio_ds_read == PNCIO_HINT_DISABLE) {
        /* if user has disabled data sieving on reads, use naive
         * approach instead.
         */
        return PNCIO_GEN_ReadStrided_naive(fd, buf, buf_view);
    }

    bufsize = buf_view.size;

    /* get max_bufsize from the info object. */
    value = (char *) NCI_Malloc((MPI_MAX_INFO_VAL + 1) * sizeof(char));
    MPI_Info_get(fd->info, "ind_rd_buffer_size", MPI_MAX_INFO_VAL, value, &info_flag);
    max_bufsize = atoi(value);
    NCI_Free(value);

    if (buf_view.count > 1 && fd->file_view.count <= 1) {
        /* noncontiguous in memory, contiguous in file. */

        off = fd->file_view.off[0];

        start_off = off;
        end_offset = start_off + bufsize - 1;
        readbuf_off = start_off;
        readbuf = (char *) NCI_Malloc(max_bufsize);
        readbuf_len = MIN(max_bufsize, end_offset - readbuf_off + 1);

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        r_len = PNCIO_ReadContig(fd, readbuf, readbuf_len, readbuf_off);
        if (r_len < 0) return r_len;

        for (i = 0; i < buf_view.count; i++) {
            userbuf_off = buf_view.off[i];
            req_off = off;
            req_len = buf_view.len[i];
            BUFFERED_READ
            off += buf_view.len[i];
        }

        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        NCI_Free(readbuf);
    }

    else { /* noncontiguous in file */
        MPI_Offset offset=0;
        MPI_Offset sum=0;

        for (i = 0; i < fd->file_view.count; i++) {
            sum += fd->file_view.len[i];
            if (sum > 0) {
                st_index = i;
                frd_size = sum;
                /* abs. offset in bytes in the file */
                offset = fd->file_view.off[i] - (sum - fd->file_view.len[i]);
                break;
            }
        }

        start_off = offset;

        /* Wei-keng Liao: read request is within a single file_view contig
         * block e.g. with subarray types that actually describe the whole
         * array */
        if (buf_view.count <= 1 && bufsize <= frd_size) {
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
        frd_size = MIN(st_frd_size, bufsize);
        while (i_offset < bufsize) {
            i_offset += frd_size;
            end_offset = off + frd_size - 1;

if (i_offset >= bufsize) break;
            j++;
            off = fd->file_view.off[j];
            frd_size = MIN(fd->file_view.len[j], bufsize - i_offset);
        }

        /* if atomicity is true, lock (exclusive) the region to be accessed */
        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        readbuf_off = 0;
        readbuf_len = 0;
        readbuf = (char *) NCI_Malloc(max_bufsize);

        if (buf_view.count <= 1 && fd->file_view.count > 1) {
            /* contiguous in memory, noncontiguous in file should be the most
             * common case.
             */
            i_offset = 0;
            j = st_index;
            off = offset;
            frd_size = MIN(st_frd_size, bufsize);
            while (i_offset < bufsize) {
                if (frd_size) {
                    req_off = off;
                    req_len = frd_size;
                    userbuf_off = i_offset;
                    BUFFERED_READ
                }

                i_offset += frd_size;
                if (i_offset >= bufsize) break;

                if (off + frd_size < fd->file_view.off[j] + fd->file_view.len[j])
                    off += frd_size; /* off is incremented by frd_size */
                else {
                    j++;
assert(j < fd->file_view.count);
                    off = fd->file_view.off[j];
                    frd_size = MIN(fd->file_view.len[j],
                                       bufsize - i_offset);
                }
            }
        } else {
            /* noncontiguous in memory as well as in file */
            k = num = 0;
            i_offset = buf_view.off[0];
            j = st_index;
            off = offset;
            frd_size = st_frd_size;
            brd_size = buf_view.len[0];

            while (num < bufsize) {
                size = MIN(frd_size, brd_size);
                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;
                    BUFFERED_READ
                }

                num += size;
                if (num >= bufsize) break;

                new_frd_size = frd_size;
                new_brd_size = brd_size;

                if (size == frd_size) {
                    /* reached end of contiguous block in file */
                    j++;
assert(j < fd->file_view.count);
                    off = fd->file_view.off[j];
                    new_frd_size = fd->file_view.len[j];
                    if (size != brd_size) {
                        i_offset += size;
                        new_brd_size -= size;
                    }
                }

                if (size == brd_size) {
                    /* reached end of contiguous block in memory */
                    k++;
assert(k < buf_view.count);
                    i_offset = buf_view.off[k];
                    new_brd_size = buf_view.len[k];
                    if (size != frd_size) {
                        off += size;
                        new_frd_size -= size;
                    }
                }
                frd_size = new_frd_size;
                brd_size = new_brd_size;
            }
        }

        if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        NCI_Free(readbuf);    /* malloced in the buffered_read macro */
    }

    assert(total_r_len >= buf_view.size);

    return buf_view.size;
}
