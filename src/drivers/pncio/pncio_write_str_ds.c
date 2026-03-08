/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

#define BUFFERED_WRITE {                                                      \
    if (req_off >= writebuf_off + writebuf_len) {                             \
        if (writebuf_len) {                                                   \
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,             \
                                      writebuf_off);                          \
            if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)  \
                    PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);   \
            if (w_len < 0) goto fn_exit;                                      \
            total_w_len += w_len;                                             \
        }                                                                     \
        writebuf_off = req_off;                                               \
        writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);            \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) goto fn_exit;                                          \
    }                                                                         \
    write_sz = (MPI_Aint)MIN(req_len, writebuf_off+writebuf_len-req_off);     \
    memcpy(writebuf+req_off-writebuf_off, (char*)buf +userbuf_off, write_sz); \
    while (write_sz != req_len) {                                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);           \
        if (w_len < 0) goto fn_exit;                                          \
        total_w_len += w_len;                                                 \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);            \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) goto fn_exit;                                          \
        write_sz = MIN(req_len, writebuf_len);                                \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}

/*----< PNCIO_GEN_ReadStrided_ds() >-----------------------------------------*/
/* This subroutine implements independent write when data sieving is NOT
 * disabled.
 */
MPI_Offset PNCIO_GEN_WriteStrided_ds(PNCIO_File *fd,
                                     const void *buf,
                                     PNCIO_View  buf_view)
{
    char *writebuf = NULL;
    int i, j, k, st_index = 0;
    MPI_Aint writebuf_len, max_bufsize, write_sz, bufsize;
    MPI_Offset i_offset, num, size;
    MPI_Offset userbuf_off, off, req_off, end_offset=0;
    MPI_Offset writebuf_off, start_off, new_bwr_size, new_fwr_size;
    MPI_Offset st_fwr_size, fwr_size = 0, bwr_size, req_len;
    MPI_Offset r_len, w_len, total_w_len=0;

#ifdef PNETCDF_DEBUG
    /* When both file_view and buf_view are contiguous, file_write() calls
     * PNCIO_WriteContig().
     */
    assert(!(buf_view.count <= 1 && fd->file_view.count <= 1));

    /* When data sieving write is disabled, file_write() calls
     * PNCIO_GEN_WriteStrided_nods() */
    assert(fd->hints->romio_ds_write != PNCIO_HINT_DISABLE);
#endif

#ifdef PNETCDF_DEBUG
assert(fd->file_view.size == buf_view.size);
#endif

    bufsize = buf_view.size;

    /* get max_bufsize from the info object. */
    max_bufsize = fd->hints->ind_wr_buffer_size;

    if (buf_view.count > 1 && fd->file_view.count <= 1) {
        /* noncontiguous in memory, contiguous in file. */

        off = fd->file_view.off[0];

        start_off = off;
        end_offset = start_off + bufsize - 1;
        writebuf_off = start_off;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        writebuf_len = MIN(max_bufsize, end_offset - writebuf_off + 1);

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        for (i = 0; i < buf_view.count; i++) {
            userbuf_off = buf_view.off[i];
            req_off = off;
            req_len = buf_view.len[i];

            /* BUFFERED_WRITE_WITHOUT_READ does neither read-modify-write nor
             * file lock
             */
            if (req_off >= writebuf_off + writebuf_len) {
                w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,
                                          writebuf_off);
                if (w_len < 0) goto fn_exit;
                total_w_len += w_len;
                writebuf_off = req_off;
                writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);
            }
            write_sz = MIN(req_len, writebuf_off + writebuf_len - req_off);
            memcpy(writebuf+req_off-writebuf_off, (char*)buf +userbuf_off,
                   write_sz);
            while (write_sz != req_len) {
                w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,
                                          writebuf_off);
                if (w_len < 0) goto fn_exit;
                total_w_len += w_len;
                req_len -= write_sz;
                userbuf_off += write_sz;
                writebuf_off += writebuf_len;
                writebuf_len = MIN(max_bufsize,end_offset-writebuf_off+1);
                write_sz = MIN(req_len, writebuf_len);
                memcpy(writebuf, (char *)buf + userbuf_off, write_sz);
            }

            off += buf_view.len[i];
        }

        /* write the buffer out finally */
        if (writebuf_len) {
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (w_len >= 0) total_w_len += w_len;
        }
        else
            w_len = 0;

        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        if (w_len < 0)
            goto fn_exit;
    }
    else { /* noncontiguous in file */
        MPI_Offset offset=0;
        MPI_Offset sum=0;

        for (i = 0; i < fd->file_view.count; i++) {
            sum += fd->file_view.len[i];
            if (sum > 0) {
                st_index = i;
                fwr_size = sum;
                /* abs. offset in bytes in the file */
                offset = fd->file_view.off[i] - (sum - fd->file_view.len[i]);
                break;
            }
        }

        start_off = offset;

        /* Write request is within single file_view contig block. This could
         * happen, for example, with subarray types that are actually fairly
         * contiguous.
         */
        if (buf_view.count <= 1 && bufsize <= fwr_size) {
            /* though MPI api has an integer 'count' parameter, derived
             * datatypes might describe more bytes than can fit into an integer.
             * if we've made it this far, we can pass a count of original
             * datatypes, instead of a count of bytes (which might overflow)
             * Other WriteContig calls in this path are operating on data
             * sieving buffer */
            PNCIO_WRITE_LOCK(fd, offset, SEEK_SET, bufsize);
            w_len = PNCIO_WriteContig(fd, buf, buf_view.size, offset);
            if (w_len > 0) total_w_len += w_len;
            PNCIO_UNLOCK(fd, offset, SEEK_SET, bufsize);

            goto fn_exit;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be write, end_offset=99 */

        st_fwr_size = fwr_size;
        j = st_index;
        fwr_size = MIN(fwr_size, bufsize);
        i_offset = fwr_size;
        end_offset = offset + fwr_size - 1;
        while (i_offset < bufsize) {
            j++;
            fwr_size = MIN(fd->file_view.len[j], bufsize - i_offset);
            i_offset += fwr_size;
            end_offset = fd->file_view.off[j] + fwr_size - 1;
        }

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        memset(writebuf, -1, max_bufsize);

        if (buf_view.count <= 1 && fd->file_view.count > 1) {
            /* contiguous in memory, noncontiguous in file should be the most
             * common case.
             */
            i_offset = 0;
            j = st_index;
            off = offset;
            fwr_size = MIN(st_fwr_size, bufsize);
            while (i_offset < bufsize) {
                if (fwr_size) {
                    req_off = off;
                    req_len = fwr_size;
                    userbuf_off = i_offset;
                    BUFFERED_WRITE;
                }

                i_offset += fwr_size;
                if (i_offset >= bufsize) break;

                if (off + fwr_size < fd->file_view.off[j] +
                                     fd->file_view.len[j])
                    off += fwr_size; /* off is incremented by fwr_size. */
                else {
                    j++;
#ifdef PNETCDF_DEBUG
assert(j < fd->file_view.count);
#endif
                    off = fd->file_view.off[j];
                    fwr_size = MIN(fd->file_view.len[j],
                                       bufsize - i_offset);
                }
            }
        } else {
            /* noncontiguous in memory as well as in file */
            k = num = 0;
            i_offset = buf_view.off[0];
            j = st_index;
            off = offset;
            fwr_size = st_fwr_size;
            bwr_size = buf_view.len[0];

            while (num < bufsize) {
                size = MIN(fwr_size, bwr_size);
                if (size) {
                    req_off = off;
                    req_len = size;
                    userbuf_off = i_offset;
                    BUFFERED_WRITE;
                }

                num += size;
                if (num >= bufsize) break;

                new_fwr_size = fwr_size;
                new_bwr_size = bwr_size;

                if (size == fwr_size) {
                    j++;
#ifdef PNETCDF_DEBUG
assert(j < fd->file_view.count);
#endif
                    off = fd->file_view.off[j];
                    new_fwr_size = fd->file_view.len[j];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
                    /* reached end of contiguous block in memory */

                    k++;
#ifdef PNETCDF_DEBUG
assert(k < buf_view.count);
#endif
                    i_offset = buf_view.off[k];
                    new_bwr_size = buf_view.len[k];
                    if (size != fwr_size) {
                        off += size;
                        new_fwr_size -= size;
                    }
                }
                fwr_size = new_fwr_size;
                bwr_size = new_bwr_size;
            }
        }

        /* write the buffer out finally */
        if (writebuf_len) {
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);
            if (w_len < 0) goto fn_exit;
            total_w_len += w_len;
        }
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
    }

fn_exit:
    if (writebuf != NULL)
        NCI_Free(writebuf);

    return total_w_len;
}
