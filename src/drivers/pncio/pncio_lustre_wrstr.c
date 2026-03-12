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
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
            if (w_len < 0) {                                                  \
                NCI_Free(writebuf);                                           \
                return w_len;                                                 \
            }                                                                 \
            total_w_len += w_len;                                             \
            writebuf_off = req_off;                                           \
        }                                                                     \
        writebuf_off = req_off;                                               \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return r_len;                                                     \
        }                                                                     \
    }                                                                         \
    write_sz = (MIN(req_len, writebuf_off + writebuf_len - req_off));         \
    memcpy(writebuf + req_off - writebuf_off, (char *)buf + userbuf_off,      \
           write_sz);                                                         \
    while (write_sz != req_len) {                                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);           \
        if (w_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return w_len;                                                     \
        }                                                                     \
        total_w_len += w_len;                                                 \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
        if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)      \
            PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);       \
        r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len, writebuf_off);   \
        if (r_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return r_len;                                                     \
        }                                                                     \
        write_sz = MIN(req_len, writebuf_len);                                \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}

/* This macro is used when file_view is contiguous and buf_view is not
 * contiguous. It does not do a read-modify-write and does not lock.
 */
#define BUFFERED_WRITE_WITHOUT_READ {                                         \
    if (req_off >= writebuf_off + writebuf_len) {                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (w_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return w_len;                                                     \
        }                                                                     \
        total_w_len += w_len;                                                 \
        writebuf_off = req_off;                                               \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
    }                                                                         \
    write_sz = MIN(req_len, writebuf_off + writebuf_len - req_off);           \
    memcpy(writebuf + req_off - writebuf_off,                                 \
           (char *)buf + userbuf_off, write_sz);                              \
    while (write_sz != req_len) {                                             \
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);  \
        if (w_len < 0) {                                                      \
            NCI_Free(writebuf);                                               \
            return w_len;                                                     \
        }                                                                     \
        total_w_len += w_len;                                                 \
        req_len -= write_sz;                                                  \
        userbuf_off += write_sz;                                              \
        writebuf_off += writebuf_len;                                         \
        /* stripe_size alignment */                                           \
        writebuf_len = MIN(end_offset - writebuf_off + 1,                     \
                           (writebuf_off / stripe_size + 1) * stripe_size     \
                           - writebuf_off);                                   \
        write_sz = MIN(req_len, writebuf_len);                                \
        memcpy(writebuf, (char *)buf + userbuf_off, write_sz);                \
    }                                                                         \
}

MPI_Offset PNCIO_LUSTRE_WriteStrided(PNCIO_File *fd,
                                     const void *buf,
                                     PNCIO_View  buf_view)
{
    char *writebuf;
    int i, j, k, st_index=0, stripe_size;
    MPI_Offset i_offset, num, size, off;
    MPI_Offset userbuf_off, req_off, end_offset=0, writebuf_off, start_off;
    MPI_Offset new_bwr_size, new_fwr_size, st_fwr_size, fwr_size=0, bwr_size;
    MPI_Offset req_len, r_len, w_len, total_w_len=0;
    MPI_Count bufsize, writebuf_len, write_sz;

    /* The case of both buf_view and file_view being contiguous has gone to
     * PNCIO_WriteContig().
     */

// printf("%s at %d:\n",__func__,__LINE__);

#ifdef PNETCDF_DEBUG
    /* When data sieveing is disabled, PNCIO_GEN_Write_indep() should be called
     * instead. */
    assert(fd->hints->romio_ds_write != PNCIO_HINT_DISABLE);

    assert(fd->file_view.size == buf_view.size);
#endif

    bufsize = buf_view.size;

    /* get striping info */
    stripe_size = fd->hints->striping_unit;

    if (buf_view.count > 1 && fd->file_view.count <= 1) {
        /* noncontiguous in write buffer, contiguous in file. */

        off = fd->file_view.off[0];

        start_off = off;
        end_offset = start_off + bufsize - 1;

        /* write stripe size buffer each time */
        writebuf = (char *) NCI_Malloc(MIN(bufsize, stripe_size));
        writebuf_off = 0;
        writebuf_len = 0;

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed
         */
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, bufsize);

        for (i = 0; i < buf_view.count; i++) {
            userbuf_off = buf_view.off[i];
            req_off = off;
            req_len = buf_view.len[i];
            BUFFERED_WRITE_WITHOUT_READ;
            off += buf_view.len[i];
        }

        /* write the buffer out the last round */
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);

        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, bufsize);

        NCI_Free(writebuf);

        if (w_len < 0) return w_len;
        total_w_len += w_len;

    } else { /* contiguous buffer and non-contiguous in file */
        /* find the starting index in fd->file_view offset-length pairs */
        MPI_Offset offset = 0;
        MPI_Offset sum = 0;

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
            req_off = start_off;
            req_len = bufsize;
            end_offset = start_off + bufsize - 1;
            writebuf = (char *) NCI_Malloc(MIN(bufsize, stripe_size));
            memset(writebuf, -1, (size_t)MIN(bufsize, stripe_size));
            writebuf_off = 0;
            writebuf_len = 0;
            userbuf_off = 0;
            BUFFERED_WRITE_WITHOUT_READ;

            /* write the buffer out the last round */
            if (fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);

            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (w_len > 0) total_w_len += w_len;

            if (fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);

            NCI_Free(writebuf);

            return total_w_len;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be write, end_offset=99 */

        st_fwr_size = fwr_size;
        j = st_index;
        i_offset = fwr_size = MIN(st_fwr_size, bufsize);
        end_offset = offset + fwr_size - 1;
        while (i_offset < bufsize) {
            j++;
assert(j < fd->file_view.count);
            off = fd->file_view.off[j];
            fwr_size = MIN(fd->file_view.len[j], bufsize - i_offset);
            i_offset += fwr_size;
            end_offset = off + fwr_size - 1;
        }

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset-start_off+1);

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) NCI_Malloc(stripe_size);
        memset(writebuf, -1, stripe_size);

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

                if (off + fwr_size < fd->file_view.off[j] + fd->file_view.len[j])
                    off += fwr_size;
                    /* no more I/O needed. off is incremented by fwr_size. */
                else {
                    j++;
assert(j < fd->file_view.count);
                    off = fd->file_view.off[j];
                    fwr_size = MIN(fd->file_view.len[j], bufsize - i_offset);
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
                    /* reached end of contiguous block in file */
                    j++;
assert(j < fd->file_view.count);
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
assert(k < buf_view.count);
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

        /* write the buffer out the last round */
        if (writebuf_len) {
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (!fd->atomicity && fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);
            if (w_len < 0) return w_len;
            total_w_len += w_len;
        }
        if (fd->atomicity || fd->hints->romio_ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        NCI_Free(writebuf);
    }

    return buf_view.size;
}
