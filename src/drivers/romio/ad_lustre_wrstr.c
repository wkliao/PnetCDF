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
                w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,   \
                                          writebuf_off);                \
                if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                    PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
                if (w_len < 0) {                                        \
                    NCI_Free(writebuf);                                 \
                    return w_len;                                       \
                }                                                       \
                total_w_len += w_len;                                   \
                writebuf_off = req_off;                                 \
            }                                                           \
            writebuf_off = req_off;                                     \
            /* stripe_size alignment */                                 \
            writebuf_len = MPL_MIN(end_offset - writebuf_off + 1,       \
                                   (writebuf_off / stripe_size + 1) *   \
                                   stripe_size - writebuf_off);         \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len,        \
                                     writebuf_off);                     \
            if (r_len < 0) {                                            \
                NCI_Free(writebuf);                                     \
                return r_len;                                           \
            }                                                           \
        }                                                               \
        write_sz = (MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        memcpy(writebuf + req_off - writebuf_off, (char *)buf +userbuf_off, write_sz); \
        while (write_sz != req_len) {                                   \
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,       \
                                      writebuf_off);                    \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            if (w_len < 0) {                                            \
                NCI_Free(writebuf);                                     \
                return w_len;                                           \
            }                                                           \
            total_w_len += w_len;                                       \
            req_len -= write_sz;                                        \
            userbuf_off += write_sz;                                    \
            writebuf_off += writebuf_len;                               \
            /* stripe_size alignment */                                 \
            writebuf_len = MPL_MIN(end_offset - writebuf_off + 1,       \
                                   (writebuf_off / stripe_size + 1) *   \
                                   stripe_size - writebuf_off);         \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len,        \
                                     writebuf_off);                     \
            if (r_len < 0) {                                            \
                NCI_Free(writebuf);                                     \
                return r_len;                                           \
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
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,       \
                                      writebuf_off);                    \
            if (w_len < 0) {                                            \
                NCI_Free(writebuf);                                     \
                return w_len;                                           \
            }                                                           \
            total_w_len += w_len;                                       \
            writebuf_off = req_off;                                     \
            /* stripe_size alignment */                                 \
            writebuf_len = MPL_MIN(end_offset - writebuf_off + 1,       \
                                   (writebuf_off / stripe_size + 1) *   \
                                   stripe_size - writebuf_off);         \
        }                                                               \
        write_sz = MPL_MIN(req_len, writebuf_off + writebuf_len - req_off); \
        memcpy(writebuf + req_off - writebuf_off,                       \
               (char *)buf + userbuf_off, write_sz);                    \
        while (write_sz != req_len) {                                   \
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,       \
                                      writebuf_off);                    \
            if (w_len < 0) {                                            \
                NCI_Free(writebuf);                                     \
                return w_len;                                           \
            }                                                           \
            total_w_len += w_len;                                       \
            req_len -= write_sz;                                        \
            userbuf_off += write_sz;                                    \
            writebuf_off += writebuf_len;                               \
            /* stripe_size alignment */                                 \
            writebuf_len = MPL_MIN(end_offset - writebuf_off + 1,       \
                                   (writebuf_off / stripe_size + 1) *   \
                                   stripe_size - writebuf_off);         \
            write_sz = MPL_MIN(req_len, writebuf_len);                  \
            memcpy(writebuf, (char *)buf + userbuf_off, write_sz);      \
        }                                                               \
    }

MPI_Offset PNCIO_LUSTRE_WriteStrided(PNCIO_File *fd,
                                     const void *buf,
                                     PNCIO_Flat_list buf_view,
                                     MPI_Offset offset)
{
    /* offset is in units of etype relative to the filetype. */
    PNCIO_Flatlist_node *flat_buf;
    MPI_Offset i_offset, sum;
    int i, j, k, st_index = 0;
    MPI_Offset num, size;
    MPI_Offset abs_off_in_filetype = 0;
    int buf_count, filetype_is_contig;
    MPI_Offset userbuf_off;
    MPI_Offset off, req_off, disp, end_offset = 0, writebuf_off, start_off;
    char *writebuf;
    MPI_Count bufsize, writebuf_len, write_sz;
    MPI_Offset new_bwr_size, new_fwr_size, st_fwr_size, fwr_size = 0, bwr_size, req_len;
    int stripe_size;
    MPI_Offset r_len, w_len, total_w_len=0;

    /* The case of both buftype and filetype being contiguous has gone to
     * PNCIO_WriteContig().
     */

// printf("%s at %d:\n",__func__,__LINE__);

    if (fd->hints->ds_write == PNCIO_HINT_DISABLE) {
        /* if user has disabled data sieving on writes, use naive
         * approach instead.
         */
        return PNCIO_GEN_WriteStrided_naive(fd, buf, buf_view, offset);
    }


assert(fd->filetype == MPI_BYTE);
assert(fd->flat_file.size == buf_view.size);

    filetype_is_contig = (fd->flat_file.count <= 1);

if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

    bufsize = buf_view.size;

    /* get striping info */
    stripe_size = fd->hints->striping_unit;

// printf("%s at %d: filetype_is_contig=%d buf_view.is_contig=%d offset=%lld\n",__func__,__LINE__, filetype_is_contig,buf_view.is_contig,offset);


PNCIO_Flatlist_node tmp_buf;
    flat_buf = &tmp_buf;
    flat_buf->count = buf_view.count;
    flat_buf->indices = buf_view.off;
    flat_buf->blocklens = buf_view.len;

    /* Different buftype to different filetype */
    if (!buf_view.is_contig && filetype_is_contig) {
        /* noncontiguous in memory, contiguous in file. */


        off = fd->disp + offset;

        start_off = off;
        end_offset = start_off + bufsize - 1;
        /* write stripe size buffer each time */
        writebuf = (char *) NCI_Malloc(MPL_MIN(bufsize, stripe_size));
        writebuf_off = 0;
        writebuf_len = 0;

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, bufsize);

        for (i = 0; i < flat_buf->count; i++) {
            userbuf_off = flat_buf->indices[i];
            req_off = off;
            req_len = flat_buf->blocklens[i];
            BUFFERED_WRITE_WITHOUT_READ;
            off += flat_buf->blocklens[i];
        }

        /* write the buffer out finally */
        w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);

        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, bufsize);

        NCI_Free(writebuf);

        if (w_len < 0) return w_len;
        total_w_len += w_len;

    } else { /* contiguous buffer and non-contiguous in file */
        MPI_Offset size_in_filetype = offset;

        disp = fd->disp;

        sum = 0;
        for (i = 0; i < fd->flat_file.count; i++) {
// printf("%s at %d: disp=%lld flat_file i=%d indices=%lld blocklens=%lld\n",__func__,__LINE__, disp,i,fd->flat_file.indices[0],fd->flat_file.blocklens[0]);
            sum += fd->flat_file.blocklens[i];
            if (sum > size_in_filetype) {
                st_index = i;
                fwr_size = sum - size_in_filetype;
                abs_off_in_filetype = fd->flat_file.indices[i] +
                    size_in_filetype - (sum - fd->flat_file.blocklens[i]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        offset = disp + abs_off_in_filetype;

        start_off = offset;

        /* Wei-keng Liao:write request is within single flat_file
         * contig block*/
        /* this could happen, for example, with subarray types that are
         * actually fairly contiguous */
        if (buf_view.is_contig && bufsize <= fwr_size) {
            req_off = start_off;
            req_len = bufsize;
            end_offset = start_off + bufsize - 1;
            writebuf = (char *) NCI_Malloc(MPL_MIN(bufsize, stripe_size));
            memset(writebuf, -1, MPL_MIN(bufsize, stripe_size));
            writebuf_off = 0;
            writebuf_len = 0;
            userbuf_off = 0;
            BUFFERED_WRITE_WITHOUT_READ;

            /* write the buffer out finally */
            if (fd->hints->ds_write != PNCIO_HINT_DISABLE)
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len);

            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (w_len > 0) total_w_len += w_len;

            if (fd->hints->ds_write != PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);

            NCI_Free(writebuf);

            return total_w_len;
        }

        /* Calculate end_offset, the last byte-offset that will be accessed.
         * e.g., if start_offset=0 and 100 bytes to be write, end_offset=99 */

        st_fwr_size = fwr_size;
        i_offset = 0;
        j = st_index;
        off = offset;
        fwr_size = MPL_MIN(st_fwr_size, bufsize);
        while (i_offset < bufsize) {
            i_offset += fwr_size;
            end_offset = off + fwr_size - 1;

            j = (j + 1) % fd->flat_file.count;
            while (fd->flat_file.blocklens[j] == 0) {
                j = (j + 1) % fd->flat_file.count;
            }

            off = disp + fd->flat_file.indices[j];
            fwr_size = MPL_MIN(fd->flat_file.blocklens[j], bufsize - i_offset);
        }

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) NCI_Malloc(stripe_size);
        memset(writebuf, -1, stripe_size);

        if (buf_view.is_contig && !filetype_is_contig) {

/* contiguous in memory, noncontiguous in file. should be the most
           common case. */

            i_offset = 0;
            j = st_index;
            off = offset;
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

                if (off + fwr_size < disp + fd->flat_file.indices[j] +
                    fd->flat_file.blocklens[j])
                    off += fwr_size;
                /* did not reach end of contiguous block in filetype.
                 * no more I/O needed. off is incremented by fwr_size. */
                else {
                    j = (j + 1) % fd->flat_file.count;
                    while (fd->flat_file.blocklens[j] == 0) {
                        j = (j + 1) % fd->flat_file.count;
                    }
                    off = disp + fd->flat_file.indices[j];
                    fwr_size = MPL_MIN(fd->flat_file.blocklens[j], bufsize - i_offset);
                }
            }
        } else {
/* noncontiguous in memory as well as in file */
            // flat_buf = PNCIO_Flatten_and_find(buftype);

            k = num = buf_count = 0;
            i_offset = flat_buf->indices[0];
            j = st_index;
            off = offset;
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
                    j = (j + 1) % fd->flat_file.count;
                    while (fd->flat_file.blocklens[j] == 0) {
                        j = (j + 1) % fd->flat_file.count;
                    }

                    off = disp + fd->flat_file.indices[j];

                    new_fwr_size = fd->flat_file.blocklens[j];
                    if (size != bwr_size) {
                        i_offset += size;
                        new_bwr_size -= size;
                    }
                }

                if (size == bwr_size) {
/* reached end of contiguous block in memory */

                    k = (k + 1) % flat_buf->count;
                    buf_count++;
                    i_offset = flat_buf->indices[k];
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
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE)
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len);
            if (w_len < 0) return w_len;
            total_w_len += w_len;
        }
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        NCI_Free(writebuf);
    }

    return buf_view.size;
}
