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
                                  writebuf_off);                        \
                if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                    PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
                if (w_len < 0) goto fn_exit;                            \
                total_w_len += w_len;                                   \
            }                                                           \
            writebuf_off = req_off;                                     \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len,        \
                                     writebuf_off);                     \
            if (r_len < 0) goto fn_exit;                                \
        }                                                               \
        write_sz = (MPI_Aint) (MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        assert((MPI_Offset)write_sz == MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        memcpy(writebuf+req_off-writebuf_off, (char *)buf +userbuf_off, write_sz); \
        while (write_sz != req_len) {                                   \
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,       \
                                      writebuf_off);                    \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_UNLOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            if (w_len < 0) goto fn_exit;                                \
            total_w_len += w_len;                                       \
            req_len -= write_sz;                                        \
            userbuf_off += write_sz;                                    \
            writebuf_off += writebuf_len;                               \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
            if (!fd->atomicity && fd->hints->ds_write == PNCIO_HINT_DISABLE) \
                PNCIO_WRITE_LOCK(fd, writebuf_off, SEEK_SET, writebuf_len); \
            r_len = PNCIO_ReadContig(fd, writebuf, writebuf_len,        \
                                     writebuf_off);                     \
            if (r_len < 0) goto fn_exit;                                \
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
            if (w_len < 0) goto fn_exit;                                \
            total_w_len += w_len;                                       \
            writebuf_off = req_off;                                     \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
        }                                                               \
        write_sz = (MPI_Aint) (MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        assert((MPI_Offset)write_sz == MPL_MIN(req_len, writebuf_off + writebuf_len - req_off)); \
        memcpy(writebuf+req_off-writebuf_off, (char *)buf +userbuf_off, write_sz); \
        while (write_sz != req_len) {                                   \
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len,       \
                             writebuf_off);                             \
            if (w_len < 0) goto fn_exit;                                \
            total_w_len += w_len;                                       \
            req_len -= write_sz;                                        \
            userbuf_off += write_sz;                                    \
            writebuf_off += writebuf_len;                               \
            writebuf_len = MPL_MIN(max_bufsize,end_offset-writebuf_off+1); \
            write_sz = MPL_MIN(req_len, writebuf_len);                  \
            memcpy(writebuf, (char *)buf + userbuf_off, write_sz);      \
        }                                                               \
    }

MPI_Offset PNCIO_GEN_WriteStrided(PNCIO_File *fd,
                                  const void *buf,
                                  PNCIO_Flat_list buf_view,
                                  MPI_Offset offset)
{

/* offset is in units of etype relative to the filetype. */

    PNCIO_Flatlist_node *flat_buf;
    MPI_Offset i_offset, sum, size_in_filetype;
    int i, j, k, st_index = 0;
    MPI_Offset num, size, n_filetypes, etype_in_filetype, st_n_filetypes;
    MPI_Offset n_etypes_in_filetype, abs_off_in_filetype = 0;
    MPI_Count filetype_size;
    MPI_Aint lb, filetype_extent;
    int buf_count, filetype_is_contig;
    MPI_Offset userbuf_off;
    MPI_Offset off, req_off, disp, end_offset = 0, writebuf_off, start_off;
    char *writebuf = NULL;
    MPI_Aint writebuf_len, max_bufsize, write_sz;
    MPI_Aint bufsize;
    MPI_Offset new_bwr_size, new_fwr_size, st_fwr_size, fwr_size = 0, bwr_size, req_len;
    MPI_Offset r_len, w_len, total_w_len=0;

    if (fd->hints->ds_write == PNCIO_HINT_DISABLE) {
        /* if user has disabled data sieving on reads, use naive
         * approach instead.
         */
        return PNCIO_GEN_WriteStrided_naive(fd, buf, buf_view, offset);
    }

assert(fd->filetype == MPI_BYTE);
assert(fd->flat_file.size == buf_view.size);

    filetype_size = fd->flat_file.size;
    filetype_is_contig = (fd->flat_file.count <= 1);

// printf("%s at %d: offset=%lld\n",__func__,__LINE__, offset);
if (fd->flat_file.count > 0) assert(offset == 0); /* not whole file visible */

// printf("%s at %d: offset=%lld filetype_size=%lld\n",__func__,__LINE__, offset,filetype_size);

// TODO: remove use of filetype_extent
    if (fd->flat_file.count == 0)
        filetype_extent = filetype_size;
    else
        filetype_extent = fd->flat_file.indices[fd->flat_file.count-1]
                        + fd->flat_file.blocklens[fd->flat_file.count-1]
                        - fd->flat_file.indices[0];
#if 0
    else if (fd->filetype == MPI_BYTE) {
        filetype_is_contig = 1;
        filetype_size = 1;
        filetype_extent = 1;
    }
    else {
        // PNCIO_Datatype_iscontig(fd->filetype, &filetype_is_contig);
        filetype_is_contig = (fd->flat_file.count <= 1);
        MPI_Type_size_x(fd->filetype, &filetype_size);
        if (filetype_size == 0)
            return NC_NOERR;
        MPI_Type_get_extent(fd->filetype, &lb, &filetype_extent);
    }
#endif

    bufsize = buf_view.size;

/* get max_bufsize from the info object. */

    max_bufsize = fd->hints->ind_wr_buffer_size;

    /* Contiguous both in buftype and filetype should have been handled in a
     * call to PNCIO_WriteContig() earlier.
     */
    assert(!(buf_view.is_contig && filetype_is_contig));

PNCIO_Flatlist_node tmp_buf;
    flat_buf = &tmp_buf;
    flat_buf->count = buf_view.count;
    flat_buf->indices = buf_view.off;
    flat_buf->blocklens = buf_view.len;

    if (!buf_view.is_contig && filetype_is_contig) {

/* noncontiguous in memory, contiguous in file. */

        off = fd->disp + offset;
assert(fd->disp == 0);

        start_off = off;
        end_offset = off + bufsize - 1;
        writebuf_off = off;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        writebuf_len = MPL_MIN(max_bufsize, end_offset - writebuf_off + 1);

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        for (i = 0; i < flat_buf->count; i++) {
            userbuf_off = flat_buf->indices[i];
            req_off = off;
            req_len = flat_buf->blocklens[i];
            BUFFERED_WRITE_WITHOUT_READ;
            off += flat_buf->blocklens[i];
        }

        /* write the buffer out finally */
        if (writebuf_len) {
            w_len = PNCIO_WriteContig(fd, writebuf, writebuf_len, writebuf_off);
            if (w_len >= 0) total_w_len += w_len;
        }
        else
            w_len = 0;

        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        if (w_len < 0)
            goto fn_exit;
    }

    else { /* noncontiguous in file */

        disp = fd->disp;
assert(fd->disp == 0);

        n_etypes_in_filetype = filetype_size;
        n_filetypes = offset / n_etypes_in_filetype;
        etype_in_filetype = offset % n_etypes_in_filetype;
        size_in_filetype = etype_in_filetype;

        sum = 0;
        for (i = 0; i < fd->flat_file.count; i++) {
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
        offset = disp + (MPI_Offset) n_filetypes *filetype_extent + abs_off_in_filetype;

        start_off = offset;
assert(offset == abs_off_in_filetype);

// printf("%s at %d: start_off=%lld n_filetypes=%lld abs_off_in_filetype=%lld\n",__func__,__LINE__,start_off,n_filetypes,abs_off_in_filetype);

        /* Wei-keng Liao:write request is within single flat_file contig block */
        /* this could happen, for example, with subarray types that are
         * actually fairly contiguous */
        if (buf_view.is_contig && bufsize <= fwr_size) {
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
        st_n_filetypes = n_filetypes;
        i_offset = 0;
        j = st_index;
        off = offset;
        fwr_size = MPL_MIN(st_fwr_size, bufsize);
        while (i_offset < bufsize) {
            i_offset += fwr_size;
            end_offset = off + fwr_size - 1;

            j = (j + 1) % fd->flat_file.count;
            n_filetypes += (j == 0) ? 1 : 0;
            while (fd->flat_file.blocklens[j] == 0) {
                j = (j + 1) % fd->flat_file.count;
                n_filetypes += (j == 0) ? 1 : 0;
            }

            off = disp + fd->flat_file.indices[j] + n_filetypes * (MPI_Offset) filetype_extent;
            fwr_size = MPL_MIN(fd->flat_file.blocklens[j], bufsize - i_offset);
        }

        /* if atomicity is true or data sieving is not disable, lock the region
         * to be accessed */
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_WRITE_LOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);

        writebuf_off = 0;
        writebuf_len = 0;
        writebuf = (char *) NCI_Malloc(max_bufsize);
        memset(writebuf, -1, max_bufsize);

        if (buf_view.is_contig && !filetype_is_contig) {

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

                if (off + fwr_size < disp + fd->flat_file.indices[j] +
                    fd->flat_file.blocklens[j] + n_filetypes * (MPI_Offset) filetype_extent)
                    off += fwr_size;
                /* did not reach end of contiguous block in filetype.
                 * no more I/O needed. off is incremented by fwr_size. */
                else {
                    j = (j + 1) % fd->flat_file.count;
                    n_filetypes += (j == 0) ? 1 : 0;
                    while (fd->flat_file.blocklens[j] == 0) {
                        j = (j + 1) % fd->flat_file.count;
                        n_filetypes += (j == 0) ? 1 : 0;
                    }
                    off = disp + fd->flat_file.indices[j] +
                        n_filetypes * (MPI_Offset) filetype_extent;
                    fwr_size = MPL_MIN(fd->flat_file.blocklens[j], bufsize - i_offset);
                }
            }
        } else {
/* noncontiguous in memory as well as in file */

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
                    j = (j + 1) % fd->flat_file.count;
                    n_filetypes += (j == 0) ? 1 : 0;
                    while (fd->flat_file.blocklens[j] == 0) {
                        j = (j + 1) % fd->flat_file.count;
                        n_filetypes += (j == 0) ? 1 : 0;
                    }

                    off = disp + fd->flat_file.indices[j] +
                        n_filetypes * (MPI_Offset) filetype_extent;

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
            if (w_len < 0) goto fn_exit;
            total_w_len += w_len;
        }
        if (fd->atomicity || fd->hints->ds_write != PNCIO_HINT_DISABLE)
            PNCIO_UNLOCK(fd, start_off, SEEK_SET, end_offset - start_off + 1);
    }

fn_exit:
    if (writebuf != NULL)
        NCI_Free(writebuf);

    return total_w_len;
}
