/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

/*----< PNCIO_GEN_WriteStrided_nods() >--------------------------------------*/
/* This subroutine implements independent write when data sieving is disabled.
 * Note in PnetCDF, the file_view and buf_view are never used for more than
 * one round, which greatly simplifies the implementation.
 */
MPI_Offset PNCIO_GEN_WriteStrided_nods(PNCIO_File *fd,
                                       const void *buf,
                                       PNCIO_View  buf_view)
{
    char *ptr, *cpy_ptr, *tmp_buf;
    MPI_Count i, j, k, ntimes;
    MPI_Offset lock_off, lock_len, len, total_len=0;

    MPI_Offset tmp_buf_size, file_off, zero=0;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset file_rem, buf_rem, buf_size;
#else
    int file_rem, buf_rem, buf_size;
#endif

// printf("%s at %d: buf_view count %lld file_view count %lld\n",__func__,__LINE__, buf_view.count,fd->file_view.count);

#ifdef PNETCDF_DEBUG
    /* fd->file_view.off and fd->file_view.len should not be NULL */
    assert(fd->file_view.count > 0);
    assert(fd->file_view.off != NULL);
    assert(fd->file_view.len != NULL);

    if (fd->file_view.count == 1)
        assert(fd->file_view.size == fd->file_view.len[0]);

    if (buf_view.count == 1)
        assert(buf_view.size == buf_view.len[0]);

    /* Contiguous both in buf_view and file_view should have already been
     * handled earlier in a call to PNCIO_WriteContig().
     */
    assert(!(buf_view.count <= 1 && fd->file_view.count <= 1));

    /* In PnetCDF, fd->file_view.size always == buf_view.size. */
    assert(fd->file_view.size == buf_view.size);
#endif

    if (fd->file_view.size == 0) /* zero-sized request */
        return 0;

    if (buf_view.count == 0) {
        /* In this case, the buffer view is contiguous, adjust buf_view.count,
         * buf_view.off, and buf_view.len, so they can be used in the while
         * loop below.
         */
        buf_size = buf_view.size;
        buf_view.count = 1;
        buf_view.off = &zero;
        buf_view.len = &buf_size;
    }

    lock_off = fd->file_view.off[0];
    if (fd->file_view.count > 1)
        lock_len = fd->file_view.off[fd->file_view.count-1]
                 + fd->file_view.len[fd->file_view.count-1]
                 - lock_off;
    else
        lock_len = fd->file_view.size;

    /* if atomicity is true, lock (exclusive) the region to be accessed */
    if (fd->atomicity)
        PNCIO_WRITE_LOCK(fd, lock_off, SEEK_SET, lock_len);

    if (buf_view.count <= 1) { /* use buf directly to write */
        tmp_buf = (char*)buf;
        buf_rem = buf_view.size;
        ntimes = 1;
        tmp_buf_size = buf_view.size;
    }
    else { /* buf is noncontiguous */
        tmp_buf_size = MIN(buf_view.size, fd->hints->ind_wr_buffer_size);
        tmp_buf = (char*) NCI_Malloc(tmp_buf_size);
        buf_rem = buf_view.len[0];
        ntimes = buf_view.size / tmp_buf_size;
        if (buf_view.size % tmp_buf_size)
            ntimes++;
    }

    file_off = fd->file_view.off[0];
    file_rem = fd->file_view.len[0];

    /* pointer to buf, starting location to copy over to tmp_buf */
    cpy_ptr = (char*)buf;

    k = (buf_view.count <= 1) ? 1 : 0; /* whether to skip while loop k */
    j = 0;
    for (i=0; i<ntimes; i++) { /* perform write in ntimes rounds */
        MPI_Offset req_len, tmp_buf_rem;

        /* copy data from buf to tmp_buf */
        tmp_buf_rem = tmp_buf_size;
        ptr = tmp_buf;
        while (k < buf_view.count) {
            req_len = MIN(tmp_buf_rem, buf_rem);
            memcpy(ptr, cpy_ptr, req_len);

            ptr += req_len;
            tmp_buf_rem -= req_len;
            if (tmp_buf_rem == 0) break;

            if (buf_rem == req_len) { /* done with pair k */
                k++;
                cpy_ptr = (char*)buf + buf_view.off[k];
                buf_rem = buf_view.len[k];
            }
            else { /* there is still data remained in pair k */
                cpy_ptr += req_len;
                buf_rem -= req_len;
            }
        }

        /* using tmp_buf to write to the file */
        tmp_buf_rem = tmp_buf_size;
        ptr = tmp_buf;
        while (j < fd->file_view.count) {
            req_len = MIN(tmp_buf_rem, file_rem);
            /* write from offset file_off of length req_len */
            len = PNCIO_WriteContig(fd, ptr, req_len, file_off);
            if (len < 0) return len;
            total_len += len;

            ptr += req_len;
            tmp_buf_rem -= req_len;
            if (tmp_buf_rem == 0) break;

            if (file_rem == req_len) { /* done with pair j */
                j++;
                file_off = fd->file_view.off[j];
                file_rem = fd->file_view.len[j];
            }
            else { /* there is still data remained in pair j */
                file_off += req_len;
                file_rem -= req_len;
            }
        }
    }

    /* free tmp_buf if allocated */
    if (tmp_buf != buf) NCI_Free(tmp_buf);

    if (fd->atomicity)
        PNCIO_UNLOCK(fd, lock_off, SEEK_SET, lock_len);

    return total_len;
}

