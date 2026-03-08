/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

/*----< PNCIO_GEN_ReadStrided_nods() >---------------------------------------*/
/* This subroutine implements independent read when data sieving is disabled.
 * Note in PnetCDF, the file_view and buf_view are never used for more than
 * one round, which greatly simplifies the implementation.
 */
MPI_Offset PNCIO_GEN_ReadStrided_nods(PNCIO_File *fd,
                                      void       *buf,
                                      PNCIO_View  buf_view)
{
    char *ptr;
    MPI_Count j, k;
    MPI_Offset lock_off, lock_len, len, total_len=0;

    MPI_Offset file_off, zero=0;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset file_rem, buf_rem, buf_size;
#else
    int file_rem, buf_rem, buf_size;
#endif

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
     * handled earlier in a call to PNCIO_ReadContig().
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
    if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
        PNCIO_WRITE_LOCK(fd, lock_off, SEEK_SET, lock_len);

    file_off = fd->file_view.off[0];
    file_rem = fd->file_view.len[0];
    buf_rem  = buf_view.len[0];
    ptr = (char*)buf + buf_view.off[0];

    j = 0;
    k = 0;
    while (1) {
        /* whichever is shorter */
        MPI_Offset req_len = MIN(file_rem, buf_rem);

        /* read from offset file_off of length req_len */
        len = PNCIO_ReadContig(fd, ptr, req_len, file_off);
        if (len < 0) return len;
        total_len += len;

        if (req_len == file_rem) { /* done with pair j */
            j++;
            if (j == fd->file_view.count) break;
            file_off = fd->file_view.off[j];
            file_rem = fd->file_view.len[j];
        }
        else { /* req_len < file_rem, remains in pair j */
            file_off += req_len;
            file_rem -= req_len;
        }

        if (req_len == buf_rem) { /* done with pair k */
            k++;
            if (k == buf_view.count) break;
            buf_rem = buf_view.len[k];
            ptr = (char*)buf + buf_view.off[k];
        }
        else { /* req_len < buf_rem, remains in pair k */
            buf_rem -= req_len;
            ptr += req_len;
        }
    }

    if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
        PNCIO_UNLOCK(fd, lock_off, SEEK_SET, lock_len);

    return total_len;
}

