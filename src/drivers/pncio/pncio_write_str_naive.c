/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <pncio.h>

/*----< PNCIO_GEN_WriteStrided_naive() >-------------------------------------*/
/* This subroutine implements independent write when data sieving is disabled.
 */
MPI_Offset PNCIO_GEN_WriteStrided_naive(PNCIO_File *fd,
                                        const void *buf,
                                        PNCIO_View  buf_view)
{
    MPI_Count j, k;
    MPI_Offset lock_off, lock_len, w_len, total_w_len=0;

#ifdef PNETCDF_DEBUG
    /* Contiguous both in buf_view and file_view should have already been
     * handled earlier in a call to PNCIO_WriteContig().
     */
    assert(!(buf_view.count <= 1 && fd->file_view.count <= 1));

    /* In PnetCDF, fd->file_view.size always == buf_view.size, i.e.
     * file_view and buf_view are never used for more than one round.
     */
    assert(fd->file_view.size == buf_view.size);
#endif

    lock_off = fd->file_view.off[0];
    lock_len = fd->file_view.size;
    if (fd->file_view.count > 1)
        lock_len = fd->file_view.off[fd->file_view.count-1]
                 + fd->file_view.len[fd->file_view.count-1];

    /* if atomicity is true, lock (exclusive) the region to be accessed */
    if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
        PNCIO_WRITE_LOCK(fd, lock_off, SEEK_SET, lock_len);

    if (buf_view.count > 1 && fd->file_view.count <= 1) {
        /* noncontiguous buffer view, contiguous file view */

        MPI_Offset off = fd->file_view.off[0];
        for (j=0; j<buf_view.count; j++) {
            /* write one buf_view's offset-length pair at a time */
            w_len = PNCIO_WriteContig(fd, (char*)buf + buf_view.off[j],
                                      buf_view.len[j], off);
            if (w_len < 0) return w_len;
            total_w_len += w_len;
            off += buf_view.len[j];
        }
    }
    else { /* noncontiguous file view */
        char *ptr;

        if (buf_view.count <= 1 && fd->file_view.count > 1) {
            /* contiguous buffer view, noncontiguous file view */

            ptr = (char*)buf;
            for (j=0; j<fd->file_view.count; j++) {
                /* write one file_view's offset-length pair at a time */
                w_len = PNCIO_WriteContig(fd, ptr, fd->file_view.len[j],
                                          fd->file_view.off[j]);
                if (w_len < 0) return w_len;
                total_w_len += w_len;
                ptr += fd->file_view.len[j];
            }
        }
        else {
            /* Both buffer view and file view are noncontiguous. */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Offset file_off, buf_off;
            MPI_Offset file_rem, buf_rem;
#else
            MPI_Offset file_off, buf_off;
            int        file_rem, buf_rem;
#endif

            file_off = fd->file_view.off[0];
            file_rem = fd->file_view.len[0];
            buf_off  = buf_view.off[0];
            buf_rem  = buf_view.len[0];
            ptr = (char*)buf + buf_off;

            j = 0;
            k = 0;
            while (j < fd->file_view.count) {

                while (k < buf_view.count) {

                    /* whichever is shorter */
                    MPI_Offset req_len = MIN(file_rem, buf_rem);

                    /* write from offset file_off for length of req_len */
                    w_len = PNCIO_WriteContig(fd, ptr, req_len, file_off);
                    if (w_len < 0) return w_len;
                    total_w_len += w_len;

                    if (req_len == file_rem) { /* done with pair j */
                        j++;
                        file_off = fd->file_view.off[j];
                        file_rem = fd->file_view.len[j];
                    }
                    else { /* req_len < file_rem, remains in pair j */
                        file_off += req_len;
                        file_rem -= req_len;
                    }

                    if (req_len == buf_rem) { /* done with pair k */
                        k++;
                        buf_off = buf_view.off[k];
                        buf_rem = buf_view.len[k];
                        ptr = (char*)buf + buf_off;
                    }
                    else { /* req_len < buf_rem, remains in pair k */
                        buf_off += req_len;
                        buf_rem -= req_len;
                        ptr += req_len;
                    }
                } /* while loop k */
            } /* while loop j */
        }
    }

    if ((fd->atomicity) && PNCIO_Feature(fd, PNCIO_LOCKS))
        PNCIO_UNLOCK(fd, lock_off, SEEK_SET, lock_len);

    return total_w_len;
}
