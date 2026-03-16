/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdbool.h> /* type bool */

#include <pncio.h>

/*----< PNCIO_File_read_at_all() >-------------------------------------------*/
/* This is a collective call. */
MPI_Offset PNCIO_File_read_at_all(PNCIO_File *fh,
                                  MPI_Offset  offset,
                                  void       *buf,
                                  PNCIO_View  buf_view)
{
    int err=NC_NOERR;
    MPI_Offset r_len;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Offset one_len = (MPI_Offset)buf_view.size;
#else
    int one_len = (int)buf_view.size;
#endif

#ifdef PNETCDF_DEBUG
    assert(fh != NULL);

    if (offset > 0)
        assert(fh->file_view.off == NULL &&
               fh->file_view.len == NULL &&
               fh->file_view.count == 0);
#endif

    if (buf_view.size < 0) err = NC_ENEGATIVECNT;

    if (fh->file_view.off == NULL) {
        /* This is when calling this subroutione without set fileview first.
         * We store offset into fh->file_view.off.
         */
        fh->file_view.off   = &offset;
        fh->file_view.count = 1;
        fh->file_view.size  = buf_view.size;
        fh->file_view.len   = &one_len;
    }

    r_len = PNCIO_UFS_read_coll(fh, buf, buf_view);

    /* reset fileview, as PnetCDF never reuses a fileview */
    fh->file_view.off = NULL;
    fh->file_view.len = NULL;
    fh->file_view.size = 0;
    fh->file_view.count = 0;

    return (err == NC_NOERR) ? r_len : err;
}

