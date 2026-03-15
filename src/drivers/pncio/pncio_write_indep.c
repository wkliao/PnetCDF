/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/errno.h>
#include <unistd.h>   /* pwrite() */

#include <mpi.h>

#include "pncio.h"

/*----< PNCIO_File_write_at() >----------------------------------------------*/
/* This is an independent call. */
MPI_Offset PNCIO_File_write_at(PNCIO_File *fh,
                               MPI_Offset  offset,
                               const void *buf,
                               PNCIO_View  buf_view)
{
    MPI_Offset w_len;

#ifdef PNETCDF_DEBUG
    assert(fh != NULL);

    if (offset > 0)
        assert(fh->file_view.off == NULL &&
               fh->file_view.len == NULL &&
               fh->file_view.count == 0);
#endif

    if (buf_view.size == 0) /* zero-sized request */
        return NC_NOERR;

    if (buf_view.size < 0)
        return NC_ENEGATIVECNT;

    if (fh->file_view.off == NULL)
        /* This is when calling this subroutione without set fileview first.
         * We store offset into fh->file_view.off.
         */
        fh->file_view.off = &offset;

    if (buf_view.count <= 1 && fh->file_view.count <= 1)
        w_len = PNCIO_UFS_write_contig(fh, buf, buf_view.size, fh->file_view.off[0], 0);
    else if (fh->file_system == PNCIO_UFS)
        w_len = PNCIO_UFS_write_indep(fh, buf, buf_view);
    else if (fh->file_system == PNCIO_LUSTRE)
        w_len = PNCIO_UFS_write_indep(fh, buf, buf_view);
    else
        w_len = NC_EFSTYPE;

    /* reset fileview, as PnetCDF never reuses a fileview */
    fh->file_view.off = NULL;
    fh->file_view.len = NULL;
    fh->file_view.size = 0;
    fh->file_view.count = 0;

    return w_len; /* when w_len < 0, it is an NetCDF error code */
}

