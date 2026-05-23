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

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "pncio.h"

/*----< PNCIO_File_set_view() >----------------------------------------------*/
/* This subroutine is an independent call.  */
int PNCIO_File_set_view(PNCIO_File   *fh,
                        MPI_Aint      npairs,
#ifdef HAVE_MPI_LARGE_COUNT
                        MPI_Count    *offsets,
                        MPI_Count    *lengths
#else
                        MPI_Offset   *offsets,
                        int          *lengths
#endif
)
{
    MPI_Aint i;

#ifdef PNETCDF_DEBUG
    if (npairs == 0) assert(offsets == NULL && lengths == NULL);
    for (i=0; i<npairs; i++) assert(lengths[i] > 0);
#endif

    if (npairs == 0 && offsets == NULL && lengths == NULL) {
        /* This is called to reset fileview to be entire file visible */
        fh->file_view.count = 0;
        fh->file_view.size = -1;
        fh->file_view.off   = NULL;
        fh->file_view.len   = NULL;
        fh->file_view.idx   = 0;
        fh->file_view.rem   = 0;
        return NC_NOERR;
    }

    fh->file_view.count = npairs;
    fh->file_view.off   = offsets;
    fh->file_view.len   = lengths;
    fh->file_view.idx   = 0;
    fh->file_view.rem   = (npairs > 0) ? lengths[0] : 0;

    /* Calculate size of fileview. Note PnetCDF has already sorted the
     * offset-length pairs into a monotonically non-decreasing order.
     *
     * fh->file_view.size is initialized to -1 and should be reset to -1 at the
     * end of PNCIO_File_read/write(). This allows PNCIO to tell whether a call
     * to PNCIO_File_set_view() is made prior to PNCIO_File_read/write(). If
     * PNCIO_File_set_view() is called, fh->file_view.size is always >= 0.
     */
    fh->file_view.size = 0;
    for (i=0; i<npairs; i++)
        fh->file_view.size += lengths[i];

    return NC_NOERR;
}

