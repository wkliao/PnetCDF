/*
 *  Copyright (C) 2025, Northwestern University
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
int PNCIO_File_set_view(PNCIO_File   *fd,
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

    fd->file_view.count = npairs;
    fd->file_view.off   = offsets;
    fd->file_view.len   = lengths;
    fd->file_view.idx   = 0;
    fd->file_view.rem   = (npairs > 0) ? lengths[0] : 0;

    /* Calculate size of fileview. Note PnetCDF has already sorted the
     * offset-length pairs into a monotonically non-decreasing order.
     */
    fd->file_view.size = 0;
    for (i=0; i<npairs; i++)
        fd->file_view.size += lengths[i];

    return NC_NOERR;
}

