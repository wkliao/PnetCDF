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

#include "pncio.h"

/*----< PNCIO_File_close() >--------------------------------------------------*/
int PNCIO_File_close(PNCIO_File *fh)
{
    int err = NC_NOERR;

    if (fh->is_open) {
        err = close(fh->fd_sys);
        if (err != 0)
            err = ncmpii_error_posix2nc("close");
    }

    if (fh->hints != NULL) {
        if (fh->hints->aggr_ranks != NULL)
            NCI_Free(fh->hints->aggr_ranks);
        NCI_Free(fh->hints);
    }
    if (fh->info != MPI_INFO_NULL)
        MPI_Info_free(&(fh->info));
    if (fh->io_buf != NULL)
        NCI_Free(fh->io_buf);

    return err;
}

