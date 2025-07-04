/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* strdup() */
#include <assert.h>
#include <sys/errno.h>

#include <mpi.h>

#include "adio.h"

/*----< ADIO_File_close() >--------------------------------------------------*/
int ADIO_File_close(ADIO_File *fh)
{
    int err = NC_NOERR;

    err = close((*fh)->fd_sys);
    if (err != 0)
        err = ncmpii_error_posix2nc("close");

    if ((*fh)->hints->ranklist != NULL)
        ADIOI_Free((*fh)->hints->ranklist);
    if ((*fh)->hints->cb_config_list != NULL)
        ADIOI_Free((*fh)->hints->cb_config_list);
    if ((*fh)->hints != NULL)
        ADIOI_Free((*fh)->hints);
    if ((*fh)->info != MPI_INFO_NULL)
        MPI_Info_free(&((*fh)->info));
    if ((*fh)->io_buf != NULL)
        ADIOI_Free((*fh)->io_buf);
    ADIOI_Type_dispose(&(*fh)->filetype);

    if (ADIOI_Flattened_type_keyval != MPI_KEYVAL_INVALID) {
        MPI_Type_free_keyval(&ADIOI_Flattened_type_keyval);
        ADIOI_Flattened_type_keyval = MPI_KEYVAL_INVALID;
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    int i, ntimers, rank;
    double tt[16], max_t[16];

    /* print two-phase I/O timing breakdown */
    MPI_Comm_rank((*fh)->comm, &rank);
    ntimers = 12;
    for (i=0; i<ntimers; i++) tt[i] = (*fh)->coll_write[i];
    MPI_Reduce(tt, max_t, ntimers, MPI_DOUBLE, MPI_MAX, 0, (*fh)->comm);
    if (rank == 0 && max_t[11] > 0)
        printf("%s: TWO-PHASE write init %5.2f pwrite %5.2f post %5.2f comm %5.2f collw %5.2f ntimes %d\n",
        __func__, max_t[1], max_t[2], max_t[4], max_t[3], max_t[0], (int)(max_t[11]));

    ntimers = 12;
    for (i=0; i<ntimers; i++) tt[i] = (*fh)->coll_write[i];
    MPI_Reduce(tt, max_t, ntimers, MPI_DOUBLE, MPI_MAX, 0, (*fh)->comm);
    if (rank == 0 && max_t[11] > 0)
        printf("%s: TWO-PHASE read  init %5.2f pread  %5.2f post %5.2f wait %5.2f collr %5.2f ntimes %d\n",
        __func__, max_t[1], max_t[2], max_t[4], max_t[3], max_t[0], (int)(max_t[11]));
#endif

    return err;
}
