/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>      /* open(), O_CREAT */
#include <sys/types.h>  /* open(), umask(), fstat() */

#if defined(HAVE_SYS_STAT_H) && HAVE_SYS_STAT_H == 1
#include <sys/stat.h>   /* fstat() */
#endif
#include <unistd.h>     /* fstat() */

#include <assert.h>
#include <sys/errno.h>

#include <mpi.h>

#include "pncio.h"

/*----< PNCIO_File_open() >--------------------------------------------------*/
/* This is a collective call.
 * This subroutine is called as part of ncmpi_create() or ncmpi_open().
 *
 * There are two cases that reach to this subroutine.
 * 1. When PnetCDF INA is disabled, all processes participate to call this
 *    subroutine which returns a non-NULL file handle.
 * 2. When PnetCDF INA is enabled, only the INA aggregators call this
 *    subroutine. Non-INA aggregators will have to set the contents of its file
 *    handler separately/independently, such as setting up the I/O hints.
 *
 * Note this subroutine will not be called as part of ncmpi_begin_indep_data().
 * When PNCIO driver is used, ncmpi_begin_indep_data() actually does nothing.
 * When an independent I/O is performed by a process, including non-INA
 * aggregators, PNCIO driver will first check fh->is_open, which indicates
 * whether this process has opened the file and get its fd->fs_sys set to >=0.
 * Check such file-open-on-demand implementation in PNCIO_UFS_write_indep() and
 * PNCIO_UFS_read_indep(), which does not call this subroutine.
 */
int PNCIO_File_open(MPI_Comm    comm,
                    const char *filename,
                    int         amode, /* O_CREAT|O_RDWR, O_RDWR, or O_RDONLY */
                    MPI_Info    info,
                    PNCIO_File *fh)
{
    /* Before arriving at this subroutine, PNCIO_FileSysType() should have been
     * called to check the file system type.
     */
    int err, min_err, status=NC_NOERR;

    fh->comm      = comm;
    fh->filename  = filename; /* without file system type name prefix */
    fh->fd_sys    = -1;       /* file has not yet been opened */
    fh->atomicity = 0;
    fh->is_open   = 0;    /* this rank has opened the file */
    fh->is_agg    = 0;    /* whether this rank is an I/O aggregator */
    fh->amode     = amode;
    fh->io_buf    = NULL; /* collective buffer used by aggregators only */

    fh->file_view.count = 0; /* flattened fileview in offset-length pairs */
    fh->file_view.size  = -1;
    fh->file_view.off   = NULL;
    fh->file_view.len   = NULL;

    /* allocate and initialize info object */
    fh->hints = (PNCIO_Hints*) NCI_Calloc(1, sizeof(PNCIO_Hints));
    status = PNCIO_File_set_info(fh, info);
    if (status != NC_NOERR && status != NC_EMULTIDEFINE_HINTS) {
        /* Inconsistent I/O hints is not a fatal error.
         * In PNCIO_File_set_info(), root's hints overwrite local's.
         */
        goto err_out;
    }

    /* Now, create/open the file. Note fh->is_agg, indicating whether this rank
     * is an I/O aggregator,  will be set at the end of create/open calls.
     */
    if (fh->fstype == PNCIO_FS_LUSTRE) {
        if (amode & O_CREAT)
            err = PNCIO_Lustre_create(fh);
        else
            err = PNCIO_Lustre_open(fh);
    }
    else if (fh->fstype == PNCIO_FS_UFS) {
        /* PNCIO_UFS_open uses fh->amode to tell if create or open */
        err = PNCIO_UFS_open(fh);
    }
    else
        err = NC_EFSTYPE;

    if (err != NC_NOERR) { /* Failer to open the file is a fatal error */
        status = err;
        goto err_out;
    }

    /* collective buffer is used only by I/O aggregators only */
    if (fh->is_agg) {
        fh->io_buf = NCI_Calloc(1, fh->hints->cb_buffer_size);
        if (fh->io_buf == NULL) /* fatal error */
            status = NC_ENOMEM;
    }

err_out:
    MPI_Allreduce(&status, &min_err, 1, MPI_INT, MPI_MIN, comm);
    /* All NC errors are < 0 */

    if (min_err != NC_NOERR) {
        if (status == NC_NOERR && fh->is_open)
            /* close file if opened successfully */
            close(fh->fd_sys);
        NCI_Free(fh->hints);
        if (fh->info != MPI_INFO_NULL)
            MPI_Info_free(&(fh->info));
        if (fh->io_buf != NULL)
            NCI_Free(fh->io_buf);
    }
    return status;
}

