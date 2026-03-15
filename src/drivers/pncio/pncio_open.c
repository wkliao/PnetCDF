/*
 *  Copyright (C) 2025, Northwestern University
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
int PNCIO_File_open(MPI_Comm    comm,
                    const char *filename,
                    int         amode,
                    MPI_Info    info,
                    PNCIO_File *fd)
{
    /* Before arriving at this subroutine, PNCIO_FileSysType() should have been
     * called to check the file system type.
     */
    char value[MPI_MAX_INFO_VAL + 1], int_str[16];
    int i, err, min_err, status=NC_NOERR;

    fd->comm        = comm;
    fd->filename    = filename;  /* without file system type name prefix */
    fd->atomicity   = 0;
    fd->is_open     = 0;
    fd->access_mode = amode;
    fd->io_buf      = NULL; /* collective buffer used by aggregators only */

    fd->file_view.count = 0; /* flattened fileview in offset-length pairs */
    fd->file_view.size = -1;
    fd->file_view.off = NULL;
    fd->file_view.len = NULL;

    /* create and initialize info object */
    fd->hints = (PNCIO_Hints*) NCI_Calloc(1, sizeof(PNCIO_Hints));
    status = PNCIO_File_set_info(fd, info);
    if (status != NC_NOERR && status != NC_EMULTIDEFINE_HINTS) {
        /* Inconsistent I/O hints is not a fatal error.
         * In PNCIO_File_set_info(), root's hints overwrite local's.
         */
        goto err_out;
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    for (i=0; i<NMEASURES; i++) {
        fd->write_timing[i]  = fd->read_timing[i]  = 0;
        fd->write_counter[i] = fd->read_counter[i] = 0;
    }
#endif

    assert(fd->file_system != PNCIO_FSTYPE_MPIIO);

    /* When hint romio_no_indep_rw hint is set to true, only aggregators open
     * the file. Note fd->is_agg will be set at the end of create/open calls
     * below.
     */
    if (fd->file_system == PNCIO_LUSTRE) {
        if (amode & MPI_MODE_CREATE)
            err = PNCIO_Lustre_create(fd, amode);
        else
            err = PNCIO_Lustre_open(fd);
    }
    else {
        if (amode & MPI_MODE_CREATE)
            err = PNCIO_UFS_create(fd, amode);
        else
            err = PNCIO_UFS_open(fd);
    }
    if (err != NC_NOERR) { /* fatal error */
        status = err;
        goto err_out;
    }

    /* set file striping hints */
    snprintf(int_str, 16, "%d", fd->hints->striping_unit);
    MPI_Info_set(fd->info, "striping_unit", int_str);

    snprintf(int_str, 16, "%d", fd->hints->striping_factor);
    MPI_Info_set(fd->info, "striping_factor", int_str);

    snprintf(int_str, 16, "%d", fd->hints->start_iodevice);
    MPI_Info_set(fd->info, "start_iodevice", int_str);

    /* set file striping hints */
    snprintf(int_str, 16, "%d", fd->hints->cb_nodes);
    MPI_Info_set(fd->info, "cb_nodes", int_str);

    /* add hint "cb_node_list", list of aggregators' rank IDs */
    snprintf(value, 16, "%d", fd->hints->ranklist[0]);
    for (i=1; i<fd->hints->cb_nodes; i++) {
        snprintf(int_str, 16, " %d", fd->hints->ranklist[i]);
        if (strlen(value) + strlen(int_str) >= MPI_MAX_INFO_VAL-5) {
            strcat(value, " ...");
            break;
        }
        strcat(value, int_str);
    }
    MPI_Info_set(fd->info, "cb_node_list", value);

    /* collective buffer size must be at least file striping size */
    if (fd->hints->cb_buffer_size < fd->hints->striping_unit) {
        fd->hints->cb_buffer_size = fd->hints->striping_unit;
        snprintf(int_str, 16, " %d", fd->hints->cb_buffer_size);
        MPI_Info_set(fd->info, "cb_buffer_size", int_str);
    }

    /* collective buffer is used only by I/O aggregators only */
    if (fd->is_agg) {
        fd->io_buf = NCI_Calloc(1, fd->hints->cb_buffer_size);
        if (fd->io_buf == NULL) /* fatal error */
            status = NC_ENOMEM;
    }

err_out:
    MPI_Allreduce(&status, &min_err, 1, MPI_INT, MPI_MIN, comm);
    /* All NC errors are < 0 */

    if (min_err != NC_NOERR) {
        if (status == NC_NOERR && fd->is_open)
            /* close file if opened successfully */
            close(fd->fd_sys);
        NCI_Free(fd->hints);
        if (fd->info != MPI_INFO_NULL)
            MPI_Info_free(&(fd->info));
        if (fd->io_buf != NULL)
            NCI_Free(fd->io_buf);
    }
    return status;
}

