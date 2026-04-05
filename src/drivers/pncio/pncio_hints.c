/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "pncio.h"

#define GET_INFO_INT(key) { \
    MPI_Info_get(users_info, #key, MPI_MAX_INFO_VAL, value, &flag); \
    if (flag) { \
        MPI_Info_set(fh->info, #key, value); \
        fh->hints->key = atoi(value); \
    } \
}

#define GET_INFO_STR(key) { \
    MPI_Info_get(users_info, #key, MPI_MAX_INFO_VAL, value, &flag); \
    if (flag) { \
        MPI_Info_set(fh->info, #key, value); \
        if (!strcasecmp(value, "true")) \
            fh->hints->key = PNCIO_HINT_ENABLE; \
        else if (!strcasecmp(value, "false")) \
            fh->hints->key = PNCIO_HINT_DISABLE; \
        else if (!strcasecmp(value, "automatic")) \
            fh->hints->key = PNCIO_HINT_AUTO; \
        else if (!strcasecmp(value, "enable")) \
            fh->hints->key = PNCIO_HINT_ENABLE; \
        else if (!strcasecmp(value, "disable")) \
            fh->hints->key = PNCIO_HINT_DISABLE; \
        else if (!strcasecmp(value, "inherit")) \
            fh->hints->key = PNCIO_STRIPING_INHERIT; \
    } \
}

#ifdef PNETCDF_DEBUG
#define CHECK_HINT(hint) { \
    if (fh->hints->hint != root_hints->hint) { \
        char int_str[16]; \
        fprintf(stderr, "Error: inconsistent I/O hint %s (%d at rank %d, %d at root)\n", \
                #hint, fh->hints->hint, rank, root_hints->hint); \
        /* overwrite local's hint with root's */ \
        snprintf(int_str, 16, "%d", root_hints->hint); \
        MPI_Info_set(fh->info, #hint, int_str); \
        err = NC_EMULTIDEFINE_HINTS; \
    } \
}
#else
#define CHECK_HINT(hint) { \
    if (fh->hints->hint != root_hints->hint) { \
        /* overwrite local's hint with root's */ \
        char int_str[16]; \
        snprintf(int_str, 16, "%d", root_hints->hint); \
        MPI_Info_set(fh->info, #hint, int_str); \
        err = NC_EMULTIDEFINE_HINTS; \
    } \
}
#endif

/*----< hint_consistency_check() >-------------------------------------------*/
static
int hint_consistency_check(PNCIO_File *fh)
{
    int err, rank;

    MPI_Comm_rank(fh->comm, &rank);

    err = NC_NOERR;

    if (rank == 0)
        /* broadcast root's hints */
        MPI_Bcast(fh->hints, sizeof(PNCIO_Hints), MPI_BYTE, 0, fh->comm);
    else {
        PNCIO_Hints *root_hints;
        root_hints = (PNCIO_Hints*) NCI_Malloc(sizeof(PNCIO_Hints));

        /* broadcast root's hints */
        MPI_Bcast(root_hints, sizeof(PNCIO_Hints), MPI_BYTE, 0, fh->comm);

        /* check hints individually against root's */
        CHECK_HINT(nc_striping);
        CHECK_HINT(striping_factor);
        CHECK_HINT(striping_unit);
        CHECK_HINT(start_iodevice);
        CHECK_HINT(cb_nodes);
        CHECK_HINT(cb_buffer_size);
        CHECK_HINT(ind_rd_buffer_size);
        CHECK_HINT(ind_wr_buffer_size);

        CHECK_HINT(romio_cb_read);
        CHECK_HINT(romio_cb_write);
        CHECK_HINT(romio_ds_read);
        CHECK_HINT(romio_ds_write);

        CHECK_HINT(lustre_overstriping_ratio);

        NCI_Free(root_hints);
    }

    /* All NetCDF erro codes are negative */
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MIN, fh->comm);

    return err;
}

/*----< PNCIO_File_set_info() >----------------------------------------------*/
/* For PnetCDF, a file info object can only be passed to PnetCDF at file create
 * or open call, i.e. I/O hints cannot be changed after file create/open.
 *
 * This subroutine is a collective call, because it checks consistency of all
 * hints among all processes.
 */
int
PNCIO_File_set_info(PNCIO_File *fh,
                    MPI_Info    users_info)
{
    int err=NC_NOERR, flag, nprocs;
    char value[MPI_MAX_INFO_VAL + 1];

    if (users_info == MPI_INFO_NULL)
        MPI_Info_create(&fh->info);
    else
        MPI_Info_dup(users_info, &fh->info);

    MPI_Comm_size(fh->comm, &nprocs);

    /* initialize fh->info and hints to default values */

    /* buffer size for collective I/O */
    MPI_Info_set(fh->info, "cb_buffer_size", PNCIO_CB_BUFFER_SIZE_DFLT);
    fh->hints->cb_buffer_size = atoi(PNCIO_CB_BUFFER_SIZE_DFLT);

    /* default is to let pncio automatically decide whether or not to use
     * collective buffering
     */
    MPI_Info_set(fh->info, "romio_cb_read", "automatic");
    fh->hints->romio_cb_read = PNCIO_HINT_AUTO;
    MPI_Info_set(fh->info, "romio_cb_write", "automatic");
    fh->hints->romio_cb_write = PNCIO_HINT_AUTO;

    /* cb_nodes may be set later right after file open call */
    fh->hints->cb_nodes = 0;

    /* buffer size for data sieving in independent reads */
    MPI_Info_set(fh->info, "ind_rd_buffer_size", PNCIO_IND_RD_BUFFER_SIZE_DFLT);
    fh->hints->ind_rd_buffer_size = atoi(PNCIO_IND_RD_BUFFER_SIZE_DFLT);

    /* buffer size for data sieving in independent writes */
    MPI_Info_set(fh->info, "ind_wr_buffer_size", PNCIO_IND_WR_BUFFER_SIZE_DFLT);
    fh->hints->ind_wr_buffer_size = atoi(PNCIO_IND_WR_BUFFER_SIZE_DFLT);

    /* default is to let romio automatically decide when to use data
     * sieving
     */
    MPI_Info_set(fh->info, "romio_ds_read", "automatic");
    fh->hints->romio_ds_read = PNCIO_HINT_AUTO;
    MPI_Info_set(fh->info, "romio_ds_write", "automatic");
    fh->hints->romio_ds_write = PNCIO_HINT_AUTO;

    /* File striping parameters will be retrieved from the file system set,
     * once the file is opened. These parameters can also be customized by
     * a user's info. Thus, default values used below are to indicate
     * whether or not they have been customized by the users.
     */
    fh->hints->nc_striping = PNCIO_STRIPING_AUTO;
    fh->hints->striping_unit = 0;
    fh->hints->striping_factor = 0;
    fh->hints->start_iodevice = -1;
    /* Lustre overstriping ratio. 0 or 1 means disabled */
    fh->hints->lustre_overstriping_ratio = 1;

    /* add in user's info --------------------------------------------------*/

    if (users_info == MPI_INFO_NULL) goto err_out;

    /* size of internal buffer to be used in collective reads and writes */
    GET_INFO_INT(cb_buffer_size);

    /* enable/disable collective buffering */
    GET_INFO_STR(romio_cb_read);

    GET_INFO_STR(romio_cb_write);

    /* enable/disable data sieving */
    GET_INFO_STR(romio_ds_read);
    GET_INFO_STR(romio_ds_write);

    /* number of I/O aggregators */
    GET_INFO_INT(cb_nodes);
    /* check ill value */
    if (fh->hints->cb_nodes > 0 && fh->hints->cb_nodes <= nprocs) {
        snprintf(value, MPI_MAX_INFO_VAL + 1, "%d", fh->hints->cb_nodes);
        MPI_Info_set(fh->info, "cb_nodes", value);
    }
    else {
        fh->hints->cb_nodes = 0;
        MPI_Info_set(fh->info, "cb_nodes", "0");
    }

    GET_INFO_INT(ind_wr_buffer_size);
    GET_INFO_INT(ind_rd_buffer_size);

    /* file striping configuration */
    GET_INFO_STR(nc_striping);
    GET_INFO_INT(striping_unit);
    GET_INFO_INT(striping_factor);
    GET_INFO_INT(start_iodevice);

    /* Lustre overstriping ratio. 0 or 1 means disabled */
    GET_INFO_INT(lustre_overstriping_ratio);

    /* Check hint consistency among all processes */
err_out:
    if (nprocs > 1)
        err = hint_consistency_check(fh);

    /* PnetCDF ignores the following hints.
     *    cb_config_list
     *    deferred_open
     */

    return err;
}

/*----< PNCIO_File_get_info() >-----------------------------------------------*/
int PNCIO_File_get_info(PNCIO_File *fh,
                        MPI_Info   *info_used)
{
    int err;

    err = MPI_Info_dup(fh->info, info_used);
    if (err == MPI_SUCCESS)
        err = NC_NOERR;
    else
        err = ncmpii_error_mpi2nc(err, "MPI_Info_dup");

    return err;
}

