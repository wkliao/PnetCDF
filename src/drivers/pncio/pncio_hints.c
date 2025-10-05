/*
 *  Copyright (C) 2025, Northwestern University
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

/*----< PNCIO_File_get_info() >-----------------------------------------------*/
int PNCIO_File_get_info(PNCIO_File *fd,
                        MPI_Info   *info_used)
{
    int err;

    err = MPI_Info_dup(fd->info, info_used);
    if (err == MPI_SUCCESS)
        err = NC_NOERR;
    else
        err = ncmpii_error_mpi2nc(err, "MPI_Info_dup");

    return err;
}

/*----< Info_check_and_install_int() >---------------------------------------*/
static
int Info_check_and_install_int(PNCIO_File *fd,
                               MPI_Info    info,
                               const char *key,
                               int        *local_cache)
{
    int intval, tmp_val, flag, ret = 0;
    char value[MPI_MAX_INFO_VAL + 1];

    MPI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
    if (flag) {
        intval = atoi(value);
        tmp_val = intval;

        MPI_Bcast(&tmp_val, 1, MPI_INT, 0, fd->comm);
        /* --BEGIN ERROR HANDLING-- */
        if (tmp_val != intval) {
            ret = ncmpii_error_mpi2nc(MPI_ERR_NOT_SAME, __func__);
            goto fn_exit;
        }
        /* --END ERROR HANDLING-- */

        MPI_Info_set(fd->info, key, value);
        /* some file systems do not cache hints in the fd struct */
        if (local_cache != NULL)
            *local_cache = intval;
    }
fn_exit:
    return ret;
}

/*----< Info_check_and_install_enabled() >-----------------------------------*/
static
int Info_check_and_install_enabled(PNCIO_File *fd,
                                   MPI_Info    info,
                                   const char *key,
                                   int        *local_cache)
{
    int tmp_val, flag, ret = 0;
    char value[MPI_MAX_INFO_VAL + 1];

    MPI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
    if (flag) {
        if (!strcmp(value, "enable") || !strcmp(value, "ENABLE")) {
            MPI_Info_set(fd->info, key, value);
            *local_cache = PNCIO_HINT_ENABLE;
        } else if (!strcmp(value, "disable") || !strcmp(value, "DISABLE")) {
            MPI_Info_set(fd->info, key, value);
            *local_cache = PNCIO_HINT_DISABLE;
        } else if (!strcmp(value, "automatic") || !strcmp(value, "AUTOMATIC")) {
            MPI_Info_set(fd->info, key, value);
            *local_cache = PNCIO_HINT_AUTO;
            /* treat the user-provided string like "enabled":  either it is a
             * hint ROMIO knows about and can support it, or ROMIO will not
             * return the hint at all in the MPI_File_get_info info object
             */
        } else if (!strcmp(value, "requested") || !strcmp(value, "REQUESTED")) {
            MPI_Info_set(fd->info, key, "enable");
            *local_cache = PNCIO_HINT_ENABLE;
        }

        tmp_val = *local_cache;

        MPI_Bcast(&tmp_val, 1, MPI_INT, 0, fd->comm);
        /* --BEGIN ERROR HANDLING-- */
        if (tmp_val != *local_cache) {
            ret = ncmpii_error_mpi2nc(MPI_ERR_NOT_SAME, __func__);
            goto fn_exit;
        }
        /* --END ERROR HANDLING-- */
    }
fn_exit:
    return ret;
}

/*----< Info_check_and_install_true() >--------------------------------------*/
static
int Info_check_and_install_true(PNCIO_File *fd,
                                MPI_Info    info,
                                const char *key,
                                int        *local_cache)
{
    int flag, tmp_val, ret = 0;
    char value[MPI_MAX_INFO_VAL + 1];

    MPI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
    if (flag) {
        if (!strcmp(value, "true") || !strcmp(value, "TRUE")) {
            MPI_Info_set(fd->info, key, value);
            *local_cache = 1;
        } else if (!strcmp(value, "false") || !strcmp(value, "FALSE")) {
            MPI_Info_set(fd->info, key, value);
            *local_cache = 0;
        }
        tmp_val = *local_cache;

        MPI_Bcast(&tmp_val, 1, MPI_INT, 0, fd->comm);
        /* --BEGIN ERROR HANDLING-- */
        if (tmp_val != *local_cache) {
            ret = ncmpii_error_mpi2nc(MPI_ERR_NOT_SAME, __func__);
            goto fn_exit;
        }
        /* --END ERROR HANDLING-- */
    }
fn_exit:
    return ret;
}

#if 0
/*----< Info_check_and_install_str() >---------------------------------------*/
static
int Info_check_and_install_str(PNCIO_File   *fd,
                               MPI_Info     info,
                               const char  *key,
                               char       **local_cache)
{
    int flag, ret = 0;
    size_t len;
    char value[MPI_MAX_INFO_VAL + 1];

    MPI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
    if (flag) {
        MPI_Info_set(fd->info, key, value);
        len = (strlen(value) + 1) * sizeof(char);
        *local_cache = NCI_Malloc(len);
        if (*local_cache == NULL) {
            ret = NC_ENOMEM;
            goto fn_exit;
        }
        strncpy(*local_cache, value, len);
    }
fn_exit:
    return ret;
}
#endif

/*----< PNCIO_File_SetInfo() >------------------------------------------------*/
/* For PnetCDF, a file info object can only be passed to PnetCDF at file create
 * or open call, i.e. I/O hints cannot be changed after file create/open.
 *
 * When users_info == MPI_INFO_NULL, this subroutine is an independent call.
 * When users_info != MPI_INFO_NULL, this subroutine is a collective call,
 * because it calls Info_check_and_install_xxx(), which checks the consistency
 * of all hints values set in user's info object.
 */
int
PNCIO_File_SetInfo(PNCIO_File *fd,
                   MPI_Info    users_info)
{
    int nprocs = 0;
    char value[MPI_MAX_INFO_VAL + 1];

    if (fd->hints->initialized && users_info == MPI_INFO_NULL)
        return NC_NOERR;

    /* set other hints */
    MPI_Comm_size(fd->comm, &nprocs);

    /* initialize info and hints to default values if they haven't been
     * previously initialized
     */
    if (!fd->hints->initialized) {
        /* Hints here are the ones that can only be set at MPI_File_open */
        MPI_Info info;
        MPI_Info_create(&(fd->info));

        info = fd->info;

        /* buffer size for collective I/O */
        MPI_Info_set(info, "cb_buffer_size", PNCIO_CB_BUFFER_SIZE_DFLT);
        fd->hints->cb_buffer_size = atoi(PNCIO_CB_BUFFER_SIZE_DFLT);

        /* default is to let romio automatically decide when to use
         * collective buffering
         */
        MPI_Info_set(info, "romio_cb_read", "automatic");
        fd->hints->cb_read = PNCIO_HINT_AUTO;
        MPI_Info_set(info, "romio_cb_write", "automatic");
        fd->hints->cb_write = PNCIO_HINT_AUTO;

        /* cb_nodes may be set later right after file open call */
        fd->hints->cb_nodes = 0;

        /* hint indicating that no indep. I/O will be performed on this file */
        MPI_Info_set(info, "romio_no_indep_rw", "false");
        fd->hints->no_indep_rw = 0;

        /* buffer size for data sieving in independent reads */
        MPI_Info_set(info, "ind_rd_buffer_size",
                       PNCIO_IND_RD_BUFFER_SIZE_DFLT);
        fd->hints->ind_rd_buffer_size = atoi(PNCIO_IND_RD_BUFFER_SIZE_DFLT);

        /* buffer size for data sieving in independent writes */
        MPI_Info_set(info, "ind_wr_buffer_size",
                       PNCIO_IND_WR_BUFFER_SIZE_DFLT);
        fd->hints->ind_wr_buffer_size = atoi(PNCIO_IND_WR_BUFFER_SIZE_DFLT);

        /* default is to let romio automatically decide when to use data
         * sieving
         */
        MPI_Info_set(info, "romio_ds_read", "automatic");
        fd->hints->ds_read = PNCIO_HINT_AUTO;
        MPI_Info_set(info, "romio_ds_write", "automatic");
        fd->hints->ds_write = PNCIO_HINT_AUTO;

        /* File striping parameters will be retrieved from the file system set,
         * once the file is opened. These parameters can also be customized by
         * a user's info. Thus, default values used below are to indicate
         * whether or not they have been customized by the users.
         */
        fd->hints->striping_unit = 0;
        fd->hints->striping_factor = 0;
        fd->hints->start_iodevice = -1;
        /* Lustre overstriping ratio. 0 or 1 means disabled */
        fd->hints->fs_hints.lustre.overstriping_ratio = 1;

        fd->hints->initialized = 1;
    }

    /* add in user's info if supplied */
    if (users_info != MPI_INFO_NULL) {
        Info_check_and_install_int(fd, users_info, "cb_buffer_size",
                                   &(fd->hints->cb_buffer_size));

        /* new hints for enabling/disabling coll. buffering on
         * reads/writes
         */
        Info_check_and_install_enabled(fd, users_info, "romio_cb_read",
                                       &(fd->hints->cb_read));
        if (fd->hints->cb_read == PNCIO_HINT_DISABLE) {
            /* romio_cb_read overrides no_indep_rw */
            MPI_Info_set(fd->info, "romio_no_indep_rw", "false");
            fd->hints->no_indep_rw = PNCIO_HINT_DISABLE;
        }

        Info_check_and_install_enabled(fd, users_info, "romio_cb_write",
                                       &(fd->hints->cb_write));
        if (fd->hints->cb_write == PNCIO_HINT_DISABLE) {
            /* romio_cb_write overrides no_indep_rw */
            MPI_Info_set(fd->info, "romio_no_indep_rw", "false");
            fd->hints->no_indep_rw = PNCIO_HINT_DISABLE;
        }

        /* Has the user indicated all I/O will be done collectively? */
        Info_check_and_install_true(fd, users_info, "romio_no_indep_rw",
                                    &(fd->hints->no_indep_rw));
        if (fd->hints->no_indep_rw == 1) {
            /* if 'no_indep_rw' set, also hint that we will do
             * collective buffering: if we aren't doing independent io,
             * then we have to do collective  */
            MPI_Info_set(fd->info, "romio_cb_write", "enable");
            MPI_Info_set(fd->info, "romio_cb_read", "enable");
            fd->hints->cb_read = PNCIO_HINT_ENABLE;
            fd->hints->cb_write = PNCIO_HINT_ENABLE;
        }
        /* new hints for enabling/disabling data sieving on
         * reads/writes
         */
        Info_check_and_install_enabled(fd, users_info, "romio_ds_read",
                                       &(fd->hints->ds_read));
        Info_check_and_install_enabled(fd, users_info, "romio_ds_write",
                                       &(fd->hints->ds_write));

        if (fd->hints->cb_nodes == 0) { /* never set before */
            /* MPI_File_open path sets up some data structures that don't
             * get resized in the MPI_File_set_view path, so ignore
             * cb_nodes in the set_view case */
            Info_check_and_install_int(fd, users_info, "cb_nodes",
                                       &(fd->hints->cb_nodes));
            /* check ill value */
            if (fd->hints->cb_nodes > 0 && fd->hints->cb_nodes <= nprocs) {
                snprintf(value, MPI_MAX_INFO_VAL + 1, "%d", fd->hints->cb_nodes);
                MPI_Info_set(fd->info, "cb_nodes", value);
            }
            else {
                fd->hints->cb_nodes = 0;
                MPI_Info_set(fd->info, "cb_nodes", "0");
            }
        }

        Info_check_and_install_int(fd, users_info, "ind_wr_buffer_size",
                                   &(fd->hints->ind_wr_buffer_size));
        Info_check_and_install_int(fd, users_info, "ind_rd_buffer_size",
                                   &(fd->hints->ind_rd_buffer_size));

        /* Now we use striping information in common code so we should
         * process hints for it. */
        Info_check_and_install_int(fd, users_info, "striping_unit",
                                   &(fd->hints->striping_unit));

        Info_check_and_install_int(fd, users_info, "striping_factor",
                                   &(fd->hints->striping_factor));

        Info_check_and_install_int(fd, users_info, "start_iodevice",
                                   &(fd->hints->start_iodevice));

        /* Lustre overstriping ratio. 0 or 1 means disabled */
        Info_check_and_install_int(fd, users_info, "lustre_overstriping_ratio",
                         &(fd->hints->fs_hints.lustre.overstriping_ratio));
    }

    /* PnetCDF ignores the following hints.
     *    cb_config_list
     *    deferred_open
     */

    return NC_NOERR;
}

