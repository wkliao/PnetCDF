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
#include "adio.h"

/*----< PNCIO_File_get_info() >-----------------------------------------------*/
int PNCIO_File_get_info(PNCIO_File *fd,
                        MPI_Info  *info_used)
{
    int err;

    err = MPI_Info_dup(fd->info, info_used);
    if (err == MPI_SUCCESS)
        err = NC_NOERR;
    else
        err = ncmpii_error_mpi2nc(err, "MPI_Info_dup");

    return err;
}

/*----< ADIOI_Info_check_and_install_int() >---------------------------------*/
static
int ADIOI_Info_check_and_install_int(PNCIO_File  *fd,
                                     MPI_Info    info,
                                     const char *key,
                                     int        *local_cache)
{
    int intval, tmp_val, flag, ret = 0;
    char value[MPI_MAX_INFO_VAL + 1];

    ADIOI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
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

        ADIOI_Info_set(fd->info, key, value);
        /* some file systems do not cache hints in the fd struct */
        if (local_cache != NULL)
            *local_cache = intval;
    }
fn_exit:
    return ret;
}

/*----< ADIOI_Info_check_and_install_enabled() >-----------------------------*/
static
int ADIOI_Info_check_and_install_enabled(PNCIO_File  *fd,
                                         MPI_Info    info,
                                         const char *key,
                                         int        *local_cache)
{
    int tmp_val, flag, ret = 0;
    char value[MPI_MAX_INFO_VAL + 1];

    ADIOI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
    if (flag) {
        if (!strcmp(value, "enable") || !strcmp(value, "ENABLE")) {
            ADIOI_Info_set(fd->info, key, value);
            *local_cache = ADIOI_HINT_ENABLE;
        } else if (!strcmp(value, "disable") || !strcmp(value, "DISABLE")) {
            ADIOI_Info_set(fd->info, key, value);
            *local_cache = ADIOI_HINT_DISABLE;
        } else if (!strcmp(value, "automatic") || !strcmp(value, "AUTOMATIC")) {
            ADIOI_Info_set(fd->info, key, value);
            *local_cache = ADIOI_HINT_AUTO;
            /* treat the user-provided string like "enabled":  either it is a
             * hint ROMIO knows about and can support it, or ROMIO will not
             * return the hint at all in the MPI_File_get_info info object
             */
        } else if (!strcmp(value, "requested") || !strcmp(value, "REQUESTED")) {
            ADIOI_Info_set(fd->info, key, "enable");
            *local_cache = ADIOI_HINT_ENABLE;
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

/*----< ADIOI_Info_check_and_install_true() >--------------------------------*/
static
int ADIOI_Info_check_and_install_true(PNCIO_File  *fd,
                                      MPI_Info    info,
                                      const char *key,
                                      int        *local_cache)
{
    int flag, tmp_val, ret = 0;
    char value[MPI_MAX_INFO_VAL + 1];

    ADIOI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
    if (flag) {
        if (!strcmp(value, "true") || !strcmp(value, "TRUE")) {
            ADIOI_Info_set(fd->info, key, value);
            *local_cache = 1;
        } else if (!strcmp(value, "false") || !strcmp(value, "FALSE")) {
            ADIOI_Info_set(fd->info, key, value);
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

/*----< ADIOI_Info_check_and_install_str() >---------------------------------*/
static
int ADIOI_Info_check_and_install_str(PNCIO_File   *fd,
                                     MPI_Info     info,
                                     const char  *key,
                                     char       **local_cache)
{
    int flag, ret = 0;
    size_t len;
    char value[MPI_MAX_INFO_VAL + 1];

    ADIOI_Info_get(info, key, MPI_MAX_INFO_VAL, value, &flag);
    if (flag) {
        ADIOI_Info_set(fd->info, key, value);
        len = (strlen(value) + 1) * sizeof(char);
        *local_cache = ADIOI_Malloc(len);
        if (*local_cache == NULL) {
            ret = NC_ENOMEM;
            goto fn_exit;
        }
        ADIOI_Strncpy(*local_cache, value, len);
    }
    /* if it has been set already, we ignore it the second time.
     * otherwise we would get an error if someone used the same
     * info value with a cb_config_list value in it in a couple
     * of calls, which would be irritating. */
fn_exit:
    return ret;
}

/*----< PNCIO_File_SetInfo() >------------------------------------------------*/
/* When users_info == MPI_INFO_NULL, this subroutine is an independent call.
 * When users_info != MPI_INFO_NULL, this subroutine is a collective call.
 * because it calls ADIOI_Info_check_and_install_xxx(), which checks the
 * consistency of all hints values set in user's info object.
 */
int
PNCIO_File_SetInfo(PNCIO_File *fd,
                   MPI_Info   users_info)
{
    int flag, nprocs = 0, len;
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

        /* has user specified striping or server buffering parameters
         * and do they have the same value on all processes? */
        if (users_info != MPI_INFO_NULL) {

            /* striping information */
            ADIOI_Info_get(users_info, "striping_unit", MPI_MAX_INFO_VAL,
                           value, &flag);
            if (flag)
                ADIOI_Info_set(info, "striping_unit", value);

            ADIOI_Info_get(users_info, "striping_factor", MPI_MAX_INFO_VAL,
                           value, &flag);
            if (flag)
                ADIOI_Info_set(info, "striping_factor", value);

            ADIOI_Info_get(users_info, "start_iodevice", MPI_MAX_INFO_VAL,
                           value, &flag);
            if (flag)
                ADIOI_Info_set(info, "start_iodevice", value);

            /* Lustre overstriping ratio. 0 or 1 means disabled */
            ADIOI_Info_get(users_info, "lustre_overstriping_ratio",
                           MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                ADIOI_Info_set(info, "lustre_overstriping_ratio", value);
        }

        /* buffer size for collective I/O */
        ADIOI_Info_set(info, "cb_buffer_size", ADIOI_CB_BUFFER_SIZE_DFLT);
        fd->hints->cb_buffer_size = atoi(ADIOI_CB_BUFFER_SIZE_DFLT);

        /* default is to let romio automatically decide when to use
         * collective buffering
         */
        ADIOI_Info_set(info, "romio_cb_read", "automatic");
        fd->hints->cb_read = ADIOI_HINT_AUTO;
        ADIOI_Info_set(info, "romio_cb_write", "automatic");
        fd->hints->cb_write = ADIOI_HINT_AUTO;

        fd->hints->cb_config_list = NULL;

        /* cb_nodes may be set later right after file open call */
        fd->hints->cb_nodes = 0;

        /* hint indicating that no indep. I/O will be performed on this file */
        ADIOI_Info_set(info, "romio_no_indep_rw", "false");
        fd->hints->no_indep_rw = 0;

        /* hint to set a threshold percentage for a datatype's size/extent at
         * which data sieving should be done in collective I/O */
        ADIOI_Info_set(info, "romio_cb_ds_threshold", "0");
        fd->hints->cb_ds_threshold = 0;

        /* buffer size for data sieving in independent reads */
        ADIOI_Info_set(info, "ind_rd_buffer_size",
                       ADIOI_IND_RD_BUFFER_SIZE_DFLT);
        fd->hints->ind_rd_buffer_size = atoi(ADIOI_IND_RD_BUFFER_SIZE_DFLT);

        /* buffer size for data sieving in independent writes */
        ADIOI_Info_set(info, "ind_wr_buffer_size",
                       ADIOI_IND_WR_BUFFER_SIZE_DFLT);
        fd->hints->ind_wr_buffer_size = atoi(ADIOI_IND_WR_BUFFER_SIZE_DFLT);

        /* default is to let romio automatically decide when to use data
         * sieving
         */
        ADIOI_Info_set(info, "romio_ds_read", "automatic");
        fd->hints->ds_read = ADIOI_HINT_AUTO;
        ADIOI_Info_set(info, "romio_ds_write", "automatic");
        fd->hints->ds_write = ADIOI_HINT_AUTO;

        /* still to do: tune this a bit for a variety of file systems. there's
         * no good default value so just leave it unset */
        fd->hints->striping_unit = 0;
        fd->hints->striping_factor = 0;
        fd->hints->start_iodevice = -1;
        /* Lustre overstriping ratio. 0 or 1 means disabled */
        fd->hints->fs_hints.lustre.overstriping_ratio = 1;

        fd->hints->initialized = 1;
    }

    /* add in user's info if supplied */
    if (users_info != MPI_INFO_NULL) {
        ADIOI_Info_check_and_install_int(fd, users_info, "cb_buffer_size",
                                         &(fd->hints->cb_buffer_size));

        /* for collective I/O, try to be smarter about when to do data sieving
         * using a specific threshold for the datatype size/extent
         * (percentage 0-100%) */
        ADIOI_Info_check_and_install_int(fd, users_info,
                                         "romio_cb_ds_threshold",
                                         &(fd->hints->cb_ds_threshold));

        /* new hints for enabling/disabling coll. buffering on
         * reads/writes
         */
        ADIOI_Info_check_and_install_enabled(fd, users_info, "romio_cb_read",
                                             &(fd->hints->cb_read));
        if (fd->hints->cb_read == ADIOI_HINT_DISABLE) {
            /* romio_cb_read overrides no_indep_rw */
            ADIOI_Info_set(fd->info, "romio_no_indep_rw", "false");
            fd->hints->no_indep_rw = ADIOI_HINT_DISABLE;
        }

        ADIOI_Info_check_and_install_enabled(fd, users_info, "romio_cb_write",
                                             &(fd->hints->cb_write));
        if (fd->hints->cb_write == ADIOI_HINT_DISABLE) {
            /* romio_cb_write overrides no_indep_rw */
            ADIOI_Info_set(fd->info, "romio_no_indep_rw", "false");
            fd->hints->no_indep_rw = ADIOI_HINT_DISABLE;
        }

        /* Has the user indicated all I/O will be done collectively? */
        ADIOI_Info_check_and_install_true(fd, users_info, "romio_no_indep_rw",
                                          &(fd->hints->no_indep_rw));
        if (fd->hints->no_indep_rw == 1) {
            /* if 'no_indep_rw' set, also hint that we will do
             * collective buffering: if we aren't doing independent io,
             * then we have to do collective  */
            ADIOI_Info_set(fd->info, "romio_cb_write", "enable");
            ADIOI_Info_set(fd->info, "romio_cb_read", "enable");
            fd->hints->cb_read = ADIOI_HINT_ENABLE;
            fd->hints->cb_write = ADIOI_HINT_ENABLE;
        }
        /* new hints for enabling/disabling data sieving on
         * reads/writes
         */
        ADIOI_Info_check_and_install_enabled(fd, users_info, "romio_ds_read",
                                             &(fd->hints->ds_read));
        ADIOI_Info_check_and_install_enabled(fd, users_info, "romio_ds_write",
                                             &(fd->hints->ds_write));

        if (fd->hints->cb_nodes == 0) { /* never set before */
            /* MPI_File_open path sets up some data structrues that don't
             * get resized in the MPI_File_set_view path, so ignore
             * cb_nodes in the set_view case */
            ADIOI_Info_check_and_install_int(fd, users_info, "cb_nodes",
                                             &(fd->hints->cb_nodes));
            /* check ill value */
            if (fd->hints->cb_nodes > 0 && fd->hints->cb_nodes <= nprocs) {
                snprintf(value, MPI_MAX_INFO_VAL + 1, "%d", fd->hints->cb_nodes);
                ADIOI_Info_set(fd->info, "cb_nodes", value);
            }
            else {
                fd->hints->cb_nodes = 0;
                ADIOI_Info_set(fd->info, "cb_nodes", "0");
            }
        }

        ADIOI_Info_check_and_install_int(fd, users_info, "ind_wr_buffer_size",
                                         &(fd->hints->ind_wr_buffer_size));
        ADIOI_Info_check_and_install_int(fd, users_info, "ind_rd_buffer_size",
                                         &(fd->hints->ind_rd_buffer_size));

        if (fd->hints->cb_config_list == NULL) {
            /* only set cb_config_list if it isn't already set.  Note that
             * since we set it below, this ensures that the cb_config_list hint
             * will be set at file open time either by the user or to the
             * default */
            /* if it has been set already, we ignore it the second time.
             * otherwise we would get an error if someone used the same info
             * value with a cb_config_list value in it in a couple of calls,
             * which would be irritating. */
            ADIOI_Info_check_and_install_str(fd, users_info, "cb_config_list",
                                             &(fd->hints->cb_config_list));

        }
        /* Now we use striping information in common code so we should
         * process hints for it. */
        ADIOI_Info_check_and_install_int(fd, users_info, "striping_unit",
                                         &(fd->hints->striping_unit));

        ADIOI_Info_check_and_install_int(fd, users_info, "striping_factor",
                                         &(fd->hints->striping_factor));

        ADIOI_Info_check_and_install_int(fd, users_info, "start_iodevice",
                                         &(fd->hints->start_iodevice));

        /* Lustre overstriping ratio. 0 or 1 means disabled */
        ADIOI_Info_check_and_install_int(fd, users_info,
                         "lustre_overstriping_ratio",
                         &(fd->hints->fs_hints.lustre.overstriping_ratio));
    }

    /* Begin hint post-processig: some hints take precedence over or conflict
     * with others, or aren't supported by some file systems */

    /* handle cb_config_list default value here; avoids an extra
     * free/alloc and insures it is always set
     */
    if (fd->hints->cb_config_list == NULL) {
        ADIOI_Info_set(fd->info, "cb_config_list", ADIOI_CB_CONFIG_LIST_DFLT);
        len = (strlen(ADIOI_CB_CONFIG_LIST_DFLT) + 1) * sizeof(char);
        fd->hints->cb_config_list = ADIOI_Malloc(len);
        if (fd->hints->cb_config_list == NULL)
            return NC_ENOMEM;

        ADIOI_Strncpy(fd->hints->cb_config_list, ADIOI_CB_CONFIG_LIST_DFLT, len);
    }
    /* deferred_open won't be set by callers, but if the user doesn't
     * explicitly disable collecitve buffering (two-phase) and does hint that
     * io w/o independent io is going on, we'll set this internal hint as a
     * convenience */
    if (((fd->hints->cb_read != ADIOI_HINT_DISABLE)
         && (fd->hints->cb_write != ADIOI_HINT_DISABLE)
         && fd->hints->no_indep_rw)) {
        fd->hints->deferred_open = 1;
    } else {
        /* setting romio_no_indep_rw enable and romio_cb_{read,write}
         * disable at the same time doesn't make sense. honor
         * romio_cb_{read,write} and force the no_indep_rw hint to
         * 'disable' */
        ADIOI_Info_set(fd->info, "romio_no_indep_rw", "false");
        fd->hints->no_indep_rw = 0;
        fd->hints->deferred_open = 0;
    }

    return NC_NOERR;
}

