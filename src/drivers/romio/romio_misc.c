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
#ifdef HAVE_UNISTD_H
#include <unistd.h> /* fsync(), unlink(), ftruncate(), lseek() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "adio.h"



/*----< ADIO_File_sync() >---------------------------------------------------*/
int ADIO_File_sync(ADIO_File fd)
{
    int err = NC_NOERR;

    if (fd->is_open > 0) {
        err = fsync(fd->fd_sys);
        if (err != 0)
            err = ncmpii_error_posix2nc("fsync");
    }

    return err;
}

/*----< ADIO_File_delete() >-------------------------------------------------*/
int ADIO_File_delete(const char *filename)
{
    int err = NC_NOERR;
    char *path = ncmpii_remove_file_system_type_prefix(filename);

    err = unlink(path);
    if (err != 0)
        err = ncmpii_error_posix2nc("unlink");

    return err;
}

/*----< ADIO_File_set_size() >-----------------------------------------------*/
int ADIO_File_set_size(ADIO_File  fd,
                      MPI_Offset size)
{
    int err = NC_NOERR, rank;

    MPI_Comm_rank(fd->comm, &rank);

    if (rank == 0) {
        err = ftruncate(fd->fd_sys, (off_t) size);
        if (err != 0)
            err = ncmpii_error_posix2nc("ftruncate");
    }

    MPI_Bcast(&err, 1, MPI_INT, 0, fd->comm);

    return err;
}

/*----< ADIO_File_get_size() >-----------------------------------------------*/
int ADIO_File_get_size(ADIO_File   fd,
                      MPI_Offset *size)
{
    int err = NC_NOERR, rank;
    MPI_Offset msg[2];

    MPI_Comm_rank(fd->comm, &rank);

    if (rank == 0) {
        *size = lseek(fd->fd_sys, 0, SEEK_END);
        if (*size == -1)
            err = ncmpii_error_posix2nc("lseek");
        msg[0] = err;
        msg[1] = *size;
    }

    MPI_Bcast(msg, 2, MPI_OFFSET, 0, fd->comm);
    err = (int)msg[0];
    *size = msg[1];

    return err;
}

/*----< ADIO_File_get_info() >-----------------------------------------------*/
int ADIO_File_get_info(ADIO_File fd,
                      MPI_Info *info_used)
{
    int err;

    err = MPI_Info_dup(fd->info, info_used);
    if (err == MPI_SUCCESS)
        err = NC_NOERR;
    else
        err = ncmpii_error_mpi2nc(err, "MPI_Info_dup");

    return err;
}

int ADIOI_Info_check_and_install_int(ADIO_File fd, MPI_Info info, const char *key,
                                     int *local_cache)
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

int ADIOI_Info_check_and_install_enabled(ADIO_File fd, MPI_Info info, const char *key,
                                         int *local_cache)
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
            /* treat the user-provided string like "enabled":  either it is a hint
             * ROMIO knows about and can support it, or ROMIO will not return the
             * hint at all in the MPI_File_get_info info object */
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

int ADIOI_Info_check_and_install_true(ADIO_File fd, MPI_Info info, const char *key,
                                      int *local_cache)
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

int ADIOI_Info_check_and_install_str(ADIO_File fd, MPI_Info info, const char *key,
                                     char **local_cache)
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

/*----< ADIO_File_SetInfo() >------------------------------------------------*/
/* When users_info == MPI_INFO_NULL, this subroutine is an independent call.
 * When users_info != MPI_INFO_NULL, this subroutine is a collective call.
 * because it calls ADIOI_Info_check_and_install_xxx(), which checks the
 * consistency of all hints values set in user's info object.
 */
int
ADIO_File_SetInfo(ADIO_File fd,
                  MPI_Info users_info)
{
    int flag, nprocs = 0, len, ok_to_override_cb_nodes = 0;
    char value[MPI_MAX_INFO_VAL + 1];
    MPI_Info info;

    if (fd->hints->initialized && users_info == MPI_INFO_NULL)
        return NC_NOERR;

    /* set other hints */
    MPI_Comm_size(fd->comm, &nprocs);

    /* initialize info and hints to default values if they haven't been
     * previously initialized
     */
    if (!fd->hints->initialized) {
        /* Hints here are the ones that can only be set at MPI_File_open */
        MPI_Info_create(&(fd->info));

        info = fd->info;

        /* has user specified striping or server buffering parameters
         * and do they have the same value on all processes? */
        if (users_info != MPI_INFO_NULL) {

            /* striping information */
            ADIOI_Info_get(users_info, "striping_unit", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                ADIOI_Info_set(fd->info, "striping_unit", value);

            ADIOI_Info_get(users_info, "striping_factor", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                ADIOI_Info_set(fd->info, "striping_factor", value);

            ADIOI_Info_get(users_info, "start_iodevice", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                ADIOI_Info_set(fd->info, "start_iodevice", value);
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

        /* for Lustre, default is set in construct_aggr_list() */
        if (fd->file_system != ADIO_LUSTRE) {
            /* number of processes that perform I/O in collective I/O */
            snprintf(value, MPI_MAX_INFO_VAL + 1, "%d", nprocs);
            ADIOI_Info_set(info, "cb_nodes", value);
            fd->hints->cb_nodes = nprocs;
        }

        /* hint indicating that no indep. I/O will be performed on this file */
        ADIOI_Info_set(info, "romio_no_indep_rw", "false");
        fd->hints->no_indep_rw = 0;

        /* hint to set a threshold percentage for a datatype's size/extent at
         * which data sieving should be done in collective I/O */
        ADIOI_Info_set(info, "romio_cb_ds_threshold", "0");
        fd->hints->cb_ds_threshold = 0;

        /* buffer size for data sieving in independent reads */
        ADIOI_Info_set(info, "ind_rd_buffer_size", ADIOI_IND_RD_BUFFER_SIZE_DFLT);
        fd->hints->ind_rd_buffer_size = atoi(ADIOI_IND_RD_BUFFER_SIZE_DFLT);

        /* buffer size for data sieving in independent writes */
        ADIOI_Info_set(info, "ind_wr_buffer_size", ADIOI_IND_WR_BUFFER_SIZE_DFLT);
        fd->hints->ind_wr_buffer_size = atoi(ADIOI_IND_WR_BUFFER_SIZE_DFLT);

        /* default is to let romio automatically decide when to use data
         * sieving
         */
        ADIOI_Info_set(info, "romio_ds_read", "automatic");
        fd->hints->ds_read = ADIOI_HINT_AUTO;
        ADIOI_Info_set(info, "romio_ds_write", "automatic");
        fd->hints->ds_write = ADIOI_HINT_AUTO;
        ADIOI_Info_set(info, "romio_ds_wr_lb", ADIOI_DS_WR_LB);
        fd->hints->ds_wr_lb = atoi(ADIOI_DS_WR_LB);

        /* still to do: tune this a bit for a variety of file systems. there's
         * no good default value so just leave it unset */
        fd->hints->striping_unit = 0;

        fd->hints->initialized = 1;

        /* ADIO_Open sets up collective buffering arrays.  If we are in this
         * path from say set_file_view, then we've don't want to adjust the
         * array: we'll get a segfault during collective i/o.  We only want to
         * look at the users cb_nodes if it's open time  */
        ok_to_override_cb_nodes = 1;

    }

    /* add in user's info if supplied */
    if (users_info != MPI_INFO_NULL) {
        ADIOI_Info_check_and_install_int(fd, users_info, "cb_buffer_size",
                                         &(fd->hints->cb_buffer_size));

        /* for collective I/O, try to be smarter about when to do data sieving
         * using a specific threshold for the datatype size/extent
         * (percentage 0-100%) */
        ADIOI_Info_check_and_install_int(fd, users_info, "romio_cb_ds_threshold",
                                         &(fd->hints->cb_ds_threshold));

        /* new hints for enabling/disabling coll. buffering on
         * reads/writes
         */
        ADIOI_Info_check_and_install_enabled(fd, users_info, "romio_cb_read",
                                             &(fd->hints->cb_read));
        if (fd->hints->cb_read == ADIOI_HINT_DISABLE) {
            /* romio_cb_read overrides no_indep_rw */
            ADIOI_Info_set(info, "romio_no_indep_rw", "false");
            fd->hints->no_indep_rw = ADIOI_HINT_DISABLE;
        }

        ADIOI_Info_check_and_install_enabled(fd, users_info, "romio_cb_write",
                                             &(fd->hints->cb_write));
        if (fd->hints->cb_write == ADIOI_HINT_DISABLE) {
            /* romio_cb_write overrides no_indep_rw */
            ADIOI_Info_set(info, "romio_no_indep_rw", "false");
            fd->hints->no_indep_rw = ADIOI_HINT_DISABLE;
        }

        /* Has the user indicated all I/O will be done collectively? */
        ADIOI_Info_check_and_install_true(fd, users_info, "romio_no_indep_rw",
                                          &(fd->hints->no_indep_rw));
        if (fd->hints->no_indep_rw == 1) {
            /* if 'no_indep_rw' set, also hint that we will do
             * collective buffering: if we aren't doing independent io,
             * then we have to do collective  */
            ADIOI_Info_set(info, "romio_cb_write", "enable");
            ADIOI_Info_set(info, "romio_cb_read", "enable");
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
        ADIOI_Info_check_and_install_int(fd, users_info, "romio_ds_wr_lb",
                                         &(fd->hints->ds_wr_lb));

        if (ok_to_override_cb_nodes) {
            /* MPI_File_open path sets up some data structrues that don't
             * get resized in the MPI_File_set_view path, so ignore
             * cb_nodes in the set_view case */
            ADIOI_Info_check_and_install_int(fd, users_info, "cb_nodes",
                                             &(fd->hints->cb_nodes));
            /* for Lustre, default is set in construct_aggr_list() */
            if (fd->file_system != ADIO_LUSTRE &&
                (fd->hints->cb_nodes <= 0 || fd->hints->cb_nodes > nprocs)) {
                /* can't ask for more aggregators than mpi processes, though it
                 * might be interesting to think what such oversubscription
                 * might mean... someday */
                snprintf(value, MPI_MAX_INFO_VAL + 1, "%d", nprocs);
                ADIOI_Info_set(info, "cb_nodes", value);
                fd->hints->cb_nodes = nprocs;
            }
        }
        /* if (ok_to_override_cb_nodes) */
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
    }

    /* Begin hint post-processig: some hints take precedence over or conflict
     * with others, or aren't supported by some file systems */

    /* handle cb_config_list default value here; avoids an extra
     * free/alloc and insures it is always set
     */
    if (fd->hints->cb_config_list == NULL) {
        ADIOI_Info_set(info, "cb_config_list", ADIOI_CB_CONFIG_LIST_DFLT);
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
        ADIOI_Info_set(info, "romio_no_indep_rw", "false");
        fd->hints->no_indep_rw = 0;
        fd->hints->deferred_open = 0;
    }

    return NC_NOERR;
}

static
int check_type(ADIOI_Flatlist_node *flat_type,
               int access_mode,
               const char *type_kind)
{
    char err_msg[128];

    err_msg[0] = '\0';

    /* MPI standard requires the displacements of etype and filetype be
     * non-negative */
    if (flat_type->flag & ADIOI_TYPE_NEGATIVE) {
        sprintf(err_msg, "displacements of %s must be non-negative", type_kind);
        goto err_check;
    }

    /* MPI standard requires the displacements of etype and filetype be in a
     * monotonically nondecreasing order */
    if (flat_type->flag & ADIOI_TYPE_DECREASE) {
        sprintf(err_msg, "displacements of %s must be in a monotonically nondecreasing order",
                type_kind);
        goto err_check;
    }

    /* If the file is opened for writing, neither the etype nor the
     * filetype is permitted to contain overlapping regions.
     */
    if (((access_mode & MPI_MODE_WRONLY) || (access_mode & MPI_MODE_RDWR)) &&
        (flat_type->flag & ADIOI_TYPE_OVERLAP)) {
        sprintf(err_msg, "%s is not permitted to contain overlapping regions", type_kind);
        goto err_check;
    }

    return NC_NOERR;

err_check:
    return ncmpii_error_mpi2nc(MPI_ERR_IO, err_msg);
}

/*----< ADIO_File_set_view() >-----------------------------------------------*/
/* For PnetCDF, this subroutine is always an independent call.
 *
 * Note PnetCDF calls MPI_File_set_view() only using the followings.
 * Argument etype is always MPI_BYTE.
 * Argument datarep is always "native".
 * Argument info is always MPI_INFO_NULL. When info is MPI_INFO_NULL, this
 * subroutine is an independent call, because there is no need to check hint
 * consistency among all processes.
 */
int ADIO_File_set_view(ADIO_File     fd,
                       MPI_Offset    disp,
                       MPI_Datatype  filetype,
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
    int is_predef, i, err=NC_NOERR, filetype_is_contig;
    MPI_Datatype copy_filetype;

// printf("%s line %d: filetype = %s\n",__func__,__LINE__,(filetype == MPI_DATATYPE_NULL)?"NULL":"NOT NULL");

    if (filetype == MPI_DATATYPE_NULL) { /* called from intra_node_aggregation() */
        fd->flat_file = NCI_Malloc(sizeof(ADIOI_Flatlist_node));
        fd->flat_file->indices   = offsets;
        fd->flat_file->blocklens = lengths;
        fd->flat_file->count     = npairs;
        fd->flat_file->lb_idx    = -1;
        fd->flat_file->ub_idx    = -1;
        fd->flat_file->flag      = 0;
        fd->flat_file->refct     = 1;

        /* mark this request comes from intra_node_aggregation() */
        fd->filetype = MPI_DATATYPE_NULL;

        /* Note: PnetCDF uses only ADIO_EXPLICIT_OFFSET. In addition, the
         * passed-in offsets and lengths are not based on MPI-IO fileview. They
         * are flattened byte offsets and sizes.
         */
        fd->disp = 0;
        fd->fp_ind = 0;
        return NC_NOERR;
    }

    if ((disp < 0) && (disp != MPI_DISPLACEMENT_CURRENT))
        return ncmpii_error_mpi2nc(MPI_ERR_ARG, "ADIO_File_set_view, disp");

#if 0
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count lb;
    MPI_Type_size_c(filetype, &fd->ftype_size);
    MPI_Type_get_extent_c(filetype, &lb, &fd->ftype_extent);
#else
    MPI_Aint lb;
    MPI_Type_size(filetype, &fd->ftype_size);
    MPI_Type_get_extent(filetype, &lb, &fd->ftype_extent);
#endif
#endif

    /* When info is MPI_INFO_NULL, ADIO_File_SetInfo() is an independent call.
     * Otherwise, it is collective, because it checks hint consistency.
     * PnetCDF always uses MPI_INFO_NULL when setting file view.
    err = ADIO_File_SetInfo(fd, info);
    if (err != NC_NOERR)
        return err;
     */

    /* free fileview if set previously */
    ADIOI_Type_dispose(&fd->filetype);

    /* fd->flat_file should have need free by the callback of MPI_Type_free() */
    fd->flat_file = NULL;

    ADIOI_Type_ispredef(filetype, &is_predef);
    if (is_predef) {
        fd->filetype = filetype;
        filetype_is_contig = 1;
    } else {
        MPI_Type_dup(filetype, &copy_filetype);
        MPI_Type_commit(&copy_filetype);

        fd->filetype = copy_filetype;
        ADIOI_Datatype_iscontig(fd->filetype, &filetype_is_contig);

        /* check filetype only if it is not a predefined MPI datatype */
        fd->flat_file = ADIOI_Flatten_and_find(fd->filetype);
        err = check_type(fd->flat_file, fd->access_mode, "filetype");
        if (err != NC_NOERR)
            return err;
    }

    /* file displacement is an absolute byte position relative to the beginning
     * of a file. The displacement defines the location where a view begins.
     */
    fd->disp = disp;

    /* reset MPI-IO file pointer to point to the first byte that can
     * be accessed in this view. */

// printf("%s line %d: filetype_is_contig=%d\n",__func__,__LINE__,(filetype_is_contig!=0));
    if (filetype_is_contig)
        fd->fp_ind = disp;
    else {
// for (i = 0; i < fd->flat_file->count; i++) printf("%s line %d: fd->flat_file count=%lld [%d] indices=%lld blocklens=%lld\n",__func__,__LINE__,fd->flat_file->count,i,fd->flat_file->indices[i],fd->flat_file->blocklens[i]);

        for (i = 0; i < fd->flat_file->count; i++) {
            if (fd->flat_file->blocklens[i]) {
                fd->fp_ind = disp + fd->flat_file->indices[i];
                break;
            }
        }
    }

    return err;
}

