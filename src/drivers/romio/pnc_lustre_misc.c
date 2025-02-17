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
#include "ncmpio_NC.h"



/*----< PNC_File_sync() >----------------------------------------------------*/
int PNC_File_sync(PNC_File fd)
{
    int err = NC_NOERR;

    if (fd->is_open > 0) {
        err = fsync(fd->fd_sys);
        if (err != 0)
            err = ncmpii_error_posix2nc("fsync");
    }

    return err;
}

/*----< PNC_File_delete() >--------------------------------------------------*/
int PNC_File_delete(const char *filename)
{
    int err = NC_NOERR;

    err = unlink(filename);
    if (err != 0)
        err = ncmpii_error_posix2nc("unlink");

    return err;
}

/*----< PNC_File_set_size() >------------------------------------------------*/
int PNC_File_set_size(PNC_File   fd,
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

/*----< PNC_File_get_size() >------------------------------------------------*/
int PNC_File_get_size(PNC_File    fd,
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

#if 0
File seek is not used in PnetCDF !

/*----< PNC_File_seek() >----------------------------------------------------*/
int PNC_File_seek(PNC_File   fd,
                  MPI_Offset offset,
                  int        whence)
{
    int err = NC_NOERR, rank, posix_whence;
    off_t file_off;

file_off = ADIOI_GEN_SeekIndividual(fd, offset, whence, &err);
/*
    switch (whence) {
        case MPI_SEEK_SET: posix_whence = SEEK_SET; break;
        case MPI_SEEK_CUR: posix_whence = SEEK_CUR; break;
        case MPI_SEEK_END: posix_whence = SEEK_END; break;
        default:
            MPI_Comm_rank(fd->comm, &rank);
            fprintf(stderr, "%s line %d: rank %d invalid whence(%d)\n",
                    __func__,__LINE__,rank,whence);
            return NC_EINVAL;
    }
    file_off = lseek(fd->fd_sys, (off_t)offset, posix_whence);
*/
    if (file_off == -1)
        err = ncmpii_error_posix2nc("lseek");

    return err;
}
#endif

/*----< PNC_File_get_info() >------------------------------------------------*/
int PNC_File_get_info(PNC_File  fd,
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

#if 0
static
int PNC_Type_get_combiner(MPI_Datatype datatype, int *combiner)
{
    int ret;
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count ni, na, nc, nt;
    ret = MPI_Type_get_envelope_c(datatype, &ni, &na, &nc, &nt, combiner);
#else
    int ni, na, nt;
    ret = MPI_Type_get_envelope(datatype, &ni, &na, &nt, combiner);
#endif
    return ret;
}

int PNC_Type_ispredef(MPI_Datatype datatype, int *flag)
{
    int ret, combiner;
    ret = PNC_Type_get_combiner(datatype, &combiner);
    switch (combiner) {
        case MPI_COMBINER_NAMED:
        case MPI_COMBINER_F90_INTEGER:
        case MPI_COMBINER_F90_REAL:
        case MPI_COMBINER_F90_COMPLEX:
            *flag = 1;
            break;
        default:
            *flag = 0;
            break;
    }
    return ret;
}

/* utility function for freeing user-defined datatypes,
 * MPI_DATATYPE_NULL and predefined datatypes are ignored,
 * datatype is set to MPI_DATATYPE_NULL upon return */
int PNC_Type_dispose(MPI_Datatype * datatype)
{
    int ret, flag;
    if (*datatype == MPI_DATATYPE_NULL)
        return MPI_SUCCESS;
    ret = PNC_Type_ispredef(*datatype, &flag);
    if (ret == MPI_SUCCESS && !flag)
        ret = MPI_Type_free(datatype);
    *datatype = MPI_DATATYPE_NULL;
    return ret;
}
#endif

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

int
PNC_File_SetInfo(PNC_File fd,
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

#if 0
void PNC_Datatype_iscontig(MPI_Datatype datatype, int *flag)
{
    int combiner;

    PNC_Type_get_combiner(datatype, &combiner);

    switch (combiner) {
        case MPI_COMBINER_NAMED:
            *flag = 1;
            break;
#ifdef MPIIMPL_HAVE_MPI_COMBINER_DUP
        case MPI_COMBINER_DUP:
#endif
        case MPI_COMBINER_CONTIGUOUS:
            {
                int *ints;
                MPI_Aint *adds;
                MPI_Count *cnts;
                MPI_Datatype *types;
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count nints, nadds, ncnts, ntypes;
                MPI_Type_get_envelope_c(datatype, &nints, &nadds, &ncnts, &ntypes, &combiner);
#else
                int nints, nadds, ncnts = 0, ntypes;
                MPI_Type_get_envelope(datatype, &nints, &nadds, &ntypes, &combiner);
#endif
                ints = (int *) ADIOI_Malloc((nints + 1) * sizeof(int));
                adds = (MPI_Aint *) ADIOI_Malloc((nadds + 1) * sizeof(MPI_Aint));
                cnts = (MPI_Count *) ADIOI_Malloc((ncnts + 1) * sizeof(MPI_Count));
                types = (MPI_Datatype *) ADIOI_Malloc((ntypes + 1) * sizeof(MPI_Datatype));
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Type_get_contents_c(datatype, nints, nadds, ncnts, ntypes, ints, adds, cnts,
                                        types);
#else
                MPI_Type_get_contents(datatype, nints, nadds, ntypes, ints, adds, types);
#endif
                PNC_Datatype_iscontig(types[0], flag);

                PNC_Type_dispose(types);
                ADIOI_Free(ints);
                ADIOI_Free(adds);
                ADIOI_Free(cnts);
                ADIOI_Free(types);
            }
            break;
        case MPI_COMBINER_F90_INTEGER:
        case MPI_COMBINER_F90_REAL:
        case MPI_COMBINER_F90_COMPLEX:
            *flag = 1;
            break;
        default:
            *flag = 0;
            break;
    }

    /* This function needs more work. It should check for contiguity
     * in other cases as well. */
}
#endif

/*----< PNC_File_set_view() >------------------------------------------------*/
int PNC_File_set_view(PNC_File      fd,
                      MPI_Offset    disp,
                      MPI_Datatype  etype,
                      MPI_Datatype  filetype,
                      char         *datarep,
                      MPI_Info      info)
{
    int is_predef, i, err, filetype_is_contig;
    MPI_Datatype copy_filetype;
    ADIOI_Flatlist_node *flat_file;

    /* rudimentary checks for incorrect etype/filetype. */
    if (etype == MPI_DATATYPE_NULL)
        return ncmpii_error_mpi2nc(MPI_ERR_ARG, "PNC_File_set_view, etype");

    if (filetype == MPI_DATATYPE_NULL)
        return ncmpii_error_mpi2nc(MPI_ERR_ARG, "PNC_File_set_view, filetype");

    if ((disp < 0) && (disp != MPI_DISPLACEMENT_CURRENT))
        return ncmpii_error_mpi2nc(MPI_ERR_ARG, "PNC_File_set_view, disp");

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count lb;
    MPI_Type_size_c(filetype, &fd->ftype_size);
    MPI_Type_size_c(etype, &fd->etype_size);
    MPI_Type_get_extent_c(filetype, &lb, &fd->ftype_extent);
#else
    MPI_Aint lb;
    MPI_Type_size(filetype, &fd->ftype_size);
    MPI_Type_size(etype, &fd->etype_size);
    MPI_Type_get_extent(filetype, &lb, &fd->ftype_extent);
#endif

    if (fd->etype_size != 0 && fd->ftype_size % fd->etype_size != 0)
        return ncmpii_error_mpi2nc(MPI_ERR_ARG, "PNC_File_set_view, type size");

    err = PNC_File_SetInfo(fd, info);
    if (err != NC_NOERR)
        return err;

    /* PnetCDF only uses etype = MPI_BYTE */
    fd->etype = etype;

    PNC_Type_dispose(&fd->filetype);
    PNC_Type_ispredef(filetype, &is_predef);
    if (is_predef) {
        fd->filetype = filetype;
        filetype_is_contig = 1;
    } else {
        MPI_Type_dup(filetype, &copy_filetype);
        MPI_Type_commit(&copy_filetype);

        fd->filetype = copy_filetype;
        PNC_Datatype_iscontig(fd->filetype, &filetype_is_contig);

        /* check filetype only if it is not a predefined MPI datatype */
        flat_file = ADIOI_Flatten_and_find(fd->filetype);
        err = check_type(flat_file, fd->orig_access_mode, "filetype");
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
// for (i = 0; i < flat_file->count; i++) printf("%s line %d: flat_file count=%lld [%d] indices=%lld blocklens=%lld\n",__func__,__LINE__,flat_file->count,i,flat_file->indices[i],flat_file->blocklens[i]);

        for (i = 0; i < flat_file->count; i++) {
            if (flat_file->blocklens[i]) {
                fd->fp_ind = disp + flat_file->indices[i];
                break;
            }
        }
    }

    return err;
}

