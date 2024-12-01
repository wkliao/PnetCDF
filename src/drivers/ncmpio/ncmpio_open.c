/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_open() : dispatcher->open()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* strcpy() */
#ifdef HAVE_ACCESS
#include <unistd.h>  /* access() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"
#ifdef ENABLE_SUBFILING
#include "ncmpio_subfile.h"
#endif

/*----< ncmpio_open() >------------------------------------------------------*/
int
ncmpio_open(MPI_Comm     comm,
            const char  *path,
            int          omode,
            int          ncid,
            MPI_Info     user_info, /* user's and env info combined */
            void       **ncpp)
{
    char *filename, *env_str;
    int i, mpiomode, err, status=NC_NOERR, mpireturn, fstype;
    MPI_File fh;
    MPI_Info info_used;
    NC *ncp=NULL;
    ADIO_File adio_fh=NULL;

    *ncpp = NULL;

    /* Note path's validity and omode consistency have been checked in
     * ncmpi_open() in src/dispatchers/file.c and
     * path consistency will be done in MPI_File_open */

    /* First, check whether omode is valid or supported ---------------------*/
    /* NC_DISKLESS is not supported yet */
    if (omode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* NC_MMAP is not supported yet */
    if (omode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_OMODE)

    /* check file system type */
    fstype = ADIO_FileSysType(path);

    if (fstype != ADIO_UFS) {
        adio_fh = (ADIO_FileD*) NCI_Calloc(1,sizeof(ADIO_FileD));
        adio_fh->file_system = fstype;
    }

    /* remove the file system type prefix name if there is any.  For example,
     * when path = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * filename pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * access() below
     */
    filename = ncmpii_remove_file_system_type_prefix(path);

#if 0 && defined(HAVE_ACCESS)
    if (mpiomode == MPI_MODE_RDONLY) { /* file should already exit */
        int rank, file_exist;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            if (access(filename, F_OK) == 0) file_exist = 1;
            else                             file_exist = 0;
        }
        TRACE_COMM(MPI_Bcast)(&file_exist, 1, MPI_INT, 0, comm);
        if (!file_exist) DEBUG_RETURN_ERROR(NC_ENOENT)
    }
#endif

    /* path's validity and omode consistency have been checked in ncmpi_open()
     * in src/dispatchers/file.c */

    /* allocate buffer for header object NC */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    MPI_Comm_rank(comm, &ncp->rank);
    MPI_Comm_size(comm, &ncp->nprocs);
    ncp->comm   = comm;  /* reuse comm duplicated in dispatch layer */
    ncp->path   = path;  /* reuse path duplicated in dispatch layer */
    ncp->iomode = omode;

    /* open file collectively ---------------------------------------------- */
    mpiomode = fIsSet(omode, NC_WRITE) ? MPI_MODE_RDWR : MPI_MODE_RDONLY;
    ncp->mpiomode = mpiomode;

    if (fstype != ADIO_UFS) {
        err = ADIO_File_open(comm, filename, mpiomode, user_info, adio_fh);
        if (err != NC_NOERR) return err;
    }
    else {
        TRACE_IO(MPI_File_open)(comm, path, mpiomode, user_info, &fh);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_open");
    }

    /* Now the file has been successfully opened */
    ncp->collective_fh  = fh;
    ncp->independent_fh = (ncp->nprocs > 1) ? MPI_FILE_NULL : fh;
    ncp->fstype         = fstype;
    ncp->adio_fh        = adio_fh;

    /* get the file info used/modified by MPI-IO */
    if (fstype != ADIO_UFS) {
        err = ADIO_File_get_info(adio_fh, &info_used);
        if (err != NC_NOERR) return err;
    }
    else {
        mpireturn = MPI_File_get_info(fh, &info_used);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_get_info");
    }

    /* Now the file has been successfully opened */

    /* PnetCDF default fill mode is no fill */
    fClr(ncp->flags, NC_MODE_FILL);
    if (!fIsSet(omode, NC_WRITE)) fSet(ncp->flags, NC_MODE_RDONLY);

    ncp->ncid = ncid;

    /* chunk size for reading header (set default before check hints) */
    ncp->chunk = PNC_DEFAULT_CHUNKSIZE;

    /* buffer to pack noncontiguous user buffers when calling wait() */
    ncp->ibuf_size = PNC_DEFAULT_IBUF_SIZE;

    /* Extract PnetCDF specific I/O hints from user_info and set default hint
     * values into info_used. Note some MPI libraries, such as MPICH 3.3.1 and
     * priors fail to preserve user hints that are not recogniozed by the MPI
     * libraries.
     */
    ncmpio_set_pnetcdf_hints(ncp, user_info, info_used);

#ifdef PNETCDF_DEBUG
    /* PNETCDF_DEBUG is set at configure time, which will be overwritten by
     * the run-time environment variable PNETCDF_SAFE_MODE */
    ncp->safe_mode = 1;
#endif
    /* If environment variable PNETCDF_SAFE_MODE is set to 1, then we perform
     * a strict consistent test, i.e. arguments used in def_dim/def_var APIs
     */
    if ((env_str = getenv("PNETCDF_SAFE_MODE")) != NULL) {
        if (*env_str == '0') ncp->safe_mode = 0;
        else                 ncp->safe_mode = 1;
        /* if PNETCDF_SAFE_MODE is set but without a value, *env_str can
         * be '\0' (null character). In this case, safe_mode is enabled */
    }

    /* construct the list of compute nodes */
    ncp->node_ids = NULL;
    if (ncp->num_aggrs_per_node != 0 || fstype != ADIO_UFS) {
        err = ncmpii_construct_node_list(comm, &ncp->num_nodes, &ncp->node_ids);
        if (err != NC_NOERR) return err;
        if (adio_fh != NULL) adio_fh->num_nodes = ncp->num_nodes;
    }

    /* set cb_nodes and construct the cb_node rank list */
    if (fstype != ADIO_UFS) {
        int i;
        char value[MPI_MAX_INFO_VAL + 1];

        if (fstype == ADIO_LUSTRE) {
            ADIO_Lustre_set_aggr_list(adio_fh, ncp->num_nodes, ncp->node_ids);

            MPI_Info_set(info_used, "romio_filesystem_type", "LUSTRE:");
            sprintf(value, "%d", adio_fh->hints->num_osts);
            MPI_Info_set(info_used, "lustre_num_osts", value);
        }

        /* set file striping hints */
        sprintf(value, "%d", adio_fh->hints->cb_nodes);
        MPI_Info_set(info_used, "cb_nodes", value);

        /* add hint "aggr_list", list of aggregators' rank IDs */
        value[0] = '\0';
        for (i=0; i<adio_fh->hints->cb_nodes; i++) {
            char str[16];
            if (i == 0)
                snprintf(str, sizeof(str), "%d", adio_fh->hints->ranklist[i]);
            else
                snprintf(str, sizeof(str), " %d", adio_fh->hints->ranklist[i]);
            if (strlen(value) + strlen(str) >= MPI_MAX_INFO_VAL-5) {
                strcat(value, " ...");
                break;
            }
            strcat(value, str);
        }
        MPI_Info_set(info_used, "aggr_list", value);

#if 0
int rank; MPI_Comm_rank(comm, &rank);
if (rank == 0) {
    int  i, nkeys;
    MPI_Info_get_nkeys(info_used, &nkeys);
    printf("%s line %d: MPI File Info: nkeys = %d\n",__func__,__LINE__,nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(info_used, i, key);
        MPI_Info_get_valuelen(info_used, key, &valuelen, &flag);
        MPI_Info_get(info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}
#endif
    }
    ncp->mpiinfo = info_used; /* is not MPI_INFO_NULL */

    /* determine whether to enable intra-node aggregation and set up all
     * intra-node aggregation metadata.
     * ncp->num_aggrs_per_node = 0, or non-zero indicates whether this feature
     *     is disabled or enabled globally for all processes, respectively.
     * ncp->my_aggr = -1 or >= 0 indicates whether aggregation is effectively
     *     enabled for the aggregation group of this process.
     */
    ncp->my_aggr = -1;
    if (ncp->num_aggrs_per_node != 0) {
        err = ncmpio_intra_node_aggr_init(ncp);
        if (err != NC_NOERR) return err;
    }

    if (ncp->node_ids != NULL) NCI_Free(ncp->node_ids);

    /* read header from file into NC object pointed by ncp -------------------*/
    err = ncmpio_hdr_get_NC(ncp);
    if (err == NC_ENULLPAD) status = NC_ENULLPAD; /* non-fatal error */
    else if (err != NC_NOERR) { /* fatal error */
        ncmpio_close_files(ncp, 0);
        ncmpio_free_NC(ncp);
        return err;
    }

#ifdef ENABLE_SUBFILING
    if (ncp->subfile_mode) {
        /* check subfiling attribute */
        err = ncmpio_get_att(ncp, NC_GLOBAL, "_PnetCDF_SubFiling.num_subfiles",
                             &ncp->num_subfiles, MPI_INT);
        if (err == NC_NOERR && ncp->num_subfiles > 1) {
            int i;
            /* ignore error NC_ENOTATT if this attribute is not defined */
            for (i=0; i<ncp->vars.ndefined; i++) {
                /* variables may have different numbers of subfiles */
                err = ncmpio_get_att(ncp, i, "_PnetCDF_SubFiling.num_subfiles",
                             &ncp->vars.value[i]->num_subfiles,MPI_INT);
                if (err == NC_ENOTATT) continue;
                if (err != NC_NOERR) return err;
                if (ncp->vars.value[i]->num_subfiles > 1) {
                    /* find the orginal ndims of variable i */
                    err = ncmpio_get_att(ncp,i,"_PnetCDF_SubFiling.ndims_org",
                                 &ncp->vars.value[i]->ndims_org,MPI_INT);
                    if (err != NC_NOERR) return err;
                    ncp->vars.value[i]->dimids_org = (int*) NCI_Malloc(
                              ncp->vars.value[i]->ndims_org * SIZEOF_INT);
                    err = ncmpio_get_att(ncp,i,"_PnetCDF_SubFiling.dimids_org",
                              ncp->vars.value[i]->dimids_org, MPI_INT);
                    if (err != NC_NOERR) return err;
                }
            }
            /* open subfile */
            err = ncmpio_subfile_open(ncp);
            if (err != NC_NOERR) return err;
        }
        else ncp->num_subfiles = 0;
    }
    else
        ncp->num_subfiles = 0;
#endif

#ifndef SEARCH_NAME_LINEARLY
    /* initialize and populate name lookup tables ---------------------------*/
    ncmpio_hash_table_populate_NC_dim(&ncp->dims, ncp->dims.hash_size);
    ncmpio_hash_table_populate_NC_var(&ncp->vars, ncp->vars.hash_size);
    ncmpio_hash_table_populate_NC_attr(ncp);
    for (i=0; i<ncp->vars.ndefined; i++)
        ncp->vars.value[i]->attrs.hash_size = ncp->hash_size_attr;
#endif

    *ncpp = (void*)ncp;

    return status;
}

