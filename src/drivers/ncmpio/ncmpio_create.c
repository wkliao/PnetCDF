/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_create() : dispatcher->create()
 * ncmpi_open()   : dispatcher->open()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>    /* strcpy(), strchr() */

#if defined(HAVE_LSTAT) || defined(HAVE_ACCESS) || defined(HAVE_OPEN) || defined(HAVE_UNLINK) || defined(HAVE_CLOSE)
#include <sys/types.h> /* lstat(), open() */
#include <sys/stat.h>  /* lstat(), open() */
#include <unistd.h>    /* lstat(), access(), unlink(), open(), close() */
#include <fcntl.h>     /* open() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< ncmpio_create() >----------------------------------------------------*/
int
ncmpio_create(MPI_Comm     comm,
              const char  *path,
              int          cmode,
              int          ncid,
              MPI_Info     user_info, /* user's and env info combined */
              void       **ncpp)
{
    char *env_str, *filename, *mpi_name;
    int rank, nprocs, mpiomode, err, mpireturn, default_format, file_exist=1;
    int use_trunc=1, fstype;
    MPI_File fh;
    MPI_Info info_used;
    NC *ncp=NULL;
    ADIO_File adio_fh=NULL;

    *ncpp = NULL;

    /* Note path's validity and cmode consistency have been checked in
     * ncmpi_create() in src/dispatchers/file.c and
     * path consistency will be done in MPI_File_open */

    /* First, check whether cmode is valid or supported ---------------------*/

    /* NC_DISKLESS is not supported yet */
    if (cmode & NC_DISKLESS) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* NC_MMAP is not supported yet */
    if (cmode & NC_MMAP) DEBUG_RETURN_ERROR(NC_EINVAL_CMODE)

    /* Check cmode for other illegal flags already done in dispatcher layer */

    /* Get default format, in case cmode does not include either
     * NC_64BIT_OFFSET or NC_64BIT_DATA */
    ncmpi_inq_default_format(&default_format);

    /* Handle file clobber --------------------------------------------------*/
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    mpiomode = MPI_MODE_RDWR | MPI_MODE_CREATE;

    /* remove the file system type prefix name if there is any.  For example,
     * when path = "lustre:/home/foo/testfile.nc", remove "lustre:" to make
     * filename pointing to "/home/foo/testfile.nc", so it can be used in POSIX
     * access() below
     */
    filename = ncmpii_remove_file_system_type_prefix(path);

    /* Check if the file already exists, if lstat() or access() is available */
#ifdef HAVE_LSTAT
    /* call lstat() to check the file if exists and if is a symbolic link */
    if (rank == 0) {
        struct stat st_buf;
        st_buf.st_mode = 0;

        if (lstat(filename, &st_buf) == -1) file_exist = 0;
        errno = 0; /* reset errno */

        /* If the file is a regular file, not a symbolic link, then we can
         * delete the file first and later create it when calling
         * MPI_File_open() with MPI_MODE_CREATE. It is OK to delete and then
         * re-create the file if the file is a regular file. If there are other
         * files symbolically linked to this file, then their links will still
         * point to this file after it is re-created.
         *
         * If the file is a symbolic link, then we cannot delete the file, as
         * the link will be gone.
         */
        if (S_ISREG(st_buf.st_mode)) use_trunc = 0;
    }
#elif defined HAVE_ACCESS
    /* if access() is available, use it to check whether file already exists
     * rank 0 calls access() and broadcasts file_exist */
    if (rank == 0) {
        if (access(filename, F_OK) == -1) file_exist = 0;
        errno = 0; /* reset errno */
    }
#endif

    /* check file system type */
    fstype = ADIO_FileSysType(path);

    if (fIsSet(cmode, NC_NOCLOBBER)) {
        /* check if file exists: NC_EEXIST is returned if the file already
         * exists and NC_NOCLOBBER mode is used in ncmpi_create */
#ifdef HAVE_ACCESS
        if (nprocs > 1)
            TRACE_COMM(MPI_Bcast)(&file_exist, 1, MPI_INT, 0, comm);
        if (file_exist) DEBUG_RETURN_ERROR(NC_EEXIST)
#else
        /* add MPI_MODE_EXCL mode for MPI_File_open to check file existence */
        fSet(mpiomode, MPI_MODE_EXCL);
        errno = 0; /* reset errno, as MPI_File_open may change it */
#endif
    }
    else { /* NC_CLOBBER is the default mode in create */
        /* rank 0 truncates or deletes the file and ignores error code.
         * Note calling MPI_File_set_size is expensive as it calls truncate()
         */
        err = NC_NOERR;
        if (rank == 0 && file_exist) {
            if (!use_trunc) { /* delete the file */
#ifdef HAVE_UNLINK
                /* unlink() is likely faster then truncate(), but may be still
                 * expensive
                 */
                err = unlink(filename);
                if (err < 0 && errno != ENOENT)
                    /* ignore ENOENT: file not exist */
                    DEBUG_ASSIGN_ERROR(err, NC_EFILE) /* other error */
                else
                    err = NC_NOERR;
#else
                err = NC_NOERR;
                if (fstype != ADIO_UFS)
                    err = ADIO_File_delete(path);
                else {
                    TRACE_IO(MPI_File_delete, (path, MPI_INFO_NULL));
                    if (mpireturn != MPI_SUCCESS) {
                        int errorclass;
                        MPI_Error_class(mpireturn, &errorclass);
                        if (errorclass != MPI_ERR_NO_SUCH_FILE)
                            /* ignore file not exist */
                            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    }
                }
#endif
            }
            else { /* file is not a regular file, truncate it to zero size */
#ifdef HAVE_TRUNCATE
                err = truncate(filename, 0); /* can be expensive */
                if (err < 0 && errno != ENOENT)
                    /* ignore ENOENT: file not exist */
                    DEBUG_ASSIGN_ERROR(err, NC_EFILE) /* other error */
                else
                    err = NC_NOERR;
#elif defined HAVE_OPEN
                int fd = open(filename, O_TRUNC | O_WRONLY);
                if (fd < 0)
                    DEBUG_ASSIGN_ERROR(err, NC_EFILE)
                else {
                    err = close(fd);
                    if (err < 0)
                        DEBUG_ASSIGN_ERROR(err, NC_EFILE)
                }
#else
                /* call MPI_File_set_size() to truncate the file. Note this can
                 * be expensive.
                 */
                err = NC_NOERR;
                if (fstype != ADIO_UFS) {
                    adio_fh = (ADIO_FileD*) NCI_Calloc(1,sizeof(ADIO_FileD));
                    err = ADIO_File_open(MPI_COMM_SELF, filename, MPI_MODE_RDWR,
                                         MPI_INFO_NULL, adio_fh);
                    if (err == NC_NOERR)
                        ADIO_File_set_size(adio_fh, 0); /* can be expensive */
                    else
                        ADIO_File_close(&adio_fh);
                    NCI_Free(adio_fh);
                    adio_fh = NULL;
                }
                else {
                    TRACE_IO(MPI_File_open, (MPI_COMM_SELF, path, MPI_MODE_RDWR, MPI_INFO_NULL, &fh));
                    if (mpireturn != MPI_SUCCESS) {
                        int errorclass;
                        MPI_Error_class(mpireturn, &errorclass);
                        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    }
                    else {
                        TRACE_IO(MPI_File_set_size, (fh, 0)); /* can be expensive */
                        if (mpireturn != MPI_SUCCESS) {
                            int errorclass;
                            MPI_Error_class(mpireturn, &errorclass);
                            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                        }
                        else {
                            TRACE_IO(MPI_File_close, (&fh));
                            if (mpireturn != MPI_SUCCESS) {
                                int errorclass;
                                MPI_Error_class(mpireturn, &errorclass);
                                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                            }
                        }
                    }
                }
#endif
            }
            if (errno == ENOENT) errno = 0; /* reset errno */
        }
        /* all processes must wait here until file deletion is completed */
        if (nprocs > 1)
            TRACE_COMM(MPI_Bcast)(&err, 1, MPI_INT, 0, comm);
        if (err != NC_NOERR) return err;
    }

    if (fstype != ADIO_UFS) {
        adio_fh = (ADIO_FileD*) NCI_Calloc(1,sizeof(ADIO_FileD));
        adio_fh->file_system = fstype;
    }

    /* allocate buffer for header object NC and initialize its contents */
    ncp = (NC*) NCI_Calloc(1, sizeof(NC));
    if (ncp == NULL) DEBUG_RETURN_ERROR(NC_ENOMEM)

    /* For file create, ignore if NC_NOWRITE set in cmode by user */
    ncp->iomode   = cmode | NC_WRITE;
    ncp->comm     = comm;  /* reuse comm duplicated in dispatch layer */
    ncp->path     = path;  /* reuse path duplicated in dispatch layer */
    ncp->rank     = rank;
    ncp->nprocs   = nprocs;
    ncp->mpiomode = mpiomode;

    /* create file collectively -------------------------------------------- */
    if (fstype != ADIO_UFS) {
        err = ADIO_File_open(comm, filename, mpiomode, user_info, adio_fh);
        if (err != NC_NOERR)
            return err;
    }
    else {
        TRACE_IO(MPI_File_open, (comm, path, mpiomode, user_info, &fh));
        if (mpireturn != MPI_SUCCESS) {
#ifndef HAVE_ACCESS
            if (fIsSet(cmode, NC_NOCLOBBER)) {
                /* This is the case when NC_NOCLOBBER is used in file creation
                 * and function access() is not available. MPI_MODE_EXCL is set
                 * in open mode. When MPI_MODE_EXCL is used and the file
                 * already exists, MPI-IO should return error class
                 * MPI_ERR_FILE_EXISTS. But, some MPI-IO implementations (older
                 * ROMIO) do not correctly return this error class. In this
                 * case, we can do the followings: check errno to see if it set
                 * to EEXIST. Note usually rank 0 makes the file open call and
                 * can be the only one having errno set.
                 */
                if (nprocs > 1)
                    TRACE_COMM(MPI_Bcast)(&errno, 1, MPI_INT, 0, comm);
                if (errno == EEXIST) DEBUG_RETURN_ERROR(NC_EEXIST)
            }
#endif
            return ncmpii_error_mpi2nc(mpireturn, "MPI_File_open");
            /* for NC_NOCLOBBER, MPI_MODE_EXCL was added to mpiomode. If the
             * file already exists, MPI-IO should return error class
             * MPI_ERR_FILE_EXISTS which PnetCDF will return error code
             * NC_EEXIST. This is checked inside of ncmpii_error_mpi2nc()
             */
        }
        else
            /* reset errno, as MPI_File_open may change it, even if it returns
             * MPI_SUCCESS
             */
            errno = 0;
    }

    /* get the I/O hints used/modified by MPI-IO */
    if (fstype != ADIO_UFS) {
        err = ADIO_File_get_info(adio_fh,  &info_used);
        if (err != NC_NOERR) return err;
    }
    else {
        TRACE_IO(MPI_File_get_info, (fh, &info_used));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }

    /* Now the file has been successfully created */
    ncp->collective_fh  = fh;
    ncp->independent_fh = (nprocs > 1) ? MPI_FILE_NULL : fh;
    ncp->fstype         = fstype;
    ncp->adio_fh        = adio_fh;

    /* set the file format version based on the create mode, cmode */
         if (fIsSet(cmode, NC_64BIT_DATA))   ncp->format = 5;
    else if (fIsSet(cmode, NC_64BIT_OFFSET)) ncp->format = 2;
    else {
             if (default_format == NC_FORMAT_CDF5) ncp->format = 5;
        else if (default_format == NC_FORMAT_CDF2) ncp->format = 2;
        else                                       ncp->format = 1;
    }

    fSet(ncp->flags, NC_MODE_CREATE);
    /* create automatically enter write mode */
    fClr(ncp->flags, NC_MODE_RDONLY);
    /* create automatically enter define mode */
    fSet(ncp->flags, NC_MODE_DEF);
    /* PnetCDF default mode is no fill */
    fClr(ncp->flags, NC_MODE_FILL);

    ncp->ncid = ncid;

    /* chunk size for reading header, set to default before check hints */
    ncp->chunk = PNC_DEFAULT_CHUNKSIZE;

    /* initialize unlimited_id as no unlimited dimension yet defined */
    ncp->dims.unlimited_id = -1;

    /* buffer to pack noncontiguous user buffers when calling wait() */
    ncp->ibuf_size = PNC_DEFAULT_IBUF_SIZE;

    /* Extract PnetCDF specific I/O hints from user_info and set default hint
     * values into info_used. Note some MPI libraries, such as MPICH 3.3.1 and
     * priors fail to preserve user hints that are not recognized by the MPI
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

    *ncpp = (void*)ncp;

    return NC_NOERR;
}

