/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in src/dispatchers/file.c
 *
 * ncmpi_sync()         : dispatcher->sync()
 * ncmpi_sync_numrecs() : dispatcher->sync_numrecs()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset() */
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncx.h>
#include "ncmpio_NC.h"

/*----< ncmpio_file_sync() >-------------------------------------------------*/
/* This function must be called collectively, no matter if it is in collective
 * or independent data mode.
 */
int
ncmpio_file_sync(NC *ncp) {
    char *mpi_name;
    int mpireturn;

    if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
        if (ncp->adio_fh == NULL)
            return NC_NOERR;
        return PNCIO_File_sync(ncp->adio_fh);
    }

    /* the remaining of this subroutine are for when using MPI-IO */

    if (ncp->independent_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->independent_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }
    /* when nprocs == 1, ncp->collective_fh == ncp->independent_fh */
    if (ncp->nprocs == 1) return NC_NOERR;

    /* When intra-node aggregation is enabled, non-aggregator's
     * ncp->collective_fh is always MPI_FILE_NULL. When disabled,
     * ncp->collective_fh on all ranks is never MPI_FILE_NULL as collective
     * mode is default in PnetCDF.
     */
    if (ncp->collective_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->collective_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }

    /* Barrier is not necessary ...
    TRACE_COMM(MPI_Barrier)(ncp->comm);
     */

    return NC_NOERR;
}

#define NC_NUMRECS_OFFSET 4

/*----< ncmpio_write_numrecs() >---------------------------------------------*/
/* Only root process writes the new record number into file.
 * This function is called by:
 * 1. ncmpio_sync_numrecs
 * 2. collective nonblocking wait API, if the new number of records is bigger
 */
int
ncmpio_write_numrecs(NC         *ncp,
                     MPI_Offset  new_numrecs)
{
    char *mpi_name;
    int mpireturn, err;
    MPI_File fh;
    MPI_Status mpistatus;

    /* return now if there is no record variable defined */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* When intra-node aggregation is enabled, non-aggregators do not
     * participate any collective calls below.
     */
    if (ncp->num_aggrs_per_node > 0 && ncp->rank != ncp->my_aggr)
        return NC_NOERR;

    /* If not requiring all MPI-IO calls to be collective, non-root processes
     * can return now. This is because only root process writes numrecs to the
     * file header.
     */
    if (!fIsSet(ncp->flags, NC_HCOLL) && ncp->rank > 0)
        return NC_NOERR;

    /* select the right file handler based on the current data mode */
    fh = ncp->independent_fh;
    if (ncp->nprocs > 1 && !NC_indep(ncp))
        fh = ncp->collective_fh;

    /* If collective MPI-IO is required for all MPI-IO calls, then all non-root
     * processes participate the collective write call with zero-size requests.
     */
    if (ncp->rank > 0 && fIsSet(ncp->flags, NC_HCOLL)) {
        if (ncp->fstype != PNCIO_FSTYPE_MPIIO)
            PNCIO_File_write_at_all(ncp->adio_fh, 0, NULL, 0, MPI_BYTE, &mpistatus);
        else
            TRACE_IO(MPI_File_write_at_all, (fh, 0, NULL, 0, MPI_BYTE, &mpistatus));
        return NC_NOERR;
    }

    if (new_numrecs > ncp->numrecs || NC_ndirty(ncp)) {
        int len;
        char pos[8], *buf=pos;

        /* update ncp->numrecs */
        if (new_numrecs > ncp->numrecs) ncp->numrecs = new_numrecs;

        if (ncp->format < 5) {
            if (ncp->numrecs > NC_MAX_INT)
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
            len = X_SIZEOF_SIZE_T;
            err = ncmpix_put_uint32((void**)&buf, (uint)ncp->numrecs);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
        }
        else {
            len = X_SIZEOF_INT64;
            err = ncmpix_put_uint64((void**)&buf, (uint64)ncp->numrecs);
            if (err != NC_NOERR) DEBUG_RETURN_ERROR(err)
        }
        /* ncmpix_put_xxx advances the 1st argument with size len */

        /* explicitly initialize mpistatus object to 0. For zero-length read,
         * MPI_Get_count may report incorrect result for some MPICH version,
         * due to the uninitialized MPI_Status object passed to MPI-IO calls.
         * Thus we initialize it above to work around.
         */
        memset(&mpistatus, 0, sizeof(MPI_Status));

        if (ncp->num_aggrs_per_node > 0 && ncp->rank != ncp->my_aggr)
            /* When intra-node aggregation is enabled, non-aggregators do not
             * participate the collective call.
             */
            return NC_NOERR;

        /* root's file view always includes the entire file header */
        if (fIsSet(ncp->flags, NC_HCOLL) && ncp->nprocs > 1) {
            if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
                err = PNCIO_File_write_at_all(ncp->adio_fh, NC_NUMRECS_OFFSET, (void*)pos,
                                            len, MPI_BYTE, &mpistatus);
                if (err != NC_NOERR) return err;
            }
            else {
                TRACE_IO(MPI_File_write_at_all, (fh, NC_NUMRECS_OFFSET, (void*)pos,
                                                 len, MPI_BYTE, &mpistatus));
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    if (err == NC_EFILE) DEBUG_RETURN_ERROR(NC_EWRITE)
                }
            }
        }
        else {
            if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
                err = PNCIO_File_write_at(ncp->adio_fh, NC_NUMRECS_OFFSET, (void*)pos,
                                        len, MPI_BYTE, &mpistatus);
                if (err != NC_NOERR) return err;
            }
            else {
                TRACE_IO(MPI_File_write_at, (fh, NC_NUMRECS_OFFSET, (void*)pos,
                                             len, MPI_BYTE, &mpistatus));
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                    if (err == NC_EFILE) DEBUG_RETURN_ERROR(NC_EWRITE)
                }
                else
                    err = NC_NOERR;
            }
        }

        if (err == NC_NOERR) {
            /* update the number of bytes written since file open.
             * Because the above MPI write writes either 4 or 8 bytes,
             * calling MPI_Get_count() is sufficient. No need to call
             * MPI_Get_count_c()
             */
            int put_count;
            mpireturn = MPI_Get_count(&mpistatus, MPI_BYTE, &put_count);
            if (mpireturn != MPI_SUCCESS || put_count == MPI_UNDEFINED)
                /* partial write: in this case MPI_Get_elements() is supposed
                 * to be called to obtain the number of type map elements
                 * actually written in order to calculate the true write
                 * amount. Below skips this step and simply ignore the partial
                 * write. See an example usage of MPI_Get_count() in Example
                 * 5.12 from MPI standard document.
                 */
                ncp->put_size += len;
            else
                ncp->put_size += put_count;
        }
    }
    return NC_NOERR;
}

/*----< ncmpio_sync_numrecs() >-----------------------------------------------*/
/* Synchronize the number of records in memory among all processes and write
 * numrecs to file.
 * This function is called by:
 * 1. ncmpi_sync_numrecs(): by the user
 * 2. ncmpi_sync(): by the user
 * 3. ncmpi_end_indep_data(): exit from independent data mode
 * 4. all blocking collective put APIs when writing record variables
 * 5. ncmpi_close(): file close and currently in independent data mode
 *
 * This API is collective, but can be called in independent data mode.
 * Note numrecs is always sync-ed in memory and update in file in collective
 * data mode.
 */
int
ncmpio_sync_numrecs(void *ncdp)
{
    int status=NC_NOERR, mpireturn;
    MPI_Offset max_numrecs;
    NC *ncp=(NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

    /* check if we have defined record variables */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* check write permission */
    if (NC_readonly(ncp)) DEBUG_RETURN_ERROR(NC_EPERM)

    if (!NC_indep(ncp)) /* in collective data mode, numrecs is always sync-ed */
        return NC_NOERR;
    else /* if called in independent mode, we force sync in memory */
        set_NC_ndirty(ncp);

    /* return now if there is no record variabled defined */
    if (ncp->vars.num_rec_vars == 0) return NC_NOERR;

    /* find the max numrecs among all processes
     * Note max numrecs may be smaller than some process's ncp->numrecs
     */
    max_numrecs = ncp->numrecs;
    if (ncp->nprocs > 1) {
        TRACE_COMM(MPI_Allreduce)(&ncp->numrecs, &max_numrecs, 1, MPI_OFFSET,
                                  MPI_MAX, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Allreduce");
    }

    /* root process writes max_numrecs to file */
    status = ncmpio_write_numrecs(ncp, max_numrecs);

    if (ncp->nprocs > 1 && ncp->safe_mode == 1) {
        /* broadcast root's status, because only root writes to the file */
        int root_status = status;
        TRACE_COMM(MPI_Bcast)(&root_status, 1, MPI_INT, 0, ncp->comm);
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, "MPI_Bcast");
        /* root's write has failed, which is serious */
        if (root_status == NC_EWRITE) DEBUG_ASSIGN_ERROR(status, NC_EWRITE)
    }

    /* update numrecs in all processes's memory */
    ncp->numrecs = max_numrecs;

    /* clear numrecs dirty bit */
    fClr(ncp->flags, NC_NDIRTY);

    return status;
}

/*----< ncmpio_sync() >------------------------------------------------------*/
/* This API is a collective subroutine, and must be called in data mode, no
 * matter if it is in collective or independent data mode.
 */
int
ncmpio_sync(void *ncdp)
{
    int err;
    NC *ncp = (NC*)ncdp;

    /* cannot be in define mode */
    if (NC_indef(ncp)) DEBUG_RETURN_ERROR(NC_EINDEFINE)

#if 0
    /* In PnetCDF, header metadata is always sync-ed among all processes.
     * There is no need to re-read the header from file.
     */
    if (NC_readonly(ncp))
        /* calling sync for file opened for read only means re-read header */
        return ncmpio_read_NC(ncp);
#else
    if (NC_readonly(ncp)) return NC_NOERR;
#endif

    /* the only part of header that can be dirty is numrecs (caused only by
     * independent APIs) */
    if (ncp->vars.num_rec_vars > 0 && NC_indep(ncp)) {
        /* sync numrecs in memory among processes and in file */
        set_NC_ndirty(ncp);
        err = ncmpio_sync_numrecs(ncp);
        if (err != NC_NOERR) return err;
    }

    /* calling MPI_File_sync() on both collective and independent handlers */
    return ncmpio_file_sync(ncp);
}

