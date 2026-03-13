/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/errno.h>
#include <unistd.h>   /* pwrite() */

#include <mpi.h>

#include "pncio.h"

#ifdef PNETCDF_DEBUG
static MPI_Offset first_ost_id;
#endif

/*----< PNCIO_WriteContig() >------------------------------------------------*/
MPI_Offset PNCIO_WriteContig(PNCIO_File *fd,
                             const void *buf,
                             MPI_Offset  w_size,
                             MPI_Offset  offset)
{
    char *p;
    ssize_t err = 0;
    size_t w_count;
    MPI_Offset bytes_xfered = 0;

    if (w_size == 0) return NC_NOERR;

#ifdef PNETCDF_DEBUG
    int rank;
    MPI_Offset ost_id;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    ost_id = (offset / fd->hints->striping_unit) % fd->hints->striping_factor;
    if (first_ost_id == -1)
        first_ost_id = ost_id;
    else if (ost_id != first_ost_id)
        printf("Warning in %s at rank %d:pwrite offset=%lld w_size=%lld ost_id=%lld not same 1st ost %lld\n",__func__,rank,offset,w_size,ost_id,first_ost_id);
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
    p = (char *) buf;
    while (bytes_xfered < w_size) {
        w_count = w_size - bytes_xfered;
        err = pwrite(fd->fd_sys, p, w_count, offset + bytes_xfered);
        if (err == -1)
            goto err_out;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->write_timing[2] += MPI_Wtime() - timing;
#endif

err_out:
    if (err == -1)
        bytes_xfered = ncmpii_error_posix2nc("pwrite");

    return bytes_xfered;
}

/*----< PNCIO_File_write_at() >----------------------------------------------*/
/* This is an independent call. */
MPI_Offset PNCIO_File_write_at(PNCIO_File *fh,
                               MPI_Offset  offset,
                               const void *buf,
                               PNCIO_View  buf_view)
{
    MPI_Offset w_len;

#ifdef PNETCDF_DEBUG
    assert(fh != NULL);
#endif

    if (buf_view.size == 0) /* zero-sized request */
        return NC_NOERR;

    if (buf_view.size < 0)
        return NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        return NC_EPERM;

#ifdef PNETCDF_DEBUG
    if (offset > 0)
        assert(fh->file_view.off == NULL &&
               fh->file_view.len == NULL &&
               fh->file_view.count == 0);
#endif

    if (fh->file_view.off == NULL)
        /* This is when calling this subroutione without set fileview first.
         * We store offset into fh->file_view.off.
         */
        fh->file_view.off = &offset;

    if (buf_view.count <= 1 && fh->file_view.count <= 1)
        w_len = PNCIO_WriteContig(fh, buf, buf_view.size, fh->file_view.off[0]);
    else if (fh->file_system == PNCIO_UFS)
        w_len = PNCIO_GEN_Write_indep(fh, buf, buf_view);
    else if (fh->file_system == PNCIO_LUSTRE) {
        if (fh->hints->romio_ds_write == PNCIO_HINT_DISABLE)
            w_len = PNCIO_GEN_Write_indep(fh, buf, buf_view);
        else
            w_len = PNCIO_LUSTRE_WriteStrided(fh, buf, buf_view);
    }
    else
        w_len = NC_EFSTYPE;

    /* reset fileview, as PnetCDF never reuses a fileview */
    fh->file_view.off = NULL;
    fh->file_view.len = NULL;
    fh->file_view.size = 0;
    fh->file_view.count = 0;

    return w_len; /* when w_len < 0, it is an NetCDF error code */
}

/*----< PNCIO_File_write_at_all() >------------------------------------------*/
/* This is a collective call. */
MPI_Offset PNCIO_File_write_at_all(PNCIO_File *fh,
                                   MPI_Offset  offset,
                                   const void *buf,
                                   PNCIO_View  buf_view)
{
    int err=NC_NOERR;
    MPI_Offset w_len;

#ifdef PNETCDF_DEBUG
    assert(fh != NULL);
#endif

    if (buf_view.size < 0)
        err = NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY && err == NC_NOERR)
        err = NC_EPERM;

#ifdef PNETCDF_DEBUG
    if (offset > 0)
        assert(fh->file_view.off == NULL &&
               fh->file_view.len == NULL &&
               fh->file_view.count == 0);
#endif

    if (fh->file_view.off == NULL)
        /* This is when calling this subroutione without set fileview first.
         * We store offset into fh->file_view.off.
         */
        fh->file_view.off = &offset;

    if (fh->file_system == PNCIO_LUSTRE)
        w_len = PNCIO_LUSTRE_WriteStridedColl(fh, buf, buf_view);
    else if (fh->file_system == PNCIO_UFS)
        w_len = PNCIO_GEN_WriteStridedColl(fh, buf, buf_view);
    else
        err = NC_EFSTYPE;

    /* reset fileview, as PnetCDF never reuses a fileview */
    fh->file_view.off = NULL;
    fh->file_view.len = NULL;
    fh->file_view.size = 0;
    fh->file_view.count = 0;

    return (err == NC_NOERR) ? w_len : err;
}


