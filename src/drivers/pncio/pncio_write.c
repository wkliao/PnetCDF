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

#ifdef WKL_DEBUG
int first_ost_id;
#endif

/*----< PNCIO_WriteContig() >-------------------------------------------------*/
MPI_Offset PNCIO_WriteContig(PNCIO_File *fd,
                             const void *buf,
                             MPI_Offset  w_size,
                             MPI_Offset  offset)
{
    ssize_t err = 0;
    size_t w_count;
    MPI_Offset bytes_xfered = 0;
    char *p;

    if (w_size == 0) return NC_NOERR;

// printf("%s at %d: pwrite offset=%lld w_size=%lld\n",__func__,__LINE__,offset,w_size);
#ifdef WKL_DEBUG
int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

MPI_Offset ost_id = (offset / fd->hints->striping_unit) % fd->hints->striping_factor;
    if (first_ost_id == -1) {
        first_ost_id = ost_id;
        // printf("%2d %s file %s First pwrite offset=%lld OST %d\n",rank,__func__,fd->filename,offset,first_ost_id);
    }
    else if (ost_id != first_ost_id)
        printf("%2d Error: %s pwrite offset=%lld w_size=%lld ost_id=%lld not same 1st ost %d\n",rank,__func__,offset,w_size,ost_id,first_ost_id);

printf("%s line %d: offset=%lld count=%ld bufType_size=%d w_size=%lld\n",__func__,__LINE__,offset,count,bufType_size,w_size);

    printf("%2d %s line %d pwrite offset=%lld w_size=%lld\n",rank,__func__,__LINE__,offset,w_size);
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
    p = (char *) buf;
    while (bytes_xfered < w_size) {
        w_count = w_size - bytes_xfered;
        err = pwrite(fd->fd_sys, p, w_count, offset + bytes_xfered);
        if (err == -1)
            goto ioerr;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->write_timing[2] += MPI_Wtime() - timing;
#endif

ioerr:
    if (err == -1)
        bytes_xfered = ncmpii_error_posix2nc("pwrite");

    return bytes_xfered;
}

/*----< file_write() >-------------------------------------------------------*/
/* This is an independent call. */
static
MPI_Offset file_write(PNCIO_File *fd,
                      const void *buf,
                      PNCIO_View  buf_view)
{
    MPI_Offset w_len;

    if (buf_view.size == 0) /* zero-sized request */
        return NC_NOERR;

    if (buf_view.count <= 1 && fd->file_view.count <= 1)
        w_len = PNCIO_WriteContig(fd, buf, buf_view.size, fd->file_view.off[0]);
    else if (fd->file_system == PNCIO_UFS)
        w_len = PNCIO_GEN_Write_indep(fd, buf, buf_view);
    else if (fd->file_system == PNCIO_LUSTRE) {
        if (fd->hints->romio_ds_write == PNCIO_HINT_DISABLE)
            w_len = PNCIO_GEN_Write_indep(fd, buf, buf_view);
        else
            w_len = PNCIO_LUSTRE_WriteStrided(fd, buf, buf_view);
    }
    else
        return NC_EFSTYPE;

    return w_len; /* when w_len < 0, it is an NetCDF error code */
}

/*----< PNCIO_File_write_at() >-----------------------------------------------*/
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
if (offset > 0) assert(fh->file_view.off == NULL && fh->file_view.len == NULL && fh->file_view.count == 0);
#endif

    if (fh->file_view.off == NULL)
        /* This is when calling this subroutione without set fileview first.
         * We store offset into fh->file_view.off.
         */
        fh->file_view.off = &offset;

    w_len = file_write(fh, buf, buf_view);
    /* when w_len < 0, it is an NetCDF error code */

    /* reset fileview, as PnetCDF never reuses a fileview */
    fh->file_view.off = NULL;
    fh->file_view.len = NULL;
    fh->file_view.size = 0;
    fh->file_view.count = 0;

    return w_len;
}

/*----< PNCIO_File_write_at_all() >-------------------------------------------*/
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
if (offset > 0) assert(fh->file_view.off == NULL && fh->file_view.len == NULL && fh->file_view.count == 0);
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


