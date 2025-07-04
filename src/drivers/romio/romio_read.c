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
#include <unistd.h>   /* pread() */

#include <mpi.h>

#include "adio.h"

/*----< ADIO_ReadContig() >--------------------------------------------------*/
int ADIO_ReadContig(ADIO_File     fd,
                    void         *buf,
                    MPI_Aint      count,
                    MPI_Datatype  bufType,
                    ADIO_Offset   offset,
                    ADIO_Status  *status,
                    int          *error_code)
{
    ssize_t err = 0;
    size_t r_count;
    ADIO_Offset off, len, bytes_xfered = 0;
    char *p;

    if (error_code != NULL) *error_code = MPI_SUCCESS;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufType_size;
    MPI_Type_size_c(bufType, &bufType_size);
#else
    int bufType_size;
    MPI_Type_size(bufType, &bufType_size);
#endif
    len = count * bufType_size;

    off = offset;
    // off = fd->disp + fd->etype_size * offset;

// int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank); printf("%s at %d: pread off=%lld len=%lld\n",__func__,__LINE__,off,len);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
    p = (char *) buf;
    while (bytes_xfered < len) {
        r_count = len - bytes_xfered;
        err = pread(fd->fd_sys, p, r_count, off + bytes_xfered);
        if (err == -1)
            goto ioerr;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->coll_read[2] += MPI_Wtime() - timing;
#endif

ioerr:
    if (err == -1) {
        if (error_code != NULL) *error_code = MPI_ERR_IO;
        return ncmpii_error_posix2nc("pread");
    }

    if (status)
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
        MPI_Status_set_elements_x(status, MPI_BYTE, bytes_xfered);
#else
        MPI_Status_set_elements(status, MPI_BYTE, bytes_xfered);
#endif

    return NC_NOERR;
}

/*----< file_read() >--------------------------------------------------------*/
/* This is an independent call. */
static
int file_read(ADIO_File     fd,
              MPI_Offset    offset,
              void         *buf,
              int           count,
              MPI_Datatype  bufType,
              MPI_Status   *status)
{
    int err=NC_NOERR, buftype_is_contig, filetype_is_contig;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufType_size;
    MPI_Type_size_c(bufType, &bufType_size);
#else
    int bufType_size;
    MPI_Type_size(bufType, &bufType_size);
#endif
    if (bufType_size == 0) return NC_NOERR;

    ADIOI_Datatype_iscontig(bufType, &buftype_is_contig);

    /* when fd->filetype == MPI_DATATYPE_NULL, this is called from INA */
    if (fd->filetype == MPI_DATATYPE_NULL && fd->flat_file != NULL) {
        filetype_is_contig = (fd->flat_file->count <= 1);
        if (fd->flat_file->count > 0)
            offset = fd->flat_file->indices[0];
    }
    else if (fd->filetype == MPI_BYTE)
        filetype_is_contig = 1;
    else
        ADIOI_Datatype_iscontig(fd->filetype, &filetype_is_contig);

    if (buftype_is_contig && filetype_is_contig) {
        MPI_Aint rcount = (MPI_Aint)count * bufType_size;
        err = ADIO_ReadContig(fd, buf, rcount, MPI_BYTE, offset, status, NULL);
    }
    else {
        ADIO_ReadStrided(fd, buf, count, bufType, offset, status, &err);
        if (err != MPI_SUCCESS)
            err = ncmpii_error_mpi2nc(err, __func__);
        else
            err = NC_NOERR;
    }
    return err;
}

/*----< ADIO_File_read_at() >------------------------------------------------*/
/* This is an independent call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int ADIO_File_read_at(ADIO_File     fh,
                      MPI_Offset    offset,
                      void         *buf,
                      int           count,
                      MPI_Datatype  bufType,
                      MPI_Status   *status)
{
    int err = NC_NOERR;

    ADIOI_Assert(fh != NULL);

    if (count == 0) return NC_NOERR;

    if (count < 0) return NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY) return NC_EPERM;

    err = file_read(fh, offset, buf, count, bufType, status);

    return err;
}

/*----< ADIO_File_read() >---------------------------------------------------*/
/* This is an independent call. Note PnetCDF never calls this subroutine. */
int ADIO_File_read(ADIO_File     fh,
                   void         *buf,
                   int           count,
                   MPI_Datatype  bufType,
                   MPI_Status   *status)
{
    int err = NC_NOERR;

    ADIOI_Assert(fh != NULL);

    if (count == 0) return NC_NOERR;

    if (count < 0) return NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY) return NC_EPERM;

    err = file_read(fh, 0, buf, count, bufType, status);

    return err;
}

/*----< ADIO_File_read_at_all() >--------------------------------------------*/
/* This is a collective call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int ADIO_File_read_at_all(ADIO_File     fh,
                          MPI_Offset    offset,
                          void         *buf,
                          int           count,
                          MPI_Datatype  bufType,
                          MPI_Status   *status)
{
    int err, st=NC_NOERR;

    ADIOI_Assert(fh != NULL);

    if (count < 0) st = NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY && st == NC_NOERR) st = NC_EPERM;

    ADIO_ReadStridedColl(fh, buf, count, bufType, offset, status, &err);
    if (err != MPI_SUCCESS && st == NC_NOERR)
        st = ncmpii_error_mpi2nc(err, __func__);

    return st;
}

/*----< ADIO_File_read_all() >----------------------------------------------*/
/* This is a collective call. Note PnetCDF never calls this subroutine. */
int ADIO_File_read_all(ADIO_File     fh,
                       void         *buf,
                       int           count,
                       MPI_Datatype  bufType,
                       MPI_Status   *status)
{
    int err, st=NC_NOERR;

    ADIOI_Assert(fh != NULL);

    if (count < 0) st = NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY && st == NC_NOERR) st = NC_EPERM;

    ADIO_ReadStridedColl(fh, buf, count, bufType, 0, status, &err);
    if (err != MPI_SUCCESS && st == NC_NOERR)
        st = ncmpii_error_mpi2nc(err, __func__);

    return st;
}

