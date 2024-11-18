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

#include "pnc_lustre.h"

/*----< PNC_ReadContig() >---------------------------------------------------*/
int PNC_ReadContig(ADIO_File     fd,
                   void   *buf,
                   MPI_Aint      count,
                   MPI_Datatype  bufType,
                   int           file_ptr_type,
                   ADIO_Offset   offset,
                   ADIO_Status  *status)
{
    ssize_t err = 0;
    size_t r_count;
    ADIO_Offset off, len, bytes_xfered = 0;
    char *p;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufType_size;
    MPI_Type_size_c(bufType, &bufType_size);
#else
    int bufType_size;
    MPI_Type_size(bufType, &bufType_size);
#endif
    len = count * bufType_size;

    if (file_ptr_type == ADIO_INDIVIDUAL)
        off = fd->fp_ind; /* offset is ignored */
    else /* ADIO_EXPLICIT_OFFSET */
        off = fd->disp + fd->etype_size * offset;

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

    if (file_ptr_type == ADIO_INDIVIDUAL)
        fd->fp_ind += bytes_xfered;
    /* if ADIO_EXPLICIT_OFFSET, do not update file pointer */

ioerr:
    if (err == -1)
        return ncmpii_error_posix2nc("pread");

    if (status)
        MPI_Status_set_elements(status, MPI_BYTE, bytes_xfered);

    return NC_NOERR;
}

/*----< file_read() >--------------------------------------------------------*/
/* This is an independent call. */
static
int file_read(PNC_File      fd,
              MPI_Offset    offset,
              void         *buf,
              int           count,
              MPI_Datatype  bufType,
              int           file_ptr_type,
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

    PNC_Datatype_iscontig(bufType, &buftype_is_contig);
    PNC_Datatype_iscontig(fd->filetype, &filetype_is_contig);

    if (buftype_is_contig && filetype_is_contig) {
        MPI_Aint wcount = (MPI_Aint)count * bufType_size;
        err = PNC_ReadContig(fd, buf, wcount, MPI_BYTE, file_ptr_type,
                             offset, status);
    }
    else {
        ADIO_ReadStrided(fd, buf, count, bufType, file_ptr_type,
                         offset, status, &err);
        if (err != MPI_SUCCESS)
            err = ncmpii_error_mpi2nc(err, __func__);
        else
            err = NC_NOERR;
    }
    return err;
}

/*----< PNC_File_read_at() >-------------------------------------------------*/
/* This is an independent call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int PNC_File_read_at(PNC_File      fh,
                     MPI_Offset    offset,
                     void         *buf,
                     int           count,
                     MPI_Datatype  bufType,
                     MPI_Status   *status)
{
    int err = NC_NOERR;

    if (count == 0) return NC_NOERR;

    if (count < 0) return NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        return NC_EPERM;

    err = file_read(fh, offset, buf, count, bufType, ADIO_EXPLICIT_OFFSET,
                     status);

    return err;
}

/*----< PNC_File_read() >----------------------------------------------------*/
/* This is an independent call. */
int PNC_File_read(PNC_File      fh,
                  void         *buf,
                  int           count,
                  MPI_Datatype  bufType,
                  MPI_Status   *status)
{
    int err = NC_NOERR;

    if (count == 0) return NC_NOERR;

    if (count < 0) return NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        return NC_EPERM;

    err = file_read(fh, 0, buf, count, bufType, ADIO_INDIVIDUAL, status);

    return err;
}

/*----< PNC_File_read_at_all() >---------------------------------------------*/
/* This is a collective call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int PNC_File_read_at_all(PNC_File      fh,
                         MPI_Offset    offset,
                         void         *buf,
                         int           count,
                         MPI_Datatype  bufType,
                         MPI_Status   *status)
{
    int err, st=NC_NOERR;

    if (count < 0) st = NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY && st == NC_NOERR)
        st = NC_EPERM;

    ADIO_ReadStridedColl(fh, buf, count, bufType, ADIO_EXPLICIT_OFFSET,
                         offset, status, &err);
    if (err != MPI_SUCCESS && st == NC_NOERR)
        st = ncmpii_error_mpi2nc(err, __func__);

    return st;
}

/*----< PNC_File_read_all() >-----------------------------------------------*/
/* This is a collective call. */
int PNC_File_read_all(PNC_File      fh,
                      void         *buf,
                      int           count,
                      MPI_Datatype  bufType,
                      MPI_Status   *status)
{
    int err, st=NC_NOERR;

    if (count < 0) st = NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY && st == NC_NOERR)
        st = NC_EPERM;

    ADIO_ReadStridedColl(fh, buf, count, bufType, ADIO_INDIVIDUAL,
                         0, status, &err);
    if (err != MPI_SUCCESS && st == NC_NOERR)
        st = ncmpii_error_mpi2nc(err, __func__);

    return st;
}

