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

#include "pnc_lustre.h"

int PNC_WriteContig(ADIO_File     fd,
                    const void   *buf,
                    MPI_Aint      count,
                    MPI_Datatype  bufType,
                    int           file_ptr_type,
                    ADIO_Offset   offset,
                    ADIO_Status  *status)
{
    ssize_t err = 0;
    size_t w_count;
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

#ifdef WKL_DEBUG
{
int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
static int first_stripe_id=-1;
ADIO_Offset stripe_id = (off / fd->hints->striping_unit) % fd->hints->  cb_nodes;
    if (first_stripe_id == -1) {
        first_stripe_id = stripe_id;
       printf("%2d First: %s file %s pwrite off=%lld len=%lld stripe %d\n"  ,rank,__func__,fd->filename,off,len,first_stripe_id);
    }
    else if (stripe_id != first_stripe_id) printf("%2d Error: %s pwrite offse  t=%lld len=%lld not same stripe %d\n",rank,__func__,off,len,first_stripe_id  );
}
else
    printf("%2d %s line %d pread off=%lld len=%lld\n",rank,__func__,__LINE__,off,len);
#endif

double tt = MPI_Wtime();
    p = (char *) buf;
    while (bytes_xfered < len) {
        w_count = len - bytes_xfered;
        err = pwrite(fd->fd_sys, p, w_count, off + bytes_xfered);
        if (err == -1)
            goto ioerr;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
fd->lustre_write_metrics[1] += MPI_Wtime() - tt;

    if (file_ptr_type == ADIO_INDIVIDUAL)
        fd->fp_ind += bytes_xfered;
    /* if ADIO_EXPLICIT_OFFSET, do not update file pointer */

ioerr:
    if (err == -1)
        return ncmpii_error_posix2nc("pwrite");

    if (status)
        MPI_Status_set_elements(status, MPI_BYTE, bytes_xfered);

    return NC_NOERR;
}

/*----< file_write() >---------------------------------------------------*/
/* This is an independent call. */
static
int file_write(PNC_File      fd,
               MPI_Offset    offset,
               const void   *buf,
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
    PNC_Datatype_iscontig(fd->ftype, &filetype_is_contig);

    if (buftype_is_contig && filetype_is_contig) {
        MPI_Aint wcount = (MPI_Aint)count * bufType_size;
        err = PNC_WriteContig(fd, buf, wcount, MPI_BYTE, file_ptr_type,
                              offset, status);
    }
    else
        err = PNC_WriteStrided(fd, buf, count, bufType, file_ptr_type,
                               offset, status);
    return err;
}

/*----< PNC_File_write_at() >------------------------------------------------*/
/* This is an independent call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int PNC_File_write_at(PNC_File      fh,
                      MPI_Offset    offset,
                      const void   *buf,
                      int           count,
                      MPI_Datatype  bufType,
                      MPI_Status   *status)
{
    int err = NC_NOERR;

    if (count == 0) return NC_NOERR;

    if (count < 0) return NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        return NC_EPERM;

    err = file_write(fh, offset, buf, count, bufType, ADIO_EXPLICIT_OFFSET,
                     status);

    return err;
}

/*----< PNC_File_write() >---------------------------------------------------*/
/* This is an independent call. */
int PNC_File_write(PNC_File      fh,
                   const void   *buf,
                   int           count,
                   MPI_Datatype  bufType,
                   MPI_Status   *status)
{
    int err = NC_NOERR;

    if (count == 0) return NC_NOERR;

    if (count < 0) return NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        return NC_EPERM;

    err = file_write(fh, 0, buf, count, bufType, ADIO_INDIVIDUAL, status);

    return err;
}

/*----< PNC_File_write_at_all() >--------------------------------------------*/
/* This is a collective call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int PNC_File_write_at_all(PNC_File      fh,
                          MPI_Offset    offset,
                          const void   *buf,
                          int           count,
                          MPI_Datatype  bufType,
                          MPI_Status   *status)
{
    int err, st=NC_NOERR;

    if (count < 0) st = NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        err = NC_EPERM;
    if (st == NC_NOERR) st = err;

/*
    err = PNC_WriteStridedColl(fh, buf, count, bufType, ADIO_EXPLICIT_OFFSET,
                               offset, status);
    if (st == NC_NOERR) st = err;
 */

    return st;
}

/*----< PNC_File_write_all() >---------------------------------------------------*/
/* This is a collective call. */
int PNC_File_write_all(PNC_File      fh,
                       const void   *buf,
                       int           count,
                       MPI_Datatype  bufType,
                       MPI_Status   *status)
{
    int err, st=NC_NOERR;

    if (count < 0) st = NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        err = NC_EPERM;
    if (st == NC_NOERR) st = err;

/*
    err = PNC_WriteStridedColl(fh, buf, count, bufType, ADIO_INDIVIDUAL,
                               0, status);
    if (st == NC_NOERR) st = err;
*/

    return st;
}
