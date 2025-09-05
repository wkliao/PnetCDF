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

#include "adio.h"

#ifdef WKL_DEBUG
int first_ost_id;
#endif

/*----< PNCIO_WriteContig() >-------------------------------------------------*/
int PNCIO_WriteContig(PNCIO_File   *fd,
                      const void   *buf,
                      MPI_Offset    count,
                      MPI_Datatype  bufType,
                      MPI_Offset    offset,
                      MPI_Status   *status,
                      int          *error_code)
{
    ssize_t err = 0;
    size_t w_count;
    MPI_Offset off, len, bytes_xfered = 0;
    char *p;

    if (error_code != NULL) *error_code = MPI_SUCCESS;

    if (count == 0) return NC_NOERR;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufType_size;
    MPI_Type_size_c(bufType, &bufType_size);
#else
    int bufType_size;
    MPI_Type_size(bufType, &bufType_size);
#endif
    len = count * bufType_size;

    off = offset;

// printf("%s at %d: pwrite off=%lld len=%lld\n",__func__,__LINE__,off,len);
#ifdef WKL_DEBUG
int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);

MPI_Offset ost_id = (off / fd->hints->striping_unit) % fd->hints->striping_factor;
    if (first_ost_id == -1) {
        first_ost_id = ost_id;
        // printf("%2d %s file %s First pwrite off=%lld OST %d\n",rank,__func__,fd->filename,off,first_ost_id);
    }
    else if (ost_id != first_ost_id)
        printf("%2d Error: %s pwrite offset=%lld len=%lld ost_id=%lld not same 1st ost %d\n",rank,__func__,off,len,ost_id,first_ost_id);

printf("%s line %d: disp=%lld offset=%lld off=%lld count=%ld bufType_size=%d len=%lld\n",__func__,__LINE__,fd->disp,offset,off,count,bufType_size,len);

    printf("%2d %s line %d pwrite off=%lld len=%lld\n",rank,__func__,__LINE__,off,len);
#endif

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
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
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->write_timing[2] += MPI_Wtime() - timing;
#endif

ioerr:
    if (err == -1) {
        if (error_code != NULL) *error_code = MPI_ERR_IO;
        return ncmpii_error_posix2nc("pwrite");
    }

    if (status)
#ifdef HAVE_MPI_STATUS_SET_ELEMENTS_X
        MPI_Status_set_elements_x(status, MPI_BYTE, bytes_xfered);
#else
        MPI_Status_set_elements(status, MPI_BYTE, bytes_xfered);
#endif

    return NC_NOERR;
}

/*----< file_write() >-------------------------------------------------------*/
/* This is an independent call. */
static
int file_write(PNCIO_File    *fd,
               MPI_Offset    offset,
               const void   *buf,
               MPI_Offset    count,
               MPI_Datatype  bufType,
               MPI_Status   *status)
{
    int err=NC_NOERR, buftype_is_contig, filetype_is_contig;
    MPI_Offset off=offset;

#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count bufType_size;
    MPI_Type_size_c(bufType, &bufType_size);
#else
    int bufType_size;
    MPI_Type_size(bufType, &bufType_size);
#endif
    if (bufType_size == 0) return NC_NOERR;

    PNCIO_Datatype_iscontig(bufType, &buftype_is_contig);

/* PnetCDF always packs non-contiguous user buffer into a contiguous one in INA */
assert(buftype_is_contig == 1);

assert(fd->filetype == MPI_DATATYPE_NULL || fd->filetype == MPI_BYTE);

    /* when fd->filetype == MPI_DATATYPE_NULL, this is called from INA */
    if (fd->filetype == MPI_DATATYPE_NULL) {
        if (fd->flat_file.count == 0)
            /* the whole file is visible */
            filetype_is_contig = 1;
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count m, size;
            MPI_Type_size_c(bufType, &size);
            size *= count;
#else
            size_t m;
            int size;
            MPI_Type_size(bufType, &size);
            size *= count;
#endif
            MPI_Offset scan_sum=0;
            filetype_is_contig = 0;
            for (m=0; m<fd->flat_file.count; m++) {
                scan_sum += fd->flat_file.blocklens[m];
                if (scan_sum > offset) {
                    if (scan_sum - offset >= size) {
                        /* check if this request falls entirely in m's
                         * offset-length pair
                         */
                        filetype_is_contig = 1;
                        off = fd->flat_file.indices[m] + offset -
                              (scan_sum - fd->flat_file.blocklens[m]);
                    }
                    break;
                }
            }
// printf("%s at %d: offset=%lld size=%lld m=%lld scan_sum=%lld off=%lld filetype_is_contig=%d\n",__func__,__LINE__, offset,size,m,scan_sum,off,filetype_is_contig);
        }
#if 0
        else if (fd->flat_file.count == 1)
            filetype_is_contig = 1;
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count m;
#else
            size_t m;
#endif
            for (m=0; m<fd->flat_file.count; m++) {
                if (offset < fd->flat_file.indices[m] + fd->flat_file.blocklens[m])
                    break;
            }
            filetype_is_contig = (fd->flat_file.count - m == 1);
        }
#endif
    }
    else if (fd->filetype == MPI_BYTE)
        filetype_is_contig = 1;
    else
        PNCIO_Datatype_iscontig(fd->filetype, &filetype_is_contig);

/*
if (fd->flat_file.count == 0)
printf("%s at %d: fd->flat_file.count=0 filetype_is_contig=%d offset=%lld\n",__func__,__LINE__, filetype_is_contig,offset);
else if (fd->flat_file.count == 1)printf("%s at %d: fd->flat_file.count=%lld indices=%lld blocklens=%lld filetype_is_contig=%d offset=%lld count=%lld\n",__func__,__LINE__, fd->flat_file.count,fd->flat_file.indices[0],fd->flat_file.blocklens[0],filetype_is_contig,offset,count);
else if (fd->flat_file.count > 1)printf("%s at %d: fd->flat_file.count=%lld indices=%lld %lld blocklens=%lld %lld filetype_is_contig=%d offset=%lld count=%lld\n",__func__,__LINE__, fd->flat_file.count,fd->flat_file.indices[0],fd->flat_file.indices[1],fd->flat_file.blocklens[0],fd->flat_file.blocklens[1],filetype_is_contig,offset,count);
*/

    if (buftype_is_contig && filetype_is_contig) {
        if (fd->filetype == MPI_DATATYPE_NULL)
            offset = off;
        err = PNCIO_WriteContig(fd, buf, count * bufType_size, MPI_BYTE,
                                offset, status, NULL);
    }
    else if (fd->file_system == PNCIO_LUSTRE)
        PNCIO_LUSTRE_WriteStrided(fd, buf, count, bufType, offset, status, &err);
    else
        PNCIO_GEN_WriteStrided(fd, buf, count, bufType, offset, status, &err);
    if (err != MPI_SUCCESS)
        err = ncmpii_error_mpi2nc(err, __func__);
    else
        err = NC_NOERR;
    return err;
}

/*----< PNCIO_File_write_at() >-----------------------------------------------*/
/* This is an independent call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int PNCIO_File_write_at(PNCIO_File    *fh,
                        MPI_Offset    offset,
                        const void   *buf,
                        MPI_Offset    count,
                        MPI_Datatype  bufType,
                        MPI_Status   *status)
{
    int err = NC_NOERR;

    assert(fh != NULL);

    if (count == 0) return NC_NOERR;

    if (count < 0) return NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY)
        return NC_EPERM;

    err = file_write(fh, offset, buf, count, bufType, status);

    return err;
}

/*----< PNCIO_File_write_at_all() >-------------------------------------------*/
/* This is a collective call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
int PNCIO_File_write_at_all(PNCIO_File    *fh,
                            MPI_Offset    offset,
                            const void   *buf,
                            MPI_Offset    count,
                            MPI_Datatype  bufType,
                            MPI_Status   *status)
{
    int err, st=NC_NOERR;

    assert(fh != NULL);

    if (count < 0) st = NC_ENEGATIVECNT;

    if (fh->access_mode & MPI_MODE_RDONLY && st == NC_NOERR)
        st = NC_EPERM;

    if (fh->file_system == PNCIO_LUSTRE)
        PNCIO_LUSTRE_WriteStridedColl(fh, buf, count, bufType, offset, status,
                                      &err);
    else if (fh->file_system == PNCIO_UFS)
        PNCIO_GEN_WriteStridedColl(fh, buf, count, bufType, offset, status,
                                   &err);
    else
        return NC_EFSTYPE;

    if (err != MPI_SUCCESS && st == NC_NOERR)
        st = ncmpii_error_mpi2nc(err, "PNCIO_File_write_at_all");

    return st;
}


