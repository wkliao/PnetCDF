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

/*----< PNCIO_ReadContig() >--------------------------------------------------*/
MPI_Offset PNCIO_ReadContig(PNCIO_File *fd,
                            void       *buf,
                            MPI_Offset  r_size,
                            MPI_Offset  offset)
{
    ssize_t err = 0;
    size_t r_count;
    MPI_Offset bytes_xfered = 0;
    char *p;

// printf("%s at %d: %s pread offset=%lld r_size=%lld\n",__func__,__LINE__,fd->filename,offset,r_size);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    double timing = MPI_Wtime();
#endif
    p = (char *) buf;
    while (bytes_xfered < r_size) {
        r_count = r_size - bytes_xfered;
        err = pread(fd->fd_sys, p, r_count, offset + bytes_xfered);
        if (err == -1)
            goto ioerr;
        if (err == 0)
            break;
        bytes_xfered += err;
        p += err;
    }
#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    fd->read_timing[2] += MPI_Wtime() - timing;
#endif

ioerr:
    if (err == -1)
        bytes_xfered = ncmpii_error_posix2nc("pread");

    return bytes_xfered;
}

/*----< file_read() >--------------------------------------------------------*/
/* This is an independent call. */
static
MPI_Offset file_read(PNCIO_File      *fd,
                     MPI_Offset       offset,
                     void            *buf,
                     PNCIO_Flat_list  buf_view)
{
    int filetype_is_contig;
    MPI_Offset r_len=0;
MPI_Offset off=offset;

    if (buf_view.size == 0) /* zero-sized request */
        return NC_NOERR;

    /* when fd->filetype == MPI_DATATYPE_NULL, this is called from INA */
    if (fd->filetype == MPI_DATATYPE_NULL) {
        if (fd->flat_file.count == 0)
            /* the whole file is visible */
            filetype_is_contig = 1;
        else {
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count m;
#else
            size_t m;
#endif
            MPI_Offset scan_sum=0;
            filetype_is_contig = 0;
            for (m=0; m<fd->flat_file.count; m++) {
                scan_sum += fd->flat_file.blocklens[m];
                if (scan_sum > offset) {
                    if (scan_sum - offset >= buf_view.size) {
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
// printf("%s at %d: offset=%lld buf_view.size=%lld m=%lld scan_sum=%lld off=%lld filetype_is_contig=%d\n",__func__,__LINE__, offset,buf_view.size,m,scan_sum,off,filetype_is_contig);
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

// printf("%s at %d: flat_file.count=%lld buf_view.is_contig=%d filetype_is_contig=%d\n",__func__,__LINE__, fd->flat_file.count, buf_view.is_contig,filetype_is_contig);

    if (buf_view.is_contig && filetype_is_contig)
        r_len = PNCIO_ReadContig(fd, buf, buf_view.size, off);
    else
        r_len = PNCIO_GEN_ReadStrided(fd, buf, buf_view, offset);

    return r_len;
}

/*----< PNCIO_File_read_at() >------------------------------------------------*/
/* This is an independent call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
MPI_Offset PNCIO_File_read_at(PNCIO_File      *fh,
                              MPI_Offset       offset,
                              void            *buf,
                              PNCIO_Flat_list  buf_view)
{
    assert(fh != NULL);

    if (buf_view.size == 0) return NC_NOERR;

    if (buf_view.size < 0) return NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY) return NC_EPERM;

    return file_read(fh, offset, buf, buf_view);
}

/*----< PNCIO_File_read_at_all() >--------------------------------------------*/
/* This is a collective call.
 * offset is a position in the file relative to the current view, expressed as
 * a count of etypes.
 */
MPI_Offset PNCIO_File_read_at_all(PNCIO_File      *fh,
                                  MPI_Offset       offset,
                                  void            *buf,
                                  PNCIO_Flat_list  buf_view)
{
    int err=NC_NOERR;
    MPI_Offset r_len;

    assert(fh != NULL);

    if (buf_view.size < 0) err = NC_ENEGATIVECNT;

    /* PnetCDF has only 2 modes: read-only and read-write */
    // if (fh->access_mode & MPI_MODE_RDONLY && st == NC_NOERR) st = NC_EPERM;

    r_len = PNCIO_GEN_ReadStridedColl(fh, buf, buf_view, offset);

    return (err == NC_NOERR) ?  r_len :  err;
}

