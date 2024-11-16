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

