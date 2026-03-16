/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <fcntl.h> /* F_GETLK64, F_SETLK64, F_SETLKW64,
                      F_RDLCK, F_WRLCK, F_UNLCK */

#include <pncio.h>

static
const char *GEN_flock_cmd_to_string(int cmd)
{
    switch (cmd) {
#ifdef F_GETLK64
        case F_GETLK64:
            return "F_GETLK64";
#else
        case F_GETLK:
            return "F_GETLK";
#endif
#ifdef F_SETLK64
        case F_SETLK64:
            return "F_SETLK64";
#else
        case F_SETLK:
            return "F_SETLK";
#endif
#ifdef F_SETLKW64
        case F_SETLKW64:
            return "F_SETLKW64";
#else
        case F_SETLKW:
            return "F_SETLKW";
#endif
        default:
            return "UNEXPECTED";
    }
}

static
const char *GEN_flock_type_to_string(int type)
{
    switch (type) {
        case F_RDLCK:
            return "F_RDLCK";
        case F_WRLCK:
            return "F_WRLCK";
        case F_UNLCK:
            return "F_UNLOCK";
        default:
            return "UNEXPECTED";
    }
}

int PNCIO_GEN_SetLock(PNCIO_File *fh,
                      int         cmd,
                      int         type,
                      MPI_Offset  offset,
                      int         whence,
                      MPI_Offset  len)
{
    int err, err_count = 0, sav_errno;
    struct flock lock;
    FDTYPE fd_sys = fh->fd_sys;

    if (len == 0)
        return MPI_SUCCESS;

    /* Depending on the compiler flags and options, struct flock may not be
     * defined with types that are the same size as MPI_Offsets.
     */

    /* FIXME: This is a temporary hack until we use flock64 where available. It
     * also doesn't fix the broken Solaris header sys/types.h header file,
     * which declares off_t as a UNION ! Configure tests to see if the off64_t
     * is a union if large file support is requested; if so, it does not select
     * large file support.
     */
#ifdef NEEDS_INT_CAST_WITH_FLOCK
    lock.l_type = type;
    lock.l_start = (int) offset;
    lock.l_whence = whence;
    lock.l_len = (int) len;
#else
    lock.l_type = type;
    lock.l_whence = whence;
    lock.l_start = offset;
    lock.l_len = len;
#endif

    /* save previous errno in case we recover from retryable errors */
    sav_errno = errno;

    errno = 0;
    do {
        err = fcntl(fd_sys, cmd, &lock);
    } while (err && (errno == EINTR || (errno == EINPROGRESS &&
                                        ++err_count < 10000)));

    if (err && errno != EBADF) {
        fprintf(stderr,
                "File locking failedin %s(fd_sys %X, cmd %s/%X, type %s/%X, whence %X) with return value %X and errno %X.\n"
                "- If the file system is NFS, you need to use NFS version 3, ensure that the lockd daemon is running on all the machines, and mount the directory with the 'noac' option (no attribute caching).\n"
                "- If the file system is LUSTRE, ensure that the directory is mounted with the 'flock' option.\n",
                __func__, fd_sys, GEN_flock_cmd_to_string(cmd), cmd,
                GEN_flock_type_to_string(type), type, whence, err, errno);
        perror("PNCIO_GEN_SetLock:");
        fprintf(stderr, "PNCIO_GEN_SetLock:offset %lld, length %lld\n",
                offset, len);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!err)
        /* report fcntl failure errno's (EBADF), otherwise restore previous
         * errno in case we recovered from retryable errors
         */
        errno = sav_errno;

    return (err == 0) ? MPI_SUCCESS : MPI_ERR_UNKNOWN;
}

int PNCIO_GEN_SetLock64(PNCIO_File *fh,
                        int         cmd,
                        int         type,
                        MPI_Offset  offset,
                        int         whence,
                        MPI_Offset  len)
{
    int err;
#ifdef _LARGEFILE64_SOURCE
    struct flock64 lock;
#else
    struct flock lock;
#endif
    FDTYPE fd_sys = fh->fd_sys;

    if (len == 0)
        return MPI_SUCCESS;

    lock.l_type = type;
    lock.l_start = offset;
    lock.l_whence = whence;
    lock.l_len = len;

    do {
        err = fcntl(fd_sys, cmd, &lock);
    } while (err && errno == EINTR);

    if (err && errno != EBADF) {
        fprintf(stderr,
                "File locking failed in %s(fd_sys %X,cmd %s/%X,type %s/%X,whence %X) with return value %X and errno %X.\n"
                "If the file system is NFS, you need to use NFS version 3, ensure that the lockd daemon is running on all the machines, and mount the directory with the 'noac' option (no attribute caching).\n",
                __func__, fd_sys, GEN_flock_cmd_to_string(cmd), cmd,
                GEN_flock_type_to_string(type), type, whence, err, errno);
        perror("PNCIO_GEN_SetLock64:");
        fprintf(stderr, "PNCIO_GEN_SetLock:offset %lld, length %lld\n",
                offset, len);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    return (err == 0) ? MPI_SUCCESS : MPI_ERR_UNKNOWN;
}

