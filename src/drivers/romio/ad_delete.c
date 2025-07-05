/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/errno.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h> /* unlink() */
#endif

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "adio.h"

/*----< ADIO_File_delete() >-------------------------------------------------*/
int ADIO_File_delete(const char *filename)
{
    int err = NC_NOERR;
    char *path = ncmpii_remove_file_system_type_prefix(filename);

    err = unlink(path);
    if (err != 0)
        err = ncmpii_error_posix2nc("unlink");

    return err;
}

