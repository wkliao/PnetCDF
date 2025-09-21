/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>   /* strdup() */
#include <assert.h>
#include <sys/errno.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "adio.h"

/*----< PNCIO_File_set_view() >-----------------------------------------------*/
/* For PnetCDF, this subroutine is an independent call, because PnetCDF only
 * use the followings.
 *   Argument etype is always MPI_BYTE.
 *   Argument datarep is always "native".
 *   Argument info is always MPI_INFO_NULL.
 */
int PNCIO_File_set_view(PNCIO_File    *fd,
                        MPI_Offset    disp,
                        MPI_Datatype  filetype,
                        MPI_Aint      npairs,
#ifdef HAVE_MPI_LARGE_COUNT
                        MPI_Count    *offsets,
                        MPI_Count    *lengths
#else
                        MPI_Offset   *offsets,
                        int          *lengths
#endif
)
{
    MPI_Aint i;
/* below 4 are no longer used */
// fd->flat_file.lb_idx    = -1;
// fd->flat_file.ub_idx    = -1;
// fd->flat_file.flag      = 0;
// fd->flat_file.refct     = 1;

assert(filetype == MPI_BYTE);
assert(disp == 0);
#if 1
    fd->flat_file.count     = npairs;
    fd->flat_file.indices   = offsets;
    fd->flat_file.blocklens = lengths;

    /* Size of fileview must be calculated here, as PnetCDF may coalesce the
     * offset-length pairs in order to make offsets sorted in a monotonically
     * non-decreasing order.
     */
    fd->flat_file.size = 0;
    for (i=0; i<npairs; i++) fd->flat_file.size += lengths[i];

fd->filetype = filetype;
fd->disp = 0;
fd->flat_file.is_contig = (npairs <= 1); /* is_contig is redundant to count <= 1*/

#else
    fd->filetype = filetype;
    fd->disp = 0;

    if (filetype == MPI_BYTE) {
assert(npairs==0);
assert(offsets==NULL);
assert(lengths==NULL);
        fd->flat_file.count     = 0;
        fd->flat_file.indices   = NULL;
        fd->flat_file.blocklens = NULL;
    }
    else if (filetype == MPI_DATATYPE_NULL) {
        fd->flat_file.count     = npairs;
        fd->flat_file.indices   = offsets;
        fd->flat_file.blocklens = lengths;
    }
    else
        /* there should be no other filetype */
        assert(0);
#endif

    return NC_NOERR;
}

#if 0
/*----< check_type() >-------------------------------------------------------*/
/* check if filetype and etype are legal */
static
int check_type(PNCIO_Flatlist_node *flat_type,
               int access_mode,
               const char *type_kind)
{
    char err_msg[128];

    err_msg[0] = '\0';

    /* MPI standard requires the displacements of etype and filetype be
     * non-negative */
    if (flat_type->flag & PNCIO_TYPE_NEGATIVE) {
        sprintf(err_msg, "displacements of %s must be non-negative", type_kind);
        goto err_check;
    }

    /* MPI standard requires the displacements of etype and filetype be in a
     * monotonically nondecreasing order */
    if (flat_type->flag & PNCIO_TYPE_DECREASE) {
        sprintf(err_msg, "displacements of %s must be in a monotonically nondecreasing order",
                type_kind);
        goto err_check;
    }

    /* If the file is opened for writing, neither the etype nor the
     * filetype is permitted to contain overlapping regions.
     */
    if (((access_mode & MPI_MODE_WRONLY) || (access_mode & MPI_MODE_RDWR)) &&
        (flat_type->flag & PNCIO_TYPE_OVERLAP)) {
        sprintf(err_msg, "%s is not permitted to contain overlapping regions", type_kind);
        goto err_check;
    }

    return NC_NOERR;

err_check:
    return ncmpii_error_mpi2nc(MPI_ERR_IO, err_msg);
}

/*----< PNCIO_File_set_view() >-----------------------------------------------*/
/* For PnetCDF, this subroutine can become an independent call, because PnetCDF
 * only use the followings.
 *   Argument etype is always MPI_BYTE.
 *   Argument datarep is always "native".
 *   Argument info is always MPI_INFO_NULL.
 *
 * This subroutine can be an independent call, because there is no need to
 * check the consistency of any of the above arguments among all processes.
 */
int PNCIO_File_set_view(PNCIO_File    *fd,
                        MPI_Offset    disp,
                        MPI_Datatype  filetype,
                        MPI_Aint      npairs,
#ifdef HAVE_MPI_LARGE_COUNT
                        MPI_Count    *offsets,
                        MPI_Count    *lengths
#else
                        MPI_Offset   *offsets,
                        int          *lengths
#endif
)
{
    int is_predef, err=NC_NOERR, filetype_is_contig;
    MPI_Datatype copy_filetype;

    fd->flat_file.count = 0;

    if (filetype == MPI_DATATYPE_NULL) {
// if (npairs==0)printf("%s at %d: npairs=0\n",__func__,__LINE__);

        /* This is called from intra_node_aggregation() */
        // fd->flat_file = NCI_Malloc(sizeof(PNCIO_Flatlist_node));

        fd->flat_file.indices   = offsets;
        fd->flat_file.blocklens = lengths;
        fd->flat_file.count     = npairs;
        fd->flat_file.lb_idx    = -1;
        fd->flat_file.ub_idx    = -1;
        fd->flat_file.flag      = 0;
        fd->flat_file.refct     = 1;

// printf("%s at %d: flat_file count=%ld offsets=%lld lengths=%lld\n",__func__,__LINE__,npairs,offsets[0],lengths[0]);
        /* mark this request comes from intra_node_aggregation() */
        fd->filetype = MPI_DATATYPE_NULL;

        /* In addition, when called from INA subroutines, the passed-in offsets
         * and lengths are not based on MPI-IO fileview. They are flattened
         * byte offsets and sizes.
         */
        fd->disp = 0;
        return NC_NOERR;
    }

    if ((disp < 0) && (disp != MPI_DISPLACEMENT_CURRENT))
        return ncmpii_error_mpi2nc(MPI_ERR_ARG, "PNCIO_File_set_view, disp");

    /* When info is MPI_INFO_NULL, PNCIO_File_SetInfo() is an independent call.
     * Otherwise, it is collective, because it checks hint consistency.
     * PnetCDF always uses MPI_INFO_NULL when setting file view.
    err = PNCIO_File_SetInfo(fd, info);
    if (err != NC_NOERR)
        return err;
     */

    /* fd->flat_file should have been freed at callback of MPI_Type_free() */
    // fd->flat_file = NULL;

    /* free fileview if set previously */
    if (fd->filetype != MPI_BYTE && fd->filetype != MPI_DATATYPE_NULL)
        PNCIO_Type_dispose(&fd->filetype);

    if (filetype != MPI_BYTE && filetype != MPI_DATATYPE_NULL)
        PNCIO_Type_ispredef(filetype, &is_predef);
    else
        is_predef = 1;

    if (is_predef) {
        fd->filetype = filetype;
        filetype_is_contig = 1;
    } else {
        MPI_Type_dup(filetype, &copy_filetype);
        MPI_Type_commit(&copy_filetype);

        fd->filetype = copy_filetype;
        PNCIO_Datatype_iscontig(fd->filetype, &filetype_is_contig);

        /* check filetype only if it is not a predefined MPI datatype */
/*
        fd->flat_file = PNCIO_Flatten_and_find(fd->filetype);
*/
        PNCIO_Flatlist_node *tmp = PNCIO_Flatten_and_find(fd->filetype);
        fd->flat_file = *tmp;

        err = check_type(&fd->flat_file, fd->access_mode, "filetype");
        if (err != NC_NOERR)
            return err;
    }

    /* file displacement is an absolute byte position relative to the beginning
     * of a file. The displacement defines the location where a view begins.
     */
    fd->disp = disp;

    return err;
}
#endif

