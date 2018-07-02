/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_def_dim()    : dispatcher->def_dim()
 * ncmpi_inq_dimid()  : dispatcher->inq_dimid()
 * ncmpi_inq_dim()    : dispatcher->inq_dim()
 * ncmpi_rename_dim() : dispatcher->rename_dim()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include <ncadio_driver.h>

int
ncadio_def_dim(void       *ncdp,
              const char *name,
              MPI_Offset  size,
              int        *dimidp)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadio_inq_dimid(void       *ncdp,
                const char *name,
                int        *dimid)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* adios dim has no name */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}

int
ncadio_inq_dim(void       *ncdp,
              int         dimid,
              char       *name,
              MPI_Offset *sizep)
{
    int err;
    int i;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Artificial dimid 
     * Assume that every variable use their own unique dim
     */
    for(i = 0; i < ncadp->fp->nvars; i++){
        if (dimid < ncadp->ndims[i]){
            break;
        }
        dimid -= ncadp->ndims[i];
    }

    if (i >= ncadp->fp->nvars || dimid < 0){
        DEBUG_RETURN_ERROR(NC_EINVAL);
    }

    if (name != NULL){
        // ADIOS dim has no name
        *name = '\0';
    }

    if (sizep != NULL){
        ADIOS_VARINFO * v = adios_inq_var_byid (ncadp->fp, i);
        *sizep = (MPI_Offset)(v->dims[dimid]);
    }

    return NC_NOERR;
}

int
ncadio_rename_dim(void       *ncdp,
                 int         dimid,
                 const char *newname)
{
    int err;
    NC_ad *ncadp = (NC_ad*)ncdp;

    /* Read only driver */
    DEBUG_RETURN_ERROR(NC_ENOTSUPPORT);

    return NC_NOERR;
}
