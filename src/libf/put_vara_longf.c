/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by buildiface -infile=../lib/pnetcdf.h -deffile=defs
 * DO NOT EDIT
 */
#include "mpinetcdf_impl.h"


#ifdef F77_NAME_UPPER
#define nfmpi_put_vara_long_ NFMPI_PUT_VARA_LONG
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_put_vara_long_ nfmpi_put_vara_long__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_put_vara_long_ nfmpi_put_vara_long
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_put_vara_long_ ( int *v1, int *v2, int v3[], int v4[], long*v5 ){
    int ierr;
    ierr = ncmpi_put_vara_long( *v1, *v2, (const size_t *)(v3), (const size_t *)(v4), v5 );
    return ierr;
}
