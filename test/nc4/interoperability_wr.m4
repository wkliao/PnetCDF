dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */
dnl
include(`foreach.m4')dnl
include(`utils.m4')dnl
dnl
define(`upcase', `translit(`$*', `a-z', `A-Z')')dnl
dnl
define(`PNC_WR',dnl
`dnl
int pnc_wr_$2(int rank, int ncid, int* vid, int *did){
    int i, err, nerrs = 0;
    MPI_Offset start[2], count[2], stride[2];
    $3 buf[5];

    /* Setup buffer */
    for(i = 0; i < 5; i++){
        buf[i] = ($3)rank;
    }

    /* Define variables */
    err = ncmpi_redef(ncid); CHECK_ERR
    err = ncmpi_def_var(ncid, "$1", $1, 2, did, vid); CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    /* Write data */
    start[0] = rank;
    start[1] = 0;
    count[0] = 1;
    count[1] = 2;
    stride[0] = 1;
    stride[1] = 2;
    err = ncmpi_put_vars_$2_all(ncid, *vid, start, count, stride, buf); CHECK_ERR
    start[1] = 3;
    err = ncmpi_put_vara_$2_all(ncid, *vid, start, count, buf + 2); CHECK_ERR
    start[1] = 1;
    err = ncmpi_put_var1_$2_all(ncid, *vid, start, buf + 4); CHECK_ERR

    return err;
}
')dnl
dnl
define(`NC_RD',dnl
`dnl
int nc_rd_$2(int rank, int ncid, int vid, int *did){
    int i, err, nerrs = 0;
    size_t start[2], count[2];
    ptrdiff_t stride[2];
    $3 buf[5], ans[5];

    /* Setup ans */
    for(i = 0; i < 5; i++){
        ans[i] = ($3)rank;
    }

    /* Read data */
    start[0] = rank;
    start[1] = 0;
    count[0] = 1;
    count[1] = 2;
    stride[0] = 1;
    stride[1] = 2;
    err = nc_get_vars_$2(ncid, vid, start, count, stride, buf); CHECK_ERR
    start[1] = 3;
    err = nc_get_vara_$2(ncid, vid, start, count, buf + 2); CHECK_ERR
    start[1] = 1;
    err = nc_get_var1_$2(ncid, vid, start, buf + 4); CHECK_ERR

    /* Compare result */
    for(i = 0; i < 5; i++){
        if (buf[i] != ans[i]){
            printf("Error at line %d in %s: expecting buf[%d] = %$1 but got %$1\n", __LINE__, __FILE__, i, ans[i], buf[i]);
        }
    }

    return err;
}
')dnl
dnl
define(`CALL_PNC_WR',dnl
`dnl
    err = pnc_wr_$2(rank, ncid, vid + 0, did); CHECK_ERR
')dnl
dnl
define(`CALL_NC_RD',dnl
`dnl
    err = nc_rd_$2(rank, ncid, vid[0], did); CHECK_ERR
')dnl
dnl
/*
   This program tests the interoperability feature of PnetCDF and NetCDF
   The data is writen using NetCDF and read by PnetCDF
   $Id$
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <pnetcdf.h>
#include <testutils.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "interoperability_wr.nc"


foreach(`dt', (`(`NC_BYTE', `schar', `signed char')', dnl
               `(`NC_UBYTE', `uchar', `unsigned char')', dnl
               `(`NC_SHORT', `short', `short')', dnl
               `(`NC_USHORT', `ushort', `unsigned short')', dnl
               `(`NC_INT', `int', `int')', dnl
               `(`NC_UINT', `uint', `unsigned int')', dnl
               `(`NC_FLOAT', `float', `float')', dnl
               `(`NC_DOUBLE', `double', `double')', dnl
               `(`NC_INT64', `longlong', `long long')', dnl
               `(`NC_UINT64', `ulonglong', `unsigned long long')', dnl
               ), `PNC_WR(translit(dt, `()'))')dnl

foreach(`dt', (`(`c', `schar', `signed char')', dnl
               `(`c', `uchar', `unsigned char')', dnl
               `(`d', `short', `short')', dnl
               `(`u', `ushort', `unsigned short')', dnl
               `(`d', `int', `int')', dnl
               `(`u', `uint', `unsigned int')', dnl
               `(`f', `float', `float')', dnl
               `(`lf', `double', `double')', dnl
               `(`lld', `longlong', `long long')', dnl
               `(`llu', `ulonglong', `unsigned long long')', dnl
               ), `NC_RD(translit(dt, `()'))')dnl


int main(int argc, char **argv) {
    int i, err, nerrs = 0;
    int rank, np;
    int ncid, pncid;
    int vid[32];
    int did[2];
    char *filename = FILE_NAME;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) filename = argv[1];

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for interoperability file", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    /* Write with NetCDF */
    /* Create the file */
    err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER | NC_NETCDF4, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* Define dimensions */
    err =  ncmpi_def_dim(ncid, "x", np, did + 0); CHECK_ERR
    err =  ncmpi_def_dim(ncid, "y", 5, did + 1); CHECK_ERR

    /* Data mode */
    err =  ncmpi_enddef(ncid); CHECK_ERR

    /* Write variables */
foreach(`dt', (`(`0', `schar', `char')', dnl
               `(`1', `uchar', `unsigned char')', dnl
               `(`2', `short', `short')', dnl
               `(`3', `ushort', `unsigned short')', dnl
               `(`4', `int', `int')', dnl
               `(`5', `uint', `unsigned int')', dnl
               `(`6', `float', `float')', dnl
               `(`7', `double', `double')', dnl
               `(`8', `longlong', `long long')', dnl
               `(`9', `ulonglong', `unsigned long long')', dnl
               ), `CALL_PNC_WR(translit(dt, `()'))')dnl

    /* Close file */
    ncmpi_close(ncid);

    /* Read with PnetCDF */
    /* Open the file */
    err = nc_open_par(filename, NC_MPIIO | NC_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid);
    CHECK_ERR

    /* Read variables */
foreach(`dt', (`(`0', `schar', `char')', dnl
               `(`1', `uchar', `unsigned char')', dnl
               `(`2', `short', `short')', dnl
               `(`3', `ushort', `unsigned short')', dnl
               `(`4', `int', `int')', dnl
               `(`5', `uint', `unsigned int')', dnl
               `(`6', `float', `float')', dnl
               `(`7', `double', `double')', dnl
               `(`8', `longlong', `long long')', dnl
               `(`9', `ulonglong', `unsigned long long')', dnl
               ), `CALL_NC_RD(translit(dt, `()'))')dnl

    /* Close file */
    nc_close(ncid);

    /* check if there is any malloc residue */
    MPI_Offset malloc_size, sum_size;
    err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

   MPI_Finalize();
   return (nerrs > 0);
}
