/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/* This program tests if PnetCDF can detect file header inconsistency and
 * overwrite the inconsistent header with root's.
 * This program is designed to run on more than 2 MPI processes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h> /* basename() */
#include <mpi.h>
#include <pnetcdf.h>
#include <testutils.h>

#define EXP_ERR2(e, exp1, exp2) { \
    if (e != exp1 && e != exp2 && e != NC_EFILE) { \
        printf("Error at line %d in %s: expecting error code %s or %s but got %s\n", \
               __LINE__, __FILE__, ncmpi_strerrno(exp1), ncmpi_strerrno(exp2), ncmpi_strerrno(e)); \
        nerrs++; \
    } \
}

#define EXP_SAFE_ERR(expect) { \
    if (safe_mode) { \
        if (err != NC_EMULTIDEFINE && err != expect) { \
            printf("Error at line %d in %s: expecting error code NC_EMULTIDEFINE or %s but got %s\n", \
                   __LINE__, __FILE__, ncmpi_strerrno(expect), ncmpi_strerrno(err)); \
            nerrs++; \
        } \
    } \
    else if (rank > 0) { \
        if (err != expect) { \
            printf("Error at line %d in %s: expecting error code %s but got %s\n", \
                   __LINE__, __FILE__, ncmpi_strerrno(expect), ncmpi_strerrno(err)); \
            nerrs++; \
        } \
    } \
}

/*----< test_open_mode() >----------------------------------------------------*/
static
int test_open_mode(char *filename, int safe_mode)
{
    int err, rank, ncid, cmode, omode, nerrs=0;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);

    /* Test inconsistent cmode -----------------------------------------------*/
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;
    if (rank == 0) cmode = NC_CLOBBER;
    err = ncmpi_create(comm, filename, cmode, info, &ncid);
    if (safe_mode)
        /* all processes got the same error code */
        EXP_ERR(NC_EMULTIDEFINE_CMODE)
    else if (rank > 0)
        /* all processes except root got the same error code */
        EXP_ERR(NC_EMULTIDEFINE_CMODE)
    /* In either case, file is created with multi-defined-cmode error or not */
    err = ncmpi_close(ncid); CHECK_ERR


    /* Test inconsistent omode -----------------------------------------------*/
    omode = NC_WRITE;
    if (rank == 0) omode = NC_NOWRITE;
    err = ncmpi_open(comm, filename, omode, info, &ncid);
    if (safe_mode)
        /* all processes got the same error code */
        EXP_ERR(NC_EMULTIDEFINE_OMODE)
    else if (rank > 0)
        /* all processes except root got the same error code */
        EXP_ERR(NC_EMULTIDEFINE_OMODE)
    /* In either case, file is opened with multi-defined-omode error or not */
    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

/*----< test_dim() >----------------------------------------------------------*/
static
int test_dim(char *filename, int safe_mode)
{
    int err, rank, ncid, cmode, dimid1, dimid2, dimid3, nerrs=0;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;

    /* Test inconsistency on dimension names ---------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    if (rank == 0)
        err = ncmpi_def_dim(ncid, "y", 100, &dimid1);
    else
        err = ncmpi_def_dim(ncid, "xx", 100, &dimid1);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_DIM_NAME)
    else
        CHECK_ERR
    err = ncmpi_enddef(ncid); CHECK_ERR

    if (safe_mode) {
        /* no processes should be able to see dim "y" */
        err = ncmpi_inq_dimid(ncid, "y", &dimid2);
        EXP_SAFE_ERR(NC_EBADDIM)

        /* no process should be able to see dim "x" */
        err = ncmpi_inq_dimid(ncid, "xx", &dimid3);
        EXP_SAFE_ERR(NC_EBADDIM)
    }

    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistency on dimension size ----------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    if (rank == 0)
        err = ncmpi_def_dim(ncid, "x", 99, &dimid1);
    else
        err = ncmpi_def_dim(ncid, "x", 100, &dimid1);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_DIM_SIZE)
    else
        CHECK_ERR

    err = ncmpi_close(ncid); CHECK_ERR
    return nerrs;
}

/*----< test_attr() >---------------------------------------------------------*/
static
int test_attr(char *filename, int safe_mode)
{
    int err, rank, ncid, cmode, nerrs=0;
    char  gattr[128];
    int   int_attr;
    float flt_attr;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;

    /* Test inconsistent global attribute name -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    int_attr = 1;
    sprintf(gattr, "gattr_name.%d",rank);
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, gattr, NC_INT, 1, &int_attr);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_NAME)
    else
        CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent global attribute type -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 1, &int_attr);
    else
        err = ncmpi_put_att_float(ncid, NC_GLOBAL, "gatt", NC_FLOAT, 1, &flt_attr);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_TYPE)
    else
        CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    int intv[2]={1,2};
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 2, intv);
    else
        err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 1, intv);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_LEN)
    else
        CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    if (rank == 0) intv[1]=3;
    err = ncmpi_put_att_int(ncid, NC_GLOBAL, "gatt", NC_INT, 2, intv);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_VAL)
    else
        CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

/*----< test_var() >----------------------------------------------------------*/
static
int test_var(char *filename, int safe_mode)
{
    int err, rank, ncid, cmode, nerrs=0;
    int dimid[3], varid1, int_attr;
    float flt_attr;
    char name[128], var_attr[128];
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm comm=MPI_COMM_WORLD;

    MPI_Comm_rank(comm, &rank);
    cmode = NC_CLOBBER|NC_64BIT_OFFSET;

    /* Test inconsistent global attribute name -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); CHECK_ERR
    int_attr = 1;
    sprintf(var_attr, "var_attr_name.%d",rank);
    err = ncmpi_put_att_int(ncid, varid1, var_attr, NC_INT, 1, &int_attr);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_NAME)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent global attribute type -------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); CHECK_ERR
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 1, &int_attr);
    else
        err = ncmpi_put_att_float(ncid, varid1, "var_att", NC_FLOAT, 1, &flt_attr);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_TYPE)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); CHECK_ERR
    int intv[2]={1,2};
    if (rank == 0)
        err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 2, intv);
    else
        err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 1, intv);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_LEN)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent global attribute length -----------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    err = ncmpi_def_var(ncid, "var1", NC_INT, 1, dimid, &varid1); CHECK_ERR
    if (rank == 0) intv[1]=3;
    err = ncmpi_put_att_int(ncid, varid1, "var_att", NC_INT, 2, intv);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_ATTR_VAL)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent variable name ---------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    sprintf(name, "var.%d",rank);
    err = ncmpi_def_var(ncid, name, NC_INT, 1, dimid, &varid1);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_VAR_NAME)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent variable ndims --------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim0", 3, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", 2, &dimid[1]); CHECK_ERR
    if (rank == 0)
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid1);
    else
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 1, dimid, &varid1);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_VAR_NDIMS)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent variable type ---------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", NC_UNLIMITED, &dimid[0]); CHECK_ERR
    if (rank == 0)
        err = ncmpi_def_var(ncid, "var", NC_INT, 1, dimid, &varid1);
    else
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 1, dimid, &varid1);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_VAR_TYPE)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent variable length -------------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim0", 5, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim1", 4, &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "dim2", 3, &dimid[2]); CHECK_ERR
    if (rank == 0)
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid1);
    else
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid+1, &varid1);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_VAR_DIMIDS)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    /* Test inconsistent variable dimension IDs ------------------------------*/
    err = ncmpi_create(comm, filename, cmode, info, &ncid); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Z", 3, &dimid[0]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "Y", 3, &dimid[1]); CHECK_ERR
    err = ncmpi_def_dim(ncid, "X", 3, &dimid[2]); CHECK_ERR
    if (rank == 0)
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid+1, &varid1);
    else
        err = ncmpi_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid1);
    if (safe_mode)
        EXP_SAFE_ERR(NC_EMULTIDEFINE_VAR_DIMIDS)
    else
        CHECK_ERR
    err = ncmpi_set_fill(ncid, NC_FILL, NULL); CHECK_ERR
    err = ncmpi_close(ncid); CHECK_ERR

    return nerrs;
}

/*----< main() >--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    char *filename="testfile.nc", *mode[2] = {"0", "1"};
    int i, rank, nprocs, verbose, nerrs=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (nprocs < 2) {
        if (!rank) printf("This program is for running 2 or more processes. Exiting ...\n");
        MPI_Finalize();
        return 1;
    }

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 1;
    }
    if (argc == 2) filename = argv[1];

    if (rank == 0) {
        char *cmd_str = (char*)malloc(strlen(argv[0]) + 256);
        sprintf(cmd_str, "*** TESTING C   %s for header consistency", basename(argv[0]));
        printf("%-66s ------ ", cmd_str);
        free(cmd_str);
    }

    verbose = 1;
    for (i=verbose; i>=0; i--) {
        /* test with safe mode off and on :
         * Note even if --enable-debug is set at configure time, safe mode
         * can still be disabled by setting the environment variable
         * PNETCDF_SAFE_MODE to 0.
         */
        setenv("PNETCDF_SAFE_MODE", mode[i], 1);
        nerrs += test_open_mode(filename, i);

        nerrs += test_dim(filename, i);

        nerrs += test_attr(filename, i);

        nerrs += test_var(filename, i);
    }

    MPI_Offset malloc_size, sum_size;
    int err = ncmpi_inq_malloc_size(&malloc_size);
    if (err == NC_NOERR) {
        MPI_Reduce(&malloc_size, &sum_size, 1, MPI_OFFSET, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0 && sum_size > 0)
            printf("heap memory allocated by PnetCDF internally has %lld bytes yet to be freed\n",
                   sum_size);
        if (malloc_size > 0) ncmpi_inq_malloc_list();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nerrs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) {
        if (nerrs) printf(FAIL_STR,nerrs);
        else       printf(PASS_STR);
    }

    MPI_Finalize();
    return (nerrs > 0);
}

