/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the following PnetCDF APIs.
 *
 * ncmpi_inq_attname() : dispatcher->inq_attname()
 * ncmpi_inq_attid()   : dispatcher->inq_attid()
 * ncmpi_inq_att()     : dispatcher->inq_att()
 * ncmpi_rename_att()  : dispatcher->inq_rename_att()
 * ncmpi_copy_att()    : dispatcher->inq_copy_att()
 * ncmpi_del_att()     : dispatcher->inq_del_att()
 * ncmpi_get_att()     : dispatcher->inq_get_att()
 * ncmpi_put_att()     : dispatcher->inq_put_arr()
 *
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

nc_type ncadio_to_nc_type(enum ADIOS_DATATYPES atype){
    switch (atype) {
        case adios_unsigned_byte:
            return NC_BYTE;
        case adios_byte:
            return NC_BYTE;
        case adios_short:
            return NC_SHORT;
        case adios_unsigned_short:
            return NC_USHORT;
        case adios_integer:
            return NC_INT;
        case adios_unsigned_integer:
            return NC_UINT;
        case adios_long:
            return NC_LINT64;
        case adios_unsigned_long:
            return NC_UINT64;
        case adios_real:
            return NC_FLOAT;
        case adios_double:
            return NC_DOUBLE;
        case adios_long_double:
            return NC_DOUBLE;
        case adios_string:
            return NC_CHAR;
        case adios_string_array:
            return NC_STRING;
    }

    return NC_NAT;
}