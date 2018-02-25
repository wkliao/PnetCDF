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
define(`GETATTTYPE',dnl
`dnl
    ifelse($1, `MPI_CHAR', , `else ')if (itype == $1){
        err = nc_get_att_$2(nc4p->ncid, varid, name, ($3*) buf);
    }
')dnl
dnl
define(`PUTATTTYPE',dnl
`dnl
    ifelse($1, `MPI_CHAR', , `else ')if (itype == $1){
        err = nc_put_att_$2(nc4p->ncid, varid, name, len, ifelse($1, `MPI_CHAR', , `xtype, ')($3*) value);
    }
')dnl
dnl
define(`GETVARTYPE',dnl
`dnl
        ifelse($2, `MPI_CHAR', , `else ')if (buftype == $2){
            err = nc_get_$1_$3(nc4p->ncid, varid, ifelse($1, `var1', `(size_t*)start, ', $1, `vara', `(size_t*)start, (size_t*)count, ', $1, `vars', `(size_t*)start, (size_t*)count, (size_t*)stride, ', $1, `varm', `(size_t*)start, (size_t*)count, (size_t*)stride, (size_t*)imap, ')($4*) buf);
        }
')dnl
dnl
define(`PUTVARTYPE',dnl
`dnl
        ifelse($2, `MPI_CHAR', , `else ')if (buftype == $2){
            err = nc_put_$1_$3(nc4p->ncid, varid, ifelse($1, `var1', `(size_t*)start, ', $1, `vara', `(size_t*)start, (size_t*)count, ', $1, `vars', `(size_t*)start, (size_t*)count, (size_t*)stride, ', $1, `varm', `(size_t*)start, (size_t*)count, (size_t*)stride, (size_t*)imap, ')($4*) buf);
        }
')dnl
dnl
define(`GETVAR',dnl
`dnl
    ifelse($1, `var', , `else ')if (apikind == NC4_API_KIND_$2){
foreach(`dt', (`(`MPI_CHAR', `text', `char')', dnl
               `(`MPI_SIGNED_CHAR', `schar', `char')', dnl
               `(`MPI_UNSIGNED_CHAR', `uchar', `unsigned char')', dnl
               `(`MPI_SHORT', `short', `short')', dnl
               `(`MPI_UNSIGNED_SHORT', `ushort', `unsigned short')', dnl
               `(`MPI_INT', `int', `int')', dnl
               `(`MPI_UNSIGNED', `uint', `unsigned int')', dnl
               `(`MPI_FLOAT', `float', `float')', dnl
               `(`MPI_DOUBLE', `double', `double')', dnl
               `(`MPI_LONG_LONG_INT', `longlong', `long long')', dnl
               `(`MPI_UNSIGNED_LONG_LONG', `ulonglong', `unsigned long long')', dnl
               ), `GETVARTYPE($1, translit(dt, `()'))')dnl
    }
')dnl
dnl
define(`PUTVAR',dnl
`dnl
    ifelse($1, `var', , `else ')if (apikind == NC4_API_KIND_$2){
foreach(`dt', (`(`MPI_CHAR', `text', `char')', dnl
               `(`MPI_SIGNED_CHAR', `schar', `char')', dnl
               `(`MPI_UNSIGNED_CHAR', `uchar', `unsigned char')', dnl
               `(`MPI_SHORT', `short', `short')', dnl
               `(`MPI_UNSIGNED_SHORT', `ushort', `unsigned short')', dnl
               `(`MPI_INT', `int', `int')', dnl
               `(`MPI_UNSIGNED', `uint', `unsigned int')', dnl
               `(`MPI_FLOAT', `float', `float')', dnl
               `(`MPI_DOUBLE', `double', `double')', dnl
               `(`MPI_LONG_LONG_INT', `longlong', `long long')', dnl
               `(`MPI_UNSIGNED_LONG_LONG', `ulonglong', `unsigned long long')', dnl
               ), `PUTVARTYPE($1, translit(dt, `()'))')dnl
    }
')dnl

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

/* Note, netcdf header must come first due to conflicting constant definition */
#include <netcdf.h>

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <pnc_debug.h>
#include <common.h>
#include <nc4io_driver.h>

int
nc4io_get_att(void         *ncdp,
              int           varid,
              const char   *name,
              void         *buf,
              MPI_Datatype  itype)
{
    int err;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Call nc_del_att_<type> */
foreach(`dt', (`(`MPI_CHAR', `text', `char')', dnl
               `(`MPI_SIGNED_CHAR', `schar', `char')', dnl
               `(`MPI_UNSIGNED_CHAR', `uchar', `unsigned char')', dnl
               `(`MPI_SHORT', `short', `short')', dnl
               `(`MPI_UNSIGNED_SHORT', `ushort', `unsigned short')', dnl
               `(`MPI_INT', `int', `int')', dnl
               `(`MPI_UNSIGNED', `uint', `unsigned int')', dnl
               `(`MPI_FLOAT', `float', `float')', dnl
               `(`MPI_DOUBLE', `double', `double')', dnl
               `(`MPI_LONG_LONG_INT', `longlong', `long long')', dnl
               `(`MPI_UNSIGNED_LONG_LONG', `ulonglong', `unsigned long long')', dnl
               ), `GETATTTYPE(translit(dt, `()'))')dnl
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_put_att(void         *ncdp,
              int           varid,
              const char   *name,
              nc_type       xtype,
              MPI_Offset    nelems,
              const void    *value,
              MPI_Datatype  itype)
{
    int err;
    size_t len;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    /* Convert from MPI_Offset to size_t */
    len = (size_t)nelems;

    /* Call nc_del_att_<type> */
foreach(`dt', (`(`MPI_CHAR', `text', `char')', dnl
               `(`MPI_SIGNED_CHAR', `schar', `char')', dnl
               `(`MPI_UNSIGNED_CHAR', `uchar', `unsigned char')', dnl
               `(`MPI_SHORT', `short', `short')', dnl
               `(`MPI_UNSIGNED_SHORT', `ushort', `unsigned short')', dnl
               `(`MPI_INT', `int', `int')', dnl
               `(`MPI_UNSIGNED', `uint', `unsigned int')', dnl
               `(`MPI_FLOAT', `float', `float')', dnl
               `(`MPI_DOUBLE', `double', `double')', dnl
               `(`MPI_LONG_LONG_INT', `longlong', `long long')', dnl
               `(`MPI_UNSIGNED_LONG_LONG', `ulonglong', `unsigned long long')', dnl
               ), `PUTATTTYPE(translit(dt, `()'))')dnl
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_get_var(void             *ncdp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              void             *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int err;
    int apikind;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    if (start == NULL){
        apikind = NC4_API_KIND_VAR;
    }
    else if (count == NULL){
        apikind = NC4_API_KIND_VAR1;
    }
    else if (stride == NULL){
        apikind = NC4_API_KIND_VARA;
    }
    else if (imap == NULL){
        apikind = NC4_API_KIND_VARS;
    }
    else{
        apikind = NC4_API_KIND_VARM;
    }

foreach(`api', `(var, var1, vara, vars, varm)', `GETVAR(api, upcase(api))') dnl
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}

int
nc4io_put_var(void             *ncdp,
              int               varid,
              const MPI_Offset *start,
              const MPI_Offset *count,
              const MPI_Offset *stride,
              const MPI_Offset *imap,
              const void       *buf,
              MPI_Offset        bufcount,
              MPI_Datatype      buftype,
              int               reqMode)
{
    int err;
    int apikind;
    NC_nc4 *nc4p = (NC_nc4*)ncdp;
    
    if (start == NULL){
        apikind = NC4_API_KIND_VAR;
    }
    else if (count == NULL){
        apikind = NC4_API_KIND_VAR1;
    }
    else if (stride == NULL){
        apikind = NC4_API_KIND_VARA;
    }
    else if (imap == NULL){
        apikind = NC4_API_KIND_VARS;
    }
    else{
        apikind = NC4_API_KIND_VARM;
    }

foreach(`api', `(var, var1, vara, vars, varm)', `PUTVAR(api, upcase(api))') dnl
    if (err != NC_NOERR) DEBUG_RETURN_ERROR(err);

    return NC_NOERR;
}