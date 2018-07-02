/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <dispatch.h>
#include <ncadio_driver.h>

static PNC_driver ncadio_driver = {
    /* FILE APIs */
    ncadio_create,
    ncadio_open,
    ncadio_close,
    ncadio_enddef,
    ncadio__enddef,
    ncadio_redef,
    ncadio_sync,
    ncadio_flush,
    ncadio_abort,
    ncadio_set_fill,
    ncadio_inq,
    ncadio_inq_misc,
    ncadio_sync_numrecs,
    ncadio_begin_indep_data,
    ncadio_end_indep_data,

    /* DIMENSION APIs */
    ncadio_def_dim,
    ncadio_inq_dimid,
    ncadio_inq_dim,
    ncadio_rename_dim,

    /* ATTRIBUTE APIs */
    ncadio_inq_att,
    ncadio_inq_attid,
    ncadio_inq_attname,
    ncadio_copy_att,
    ncadio_rename_att,
    ncadio_del_att,
    ncadio_get_att,
    ncadio_put_att,

    /* VARIABLE APIs */
    ncadio_def_var,
    ncadio_def_var_fill,
    ncadio_fill_var_rec,
    ncadio_inq_var,
    ncadio_inq_varid,
    ncadio_rename_var,
    ncadio_get_var,
    ncadio_put_var,
    ncadio_get_varn,
    ncadio_put_varn,
    ncadio_get_vard,
    ncadio_put_vard,
    ncadio_iget_var,
    ncadio_iput_var,
    ncadio_bput_var,
    ncadio_iget_varn,
    ncadio_iput_varn,
    ncadio_bput_varn,

    ncadio_buffer_attach,
    ncadio_buffer_detach,
    ncadio_wait,
    ncadio_cancel
};

PNC_driver* ncadio_inq_driver(void) {
    return &ncadio_driver;
}

