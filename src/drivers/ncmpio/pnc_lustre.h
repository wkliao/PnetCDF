/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifndef H_PNC_LUSTRE
#define H_PNC_LUSTRE

#include <stddef.h>     /* size_t */
#include <sys/types.h>  /* off_t */

#include <dispatch.h>
#include "ncmpio_driver.h"

typedef struct {
    int initialized;
    int striping_factor;
    int striping_unit;
    int cb_read;
    int cb_write;
    int cb_nodes;
    int cb_buffer_size;
    int cb_ds_threshold;
    int ds_read;
    int ds_write;
    int ds_wr_lb;
    int no_indep_rw;
    int ind_rd_buffer_size;
    int ind_wr_buffer_size;
    int deferred_open;
    int start_iodevice;
    char *cb_config_list;
    int is_agg;                 /* bool: if I am an aggregator */
    int my_cb_nodes_index;      /* my index into cb_config_list. -1 if N/A */
    int *ranklist;
    int num_osts;
} PNC_Hints;

typedef struct {
    MPI_Comm comm;              /* communicator indicating who called open */
    char *filename;
    int file_system;            /* type of file system */

    int num_nodes;              /* number of unique compute nodes, MPI_Get_processor_name() */
    int access_mode;            /* Access mode (sequential, append, etc.),
                                 * possibly modified to deal with
                                 * data sieving or deferred open */
    int orig_access_mode;       /* Access mode provided by user: unmodified */
    MPI_Datatype etype;         /* reqd. for MPI-IO */
    MPI_Datatype filetype;      /* reqd. for MPI-IO */
    int atomicity;              /* true=atomic, false=nonatomic */
    char *io_buf;               /* two-phase buffer allocated out of i/o path */
    MPI_Info info;
    PNC_Hints *hints;         /* structure containing fs-indep. info values */
double lustre_write_metrics[3];
} PNC_File;



#endif
