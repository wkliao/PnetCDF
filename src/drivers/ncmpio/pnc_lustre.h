/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifndef H_PNC_LUSTRE
#define H_PNC_LUSTRE

#include <stddef.h>     /* size_t */
#include <sys/types.h>  /* off_t */
#include <assert.h>

#include <dispatch.h>
#include "ncmpio_driver.h"
#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

#define ADIO_LUSTRE 163
#define ADIOI_Strdup strdup
#define ADIOI_Malloc NCI_Malloc
#define ADIOI_Calloc NCI_Calloc
#define ADIOI_Free NCI_Free
#define ADIOI_Strncpy strncpy
#define ADIOI_Info_get MPI_Info_get
#define ADIOI_Info_set MPI_Info_set
#define ADIOI_Assert assert

#define ADIOI_CB_BUFFER_SIZE_DFLT "16777216"
#define ADIOI_IND_RD_BUFFER_SIZE_DFLT     "4194304"
#define ADIOI_IND_WR_BUFFER_SIZE_DFLT     "524288"
#define ADIOI_CB_CONFIG_LIST_DFLT "*:1"

    /* ADIOI_DS_WR_LB is the lower bound of the number of noncontiguous
     * offset-length pairs to trigger data sieving write. When hint ds_write is
     * set to 'auto' and the number of offset-length pairs is more than hint
     * romio_cb_wr_lb, then data sieving is activated and checking holes in the
     * file domains is skipped. When the number of holes is large, checking
     * holes can be expensive, because it requires to a merge-sort of all the
     * offset-length pairs.
     */
#define ADIOI_DS_WR_LB "8192"

#define ADIOI_TYPE_DECREASE 0x00000001  /* if not monotonic nondecreasing */
#define ADIOI_TYPE_OVERLAP  0x00000002  /* if contains overlapping regions */
#define ADIOI_TYPE_NEGATIVE 0x00000004  /* if one of displacements is negative */

typedef MPI_Offset ADIO_Offset;

enum {
    ADIOI_HINT_AUTO = 0,
    ADIOI_HINT_ENABLE = 1,
    ADIOI_HINT_DISABLE = 2
};

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
    int *ranklist;
    int num_osts;
} PNC_Hints;

typedef struct {
    MPI_Comm comm;              /* communicator indicating who called open */
    char *filename;
    int file_system;            /* type of file system */
    int fd_sys;                 /* system file descriptor */
    int num_nodes;              /* number of unique compute nodes, MPI_Get_processor_name() */
    int access_mode;            /* Access mode (sequential, append, etc.),
                                 * possibly modified to deal with
                                 * data sieving or deferred open */
    int orig_access_mode;       /* Access mode provided by user: unmodified */
    int is_open;                /* no_indep_rw, 0: not open yet 1: is open */
    ADIO_Offset disp;
    ADIO_Offset fp_ind;         /* individual file pointer in MPI-IO (in bytes) */
    MPI_Datatype etype;         /* reqd. for MPI-IO */
    MPI_Datatype filetype;      /* reqd. for MPI-IO */
    int atomicity;              /* true=atomic, false=nonatomic */
    char *io_buf;               /* two-phase buffer allocated out of i/o path */
    int is_agg;                 /* bool: if I am an aggregator */
    int my_cb_nodes_index;      /* my index into cb_config_list. -1 if N/A */
    MPI_Info info;
    PNC_Hints *hints;         /* structure containing fs-indep. info values */
double lustre_write_metrics[3];
} PNC_File;

typedef PNC_File *ADIO_File;

typedef struct ADIOI_Fl_node {
    MPI_Datatype type;
    MPI_Count count;            /* no. of contiguous blocks */
    ADIO_Offset *blocklens;     /* array of contiguous block lengths (bytes) */
    ADIO_Offset *indices;       /* array of byte offsets of each block */
    /* the type processing code in ROMIO loops through the flattened
     * representation to tile file views.  so, we cannot simply indicate a
     * lower bound and upper bound with entries here -- those are instead
     * indicated by offset-length pairs with zero length.  In order to deal
     * with repeatedly resetting the LB and UB though (as in resized of
     * resized, or struct of resized, perhaps?), indicate where in the
     * tuple-stream the lb and ub parts are kept
     * (-1 indicates "not explicitly set") */
    ADIO_Offset lb_idx;
    ADIO_Offset ub_idx;
    int refct;                  /* when storing flattened representation on a
                                 * type, attribute copy and delete routines
                                 * will manage refct */
    int flag;                   /* ADIOI_TYPE_XXX */
} ADIOI_Flatlist_node;


int PNC_Lustre_open(MPI_Comm comm, const char *filename, int amode,
                    MPI_Info info, PNC_File *fd);

int PNC_Lustre_close(PNC_File *fd);
int PNC_Lustre_sync(PNC_File *fd);
int PNC_Lustre_delete(const char *filename);
int PNC_Lustre_set_size(PNC_File *fd, MPI_Offset size);
int PNC_Lustre_get_size(PNC_File *fd, MPI_Offset *size);
int PNC_Lustre_seek(PNC_File *fd, MPI_Offset offset, int whence);
int PNC_Lustre_get_info(PNC_File *fd, MPI_Info *info_used);
int PNC_lustre_SetInfo(PNC_File *fd, MPI_Info  users_info);

void PNC_Datatype_iscontig(MPI_Datatype datatype, int *flag);
int PNC_Type_ispredef(MPI_Datatype datatype, int *flag);
int PNC_Type_dispose(MPI_Datatype *datatype);
ADIOI_Flatlist_node *PNC_Flatten_and_find(MPI_Datatype datatype);



#endif
