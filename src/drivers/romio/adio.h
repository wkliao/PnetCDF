/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifndef H_ADIO
#define H_ADIO

#include <stdio.h>
#include <stdlib.h>
#include <sys/errno.h>
#include <unistd.h>   /* pwrite() */

#include <stdbool.h>
#include <string.h>     /* memcpy() */
#include <stddef.h>     /* size_t */
#include <sys/types.h>  /* off_t */
#include <assert.h>
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif
#define FDTYPE int

#include <pnc_debug.h>
#include <common.h>

#ifndef MPL_MAX
#define MPL_MAX MAX
#endif
#ifndef MPL_MIN
#define MPL_MIN MIN
#endif
#define MPL_UNREFERENCED_ARG(a)

#define ADIO_LOCKS  300    /* file system supports fcntl()-style locking */
#define ADIO_Feature(a, b) ((b == ADIO_LOCKS) ? 1 : 0)

#if defined(F_SETLKW64)
#define ADIOI_UNLOCK(fd, offset, whence, len) \
        ADIOI_GEN_SetLock64(fd, F_SETLK, F_UNLCK, offset, whence, len)
#define ADIOI_WRITE_LOCK(fd, offset, whence, len) \
        ADIOI_GEN_SetLock64(fd, F_SETLKW, F_WRLCK, offset, whence, len)
#else
#define ADIOI_UNLOCK(fd, offset, whence, len) \
        ADIOI_GEN_SetLock(fd, F_SETLK, F_UNLCK, offset, whence, len)
#define ADIOI_WRITE_LOCK(fd, offset, whence, len) \
        ADIOI_GEN_SetLock(fd, F_SETLKW, F_WRLCK, offset, whence, len)
#endif


#define ADIO_UFS                 152    /* Unix file system */
#define ADIO_LUSTRE              163    /* Lustre */

#define ADIOI_Strdup strdup
#define ADIOI_Malloc malloc
#define ADIOI_Calloc(a,b) calloc((size_t)a,b)
#define ADIOI_Realloc realloc
#define ADIOI_Free free
#define ADIOI_Strncpy strncpy
#define ADIOI_Info_get MPI_Info_get
#define ADIOI_Info_set MPI_Info_set
#define ADIOI_Assert assert
#define ADIO_Status MPI_Status
#define ADIO_EXPLICIT_OFFSET     100
#define ADIO_INDIVIDUAL          101
#define ADIOI_COLL_TAG(rank,iter) 0
#define ADIOI_AINT_CAST_TO_OFFSET (MPI_Offset)

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

#define ADIO_OFFSET MPI_OFFSET
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

    int cb_pfr;
    int min_fdomain_size;
    union {
        struct {
            int num_osts;
            int co_ratio;
            int coll_threshold;
        } lustre;
    } fs_hints;
} ADIOI_Hints;

typedef struct {
    MPI_Comm comm;          /* communicator indicating who called open */
    char *filename;
    int file_system;        /* type of file system */
    int fd_sys;             /* system file descriptor */
    int num_nodes;          /* number of unique compute nodes from
                             * MPI_Get_processor_name() */
    int access_mode;        /* Access mode (sequential, append, etc.),
                             * possibly modified to deal with
                             * data sieving or deferred open */
    int orig_access_mode;   /* Access mode provided by user: unmodified */
    int is_open;            /* no_indep_rw, 0: not open yet 1: is open */
    ADIO_Offset fp_ind;     /* individual file pointer (in bytes) */
    ADIO_Offset fp_sys_posn;/* current location of the system file-pointer
                             * in bytes */

    ADIO_Offset disp;       /* file displacement */
    MPI_Datatype etype;     /* element datatype */
    MPI_Datatype filetype;  /* file type set in fileview */
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count etype_size;
    MPI_Count ftype_size;
    MPI_Count ftype_extent;
#else
    int etype_size;
    int ftype_size;
    MPI_Aint ftype_extent;
#endif

    int atomicity;          /* true=atomic, false=nonatomic */
    char *io_buf;           /* two-phase buffer allocated out of i/o path */
    int is_agg;             /* bool: if I am an aggregator */
    int my_cb_nodes_index;  /* my index into cb_config_list. -1 if N/A */
    ADIOI_Hints *hints;     /* structure containing fs-indep. info values */
    MPI_Info info;

double lustre_write_metrics[3];
} ADIO_FileD;

typedef ADIO_FileD *ADIO_File;
typedef ADIO_FileD *ADIO_File;

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

typedef struct {
    ADIO_Offset *offsets;       /* array of offsets */
    ADIO_Offset *lens;          /* array of lengths */
#if 1
    MPI_Count *mem_ptrs;        /* array of pointers. used in the read/write
                                 * phase to indicate where the data
                                 * is stored in memory
                                 * promoted to MPI_Count so we can construct types with _c versions */
    MPI_Count count;            /* size of above arrays */
#else
    MPI_Aint *mem_ptrs;
    size_t count;
#endif
} ADIOI_Access;

/*---- APIs -----------------------------------------------------------------*/
int PNC_Check_Lustre(const char *filename);
int PNC_File_open(MPI_Comm comm, const char *filename, int amode,
                    MPI_Info info, ADIO_File *fh);

int PNC_File_close(ADIO_File *fh);
int PNC_File_set_view(ADIO_File fh, MPI_Offset disp, MPI_Datatype etype,
                      MPI_Datatype filetype, char *datarep, MPI_Info info);
int PNC_File_sync(ADIO_File fh);
int PNC_File_delete(const char *filename);
int PNC_File_set_size(ADIO_File fh, MPI_Offset size);
int PNC_File_get_size(ADIO_File fh, MPI_Offset *size);
int PNC_File_seek(ADIO_File fh, MPI_Offset offset, int whence);
int PNC_File_get_info(ADIO_File fh, MPI_Info *info_used);
int PNC_File_SetInfo(ADIO_File fh, MPI_Info  users_info);

int PNC_File_write(ADIO_File fh, const void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status);
int PNC_File_write_at(ADIO_File fh, MPI_Offset offset, const void *buf,
                   int count, MPI_Datatype  datatype, MPI_Status *status);
int PNC_File_write_all(ADIO_File fh, const void *buf, int count,
                   MPI_Datatype datatype, MPI_Status *status);
int PNC_File_write_at_all(ADIO_File fh, MPI_Offset offset, const void *buf,
                   int count, MPI_Datatype  datatype, MPI_Status *status);

int PNC_File_read(ADIO_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
int PNC_File_read_at(ADIO_File fh, MPI_Offset offset, void *buf,
                  int count, MPI_Datatype  datatype, MPI_Status *status);
int PNC_File_read_all(ADIO_File fh, void *buf, int count,
                  MPI_Datatype datatype, MPI_Status *status);
int PNC_File_read_at_all(ADIO_File fh, MPI_Offset offset, void *buf,
                  int count, MPI_Datatype  datatype, MPI_Status *status);

void ADIOI_Datatype_iscontig(MPI_Datatype datatype, int *flag);
#define PNC_Datatype_iscontig ADIOI_Datatype_iscontig
// void PNC_Datatype_iscontig(MPI_Datatype datatype, int *flag);

int ADIOI_Type_ispredef(MPI_Datatype datatype, int *flag);
#define PNC_Type_ispredef ADIOI_Type_ispredef
// int PNC_Type_ispredef(MPI_Datatype datatype, int *flag);

int ADIOI_Type_dispose(MPI_Datatype * datatype);
#define PNC_Type_dispose ADIOI_Type_dispose
// int PNC_Type_dispose(MPI_Datatype *datatype);

ADIOI_Flatlist_node *ADIOI_Flatten_and_find(MPI_Datatype);
#define PNC_Flatten_and_find ADIOI_Flatten_and_find
// ADIOI_Flatlist_node *PNC_Flatten_and_find(MPI_Datatype datatype);

void ADIOI_LUSTRE_WriteStrided(ADIO_File fd, const void *buf, MPI_Aint count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status, int *error_code);
#define ADIO_WriteStrided ADIOI_LUSTRE_WriteStrided

/*
int PNC_WriteStrided(ADIO_File fd, const void *buf, MPI_Aint count,
                     MPI_Datatype datatype, int file_ptr_type,
                     ADIO_Offset offset, ADIO_Status *status);
*/

void ADIOI_LUSTRE_WriteStridedColl(ADIO_File fd, const void *buf,
                     MPI_Aint count, MPI_Datatype buftype, int file_ptr_type,
                     ADIO_Offset offset, ADIO_Status *status, int *error_code);


void ADIOI_GEN_ReadStrided_naive(ADIO_File fd, void *buf, MPI_Aint count,
                                 MPI_Datatype buftype, int file_ptr_type,
                                 ADIO_Offset offset, ADIO_Status * status, int
                                 *error_code);

#define ADIO_ReadStridedColl ADIOI_GEN_ReadStridedColl
void ADIOI_GEN_ReadStridedColl(ADIO_File fd, void *buf, MPI_Aint count,
                               MPI_Datatype datatype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status * status, int
                               *error_code);

void ADIOI_GEN_WriteStrided_naive(ADIO_File fd, const void *buf, MPI_Aint count,
                                  MPI_Datatype buftype, int file_ptr_type,
                                  ADIO_Offset offset, ADIO_Status * status, int
                                  *error_code);
int PNC_GEN_WriteStrided_naive(ADIO_File fd, const void *buf, MPI_Aint count,
                               MPI_Datatype buftype, int file_ptr_type,
                               ADIO_Offset offset, ADIO_Status *status);

#define ADIO_WriteContig(a,b,c,d,e,f,g,x) \
    *x = PNC_WriteContig(a,b,c,d,e,f,g);
int PNC_WriteContig(ADIO_File fd, const void *buf, MPI_Aint count,
                    MPI_Datatype bufType, int file_ptr_type,
                    ADIO_Offset offset, ADIO_Status *status);

#define ADIO_ReadContig(a,b,c,d,e,f,g,x) \
    *x = PNC_ReadContig(a,b,c,d,e,f,g);
int PNC_ReadContig(ADIO_File fd, void *buf, MPI_Aint count,
                    MPI_Datatype bufType, int file_ptr_type,
                    ADIO_Offset offset, ADIO_Status *status);


#define ADIO_ReadStrided ADIOI_GEN_ReadStrided
void ADIOI_GEN_ReadStrided(ADIO_File fd, void *buf, MPI_Aint count,
                           MPI_Datatype datatype, int file_ptr_type,
                           ADIO_Offset offset, ADIO_Status * status,
                           int *error_code);

void ADIOI_Calc_my_off_len(ADIO_File fd, MPI_Aint bufcount, MPI_Datatype
                           datatype, int file_ptr_type, ADIO_Offset
                           offset, ADIO_Offset ** offset_list_ptr, ADIO_Offset
                           ** len_list_ptr, ADIO_Offset * start_offset_ptr,
                           ADIO_Offset * end_offset_ptr, MPI_Count * contig_access_count_ptr);

void ADIOI_Calc_file_domains(ADIO_Offset * st_offsets, ADIO_Offset
                             * end_offsets, int nprocs, int nprocs_for_coll,
                             ADIO_Offset * min_st_offset_ptr,
                             ADIO_Offset ** fd_start_ptr, ADIO_Offset
                             ** fd_end_ptr, int min_fd_size,
                             ADIO_Offset * fd_size_ptr, int striping_unit);
void ADIOI_Calc_my_req(ADIO_File fd, ADIO_Offset * offset_list, ADIO_Offset * len_list,
                       MPI_Count contig_access_count,
                       ADIO_Offset min_st_offset, ADIO_Offset * fd_start,
                       ADIO_Offset * fd_end, ADIO_Offset fd_size,
                       int nprocs,
                       MPI_Count * count_my_req_procs_ptr,
                       MPI_Count ** count_my_req_per_proc_ptr,
                       ADIOI_Access ** my_req_ptr, MPI_Aint ** buf_idx_ptr);
void ADIOI_Calc_others_req(ADIO_File fd, MPI_Count count_my_req_procs,
                           MPI_Count * count_my_req_per_proc,
                           ADIOI_Access * my_req,
                           int nprocs, int myrank,
                           MPI_Count * count_others_req_procs_ptr,
                           MPI_Count ** count_others_req_per_proc_ptr,
                           ADIOI_Access ** others_req_ptr);
void ADIOI_Free_my_req(int nprocs, MPI_Count * count_my_req_per_proc,
                       ADIOI_Access * my_req, MPI_Aint * buf_idx);
void ADIOI_Free_others_req(int nprocs, MPI_Count * count_others_req_per_proc,
                           ADIOI_Access * others_req);

int ADIOI_Type_create_hindexed_x(MPI_Count count,
                                 const MPI_Count array_of_blocklengths[],
                                 const MPI_Count array_of_displacements[],
                                 MPI_Datatype oldtype, MPI_Datatype * newtype);

int ADIOI_Calc_aggregator(ADIO_File fd,
                          ADIO_Offset off,
                          ADIO_Offset min_off,
                          ADIO_Offset * len,
                          ADIO_Offset fd_size, ADIO_Offset * fd_start, ADIO_Offset * fd_end);



#define MPIO_Err_create_code(err, ...) err
#define FPRINTF fprintf
#define DBG_FPRINTF fprintf

int ADIO_Type_create_subarray(int ndims,
                              const int *array_of_sizes,
                              const int *array_of_subsizes,
                              const int *array_of_starts,
                              int order, MPI_Datatype oldtype, MPI_Datatype * newtype);
int ADIO_Type_create_darray(int size, int rank, int ndims,
                            const int *array_of_gsizes, const int *array_of_distribs,
                            const int *array_of_dargs, const int *array_of_psizes,
                            int order, MPI_Datatype oldtype, MPI_Datatype * newtype);
int ADIOI_Type_get_combiner(MPI_Datatype datatype, int *combiner);


int ADIOI_GEN_SetLock(ADIO_File fd, int cmd, int type, ADIO_Offset offset, int whence,
                      ADIO_Offset len);
int ADIOI_GEN_SetLock64(ADIO_File fd, int cmd, int type, ADIO_Offset offset, int whence,
                        ADIO_Offset len);

#endif
