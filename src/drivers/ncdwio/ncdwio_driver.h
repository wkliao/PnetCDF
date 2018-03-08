/*
 *  Copyright (C) 2017, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

/* $Id$ */

#ifndef _NCFOO_DRIVER_H
#define _NCFOO_DRIVER_H

#include <mpi.h>
#include <pnetcdf.h>
#include <dispatch.h>
#include <limits.h>
#include <unistd.h>

#define NC_LOG_TYPE_TEXT 1
#define NC_LOG_TYPE_SCHAR 2
#define NC_LOG_TYPE_UCHAR 3
#define NC_LOG_TYPE_SHORT 4
#define NC_LOG_TYPE_USHORT 5
#define NC_LOG_TYPE_INT 6
#define NC_LOG_TYPE_UINT 7
#define NC_LOG_TYPE_FLOAT 8
#define NC_LOG_TYPE_DOUBLE 9
#define NC_LOG_TYPE_LONGLONG 10
#define NC_LOG_TYPE_ULONGLONG 11
#define NC_LOG_TYPE_NATIVE 12

#define NC_LOG_API_KIND_VAR 1
#define NC_LOG_API_KIND_VAR1 2
#define NC_LOG_API_KIND_VARA 3
#define NC_LOG_API_KIND_VARS 4

#define NC_LOG_MAGIC_SIZE 8
#define NC_LOG_MAGIC "PnetCDF0"

#define NC_LOG_FORMAT_SIZE 8
#define NC_LOG_FORMAT_CDF_MAGIC "CDF0\0\0\0\0"
#define NC_LOG_FORMAT_HDF5_MAGIC "\211HDF\r\n\032\n"
#define NC_LOG_FORMAT_BP_MAGIC "BP\0\0\0\0\0\0"

#define NC_LOG_FALSE 0x00
#define NC_LOG_TRUE 0x01

#define NC_LOG_HINT_LOG_ENABLE 0x01
#define NC_LOG_HINT_DEL_ON_CLOSE 0x02
#define NC_LOG_HINT_FLUSH_ON_WAIT 0x04
#define NC_LOG_HINT_FLUSH_ON_SYNC 0x08
#define NC_LOG_HINT_FLUSH_ON_READ 0x10
#define NC_LOG_HINT_LOG_OVERWRITE 0x20
#define NC_LOG_HINT_LOG_CHECK 0x40
#define NC_LOG_HINT_LOG_SHARE 0x80

/* PATH_MAX after padding to 4 byte allignment */
#if PATH_MAX % 4 == 0
#define NC_LOG_PATH_MAX PATH_MAX
#elif PATH_MAX % 4 == 1
#define NC_LOG_PATH_MAX PATH_MAX + 3
#elif PATH_MAX % 4 == 2
#define NC_LOG_PATH_MAX PATH_MAX + 2
#elif PATH_MAX % 4 == 3
#define NC_LOG_PATH_MAX PATH_MAX + 1
#endif

/* Metadata header
 * Variable named according to the spec
 * ToDo: Replace int with 4 byte integer variable if int is not 4 byte
 */
typedef struct NC_dw_metadataheader {
    char magic[NC_LOG_MAGIC_SIZE];
    char format[NC_LOG_FORMAT_SIZE];
    int big_endian;
    int is_external;
    MPI_Offset num_ranks;
    MPI_Offset rank_id;
    MPI_Offset entry_begin;
    MPI_Offset max_ndims;
    MPI_Offset num_entries;
    int basenamelen;
    char basename[1];   /* The hack to keep basename inside the structure */
} NC_dw_metadataheader;

/* Metadata entry header
 * Variable named according to the spec
 * ToDo: Replace int with 4 byte integer variable if int is not 4 byte
 */
typedef struct NC_dw_metadataentry {
    MPI_Offset esize;
    int api_kind;
    int itype;
    int varid;
    int ndims;
    MPI_Offset data_off;
    MPI_Offset data_len;
} NC_dw_metadataentry;

typedef struct NC_dw_metadataptr {
    NC_dw_metadataentry *ptr;
    int valid;
    int reqid;
} NC_dw_metadataptr;

typedef struct NC_dw_metadataidx {
    NC_dw_metadataptr *entries;
    int nused;
    int nalloc;
} NC_dw_metadataidx;

/* Buffer structure */
typedef struct NC_dw_buffer {
    size_t nalloc;
    size_t nused;
    void *buffer;
} NC_dw_buffer;

/* Vector structure */
typedef struct NC_dw_sizevector {
    size_t nalloc;
    size_t nused;
    size_t *values;
} NC_dw_sizevector;

/* Vector structure */
typedef struct NC_dw_intvector {
    size_t nalloc;
    size_t nused;
    int *values;
} NC_dw_intvector;

/* Put_req structure */
typedef struct NC_dw_put_req {
    int valid;  // If this request object is in use (corresponding to some nonblocking request)
    int ready;  // If the corresponding log entry is flushed and the status is avaiable
    int status; // status of flushing the corresponding log entry
    int entrystart; // First log entry generated by this request
    int entryend;   // Last log entry generated by this request
} NC_dw_put_req;

/* Put_req structure */
typedef struct NC_dw_put_list {
    NC_dw_put_req *reqs;    // Array of request object
    int *ids;   // Array of request ids
    int nalloc; // Size of the pool
    int nused;  // Number of ids issued
} NC_dw_put_list;

/* Shared file object */
typedef struct NC_dw_sharedfile {
    int fd; // POSIX file descriptor
    int chanel; // Which chanel are we on (Usually according to the rank)
    int nchanel;    // How many chanel are there (how many process are sharing the file)
    size_t pos; // Logical file position within the fileview
    size_t bsize;   // Dividing blocksize
    size_t fsize;   // Current file size
} NC_dw_sharedfile;

/* File structure */
typedef struct NC_dw_bufferedfile {
    NC_dw_sharedfile *fd;    // Shared file
    size_t pos; // File position
    char *buffer;  // Buffer
    // We require buffered region always maps to an aligned block boundary
    // If we seek to an unaligned position, the region from the start of buffer to this region must be mark as unused to prevent being flushed to the file
    size_t bunused;  // Unused amount of the buffer
    size_t bused;     // Buffer used region
    size_t bsize;   // Buffer size, also write block size
    size_t fsize;   // Current file size
} NC_dw_bufferedfile;

/* Log structure */
typedef struct NC_dw {
    char metalogpath[PATH_MAX];    /* path of metadata log */
    char datalogpath[PATH_MAX];    /* path of data log */
    char logbase[PATH_MAX];        /* path of log files */
    int rank;
    int np;
    NC_dw_sharedfile *metalog_fd;    /* file handle of metadata log */
    NC_dw_bufferedfile *datalog_fd;    /* file handle of data log */
    int recdimid;
    int inited;
    int hints;
    int isindep;
    size_t datalogsize;
    NC_dw_buffer metadata; /* In memory metadata buffer that mirrors the metadata log */
    NC_dw_metadataidx metaidx;
    NC_dw_sizevector entrydatasize;    /* Array of metadata entries */
    int isflushing;   /* If log is flushing */
    MPI_Offset max_ndims;
    NC_dw_put_list putlist;
    MPI_Offset recdimsize;
    MPI_Offset flushbuffersize;
    MPI_Offset maxentrysize;
#ifdef PNETCDF_PROFILING
    /* Profiling information */
    MPI_Offset total_data;
    MPI_Offset total_meta;
    MPI_Offset max_buffer;
    double total_time;
    double create_time;
    double enddef_time;
    double put_time;
    double flush_time;
    double close_time;
    double flush_replay_time;
    double flush_data_rd_time;
    double flush_put_time;
    double flush_wait_time;
    double put_data_wr_time;
    double put_meta_wr_time;
    double put_num_wr_time;
#endif

    int                mode;        /* file _open/_create mode */
    int                flag;        /* define/data/collective/indep mode */
    int                ncid;
    char              *path;        /* path name */
    MPI_Comm           comm;        /* MPI communicator */
    MPI_Comm           logcomm;        /* MPI communicator */
    MPI_Info           info;
    void              *ncp;         /* pointer to driver's internal object */
    struct PNC_driver *ncmpio_driver;
} NC_dw;

int ncdwio_log_buffer_init(NC_dw_buffer * bp);
void ncdwio_log_buffer_free(NC_dw_buffer * bp);
char* ncdwio_log_buffer_alloc(NC_dw_buffer *bp, size_t size);
int ncdwio_log_sizearray_init(NC_dw_sizevector *sp);
void ncdwio_log_sizearray_free(NC_dw_sizevector *sp);
int ncdwio_log_sizearray_append(NC_dw_sizevector *sp, size_t size);
int log_flush(NC_dw *ncdwp);
int ncdwio_log_create(NC_dw *ncdwp, MPI_Info info);
int ncdwio_log_put_var(NC_dw *ncdwp, int varid, const MPI_Offset start[], const MPI_Offset count[], const MPI_Offset stride[], void *buf, MPI_Datatype buftype, MPI_Offset *putsize);
int ncdwio_log_close(NC_dw *ncdwp);
int ncdwio_log_flush(NC_dw *ncdwp);
int ncdwio_log_enddef(NC_dw *ncdwp);

int ncdwio_put_list_init(NC_dw *ncdwp);
int ncdwio_put_list_resize(NC_dw *ncdwp);
int ncdwio_put_list_free(NC_dw *ncdwp);
int ncdwio_put_list_add(NC_dw *ncdwp, int *id);
int ncdwio_put_list_remove(NC_dw *ncdwp, int reqid);
int ncdwio_handle_put_req(NC_dw *ncdwp, int reqid, int *stat);
int ncdwio_handle_all_put_req(NC_dw *ncdwp);
int ncdwio_cancel_put_req(NC_dw *ncdwp, int reqid, int *stat);
int ncdwio_cancel_all_put_req(NC_dw *ncdwp);
int ncdwio_metaidx_init(NC_dw *ncdwp);
int ncdwio_metaidx_add(NC_dw *ncdwp, NC_dw_metadataentry *entry);
int ncdwio_metaidx_free(NC_dw *ncdwp);
int ncdwio_log_intvector_init(NC_dw_intvector *vp);
void ncdwio_log_intvector_free(NC_dw_intvector *vp);
int ncdwio_log_intvector_append(NC_dw_intvector *vp, int size);

int ncdwio_sharedfile_open(MPI_Comm comm, char *path, int flag, MPI_Info info, NC_dw_sharedfile **fh);
int ncdwio_sharedfile_close(NC_dw_sharedfile *f);
int ncdwio_sharedfile_pwrite(NC_dw_sharedfile *f, void *buf, size_t count, off_t offset);
int ncdwio_sharedfile_write(NC_dw_sharedfile *f, void *buf, size_t count);
int ncdwio_sharedfile_pread(NC_dw_sharedfile *f, void *buf, size_t count, off_t offset);
int ncdwio_sharedfile_read(NC_dw_sharedfile *f, void *buf, size_t count);
int ncdwio_sharedfile_seek(NC_dw_sharedfile *f, off_t offset, int whence);

int ncdwio_bufferedfile_open(MPI_Comm comm, char *path, int flag, MPI_Info info, NC_dw_bufferedfile **fh);
int ncdwio_bufferedfile_close(NC_dw_bufferedfile *f);
int ncdwio_bufferedfile_pwrite(NC_dw_bufferedfile *f, void *buf, size_t count, off_t offset);
int ncdwio_bufferedfile_write(NC_dw_bufferedfile *f, void *buf, size_t count);
int ncdwio_bufferedfile_pread(NC_dw_bufferedfile *f, void *buf, size_t count, off_t offset);
int ncdwio_bufferedfile_read(NC_dw_bufferedfile *f, void *buf, size_t count);
int ncdwio_bufferedfile_seek(NC_dw_bufferedfile *f, off_t offset, int whence);

void ncdwio_extract_hint(NC_dw *ncdwp, MPI_Info info);
void ncdwio_export_hint(NC_dw *ncdwp, MPI_Info info);

extern int
ncdwio_create(MPI_Comm comm, const char *path, int cmode, int ncid, MPI_Info info, void **ncdp);

extern int
ncdwio_open(MPI_Comm comm, const char *path, int omode, int ncid, MPI_Info info, void **ncdp);

extern int
ncdwio_close(void *ncdp);

extern int
ncdwio_enddef(void *ncdp);

extern int
ncdwio__enddef(void *ncdp, MPI_Offset h_minfree, MPI_Offset v_align, MPI_Offset v_minfree, MPI_Offset r_align);

extern int
ncdwio_redef(void *ncdp);

extern int
ncdwio_sync(void *ncdp);

extern int
ncdwio_flush(void *ncdp);

extern int
ncdwio_abort(void *ncdp);

extern int
ncdwio_set_fill(void *ncdp, int fill_mode, int *old_fill_mode);

extern int
ncdwio_fill_var_rec(void *ncdp, int varid, MPI_Offset recno);

extern int
ncdwio_inq(void *ncdp, int *ndimsp, int *nvarsp, int *nattsp, int *xtendimp);

extern int
ncdwio_inq_misc(void *ncdp, int *pathlen, char *path, int *num_fix_varsp,
               int *num_rec_varsp, int *striping_size, int *striping_count,
               MPI_Offset *header_size, MPI_Offset *header_extent,
               MPI_Offset *recsize, MPI_Offset *put_size, MPI_Offset *get_size,
               MPI_Info *info_used, int *nreqs, MPI_Offset *usage,
               MPI_Offset *buf_size);

extern int
ncdwio_sync_numrecs(void *ncdp);

extern int
ncdwio_begin_indep_data(void *ncdp);

extern int
ncdwio_end_indep_data(void *ncdp);

extern int
ncdwio_def_dim(void *ncdp, const char *name, MPI_Offset size, int *dimidp);

extern int
ncdwio_inq_dimid(void *ncdp, const char *name, int *dimidp);

extern int
ncdwio_inq_dim(void *ncdp, int dimid, char *name, MPI_Offset *lengthp);

extern int
ncdwio_rename_dim(void *ncdp, int dimid, const char *newname);

extern int
ncdwio_inq_att(void *ncdp, int varid, const char *name, nc_type *xtypep, MPI_Offset *lenp);

extern int
ncdwio_inq_attid(void *ncdp, int varid, const char *name, int *idp);

extern int
ncdwio_inq_attname(void *ncdp, int varid, int attnum, char *name);

extern int
ncdwio_copy_att(void *ncdp_in, int varid_in, const char *name, void *ncdp_out, int varid_out);

extern int
ncdwio_rename_att(void *ncdp, int varid, const char *name, const char *newname);

extern int
ncdwio_del_att(void *ncdp, int varid, const char *name);

extern int
ncdwio_get_att(void *ncdp, int varid, const char *name, void *value, MPI_Datatype itype);

extern int
ncdwio_put_att(void *ncdp, int varid, const char *name, nc_type xtype, MPI_Offset nelems, const void *value, MPI_Datatype itype);

extern int
ncdwio_def_var(void *ncdp, const char *name, nc_type type, int ndims, const int *dimids, int *varidp);

extern int
ncdwio_def_var_fill(void *ncdp, int varid, int nofill, const void *fill_value);

extern int
ncdwio_inq_var(void *ncdp, int varid, char *name, nc_type *xtypep, int *ndimsp,
               int *dimids, int *nattsp, MPI_Offset *offsetp, int *no_fill, void *fill_value);

extern int
ncdwio_inq_varid(void *ncdp, const char *name, int *varid);

extern int
ncdwio_rename_var(void *ncdp, int varid, const char *newname);

extern int
ncdwio_get_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncdwio_put_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncdwio_get_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncdwio_put_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncdwio_get_vard(void *ncdp, int varid, MPI_Datatype filetype, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncdwio_put_vard(void *ncdp, int varid, MPI_Datatype filetype, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int reqMode);

extern int
ncdwio_iget_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncdwio_iput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncdwio_bput_var(void *ncdp, int varid, const MPI_Offset *start, const MPI_Offset *count, const MPI_Offset *stride, const MPI_Offset *imap, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *req, int reqMode);

extern int
ncdwio_iget_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncdwio_iput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncdwio_bput_varn(void *ncdp, int varid, int num, MPI_Offset* const *starts, MPI_Offset* const *counts, const void *buf, MPI_Offset bufcount, MPI_Datatype buftype, int *reqid, int reqMode);

extern int
ncdwio_buffer_attach(void *ncdp, MPI_Offset bufsize);

extern int
ncdwio_buffer_detach(void *ncdp);

extern int
ncdwio_wait(void *ncdp, int num_reqs, int *req_ids, int *statuses, int reqMode);

extern int
ncdwio_cancel(void *ncdp, int num_reqs, int *req_ids, int *statuses);

#endif
