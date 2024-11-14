/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#include "string.h" /*strdup() */
#include "pnc_lustre.h"

#define ADIO_LUSTRE 163
#define ADIOI_Strdup strdup
#define ADIOI_Malloc NCI_Malloc
#define ADIOI_Free NCI_Free

/*
 ADIO_FileSysType_parentdir - determines a string pathname for the
 parent directory of a given filename.

Input Parameters:
. filename - pointer to file name character array

Output Parameters:
. dirnamep - pointer to location in which to store a pointer to a string

 Note that the caller should free the memory located at the pointer returned
 after the string is no longer needed.
*/

/* In a strict ANSI environment, S_ISLNK may not be defined.  Fix that
   here.  We assume that S_ISLNK is *always* defined as a macro.  If
   that is not universally true, then add a test to the romio
   configure that tries to link a program that references S_ISLNK */
#if !defined(S_ISLNK)
#if defined(S_IFLNK)
     /* Check for the link bit */
#define S_ISLNK(mode) ((mode) & S_IFLNK)
#else
     /* no way to check if it is a link, so say false */
#define S_ISLNK(mode) 0
#endif
#endif /* !(S_ISLNK) */

/* ADIO_FileSysType_parentdir
 *
 * Returns pointer to string in dirnamep; that string is allocated with
 * strdup and must be free()'d.
 */
static void ADIO_FileSysType_parentdir(const char *filename, char **dirnamep)
{
    int err;
    char *dir = NULL, *slash;
    struct stat statbuf;

    err = lstat(filename, &statbuf);

    if (err || (!S_ISLNK(statbuf.st_mode))) {
        /* no such file, or file is not a link; these are the "normal"
         * cases where we can just return the parent directory.
         */
        dir = ADIOI_Strdup(filename);
    } else {
        /* filename is a symlink.  we've presumably already tried
         * to stat it and found it to be missing (dangling link),
         * but this code doesn't care if the target is really there
         * or not.
         */
        ssize_t namelen;
        char *linkbuf;

        linkbuf = ADIOI_Malloc(PATH_MAX + 1);
        namelen = readlink(filename, linkbuf, PATH_MAX + 1);
        if (namelen == -1) {
            /* something strange has happened between the time that
             * we determined that this was a link and the time that
             * we attempted to read it; punt and use the old name.
             */
            dir = ADIOI_Strdup(filename);
        } else {
            /* successfully read the link */
            linkbuf[namelen] = '\0';    /* readlink doesn't null terminate */
            dir = ADIOI_Strdup(linkbuf);
        }
        ADIOI_Free(linkbuf);
    }

    slash = strrchr(dir, '/');
    if (!slash)
        ADIOI_Strncpy(dir, ".", 2);
    else {
        if (slash == dir)
            *(dir + 1) = '\0';
        else
            *slash = '\0';
    }

    *dirnamep = dir;
    return;
}

#define UNKNOWN_SUPER_MAGIC (0xDEADBEEF)
#ifndef LL_SUPER_MAGIC
#define LL_SUPER_MAGIC 0x0BD00BD0
#endif

static int romio_statfs(const char *filename, int64_t * file_id)
{

    int err = 0;

#ifdef HAVE_STRUCT_STATVFS_WITH_F_BASETYPE
    /* rare: old solaris machines */
    struct statvfs vfsbuf;
#endif
#if defined(HAVE_STRUCT_STATFS_F_TYPE) || defined(HAVE_STRUCT_STATFS_F_FSTYPENAME)
    /* common fs-detection logic for any modern POSIX-compliant environment,
     * with the one wrinkle that some platforms (Darwin, BSD) give us a file
     * system as a string, not an identifier */
    struct statfs fsbuf;
#endif
#if defined (HAVE_STRUCT_STAT_ST_FSTYPE)
    struct stat sbuf;
#endif

    *file_id = UNKNOWN_SUPER_MAGIC;

#ifdef HAVE_STRUCT_STATVFS_WITH_F_BASETYPE
    err = statvfs(filename, &vfsbuf);
    if (err == 0)
        *file_id = vfsbuf.f_basetype;
#endif

/* remember above how I said 'statfs with f_type' was the common linux-y way to
 * report file system type?  Darwin (and probably the BSDs) *also* uses f_type
 * but it is "reserved" and does not give us anything meaningful.  Fine.  If
 * configure detects f_type we'll use it here and on those "reserved" platforms
 * we'll ignore that result and check the f_fstypename field  */
#ifdef HAVE_STRUCT_STATFS_F_TYPE
    err = statfs(filename, &fsbuf);
    if (err == 0)
        *file_id = fsbuf.f_type;
#endif


#if defined(HAVE_STRUCT_STATFS_F_FSTYPENAME) || defined(HAVE_STRUCT_STAT_ST_FSTYPE)
    /* these stat routines store the file system type in a string */
    char *fstype;
#ifdef HAVE_STRUCT_STATFS_F_FSTYPENAME
    err = statfs(filename, &fsbuf);
    fstype = fsbuf.f_fstypename;
#else
    err = stat(filename, &sbuf);
    fstype = sbuf.st_fstype;
#endif
    if (err == 0) {
        int i = 0;
        /* any file system type not explicitly in the fstype table (ffs, hfs)
         * will be "unknown" which ROMIO will service with ADIO_UFS */
        while (fstypes[i].fileops) {
            /* '-1' to ignore the trailing colon */
            if (!strncasecmp(fstypes[i].prefix, fstype, strlen(fstypes[i].prefix) - 1)) {
                *file_id = fstypes[i].magic;
            }
            i++;
        }
    }
#endif

    return err;
}

/* Check if file system is Lustre from file name, using a system-dependent
 * function call. This is a collective call.
 */
int PNC_Check_Lustre(MPI_Comm comm,
                     const char *filename)
{
    int err, retry_cnt;
    int64_t file_id=UNKNOWN_SUPER_MAGIC;

/* NFS can get stuck and end up returning ESTALE "forever" */
#define MAX_ESTALE_RETRY 10000

    retry_cnt = 0;
    do {
        err = romio_statfs(filename, &file_id);
    } while (err && (errno == ESTALE) && retry_cnt++ < MAX_ESTALE_RETRY);

    if (err) {
        /* ENOENT may be returned in two cases:
         * 1) no directory entry for "filename"
         * 2) "filename" is a dangling symbolic link
         *
         * ADIO_FileSysType_parentdir tries to deal with both cases.
         */
        if (errno == ENOENT) {
            char *dir;
            ADIO_FileSysType_parentdir(filename, &dir);
            err = romio_statfs(dir, &file_id);
            ADIOI_Free(dir);
        } else
            return 0;
    }

    return ((file_id == LL_SUPER_MAGIC));
}

/*
 *   1. root creates/opens the file
 *      a. root sets/obtains striping info
 *      b. root closes file
 *   2. root determines cb_nodes and ranklist
 *      a. all processes send root its proc_name
 *      b. root broadcasts cb_nodes and ranklist
 *   3. When deferred_open is true:
 *      then only aggregators open the file
 *      else only aggregators open the file
 *   4. all processes sync and set file striping info
 *      a. root bcasts striping info
 *      b. all processes set hints
 */
static int
PNC_Lustre_file_create(ADIO_File fd,
                       int       access_mode)
{
    int err=0, rank, amode, perm, old_mask;
    int stripin_info[4];

    MPI_Comm_rank(fd->comm, &rank);

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ 0666;

    /* root process creates the file first, followed by all processes open the
     * file.
     * For Lustre, we need to obtain file striping info (striping_factor,
     * striping_unit, and num_osts) in order to select the I/O aggregators
     * in fd->hints->ranklist, no matter its is open or create mode.
     */
    if (rank == 0) {
        char value[MPI_MAX_INFO_VAL+1];
        int lumlen, flag, set_layout = 0;
#ifndef MIMIC_LUSTRE
        struct lov_user_md *lum = NULL;
#endif
        ADIO_Offset str_factor = -1, str_unit = 0, start_iodev = -1;

        /* In LUSTRE_SetInfo, these hints have been validated to be consistent
         * among all processes.
         */
        if (fd->info != MPI_INFO_NULL) {
            /* get striping information from user hints */
            ADIOI_Info_get(fd->info, "striping_unit", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                str_unit = atoll(value);

            ADIOI_Info_get(fd->info, "striping_factor", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                str_factor = atoll(value);

            ADIOI_Info_get(fd->info, "start_iodevice", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                start_iodev = atoll(value);
        }
        if ((str_factor > 0) || (str_unit > 0) || (start_iodev >= 0))
            set_layout = 1;

        amode = O_CREAT;
        if (access_mode & MPI_MODE_RDWR) amode |= O_RDWR;

#ifndef MIMIC_LUSTRE
        /* if hints were set, we need to delay creation of any lustre objects.
         * However, if we open the file with O_LOV_DELAY_CREATE and don't call
         * the follow-up ioctl, subsequent writes will fail
         */
        if (set_layout)
            amode = O_CREAT | O_LOV_DELAY_CREATE;
#endif
        fd->fd_sys = open(fd->filename, amode, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d fails to create file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            err = errno;
            goto err_out;
        }

#ifdef MIMIC_LUSTRE
#define xstr(s) str(s)
#define str(s) #s
#define STRIPE_SIZE 1024
#define STRIPE_COUNT 4
        stripin_info[0] = STRIPE_SIZE;
        stripin_info[1] = STRIPE_COUNT;
        stripin_info[2] = 0;
        stripin_info[3] = STRIPE_COUNT;
#else
        /* odd length here because lov_user_md contains some fixed data and
         * then a list of 'lmm_objects' representing stripe */
        lumlen = sizeof(struct lov_user_md)
               + MAX_LOV_UUID_COUNT * sizeof(struct lov_user_ost_data);
        lum = (struct lov_user_md *) ADIOI_Calloc(1, lumlen);

        /* we can only set these hints on new files */
        if (set_layout) {
            lum->lmm_magic = LOV_USER_MAGIC;
            lum->lmm_pattern = 0;
            /* crude check for overflow of lustre internal datatypes.
             * Silently cap to large value if user provides a value
             * larger than lustre supports
             */
            if (str_unit > UINT_MAX)
                lum->lmm_stripe_size = UINT_MAX;
            else
                /* cast from 64 to 32 ok: we checked/set above */
                lum->lmm_stripe_size = (unsigned int) str_unit;

            if (str_factor > USHRT_MAX)
                lum->lmm_stripe_count = USHRT_MAX;
            else
                lum->lmm_stripe_count = str_factor;

            if (start_iodev > USHRT_MAX)
                lum->lmm_stripe_offset = USHRT_MAX;
            else
                lum->lmm_stripe_offset = start_iodev;

            /* set Lustre file stripning */
            err = ioctl(fd->fd_sys, LL_IOC_LOV_SETSTRIPE, lum);
            if (err == -1 && errno != EEXIST) {
                fprintf(stderr,"%s line %d: rank %d fails to set stripe info from file %s (%s)\n",
                        __func__,__LINE__,rank, fd->filename, strerror(errno));
                /* not a fatal error, but user might care to know */
            }
        }

        /* get Lustre file stripning, even if setting it failed */
        memset(lum, 0, lumlen);
        lum->lmm_magic = LOV_USER_MAGIC;
        err = ioctl(fd->fd_sys, LL_IOC_LOV_GETSTRIPE, (void *) lum);
        if (err == 0) {
            stripin_info[0] = lum->lmm_stripe_size;
            stripin_info[1] = lum->lmm_stripe_count;
            stripin_info[2] = lum->lmm_stripe_offset;
            /* TODO: figure out the correct way to find the number of unique
             * OSTs. This is relevant only when Lustre over-striping is used.
             * For now, set it to same as fd->hints->striping_factor.
             fd->hints->fs_hints.lustre.num_osts = num_uniq_osts(fd->filename);
             */
            stripin_info[3] = fd->hints->striping_factor;
        }
        else {
            fprintf(stderr,"%s line %d: rank %d fails to get stripe info from file %s (%s)\n",
                    __func__,__LINE__,rank, fd->filename, strerror(errno));
            err = errno;
        }
        ADIOI_Free(lum);
#endif
    }
err_out:
    if (err != 0) {
        stripin_info[0] = -1;
        stripin_info[1] = -1;
        stripin_info[2] = -1;
        stripin_info[3] = -1;
    }
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);
    if (stripin_info[0] == -1) {
        fprintf(stderr, "%s line %d: failed to create file %s\n", __func__,
                __LINE__, fd->filename);
        return err;
    }

    fd->hints->striping_unit            = stripin_info[0];
    fd->hints->striping_factor          = stripin_info[1];
    fd->hints->start_iodevice           = stripin_info[2];
    fd->hints->fs_hints.lustre.num_osts = stripin_info[3];

    if (rank > 0) { /* non-root processes */
        fd->fd_sys = open(fd->filename, O_RDWR, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d failure to open file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            return errno;
        }
    }
    return 0;
}

static int
PNC_Lustre_file_open(ADIO_File fd)
{
    int err=0, rank, amode, perm, old_mask;
    int stripin_info[4];

    MPI_Comm_rank(fd->comm, &rank);

    stripin_info[0] = -1;
    stripin_info[1] = -1;
    stripin_info[2] = -1;
    stripin_info[3] = -1;

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ 0666;

    /* All processes open the file. */
    fd->fd_sys = open(fd->filename, O_RDWR, perm);
    if (fd->fd_sys == -1) {
        fprintf(stderr, "%s line %d: rank %d failure to open file %s (%s)\n",
                __func__,__LINE__, rank, fd->filename, strerror(errno));
        err = errno;
        goto err_out;
    }

    /* Only root obtains the striping information and bcast to all other
     * processes.
     */
    if (rank == 0) {
#ifndef MIMIC_LUSTRE
        ADIO_Offset str_factor = -1, str_unit = 0, start_iodev = -1;
        struct lov_user_md *lum = NULL;

        lumlen = sizeof(struct lov_user_md)
               + MAX_LOV_UUID_COUNT * sizeof(struct lov_user_ost_data);
        lum = (struct lov_user_md *) ADIOI_Calloc(1, lumlen);

        /* get Lustre file stripning, even if setting it failed */
        lum->lmm_magic = LOV_USER_MAGIC;
        err = ioctl(fd->fd_sys, LL_IOC_LOV_GETSTRIPE, (void *) lum);
        if (!err) {
            stripin_info[0] = lum->lmm_stripe_size;
            stripin_info[1] = lum->lmm_stripe_count;
            stripin_info[2] = lum->lmm_stripe_offset;
            /* TODO: figure out the correct way to find the number of unique
             * OSTs. This is relevant only when Lustre over-striping is used.
             * For now, set it to same as fd->hints->striping_factor.
             fd->hints->fs_hints.lustre.num_osts = num_uniq_osts(fd->filename);
             */
            stripin_info[3] = fd->hints->striping_factor;
        }
        ADIOI_Free(lum);
#else
#define xstr(s) str(s)
#define str(s) #s
#define STRIPE_SIZE 1024
#define STRIPE_COUNT 4
        stripin_info[0] = STRIPE_SIZE;
        stripin_info[1] = STRIPE_COUNT;
        stripin_info[2] = 0;
        stripin_info[3] = STRIPE_COUNT;
#endif
    }

err_out:
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);
    fd->hints->striping_unit            = stripin_info[0];
    fd->hints->striping_factor          = stripin_info[1];
    fd->hints->start_iodevice           = stripin_info[2];
    fd->hints->fs_hints.lustre.num_osts = stripin_info[3];

    return err;
}

int PNC_Lustre_open(MPI_Comm    comm,
                    const char *filename,
                    int         amode,
                    MPI_Info    info,
                    PNC_File   *fd)
{
    /* Before reaching to this subroutine, PNC_Check_Lustre() should have been
     * called to verify filename is on Lustre.
     */

(MPI_Comm orig_comm,
                  MPI_Comm comm, const char *filename, int file_system,
                  ADIOI_Fns * ops,
                  int access_mode, PNC_Offset disp, MPI_Datatype etype,
                  MPI_Datatype filetype, MPI_Info info, int perm, int *error_code)
{
    MPI_File mpi_fh;
    PNC_File fd;
    int err, rank, procs;
    int max_error_code;
    MPI_Info dupinfo;
    int syshints_processed, can_skip;
    char *p;

    fd->comm = comm;    /* dup'ed in MPI_File_open */
    fd->filename = ADIOI_Strdup(filename);
    fd->atomicity = 0;
    fd->etype = MPI_BYTE;
    fd->etype_size = 1;
    fd->filetype = MPI_BYTE;
    fd->file_system = ADIO_LUSTRE;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);

    /* create and initialize info object */
    fd->hints = (PNC_Hints*) NCI_Calloc(1, sizeof(PNC_Hints));

    fd->hints->initialized = 0;
    fd->hints->cb_config_list = NULL;
    fd->hints->ranklist = NULL;
    fd->info = MPI_INFO_NULL;

    /* collective buffer */
    fd->io_buf = ADIOI_Malloc(fd->hints->cb_buffer_size);

    /* For Lustre, determining the I/O aggregators and constructing ranklist
     * requires file stripe count, which can only be obtained after file is
     * opened.
     *
     * fd->is_agg and fd->my_cb_nodes_index are set in construct_aggr_list()
     */

    ADIOI_OpenColl(fd, rank, access_mode, error_code);

    /* deferred open consideration: if an independent process lied about
     * "no_indep_rw" and opens the file later (example: HDF5 uses independent
     * i/o for metadata), that deferred open will use the access_mode provided
     * by the user.  CREATE|EXCL only makes sense here -- exclusive access in
     * the deferred open case is going to fail and surprise the user.  Turn off
     * the excl amode bit. Save user's amode for MPI_FILE_GET_AMODE */
    fd->orig_access_mode = access_mode;
    if (fd->access_mode & PNC_EXCL)
        fd->access_mode ^= PNC_EXCL;

  fn_exit:
    MPI_Allreduce(error_code, &max_error_code, 1, MPI_INT, MPI_MAX, comm);
    if (max_error_code != MPI_SUCCESS) {

        /* If the file was successfully opened, close it */
        if (*error_code == MPI_SUCCESS) {

            /* in the deferred open case, only those who have actually
             * opened the file should close it */
            if (fd->hints->deferred_open) {
                if (fd->is_agg) {
                    (*(fd->fns->ADIOI_xxx_Close)) (fd, error_code);
                }
            } else {
                (*(fd->fns->ADIOI_xxx_Close)) (fd, error_code);
            }
        }
        ADIOI_Free(fd->filename);
        if (fd->hints->ranklist != NULL)
            ADIOI_Free(fd->hints->ranklist);
        if (fd->hints->cb_config_list != NULL)
            ADIOI_Free(fd->hints->cb_config_list);
        ADIOI_Free(fd->hints);
        if (fd->info != MPI_INFO_NULL)
            MPI_Info_free(&(fd->info));
        ADIOI_Free(fd->io_buf);
        ADIOI_Free(fd);
        fd = PNC_FILE_NULL;
        if (*error_code == MPI_SUCCESS) {
            *error_code = MPIO_Err_create_code(MPI_SUCCESS,
                                               MPIR_ERR_RECOVERABLE, myname,
                                               __LINE__, MPI_ERR_IO, "**oremote_fail", 0);
        }
    }

    return fd;
}


static int construct_aggr_list(ADIO_File fd, int root, int *error_code);

            tmp_comm = fd->comm;
            fd->comm = MPI_COMM_SELF;
            (*(fd->fns->ADIOI_xxx_Open)) (fd, error_code);
            fd->comm = tmp_comm;
            MPI_Bcast(error_code, 1, MPI_INT, root, fd->comm);
            /* if no error, close the file and reopen normally below */
            if (*error_code == MPI_SUCCESS)
                (*(fd->fns->ADIOI_xxx_Close)) (fd, error_code);

            fd->access_mode = access_mode;      /* back to original */
        } else
            MPI_Bcast(error_code, 1, MPI_INT, root, fd->comm);

        if (*error_code != MPI_SUCCESS) {
            return;
        } else {
            /* turn off CREAT (and EXCL if set) for real multi-processor open */
            access_mode ^= ADIO_CREATE;
            if (access_mode & ADIO_EXCL)
                access_mode ^= ADIO_EXCL;
        }

        if (fd->file_system == ADIO_LUSTRE) {
            /* Use file striping count and number of unique OSTs to construct
             * the I/O aggregator rank list, fd->hints->ranklist[].
             */
            construct_aggr_list(fd, root, error_code);
            if (*error_code != MPI_SUCCESS)
                return;
        }
    }

    fd->blksize = 1024 * 1024 * 4;
    /* this large default value should be good for most file systems. any ROMIO
     * driver is free to stat the file and find an optimal value */

    /* add to fd->info the hint "aggr_list", list of aggregators' rank IDs */
    value[0] = '\0';
    for (i = 0; i < fd->hints->cb_nodes; i++) {
        char str[16];
        if (i == 0)
            snprintf(str, sizeof(str), "%d", fd->hints->ranklist[i]);
        else
            snprintf(str, sizeof(str), " %d", fd->hints->ranklist[i]);
        if (strlen(value) + strlen(str) >= MPI_MAX_INFO_VAL-5) {
            strcat(value, " ...");
            break;
        }
        strcat(value, str);
    }
    ADIOI_Info_set(fd->info, "aggr_list", value);

    /* if we are doing deferred open, non-aggregators should return now */
    if (fd->hints->deferred_open) {
        if (!(fd->is_agg)) {
            /* we might have turned off EXCL for the aggregators.
             * restore access_mode that non-aggregators get the right
             * value from get_amode */
            fd->access_mode = orig_amode_excl;

            /* In file-system specific open, a driver might collect some
             * information via stat().  Deferred open means not every process
             * participates in fs-specific open, but they all participate in
             * this open call.  Broadcast a bit of information in case
             * lower-level file system driver (e.g. 'bluegene') collected it
             * (not all do)*/
            stats_type = make_stats_type(fd);
            MPI_Bcast(MPI_BOTTOM, 1, stats_type, root, fd->comm);
            ADIOI_Assert(fd->blksize > 0);

            /* set file striping hints */
            snprintf(value, sizeof(value), "%d", fd->hints->striping_unit);
            ADIOI_Info_set(fd->info, "striping_unit", value);

            snprintf(value, sizeof(value), "%d", fd->hints->striping_factor);
            ADIOI_Info_set(fd->info, "striping_factor", value);

            snprintf(value, sizeof(value), "%d", fd->hints->start_iodevice);
            ADIOI_Info_set(fd->info, "start_iodevice", value);

            *error_code = MPI_SUCCESS;
            MPI_Type_free(&stats_type);
            return;
        }
    }

/* For writing with data sieving, a read-modify-write is needed. If
   the file is opened for write_only, the read will fail. Therefore,
   if write_only, open the file as read_write, but record it as write_only
   in fd, so that get_amode returns the right answer. */

    /* observation from David Knaak: file systems that do not support data
     * sieving do not need to change the mode */

    orig_amode_wronly = access_mode;
    if ((access_mode & ADIO_WRONLY) && ADIO_Feature(fd, ADIO_DATA_SIEVING_WRITES)) {
        access_mode = access_mode ^ ADIO_WRONLY;
        access_mode = access_mode | ADIO_RDWR;
    }
    fd->access_mode = access_mode;

    (*(fd->fns->ADIOI_xxx_Open)) (fd, error_code);

    /* if error, may be it was due to the change in amode above.
     * therefore, reopen with access mode provided by the user. */
    fd->access_mode = orig_amode_wronly;
    if (*error_code != MPI_SUCCESS)
        (*(fd->fns->ADIOI_xxx_Open)) (fd, error_code);

    /* if we turned off EXCL earlier, then we should turn it back on */
    if (fd->access_mode != orig_amode_excl)
        fd->access_mode = orig_amode_excl;

    /* broadcast information to all processes in communicator, not just
     * those who participated in open.
     */
    stats_type = make_stats_type(fd);
    MPI_Bcast(MPI_BOTTOM, 1, stats_type, root, fd->comm);
    MPI_Type_free(&stats_type);
    /* file domain code will get terribly confused in a hard-to-debug way
     * if gpfs blocksize not sensible */
    ADIOI_Assert(fd->blksize > 0);

    /* set file striping hints */
    snprintf(value, sizeof(value), "%d", fd->hints->striping_unit);
    ADIOI_Info_set(fd->info, "striping_unit", value);

    snprintf(value, sizeof(value), "%d", fd->hints->striping_factor);
    ADIOI_Info_set(fd->info, "striping_factor", value);

    snprintf(value, sizeof(value), "%d", fd->hints->start_iodevice);
    ADIOI_Info_set(fd->info, "start_iodevice", value);

    /* for deferred open: this process has opened the file (because if we are
     * not an aggregator and we are doing deferred open, we returned earlier)*/
    fd->is_open = 1;

    /* sync optimization: we can omit the fsync() call if we do no writes */
    fd->dirty_write = 0;
}

/*----< construct_aggr_list() >----------------------------------------------*/
/* Allocate and construct fd->hints->ranklist[].
 * Overwrite fd->hints->cb_nodes and set hint cb_nodes.
 * Set fd->is_agg, whether this rank is an I/O aggregator
 *     fd->hints->cb_nodes
 *     fd->hints->fs_hints.lustre.num_osts
 */
static int construct_aggr_list(ADIO_File fd, int root, int *error_code)
{
    int i, j, k, rank, nprocs, num_aggr, my_procname_len, num_nodes;
    int msg[3], striping_factor;
    int *all_procname_lens = NULL;
    int *nprocs_per_node, **ranks_per_node;
    char value[MPI_MAX_INFO_VAL + 1], my_procname[MPI_MAX_PROCESSOR_NAME];
    char **all_procnames = NULL;
    static char myname[] = "ADIO_OPENCOLL construct_aggr_list";

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &rank);

    /* Collect info about compute nodes in order to select I/O aggregators.
     * Note my_procname is null character terminated, but my_procname_len
     * does not include the null character.
     */
#ifdef MIMIC_LUSTRE
    /* mimic number of compute nodes = 4 */
    int node_id, np_per_node = nprocs / 4;
    if (nprocs % 4 > 0) np_per_node++;
    if (rank < np_per_node * (nprocs % 4))
        node_id = rank / np_per_node;
    else
        node_id = (rank - np_per_node * (nprocs % 4)) / (nprocs / 4) + (nprocs % 4);

    sprintf(my_procname,"compute.node.%d", node_id);
    my_procname_len = strlen(my_procname);
#else
    MPI_Get_processor_name(my_procname, &my_procname_len);
#endif
    my_procname_len++; /* to include terminate null character */

    if (rank == root) {
        /* root collects all procnames */
        all_procnames = (char **) ADIOI_Malloc(nprocs * sizeof(char *));
        if (all_procnames == NULL) {
            *error_code = MPIO_Err_create_code(*error_code,
                                               MPIR_ERR_RECOVERABLE, myname,
                                               __LINE__, MPI_ERR_OTHER,
                                               "**nomem2", 0);
            return 0;
        }

        all_procname_lens = (int *) ADIOI_Malloc(nprocs * sizeof(int));
        if (all_procname_lens == NULL) {
            ADIOI_Free(all_procnames);
            *error_code = MPIO_Err_create_code(*error_code,
                                               MPIR_ERR_RECOVERABLE, myname,
                                                __LINE__, MPI_ERR_OTHER,
                                                "**nomem2", 0);
            return 0;
        }
    }
    /* gather process name lengths from all processes first */
    MPI_Gather(&my_procname_len, 1, MPI_INT, all_procname_lens, 1, MPI_INT,
               root, fd->comm);

    if (rank == root) {
        int *disp;
        size_t alloc_size = 0;

        for (i = 0; i < nprocs; i++)
            alloc_size += all_procname_lens[i];

        all_procnames[0] = (char *) ADIOI_Malloc(alloc_size);
        if (all_procnames[0] == NULL) {
            ADIOI_Free(all_procname_lens);
            ADIOI_Free(all_procnames);
            *error_code = MPIO_Err_create_code(*error_code,
                                               MPIR_ERR_RECOVERABLE, myname,
                                               __LINE__, MPI_ERR_OTHER,
                                               "**nomem2", 0);
            return 0;
        }

        /* Construct displacement array for the MPI_Gatherv, as each process
         * may have a different length for its process name.
         */
        disp = (int *) ADIOI_Malloc(nprocs * sizeof(int));
        disp[0] = 0;
        for (i = 1; i < nprocs; i++) {
            all_procnames[i] = all_procnames[i - 1] + all_procname_lens[i - 1];
            disp[i] = disp[i - 1] + all_procname_lens[i - 1];
        }

        /* gather all process names */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    all_procnames[0], all_procname_lens, disp, MPI_CHAR,
                    root, fd->comm);

        ADIOI_Free(disp);
        ADIOI_Free(all_procname_lens);
    } else
        /* send process name to root */
        MPI_Gatherv(my_procname, my_procname_len, MPI_CHAR,
                    NULL, NULL, NULL, MPI_CHAR, root, fd->comm);

    if (rank == root) {
        /* all_procnames[] can tell us the number of nodes and number of
         * processes per node.
         */
        char **node_names;
        int last, *node_ids;

        /* number of MPI processes running on each node */
        nprocs_per_node = (int *) ADIOI_Malloc(nprocs * sizeof(int));

        /* ech MPI process's compute node ID */
        node_ids = (int *) ADIOI_Malloc(nprocs * sizeof(int));

        /* array of pointers pointing to unique host names (compute nodes) */
        node_names = (char **) ADIOI_Malloc(nprocs * sizeof(char *));

        /* calculate nprocs_per_node[] and node_ids[] */
        last = 0;
        num_nodes = 0; /* number of unique compute nodes */
        for (i = 0; i < nprocs; i++) {
            k = last;
            for (j = 0; j < num_nodes; j++) {
                /* check if [i] has already appeared in [] */
                if (!strcmp(all_procnames[i], node_names[k])) { /* found */
                    node_ids[i] = k;
                    nprocs_per_node[k]++;
                    break;
                }
                k = (k == num_nodes - 1) ? 0 : k + 1;
            }
            if (j < num_nodes)  /* found, next iteration, start with node n */
                last = k;
            else {      /* not found, j == num_nodes, add a new node */
                node_names[j] = ADIOI_Strdup(all_procnames[i]);
                nprocs_per_node[j] = 1;
                node_ids[i] = j;
                last = j;
                num_nodes++;
            }
        }
        /* num_nodes is now the number of compute nodes (unique node names) */
        fd->num_nodes = num_nodes;

        for (i = 0; i < num_nodes; i++)
            ADIOI_Free(node_names[i]);
        ADIOI_Free(node_names);
        ADIOI_Free(all_procnames[0]);
        ADIOI_Free(all_procnames);

        /* construct rank IDs of MPI processes running on each node */
        ranks_per_node = (int **) ADIOI_Malloc(num_nodes * sizeof(int *));
        ranks_per_node[0] = (int *) ADIOI_Malloc(nprocs * sizeof(int));
        for (i = 1; i < num_nodes; i++)
            ranks_per_node[i] = ranks_per_node[i - 1] + nprocs_per_node[i - 1];
        for (i = 0; i < num_nodes; i++)
            nprocs_per_node[i] = 0;

        /* Populate ranks_per_node[], list of MPI ranks running on each node.
         * Populate nprocs_per_node[], number of MPI processes on each node.
         */
        for (i = 0; i < nprocs; i++) {
            k = node_ids[i];
            ranks_per_node[k][nprocs_per_node[k]] = i;
            nprocs_per_node[k]++;
        }
        ADIOI_Free(node_ids);

        /* Given the number of nodes, num_nodes, and processes per node,
         * nprocs_per_node, we can now set num_aggr, the number of I/O
         * aggregators. At this moment, root should have obtained the file
         * striping settings.
         */
        striping_factor = fd->hints->striping_factor;

        if (striping_factor > nprocs) {
            /* When number of MPI processes is less than striping_factor, set
             * num_aggr to the max number less than nprocs that divides
             * striping_factor. An naive way is:
             *     num_aggr = nprocs;
             *     while (striping_factor % num_aggr > 0)
             *         num_aggr--;
             * Below is equivalent, but faster.
             */
            int divisor = 2;
            num_aggr = 1;
            /* try to divide */
            while (striping_factor >= divisor * divisor) {
                if ((striping_factor % divisor) == 0) {
                    if (striping_factor / divisor <= nprocs) {
                        /* The value is found ! */
                        num_aggr = striping_factor / divisor;
                        break;
                    }
                    /* if divisor is less than nprocs, divisor is a solution,
                     * but it is not sure that it is the best one
                     */
                    else if (divisor <= nprocs)
                        num_aggr = divisor;
                }
                divisor++;
            }
        }
        else { /* striping_factor <= nprocs */
            /* Select striping_factor processes to be I/O aggregators */

            if (fd->hints->cb_nodes == 0 || fd->access_mode & ADIO_RDONLY) {
                /* hint cb_nodes is not set by user and this file is opened for
                 * read only. Because collective read is using a different file
                 * domain partitioning strategy, for now we do not mess up
                 * ranklist for read operations and do not accept user hint to
                 * set ranklist for read operations
                 */
                num_aggr = striping_factor;
            }
            else if (fd->hints->cb_nodes <= striping_factor) {
                /* User has set hint cb_nodes and cb_nodes <= striping_factor.
                 * Ignore user's hint and try to set cb_nodes to be at least
                 * striping_factor.
                 */
                num_aggr = striping_factor;
            }
            else {
                /* User has set hint cb_nodes and cb_nodes > striping_factor */
                if (nprocs < fd->hints->cb_nodes)
                    num_aggr = nprocs; /* BAD cb_nodes set by users */
                else
                    num_aggr = fd->hints->cb_nodes;

                /* Number of processes per node may not be enough to be picked
                 * as aggregators. If this case, reduce num_aggr (cb_nodes).
                 * Consider the following case: number of processes = 18,
                 * number of nodes = 7, striping_factor = 8, cb_nodes = 16.
                 * cb_nodes should be reduced to 8 and the ranks of aggregators
                 * should be 0, 3, 6, 9, 12, 14, 16, 1.
                 * If the number of processes changes to 25, then cb_nodes
                 * should be 16 and the ranks of aggregators should be 0, 4, 8,
                 * 12, 16, 19, 22, 1, 2, 6, 10, 14, 18, 21, 24, 3.
                 */
                int max_nprocs_node = 0;
                for (i = 0; i < num_nodes; i++)
                    max_nprocs_node = MAX(max_nprocs_node, nprocs_per_node[i]);
                int max_naggr_node = striping_factor / num_nodes;
                if (striping_factor % num_nodes) max_naggr_node++;
                /* max_naggr_node is the max number of processes per node to be
                 * picked as aggregator in each round.
                 */
                int rounds = num_aggr / striping_factor;
                if (num_aggr % striping_factor) rounds++;
                while (max_naggr_node * rounds > max_nprocs_node) rounds--;
                num_aggr = striping_factor * rounds;
            }
        }

        /* TODO: the above setting for num_aggr is for collective writes. Reads
         * should be the number of nodes.
         */

        /* Next step is to determine the MPI rank IDs of I/O aggregators into
         * ranklist[].
         */
        fd->hints->ranklist = (int *) ADIOI_Malloc(num_aggr * sizeof(int));
        if (fd->hints->ranklist == NULL) {
            *error_code = MPIO_Err_create_code(*error_code,
                                               MPIR_ERR_RECOVERABLE, myname,
                                               __LINE__, MPI_ERR_OTHER,
                                               "**nomem2", 0);
            return 0;
        }

        if (striping_factor <= num_nodes) {
            /* When number of OSTs is less than number of compute nodes,
             * first select number of nodes equal to the number of OSTs by
             * spread the selection evenly across all compute nodes and then
             * pick processes from the selected nodes, also evenly spread
             * evenly among processes on each selected node to be aggregators.
             */
            int avg = num_aggr / striping_factor;
            int stride = num_nodes / striping_factor;
            if (num_aggr % striping_factor) avg++;
            for (i = 0; i < num_aggr; i++) {
                j = (i % striping_factor) * stride; /* to select from node j */
                k = (i / striping_factor) * (nprocs_per_node[j] / avg);
                assert(k < nprocs_per_node[j]);
                fd->hints->ranklist[i] = ranks_per_node[j][k];
            }
        }
        else { /* striping_factor > num_nodes */
            /* When number of OSTs is more than number of compute nodes, I/O
             * aggregators are selected from all nodes are selected. Within
             * each node, aggregators are spread evenly instead of the first
             * few ranks.
             */
            int *naggr_per_node, *idx_per_node, avg;
            idx_per_node = (int*) ADIOI_Calloc(num_nodes, sizeof(int));
            naggr_per_node = (int*) ADIOI_Malloc(num_nodes * sizeof(int));
            for (i = 0; i < striping_factor % num_nodes; i++)
                naggr_per_node[i] = striping_factor / num_nodes + 1;
            for (; i < num_nodes; i++)
                naggr_per_node[i] = striping_factor / num_nodes;
            avg = num_aggr / striping_factor;
            if (avg > 0)
                for (i = 0; i < num_nodes; i++)
                    naggr_per_node[i] *= avg;
            for (i = 0; i < num_nodes; i++)
                naggr_per_node[i] = MIN(naggr_per_node[i], nprocs_per_node[i]);
            /* naggr_per_node[] is the number of aggregators that can be
             * selected as I/O aggregators
             */

#ifdef WKL_DEBUG
// printf("%d: num_nodes=%d nprocs_per_node[0]=%d num_aggr=%d naggr_per_node[0]=%d\n",rank,num_nodes,nprocs_per_node[0],num_aggr,naggr_per_node[0]);
#endif
            for (i = 0; i < num_aggr; i++) {
                int stripe_i = i % striping_factor;
                j = stripe_i % num_nodes; /* to select from node j */
                k = nprocs_per_node[j] / naggr_per_node[j];
                k *= idx_per_node[j];
                idx_per_node[j]++;
                assert(k < nprocs_per_node[j]);
                fd->hints->ranklist[i] = ranks_per_node[j][k];
            }
            ADIOI_Free(naggr_per_node);
            ADIOI_Free(idx_per_node);
        }
#ifdef WKL_DEBUG
// printf("%d: num_nodes=%d nprocs_per_node[0]=%d num_aggr=%d ranklist=%d %d %d %d (%s)\n",rank,num_nodes,nprocs_per_node[0],num_aggr,fd->hints->ranklist[0],fd->hints->ranklist[1],fd->hints->ranklist[2],fd->hints->ranklist[3],fd->filename); fflush(stdout);
#endif

        /* TODO: we can keep these two arrays in case for dynamic construction
         * of fd->hints->ranklist[], such as in group-cyclic file domain
         * assignment method, used in each collective write call.
         */
        ADIOI_Free(nprocs_per_node);
        ADIOI_Free(ranks_per_node[0]);
        ADIOI_Free(ranks_per_node);

        msg[0] = num_aggr;
        msg[1] = fd->hints->fs_hints.lustre.num_osts;
        msg[2] = num_nodes;
    }

    /* bcast cb_nodes and lustre.num_osts to all processes */
    MPI_Bcast(msg, 3, MPI_INT, root, fd->comm);
    num_aggr = msg[0];

#ifdef WKL_DEBUG
    if (rank == root && fd->hints->cb_nodes > 0 && fd->hints->cb_nodes != num_aggr) {
        /* user has set hint cb_nodes and the value is not num_aggr */
        printf("Warning: %s line %d: Set cb_nodes to %d and ignore user's hint of %d\n",
               __func__, __LINE__, num_aggr, fd->hints->cb_nodes);
    }
#endif

    /* set file striping hints */
    fd->hints->cb_nodes = num_aggr;
    sprintf(value, "%d", fd->hints->cb_nodes);
    ADIOI_Info_set(fd->info, "cb_nodes", value);

    fd->hints->fs_hints.lustre.num_osts = msg[1];
    sprintf(value, "%d", fd->hints->fs_hints.lustre.num_osts);
    ADIOI_Info_set(fd->info, "lustre_num_osts", value);

    /* save the number of unique compute nodes in communicator fd->comm */
    fd->num_nodes = msg[2];

    if (rank != root) {
        /* ranklist[] contains the MPI ranks of I/O aggregators */
        fd->hints->ranklist = (int *) ADIOI_Malloc(num_aggr * sizeof(int));
        if (fd->hints->ranklist == NULL) {
            *error_code = MPIO_Err_create_code(*error_code,
                                               MPIR_ERR_RECOVERABLE, myname,
                                               __LINE__, MPI_ERR_OTHER,
                                               "**nomem2", 0);
            return 0;
        }
    }

    MPI_Bcast(fd->hints->ranklist, fd->hints->cb_nodes, MPI_INT, root, fd->comm);

    /* check whether this process is selected as an I/O aggregator */
    fd->is_agg = 0;
    fd->my_cb_nodes_index = -1;
    for (i = 0; i < num_aggr; i++) {
        if (rank == fd->hints->ranklist[i]) {
            fd->is_agg = 1;
            fd->my_cb_nodes_index = i;
            break;
        }
    }

    return 0;
}

