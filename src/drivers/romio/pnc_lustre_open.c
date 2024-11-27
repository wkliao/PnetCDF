/*
 *  Copyright (C) 2025, Northwestern University
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>   /* readlink() */
#include <string.h>   /* strdup() */
#include <assert.h>
#include <sys/errno.h>
#include <fcntl.h>    /* O_CREAT */

#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#ifndef PATH_MAX
#define PATH_MAX 65535
#endif

#ifdef HAVE_SYS_VFS_H
#include <sys/vfs.h>
#endif
#ifdef HAVE_SYS_STATVFS_H
#include <sys/statvfs.h>
#endif
#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h> /* struct statfs */
#endif
#ifdef HAVE_SYS_MOUNT_H
#include <sys/mount.h> /* struct statfs */
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h> /* fstat(), lstat(), stat() */
#endif

#ifndef MIMIC_LUSTRE
#include <lustre/lustreapi.h>

/* what is the basis for this define?
 * what happens if there are more than 1k UUIDs? */
#define MAX_LOV_UUID_COUNT      1000
#endif

#include <mpi.h>

#include "adio.h"

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
    if (err == 0 && !strncasecmp(fstype, "lustre", 6))
        *file_id = LL_SUPER_MAGIC;
#endif
    return err;
}

/* Check if file system is Lustre from file name, using a system-dependent
 * function call. This is a collective call.
 */
int PNC_Check_Lustre(const char *filename)
{
#ifdef MIMIC_LUSTRE
    return 1;
#endif
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

/*----< lustre_file_create() >-----------------------------------------------*/
/*   1. root creates the file
 *   2. root sets and obtains striping info
 *   3. root broadcasts striping info
 *   4. non-root processes receive striping info from root
 *   5. non-root processes opens the fie
 */
static int
lustre_file_create(ADIO_File fd,
                   int       access_mode)
{
    int err=NC_NOERR, rank, amode, perm, old_mask;
    int stripin_info[4] = {-1, -1, -1, -1};

    MPI_Comm_rank(fd->comm, &rank);

static int wkl=0; if (wkl == 0) {int rank; MPI_Comm_rank(fd->comm, &rank); if (rank == 0) printf("\n%s line %d: ---------%s\n",__func__,__LINE__,fd->filename); wkl++; }

    amode = O_CREAT;
    if (access_mode & MPI_MODE_RDWR)  amode |= O_RDWR;

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
#ifndef MIMIC_LUSTRE
        char value[MPI_MAX_INFO_VAL+1];
        int lumlen, flag, set_layout = 0;
        struct lov_user_md *lum = NULL;

        ADIO_Offset str_factor = -1, str_unit = 0, start_iodev = -1;

        /* In LUSTRE_SetInfo, these hints have been validated to be consistent
         * among all processes.
         */
        if (fd->info != MPI_INFO_NULL) {
            /* get striping information from user hints */
            MPI_Info_get(fd->info, "striping_unit", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                str_unit = atoll(value);

            MPI_Info_get(fd->info, "striping_factor", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                str_factor = atoll(value);

            MPI_Info_get(fd->info, "start_iodevice", MPI_MAX_INFO_VAL, value, &flag);
            if (flag)
                start_iodev = atoll(value);
        }
        if ((str_factor > 0) || (str_unit > 0) || (start_iodev >= 0))
            set_layout = 1;

#ifdef HAVE_LUSTRE
        /* if hints were set, we need to delay creation of any lustre objects.
         * However, if we open the file with O_LOV_DELAY_CREATE and don't call
         * the follow-up ioctl, subsequent writes will fail
         */
        if (set_layout)
            amode = O_CREAT | O_LOV_DELAY_CREATE;
#endif
#endif
        fd->fd_sys = open(fd->filename, amode, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d fails to create file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            err = ncmpii_error_posix2nc("open");
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
             stripin_info[3] = num_uniq_osts(fd->filename);
             */
            stripin_info[3] = fd->hints->striping_factor;
        }
        else {
            fprintf(stderr,"%s line %d: rank %d fails to get stripe info from file %s (%s)\n",
                    __func__,__LINE__,rank, fd->filename, strerror(errno));
            err = ncmpii_error_posix2nc("ioctl");
        }
        ADIOI_Free(lum);
#endif
    }
err_out:
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);
    if (stripin_info[0] == -1) {
        fprintf(stderr, "%s line %d: failed to create file %s\n", __func__,
                __LINE__, fd->filename);
        return err;
    }

    fd->hints->striping_unit   = stripin_info[0];
    fd->hints->striping_factor = stripin_info[1];
    fd->hints->start_iodevice  = stripin_info[2];
    fd->hints->num_osts        = stripin_info[3];

    if (rank > 0) { /* non-root processes */
        fd->fd_sys = open(fd->filename, O_RDWR, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d failure to open file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            return ncmpii_error_posix2nc("ioctl");
        }
    }
    return err;
}

/*----< lustre_file_open() >-------------------------------------------------*/
/*   1. all processes open the file.
 *   2. root obtains striping info and broadcasts to all others
 */
static int
lustre_file_open(ADIO_File fd)
{
    int err=NC_NOERR, rank, perm, old_mask;
    int stripin_info[4] = {-1, -1, -1, -1};

    MPI_Comm_rank(fd->comm, &rank);

static int wkl=0; if (wkl == 0) {int rank; MPI_Comm_rank(fd->comm, &rank); if (rank == 0) printf("\n%s line %d: ---------%s\n",__func__,__LINE__,fd->filename); wkl++; }

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ 0666;

    /* All processes open the file. */
    fd->fd_sys = open(fd->filename, O_RDWR, perm);
    if (fd->fd_sys == -1) {
        fprintf(stderr, "%s line %d: rank %d failure to open file %s (%s)\n",
                __func__,__LINE__, rank, fd->filename, strerror(errno));
        err = ncmpii_error_posix2nc("open");
        goto err_out;
    }

    /* Only root obtains the striping information and bcast to all other
     * processes.
     */
    if (rank == 0) {
#ifndef MIMIC_LUSTRE
        ADIO_Offset str_factor = -1, str_unit = 0, start_iodev = -1;
        struct lov_user_md *lum = NULL;

        int lumlen = sizeof(struct lov_user_md)
                   + MAX_LOV_UUID_COUNT * sizeof(struct lov_user_ost_data);
        lum = (struct lov_user_md *) ADIOI_Calloc(1, lumlen);

        /* get Lustre file stripning, even if setting it failed */
        lum->lmm_magic = LOV_USER_MAGIC;
        err = ioctl(fd->fd_sys, LL_IOC_LOV_GETSTRIPE, (void *) lum);
        if (err == 0) {
            stripin_info[0] = lum->lmm_stripe_size;
            stripin_info[1] = lum->lmm_stripe_count;
            stripin_info[2] = lum->lmm_stripe_offset;
            /* TODO: figure out the correct way to find the number of unique
             * OSTs. This is relevant only when Lustre over-striping is used.
             * For now, set it to same as fd->hints->striping_factor.
             stripin_info[3] = fd->hints->striping_factor;
             */
            stripin_info[3] = fd->hints->striping_factor;
        }
        else
            err =  ncmpii_error_posix2nc("ioctl");
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
    fd->hints->striping_unit   = stripin_info[0];
    fd->hints->striping_factor = stripin_info[1];
    fd->hints->start_iodevice  = stripin_info[2];
    fd->hints->num_osts        = stripin_info[3];

    return err;
}

/*----< construct_aggr_list() >----------------------------------------------*/
/* Allocate and construct fd->hints->ranklist[].
 * Overwrite fd->hints->cb_nodes and set hint cb_nodes.
 * Set fd->is_agg, whether this rank is an I/O aggregator
 *     fd->hints->cb_nodes
 *     fd->hints->num_osts
 */
static int construct_aggr_list(ADIO_File fd, int root)
{
    int i, j, k, rank, nprocs, num_aggr, my_procname_len, num_nodes;
    int msg[3], striping_factor;
    int *all_procname_lens = NULL;
    int *nprocs_per_node, **ranks_per_node;
    char value[MPI_MAX_INFO_VAL + 1], my_procname[MPI_MAX_PROCESSOR_NAME];
    char **all_procnames = NULL;

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
    my_procname_len = (int)strlen(my_procname);
#else
    MPI_Get_processor_name(my_procname, &my_procname_len);
#endif
    my_procname_len++; /* to include terminate null character */

    if (rank == root) {
        /* root collects all procnames */
        all_procnames = (char **) ADIOI_Malloc(nprocs * sizeof(char *));
        if (all_procnames == NULL)
            return NC_ENOMEM;

        all_procname_lens = (int *) ADIOI_Malloc(nprocs * sizeof(int));
        if (all_procname_lens == NULL) {
            ADIOI_Free(all_procnames);
            return NC_ENOMEM;
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
            return NC_ENOMEM;
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

            if (fd->hints->cb_nodes == 0 || fd->access_mode & MPI_MODE_RDONLY) {
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
        if (fd->hints->ranklist == NULL)
            return NC_ENOMEM;

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
        msg[1] = fd->hints->num_osts;
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
    MPI_Info_set(fd->info, "cb_nodes", value);

    fd->hints->num_osts = msg[1];
    sprintf(value, "%d", fd->hints->num_osts);
    MPI_Info_set(fd->info, "lustre_num_osts", value);

    /* save the number of unique compute nodes in communicator fd->comm */
    fd->num_nodes = msg[2];

    if (rank != root) {
        /* ranklist[] contains the MPI ranks of I/O aggregators */
        fd->hints->ranklist = (int *) ADIOI_Malloc(num_aggr * sizeof(int));
        if (fd->hints->ranklist == NULL)
            return NC_ENOMEM;
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

/*----< PNC_File_open() >----------------------------------------------------*/
int PNC_File_open(MPI_Comm    comm,
                  const char *filename,
                  int         amode,
                  MPI_Info    info,
                  ADIO_File  *fh)
{
    /* Before reaching to this subroutine, PNC_Check_Lustre() should have been
     * called to verify filename is on Lustre.
     */
    char value[MPI_MAX_INFO_VAL + 1];
    int i, err, min_err;
    ADIO_File fd = (ADIO_FileD*) NCI_Calloc(1,sizeof(ADIO_FileD));

    *fh = fd;

    fd->comm        = comm;
    fd->filename    = ADIOI_Strdup(filename);
    fd->atomicity   = 0;
    fd->etype       = MPI_BYTE;
    fd->etype_size  = 1;
    fd->filetype    = MPI_BYTE;
    fd->ftype_size  = 1;
    fd->file_system = ADIO_LUSTRE;
    fd->is_open     = 0;
    fd->access_mode = amode;
    fd->orig_access_mode = amode;

    /* create and initialize info object */
    fd->hints = (ADIOI_Hints*) NCI_Calloc(1, sizeof(ADIOI_Hints));
    if (info == MPI_INFO_NULL)
        MPI_Info_create(&fd->info);
    else
        MPI_Info_dup(info, &fd->info);

    err = PNC_File_SetInfo(fd, fd->info);
    if (err != NC_NOERR)
        return err;

    /* collective buffer */
    fd->io_buf = ADIOI_Calloc(1, fd->hints->cb_buffer_size);
    if (fd->io_buf == NULL) {
        err = NC_ENOMEM;
        goto err_out;
    }

    if (PNC_Check_Lustre(filename))
        MPI_Info_set(fd->info, "romio_filesystem_type", "LUSTRE:");

    /* For Lustre, determining the I/O aggregators and constructing ranklist
     * requires file stripe count, which can only be obtained after file is
     * opened.
     */
    if (amode & MPI_MODE_CREATE)
        err = lustre_file_create(fd, amode);
    else
        err = lustre_file_open(fd);
    if (err != NC_NOERR) goto err_out;

    /* TODO: when no_indep_rw hint is enabled, only aggregators open the file */
    fd->is_open = 1;

    /* Use file striping count and number of unique OSTs to construct
     * the I/O aggregator rank list, fd->hints->ranklist[]. Also, set
     * hint cb_nodes, fd->is_agg, and fd->my_cb_nodes_index.
     */
    err = construct_aggr_list(fd, 0);
    if (err != NC_NOERR) goto err_out;

    /* set file striping hints */
    snprintf(value, sizeof(value), "%d", fd->hints->striping_unit);
    MPI_Info_set(fd->info, "striping_unit", value);

    snprintf(value, sizeof(value), "%d", fd->hints->striping_factor);
    MPI_Info_set(fd->info, "striping_factor", value);

    snprintf(value, sizeof(value), "%d", fd->hints->start_iodevice);
    MPI_Info_set(fd->info, "start_iodevice", value);

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
    MPI_Info_set(fd->info, "aggr_list", value);

err_out:
    MPI_Allreduce(&err, &min_err, 1, MPI_INT, MPI_MIN, comm);
    /* All NC errors are < 0 */
    if (min_err < 0) {
        if (err == 0) /* close file if opened successfully */
            close(fd->fd_sys);
        ADIOI_Free(fd->filename);
        if (fd->hints->ranklist != NULL)
            ADIOI_Free(fd->hints->ranklist);
        ADIOI_Free(fd->hints);
        if (fd->info != MPI_INFO_NULL)
            MPI_Info_free(&(fd->info));
        ADIOI_Free(fd->io_buf);
    }

/*
int rank; MPI_Comm_rank(comm, &rank);
if (rank == 0) {
    int  i, nkeys;
    MPI_Info_get_nkeys(fd->info, &nkeys);
    printf("%s line %d: MPI File Info: nkeys = %d\n",__func__,__LINE__,nkeys);
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(fd->info, i, key);
        MPI_Info_get_valuelen(fd->info, key, &valuelen, &flag);
        MPI_Info_get(fd->info, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}
*/

    return err;
}

/*----< PNC_File_close() >---------------------------------------------------*/
int PNC_File_close(ADIO_File *fh)
{
    int err = NC_NOERR;

    err = close((*fh)->fd_sys);
    if (err != 0)
        err = ncmpii_error_posix2nc("close");

    ADIOI_Free((*fh)->filename);
    if ((*fh)->hints->ranklist != NULL)
        ADIOI_Free((*fh)->hints->ranklist);
    if ((*fh)->hints->cb_config_list != NULL)
        ADIOI_Free((*fh)->hints->cb_config_list);
    if ((*fh)->hints != NULL)
        NCI_Free((*fh)->hints);
    if ((*fh)->info != MPI_INFO_NULL)
        MPI_Info_free(&((*fh)->info));
    if ((*fh)->io_buf != NULL)
        ADIOI_Free((*fh)->io_buf);
    if ((*fh)->filetype != MPI_BYTE)
        PNC_Type_dispose(&(*fh)->filetype);
    NCI_Free(*fh);

    return err;
}
