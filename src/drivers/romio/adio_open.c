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
#include <fcntl.h>      /* open(), O_CREAT */
#include <sys/types.h>  /* open() */
#include <libgen.h>     /* basename() */

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
#include <sys/stat.h> /* open(), fstat(), lstat(), stat() */
#endif

#ifdef HAVE_LUSTRE
#include <lustre/lustreapi.h>

/* what is the basis for this define?
 * what happens if there are more than 1k UUIDs? */
#define MAX_LOV_UUID_COUNT      1000
#endif

#ifdef MIMIC_LUSTRE
#define xstr(s) str(s)
#define str(s) #s
#define STRIPE_SIZE 64
#define STRIPE_COUNT 4
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

/* Check if file system type from file name, using a system-dependent function
 * call.
 */
int ADIO_FileSysType(const char *filename)
{

    int err, retry_cnt;
    int64_t file_id=UNKNOWN_SUPER_MAGIC;

    char *colon = strchr(filename, ':');
    if (colon != NULL) { /* there is a prefix end with : */
        if (!strncmp(filename, "lustre", 6))
            return ADIO_LUSTRE;
        else if (!strncmp(filename, "ufs", 3))
            return ADIO_UFS;
        else
            return 0;
    }
#ifdef MIMIC_LUSTRE
    return ADIO_LUSTRE;
#endif

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

    if (file_id == LL_SUPER_MAGIC)
        return ADIO_LUSTRE;
    else
        return ADIO_UFS; /* UFS support if we don't know what else to use */
}

#ifdef HAVE_LUSTRE
#define ERR(fn) { \
    printf("Error at %s (%d) calling %s\n", __func__, __LINE__, fn); \
    return -1; \
}

static int compare(const void *a, const void *b)
{
     if (*(uint64_t*)a > *(uint64_t*)b) return (1);
     if (*(uint64_t*)a < *(uint64_t*)b) return (-1);
     return (0);
}

#include <sys/types.h>  /* open(), fgetxattr() */
#include <sys/xattr.h>  /* fgetxattr() */

/*----< num_uniq_osts() >----------------------------------------------------*/
static int num_uniq_osts(int fd_sys)
{
    int err;
    void *xattr_val;
    ssize_t xattr_size = XATTR_SIZE_MAX;
    struct llapi_layout *layout;
    uint64_t i, stripe_count, stripe_size, *osts, numOSTs;

    if ((xattr_val = ADIOI_Calloc(1, xattr_size)) == NULL)
        ERR("ADIOI_Calloc")

    xattr_size = fgetxattr(fd_sys, "lustre.lov", xattr_val, xattr_size);
    if (xattr_size == -1) {
        ADIOI_Free(xattr_val);
        ERR("fgetxattr")
    }

    layout = llapi_layout_get_by_xattr(xattr_val, xattr_size, LLAPI_LAYOUT_GET_CHECK);
    ADIOI_Free(xattr_val);
    if (layout == NULL) ERR("llapi_layout_get_by_xattr")

    /* obtain file striping count */
    err = llapi_layout_stripe_count_get(layout, &stripe_count);
    if (err != 0) ERR("llapi_layout_stripe_count_get")

    /* obtain file striping unit size */
    err = llapi_layout_stripe_size_get(layout, &stripe_size);
    if (err != 0) ERR("llapi_layout_stripe_size_get")

    /* obtain all OST IDs */
    osts = (uint64_t*) ADIOI_Malloc(sizeof(uint64_t) * stripe_count);
    for (i=0; i<stripe_count; i++) {
        uint64_t tmp_ost;
        if (llapi_layout_ost_index_get(layout, i, &tmp_ost) == -1)
            ERR("llapi_layout_ost_index_get")
        else
            osts[i] = tmp_ost;
    }

    /* if (debug) {
        printf("Input file=%s\n",path);
        printf("stripe_count=%llu\n",stripe_count);
        printf("stripe_size=%llu\n",stripe_size);
        for (i=0; i<stripe_count; i++)
            printf("\tosts[%d]=%d\n",i,osts[i]);
    } */

    /* count the number of unique OST IDs. When Lustre overstriping is
     * used, the unique OSTs may be less than stripe_count.
     */
    qsort(osts, stripe_count, sizeof(uint64_t), compare);
    numOSTs = 0;
    for (i=1; i<stripe_count; i++) {
        if (osts[i] > osts[numOSTs]) osts[++numOSTs] = osts[i];
    }
    numOSTs++;

    ADIOI_Free(osts);
    llapi_layout_free(layout);

    return numOSTs;
}
#endif

#if 1

#ifdef HAVE_LUSTRE
static
int sort_ost_ids(struct llapi_layout *layout,
                 uint64_t             stripe_count,
                 uint64_t            *osts)
{
    uint64_t i, numOSTs;

    for (i=0; i<stripe_count; i++) {
        if (llapi_layout_ost_index_get(layout, i, &osts[i]) != 0)
            printf("Error at line %d: llapi_layout_ost_index_get(%lu) (%s)\n",
                   __LINE__,i,strerror(errno));
    }
    /* count the number of unique OST IDs. When Lustre overstriping is
     * used, the unique OSTs may be less than stripe_count.
     */
    qsort(osts, stripe_count, sizeof(uint64_t), compare);
    numOSTs = 0;
    for (i=1; i<stripe_count; i++) {
        if (osts[i] > osts[numOSTs]) osts[++numOSTs] = osts[i];
    }

    return (numOSTs + 1);
}

#define ERR0(fn) { \
    printf("Error at %s (%d) calling %s\n", __func__, __LINE__, fn); \
    return 0; \
}

/*----< get_striping() >-----------------------------------------------------*/
static uint64_t get_striping(int       fd,
                             uint64_t *pattern,
                             uint64_t *stripe_count,
                             uint64_t *stripe_size,
                             uint64_t *start_iodevice)
{
    int err;
    struct llapi_layout *layout;
    uint64_t *osts, numOSTs;

    layout = llapi_layout_get_by_fd(fd, LLAPI_LAYOUT_GET_COPY);
    if (layout == NULL) ERR0("llapi_layout_get_by_fd")

    err = llapi_layout_pattern_get(layout, pattern);
    if (err != 0) ERR0("llapi_layout_pattern_get")
/*
    if (*pattern == LLAPI_LAYOUT_DEFAULT)
        printf("line %d: pattern is LLAPI_LAYOUT_DEFAULT\n",__LINE__);
    else if (*pattern == LLAPI_LAYOUT_RAID0)
        printf("line %d: pattern is LLAPI_LAYOUT_RAID0\n",__LINE__);
    else if (*pattern == LLAPI_LAYOUT_MDT)
        printf("line %d: pattern is LLAPI_LAYOUT_MDT\n",__LINE__);
    else if (*pattern == LLAPI_LAYOUT_OVERSTRIPING)
        printf("line %d: pattern is LLAPI_LAYOUT_OVERSTRIPING\n",__LINE__);
    else if (*pattern == LLAPI_LAYOUT_SPECIFIC)
        printf("line %d: pattern is LLAPI_LAYOUT_SPECIFIC\n",__LINE__);
    else
        printf("line %d: pattern is unknown\n",__LINE__);
*/

    /* obtain file striping count */
    err = llapi_layout_stripe_count_get(layout, stripe_count);
    if (err != 0) ERR0("llapi_layout_stripe_count_get")

/*
    if (*stripe_count == LLAPI_LAYOUT_DEFAULT)
        printf("%s at %d: stripe count is NOT SET %lu\n",__func__,__LINE__,*stripe_count);
    else
        printf("%s at %d: stripe count = %lu\n",__func__,__LINE__,*stripe_count);
*/

    /* obtain file striping unit size */
    err = llapi_layout_stripe_size_get(layout, stripe_size);
    if (err != 0) ERR0("llapi_layout_stripe_size_get")

/*
    if (*stripe_size == LLAPI_LAYOUT_DEFAULT)
        printf("%s at %d: stripe size is NOT SET %lu\n",__func__,__LINE__,*stripe_size);
    else
        printf("%s at %d: stripe size is %lu\n",__func__,__LINE__,*stripe_size);
*/

        /*
         * /usr/include/linux/lustre/lustre_user.h
         * The stripe size fields are shared for the extension size storage,
         * however the extension size is stored in KB, not bytes.
         *     #define SEL_UNIT_SIZE 1024llu
         * Therefore, the default stripe_size is (SEL_UNIT_SIZE * 1024)
         */

    /* obtain all OST IDs */
    osts = (uint64_t*) ADIOI_Malloc(sizeof(uint64_t) * (*stripe_count));
    if (llapi_layout_ost_index_get(layout, 0, &osts[0]) != 0) {

        /* check if is a folder */
        struct stat path_stat;
        fstat(fd, &path_stat);
        if (S_ISREG(path_stat.st_mode)) /* not a regular file */
            printf("%s at %d: is a regular file\n",__func__,__LINE__);
        else if (S_ISDIR(path_stat.st_mode))
            printf("%s at %d: is a folder\n",__func__,__LINE__);
        else
            ERR0("fstat")

        *start_iodevice = LLAPI_LAYOUT_DEFAULT;
        return *stripe_count;
    }
    // printf("%s at %d:  osts[0]=%lu\n",__func__,__LINE__, osts[0]);
    *start_iodevice = osts[0];

    numOSTs = sort_ost_ids(layout, *stripe_count, osts);
    ADIOI_Assert(numOSTs <= *stripe_count);

    ADIOI_Free(osts);
    llapi_layout_free(layout);

    return numOSTs;
}

/*----< set_striping() >-----------------------------------------------------*/
static int set_striping(const char *path,
                        uint64_t    pattern,
                        uint64_t    numOSTs,
                        uint64_t    stripe_count,
                        uint64_t    stripe_size,
                        uint64_t    start_iodevice)
{
    int fd, err;

    struct llapi_layout *layout = llapi_layout_alloc();
    if (layout == NULL) ERR("llapi_layout_alloc")

    // printf("%s at %d: stripe_count = %lu\n",__func__,__LINE__,stripe_count);
    err = llapi_layout_stripe_count_set(layout, stripe_count);
    if (err != 0) ERR("llapi_layout_stripe_count_set")

    // printf("%s at %d: stripe_size = %lu\n",__func__,__LINE__,stripe_size);
    err = llapi_layout_stripe_size_set(layout, stripe_size);
    if (err != 0) ERR("llapi_layout_stripe_size_set")

    if (pattern == LLAPI_LAYOUT_OVERSTRIPING) {
        uint64_t i, ost_id;
        if (start_iodevice == LLAPI_LAYOUT_DEFAULT)
            start_iodevice = 0;
        for (i=0; i<stripe_count; i++) {
            ost_id = start_iodevice + (i % numOSTs);
            err = llapi_layout_ost_index_set(layout, i, ost_id);
            if (err != 0) ERR("llapi_layout_ost_index_set")
        }
    }
    else {
        err = llapi_layout_ost_index_set(layout, 0, start_iodevice);
        if (err != 0) ERR("llapi_layout_ost_index_set")
    }

    err = llapi_layout_pattern_set(layout, pattern);
    if (err != 0) ERR("llapi_layout_pattern_set")

    fd = llapi_layout_file_create(path, O_CREAT|O_RDWR, 0660, layout);
    if (fd < 0) ERR("llapi_layout_file_create")

#if 0
{
    void *lmdbuf = ADIOI_Malloc(XATTR_SIZE_MAX);

    err = ioctl(fd, LL_IOC_LOV_GETSTRIPE_NEW, lmdbuf);
    if (err < 0) ERR("close")

    uint64_t count = ((struct lov_user_md *)lmdbuf)->lmm_stripe_count;
    // printf("%s at %d: count = %lu\n",__func__,__LINE__,count);

    ADIOI_Free(lmdbuf);
}
#endif

    llapi_layout_free(layout);

    return fd;
}
#endif

/*----< file_create() >------------------------------------------------------*/
/*   1. root creates the file
 *   2. root sets and obtains striping info
 *   3. root broadcasts striping info
 *   4. non-root processes receive striping info from root
 *   5. non-root processes opens the fie
 */
static int
file_create(ADIO_File fd,
            int       access_mode)
{
    int err=NC_NOERR, rank, amode, perm, old_mask;
    int stripin_info[4] = {-1, -1, -1, -1};

    MPI_Comm_rank(fd->comm, &rank);

// int world_rank; MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); printf("%s at %d: world_rank=%d local rank=%d ---- calling POSIX open()\n",__func__,__LINE__,world_rank,rank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
static int wkl=0; if (wkl == 0 && rank == 0) { printf("\nxxxx %s at %d: %s ---- %s\n",__func__,__LINE__,(fd->file_system == ADIO_LUSTRE)?"ADIO_LUSTRE":"ADIO_UFS",fd->filename); wkl++; fflush(stdout);}
#endif

    amode = O_CREAT;
    if (access_mode & MPI_MODE_RDWR)  amode |= O_RDWR;

    old_mask = umask(022);
    umask(old_mask);
    perm = old_mask ^ 0666;

    /* root process creates the file first, followed by all processes open the
     * file.
     */
    if (rank > 0) goto err_out;

    if (fd->file_system == ADIO_LUSTRE) {
        /* For Lustre, we need to obtain file striping info (striping_factor,
         * striping_unit, and num_osts) in order to select the I/O aggregators
         * in fd->hints->ranklist, no matter its is open or create mode.
         */

#ifdef HAVE_LUSTRE
        int set_user_layout = 0, overstriping_ratio;
        int str_factor, str_unit, start_iodev;

        /* In a call to ADIO_File_SetInfo() earlier, hints have been validated
         * to be consistent among all processes.
         */

        str_unit = fd->hints->striping_unit;
        str_factor = fd->hints->striping_factor;
        start_iodev = fd->hints->start_iodevice;
        overstriping_ratio = fd->hints->fs_hints.lustre.overstriping_ratio;

        /* when no file striping hint is set, their values are:
         * fd->hints->striping_unit = 0;
         * fd->hints->striping_factor = 0;
         * fd->hints->start_iodevice = -1;
         * fd->hints->fs_hints.lustre.overstriping_ratio = 1;
         */

        /* if user has set the file striping hints */
        if ((str_factor > 0) || (str_unit > 0) || (start_iodev >= 0) ||
            (overstriping_ratio > 1))
            set_user_layout = 1;

/* query the total number of available OSTs in the pool.
        char **members, *buffer, *pool="scratch.original";
        int list_size = 512;

        members = (char**)ADIOI_Malloc(sizeof(char*) * list_size);
        buffer = (char*)ADIOI_Calloc(list_size * 64, sizeof(char));

        int num_OSTs = llapi_get_poolmembers(pool, members, list_size, buffer, list_size * 64);
        printf("Lustre pool %s has %d OSTs\n",pool, num_OSTs);
        int i;
        for (i=0; i<20; i++)
            printf("member[%3d] %s\n",i,members[i]);
        printf("member[%3d] %s\n",num_OSTs-1,members[num_OSTs-1]);
        ADIOI_Free(buffer);
        ADIOI_Free(members);
*/
        uint64_t numOSTs;
        uint64_t pattern;
        uint64_t stripe_count;
        uint64_t stripe_size;
        uint64_t start_iodevice;

        if (set_user_layout) {
            /* user has set striping hints, grant wish */

            if ((str_factor == 0) || (str_unit == 0) || (start_iodev == -1) ||
                (overstriping_ratio < 0)) {
                /* For those striping hints are not set by the user, inherit
                 * striping settings from the directory containing the file.
                 */
                int dd;
                char *dirc, *dname;
                dirc = ADIOI_Strdup(fd->filename);
                dname = dirname(dirc);

                dd = open(dname, O_RDONLY, 0600);

                numOSTs = get_striping(dd, &pattern,
                                           &stripe_count,
                                           &stripe_size,
                                           &start_iodevice);
                close(dd);
                ADIOI_Free(dirc);
                if (numOSTs == 0) {
                    stripin_info[3] = numOSTs;
                    goto err_out;
                }
            }

            numOSTs = (str_factor == -1) ? numOSTs : str_factor;
            if (overstriping_ratio > 1) {
                pattern = LLAPI_LAYOUT_OVERSTRIPING;
                stripe_count = numOSTs * overstriping_ratio;
            }
            else {
                pattern = LLAPI_LAYOUT_RAID0;
                stripe_count = numOSTs;
            }
            stripe_size = (str_unit == 0) ? stripe_size : str_unit;
            start_iodevice = (start_iodev == -1) ? start_iodevice : start_iodev;
        }
        else {
            /* user has not set file striping hints, use the followings:
             * set pattern to LLAPI_LAYOUT_OVERSTRIPING
             * set overstriping_ratio to 4
             * set stripe_count to fd->num_nodes * overstriping_ratio
             * set stripe_size to 1 MiB
             */
            numOSTs = fd->num_nodes;
            if (numOSTs > 256) numOSTs = 256; /* TODO: check max available OSTs */
            overstriping_ratio = 4;
            pattern = LLAPI_LAYOUT_OVERSTRIPING;
            stripe_count = numOSTs * overstriping_ratio;
            stripe_size = 1048576;
            start_iodevice = LLAPI_LAYOUT_DEFAULT;
        }

        fd->fd_sys = set_striping(fd->filename, pattern,
                                                numOSTs,
                                                stripe_count,
                                                stripe_size,
                                                start_iodevice);
        if (fd->fd_sys < 0) {
            fprintf(stderr,"%s line %d: rank %d fails to create and set striping file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            err = ncmpii_error_posix2nc("Lustre set striping");
            goto err_out;
        }

        /* get Lustre file stripning */
        numOSTs = get_striping(fd->fd_sys, &pattern,
                                           &stripe_count,
                                           &stripe_size,
                                           &start_iodevice);

        stripin_info[0] = stripe_size;
        stripin_info[1] = stripe_count;
        stripin_info[2] = start_iodevice;
        stripin_info[3] = numOSTs;

#elif defined(MIMIC_LUSTRE)
        fd->fd_sys = open(fd->filename, amode, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d fails to create file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            err = ncmpii_error_posix2nc("open");
            goto err_out;
        }
        stripin_info[0] = STRIPE_SIZE;
        stripin_info[1] = STRIPE_COUNT;
        stripin_info[2] = 0;
        stripin_info[3] = STRIPE_COUNT;
#endif
    }
    else {
        fd->fd_sys = open(fd->filename, amode, perm);
        if (fd->fd_sys == -1) {
            fprintf(stderr,"%s line %d: rank %d fails to create file %s (%s)\n",
                    __func__,__LINE__, rank, fd->filename, strerror(errno));
            err = ncmpii_error_posix2nc("open");
            goto err_out;
        }
    }

err_out:
    MPI_Bcast(stripin_info, 4, MPI_INT, 0, fd->comm);
    if (fd->file_system == ADIO_LUSTRE &&
        (stripin_info[0] == -1 || stripin_info[3] == 0)) {
        fprintf(stderr, "%s line %d: failed to create Lustre file %s\n",
                __func__, __LINE__, fd->filename);
        return err;
    }

    fd->hints->striping_unit   = stripin_info[0];
    fd->hints->striping_factor = stripin_info[1];
    fd->hints->start_iodevice  = stripin_info[2];
    if (fd->file_system == ADIO_LUSTRE) {
        fd->hints->fs_hints.lustre.num_osts = stripin_info[3];
        fd->hints->fs_hints.lustre.overstriping_ratio = stripin_info[1] / stripin_info[3];
    }

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
#else

/*----< file_create() >------------------------------------------------------*/
/*   1. root creates the file
 *   2. root sets and obtains striping info
 *   3. root broadcasts striping info
 *   4. non-root processes receive striping info from root
 *   5. non-root processes opens the fie
 */
static int
file_create(ADIO_File fd,
            int       access_mode)
{
    int err=NC_NOERR, rank, amode, perm, old_mask;
    int stripin_info[4] = {1048576, -1, -1, -1};

    MPI_Comm_rank(fd->comm, &rank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
static int wkl=0; if (wkl == 0 && rank == 0) { printf("\nxxxx %s at %d: %s ---- %s\n",__func__,__LINE__,(fd->file_system == ADIO_LUSTRE)?"ADIO_LUSTRE":"ADIO_UFS",fd->filename); wkl++; fflush(stdout);}
#endif

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
#ifdef HAVE_LUSTRE
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
            err = ncmpii_error_posix2nc("open");
            goto err_out;
        }

#ifdef HAVE_LUSTRE
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
            /* Find the number of unique OSTs. This is relevant only when
             * Lustre over-striping is used.
             */
            stripin_info[3] = num_uniq_osts(fd->fd_sys);
        }
        else {
            fprintf(stderr,"%s line %d: rank %d fails to get stripe info from file %s (%s)\n",
                    __func__,__LINE__,rank, fd->filename, strerror(errno));
            err = ncmpii_error_posix2nc("ioctl");
        }
        ADIOI_Free(lum);
#elif defined(MIMIC_LUSTRE)
        stripin_info[0] = STRIPE_SIZE;
        stripin_info[1] = STRIPE_COUNT;
        stripin_info[2] = 0;
        stripin_info[3] = STRIPE_COUNT;
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
    fd->hints->fs_hints.lustre.num_osts = stripin_info[3];

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
#endif

/*----< file_open() >--------------------------------------------------------*/
/*   1. all processes open the file.
 *   2. root obtains striping info and broadcasts to all others
 */
static int
file_open(ADIO_File fd)
{
    int err=NC_NOERR, rank, perm, old_mask;
    int stripin_info[4] = {1048576, -1, -1, -1};

    MPI_Comm_rank(fd->comm, &rank);

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
static int wkl=0; if (wkl == 0 && rank == 0) { printf("\nxxxx %s at %d: %s ---- %s\n",__func__,__LINE__,(fd->file_system == ADIO_LUSTRE)?"ADIO_LUSTRE":"ADIO_UFS",fd->filename); wkl++; fflush(stdout);}
#endif

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
#ifdef HAVE_LUSTRE
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
            stripin_info[3] = lum->lmm_stripe_count;
        }
        else
            err = ncmpii_error_posix2nc("ioctl");
        ADIOI_Free(lum);
#elif defined(MIMIC_LUSTRE)
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
    fd->hints->fs_hints.lustre.num_osts = stripin_info[3];
    fd->hints->fs_hints.lustre.overstriping_ratio = stripin_info[1] / stripin_info[3];

    return err;
}

/*----< ADIO_GEN_set_cb_node_list() >----------------------------------------*/
/* Construct the list of I/O aggregators. It sets the followings.
 *   fd->hints->ranklist[].
 *   fd->hints->cb_nodes and set file info for hint cb_nodes.
 *   fd->is_agg: indicating whether this rank is an I/O aggregator
 *   fd->my_cb_nodes_index: index into fd->hints->ranklist[]. -1 if N/A
 */
static
int ADIO_GEN_set_cb_node_list(ADIO_File fd)
{
    int i, j, k, nprocs, rank, *nprocs_per_node, **ranks_per_node;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &rank);

    if (fd->hints->cb_nodes == 0)
        /* If hint cb_nodes is not set by user, select one rank per node to be
         * an I/O aggregator
         */
        fd->hints->cb_nodes = fd->num_nodes;
    else if (fd->hints->cb_nodes > 0)
        /* cb_nodes must be <= nprocs */
        fd->hints->cb_nodes = nprocs;

    fd->hints->ranklist = (int *) ADIOI_Malloc(sizeof(int) * fd->hints->cb_nodes);
    if (fd->hints->ranklist == NULL)
        return NC_ENOMEM;

    /* number of MPI processes running on each node */
    nprocs_per_node = (int *) ADIOI_Calloc(fd->num_nodes, sizeof(int));

    for (i=0; i<nprocs; i++) nprocs_per_node[fd->node_ids[i]]++;

    /* construct rank IDs of MPI processes running on each node */
    ranks_per_node = (int **) ADIOI_Malloc(sizeof(int*) * fd->num_nodes);
    ranks_per_node[0] = (int *) ADIOI_Malloc(sizeof(int) * nprocs);
    for (i=1; i<fd->num_nodes; i++)
        ranks_per_node[i] = ranks_per_node[i - 1] + nprocs_per_node[i - 1];

    for (i=0; i<fd->num_nodes; i++) nprocs_per_node[i] = 0;

    /* Populate ranks_per_node[], list of MPI ranks running on each node.
     * Populate nprocs_per_node[], number of MPI processes on each node.
     */
    for (i=0; i<nprocs; i++) {
        k = fd->node_ids[i];
        ranks_per_node[k][nprocs_per_node[k]] = i;
        nprocs_per_node[k]++;
    }

    /* select process ranks from nodes in a round-robin fashion to be I/O
     * aggregators
     */
    k = j = 0;
    for (i=0; i<fd->hints->cb_nodes; i++) {
        if (j >= nprocs_per_node[k]) { /* if run out of ranks in this node k */
            k++;
            if (k == fd->num_nodes) { /* round-robin to first node */
                k = 0;
                j++;
            }
        }
        /* select jth rank of node k as an I/O aggregator */
        fd->hints->ranklist[i] = ranks_per_node[k++][j];
        if (rank == fd->hints->ranklist[i]) {
            fd->is_agg = 1;
            fd->my_cb_nodes_index = i;
        }
        if (k == fd->num_nodes) { /* round-robin to first node */
            k = 0;
            j++;
        }
    }
    ADIOI_Free(ranks_per_node[0]);
    ADIOI_Free(ranks_per_node);
    ADIOI_Free(nprocs_per_node);

    return 0;
}

/*----< ADIO_Lustre_set_cb_node_list() >-------------------------------------*/
/* Construct the list of I/O aggregators. It sets the followings.
 *   fd->hints->ranklist[].
 *   fd->hints->cb_nodes and set file info for hint cb_nodes.
 *   fd->is_agg: indicating whether this rank is an I/O aggregator
 *   fd->my_cb_nodes_index: index into fd->hints->ranklist[]. -1 if N/A
 */
static
int ADIO_Lustre_set_cb_node_list(ADIO_File fd)
{
    int i, j, k, rank, nprocs, num_aggr, striping_factor;
    int *nprocs_per_node, **ranks_per_node;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &rank);

    /* number of MPI processes running on each node */
    nprocs_per_node = (int *) ADIOI_Calloc(fd->num_nodes, sizeof(int));

    for (i=0; i<nprocs; i++) nprocs_per_node[fd->node_ids[i]]++;

    /* construct rank IDs of MPI processes running on each node */
    ranks_per_node = (int **) ADIOI_Malloc(sizeof(int*) * fd->num_nodes);
    ranks_per_node[0] = (int *) ADIOI_Malloc(sizeof(int) * nprocs);
    for (i=1; i<fd->num_nodes; i++)
        ranks_per_node[i] = ranks_per_node[i - 1] + nprocs_per_node[i - 1];

    for (i=0; i<fd->num_nodes; i++) nprocs_per_node[i] = 0;

    /* Populate ranks_per_node[], list of MPI ranks running on each node.
     * Populate nprocs_per_node[], number of MPI processes on each node.
     */
    for (i=0; i<nprocs; i++) {
        k = fd->node_ids[i];
        ranks_per_node[k][nprocs_per_node[k]] = i;
        nprocs_per_node[k]++;
    }

    /* All processes run the same codes below to calculate num_aggr, number of
     * aggregators, so we can save a call to MPI_Bcast().
     * Given the number of nodes, fd->num_nodes, and processes per node,
     * nprocs_per_node, we can now set num_aggr, the number of I/O aggregators.
     * At this moment, root should have obtained the file striping settings.
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
        /* Select striping_factor processes to be I/O aggregators. Note this
         * also applies to collective reads to allow more/less aggregators. In
         * most cases, more aggregators yields better read performance.
         */
        if (fd->hints->cb_nodes <= striping_factor) {
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
            for (i = 0; i < fd->num_nodes; i++)
                max_nprocs_node = MAX(max_nprocs_node, nprocs_per_node[i]);
            int max_naggr_node = striping_factor / fd->num_nodes;
            if (striping_factor % fd->num_nodes) max_naggr_node++;
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

    /* Next step is to determine the MPI rank IDs of I/O aggregators and add
     * them into ranklist[].
     */
    fd->hints->ranklist = (int *) ADIOI_Malloc(num_aggr * sizeof(int));
    if (fd->hints->ranklist == NULL)
        return NC_ENOMEM;

    int block_assignment=0;
#ifdef TRY_AGGR_BLOCK_ASSIGNMENT
    {
        char *env_str;
        if ((env_str = getenv("PNETCDF_USE_BLOCK_ASSIGN")) != NULL)
            block_assignment = (strcasecmp(env_str, "true") == 0) ? 1 : 0;
        if (rank == 0)
            printf("%s %d: PNETCDF_USE_BLOCK_ASSIGN = %d\n",
            __func__,__LINE__,block_assignment);
    }
#endif

    if (striping_factor <= fd->num_nodes) {
        /* When number of OSTs is less than number of compute nodes, first
         * select number of nodes equal to the number of OSTs by spread the
         * selection evenly across all compute nodes (i.e. with a stride
         * between every 2 consecutive nodes).
         * Selection of MPI ranks can be done in 2 ways.
         * 1. block assignment
         *    Select ranks from a node and then move on to the next node.
         * 2. cyclic assignment
         *    Select ranks round-robin across all selected nodes.
         * Note when selecting ranks within a node, the ranks are evenly spread
         * among all processes in the node.
         */
        if (block_assignment) {
            int n=0;
            int remain = num_aggr % striping_factor;
            int node_stride = fd->num_nodes / striping_factor;
            /* walk through each node and pick aggregators */
            for (j=0; j<fd->num_nodes; j+=node_stride) {
                /* Selecting node IDs with a stride. j is the node ID */
                int nranks_per_node = num_aggr / striping_factor;
                /* front nodes may have 1 more to pick */
                if (remain > 0 && j/node_stride < remain) nranks_per_node++;
                int rank_stride = nprocs_per_node[j] / nranks_per_node;
                for (k=0; k<nranks_per_node; k++) {
                    /* Selecting rank IDs within node j with a stride */
                    fd->hints->ranklist[n] = ranks_per_node[j][k*rank_stride];
                    if (++n == num_aggr) {
                        j = fd->num_nodes; /* break loop j */
                        break; /* loop k */
                    }
                }
            }
        }
        else {
            int avg = num_aggr / striping_factor;
            int stride = fd->num_nodes / striping_factor;
            if (num_aggr % striping_factor) avg++;
            for (i = 0; i < num_aggr; i++) {
                /* j is the selected node ID. This selection is round-robin
                 * across selected nodes.
                 */
                j = (i % striping_factor) * stride;
                k = (i / striping_factor) * (nprocs_per_node[j] / avg);
                ADIOI_Assert(k < nprocs_per_node[j]);
                fd->hints->ranklist[i] = ranks_per_node[j][k];
            }
        }
    }
    else { /* striping_factor > fd->num_nodes */
        /* When number of OSTs is more than number of compute nodes, I/O
         * aggregators are selected from all nodes. Within each node,
         * aggregators are spread evenly instead of the first few ranks.
         */
        int *naggr_per_node, *idx_per_node, avg;
        idx_per_node = (int*) ADIOI_Calloc(fd->num_nodes, sizeof(int));
        naggr_per_node = (int*) ADIOI_Malloc(fd->num_nodes * sizeof(int));
        for (i = 0; i < striping_factor % fd->num_nodes; i++)
            naggr_per_node[i] = striping_factor / fd->num_nodes + 1;
        for (; i < fd->num_nodes; i++)
            naggr_per_node[i] = striping_factor / fd->num_nodes;
        avg = num_aggr / striping_factor;
        if (avg > 0)
            for (i = 0; i < fd->num_nodes; i++)
                naggr_per_node[i] *= avg;
        for (i = 0; i < fd->num_nodes; i++)
            naggr_per_node[i] = MIN(naggr_per_node[i], nprocs_per_node[i]);
        /* naggr_per_node[] is the number of aggregators that can be
         * selected as I/O aggregators
         */

        if (block_assignment) {
            int n = 0;
            for (j=0; j<fd->num_nodes; j++) {
                /* j is the node ID */
                int rank_stride = nprocs_per_node[j] / naggr_per_node[j];
                /* try stride==1 seems no effect, rank_stride = 1; */
                for (k=0; k<naggr_per_node[j]; k++) {
                    fd->hints->ranklist[n] = ranks_per_node[j][k*rank_stride];
                    if (++n == num_aggr) {
                        j = fd->num_nodes; /* break loop j */
                        break; /* loop k */
                    }
                }
            }
        }
        else {
            for (i = 0; i < num_aggr; i++) {
                int stripe_i = i % striping_factor;
                j = stripe_i % fd->num_nodes; /* to select from node j */
                k = nprocs_per_node[j] / naggr_per_node[j];
                k *= idx_per_node[j];
                /* try stride==1 seems no effect, k = idx_per_node[j]; */
                idx_per_node[j]++;
                ADIOI_Assert(k < nprocs_per_node[j]);
                fd->hints->ranklist[i] = ranks_per_node[j][k];
            }
        }
        ADIOI_Free(naggr_per_node);
        ADIOI_Free(idx_per_node);
    }

    /* TODO: we can keep these two arrays in case for dynamic construction
     * of fd->hints->ranklist[], such as in group-cyclic file domain
     * assignment method, used in each collective write call.
     */
    ADIOI_Free(nprocs_per_node);
    ADIOI_Free(ranks_per_node[0]);
    ADIOI_Free(ranks_per_node);

    /* set file striping hints */
    fd->hints->cb_nodes = num_aggr;

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

/*----< ADIO_File_open() >---------------------------------------------------*/
int ADIO_File_open(MPI_Comm    comm,
                   const char *filename,
                   int         amode,
                   MPI_Info    info,
                   ADIO_File   fd)
{
    /* Before reaching to this subroutine, ADIO_FileSysType() should have been
     * called to verify filename is on Lustre.
     */
    char value[MPI_MAX_INFO_VAL + 1], int_str[16];
    int i, err, min_err;

#ifdef WKL_DEBUG
extern int first_ost_id;
first_ost_id = -1;
#endif

    fd->comm        = comm;
    fd->filename    = filename;  /* without file system type name prefix */
    fd->atomicity   = 0;
    fd->filetype    = MPI_BYTE;
    // fd->ftype_size  = 1;
    fd->is_open     = 0;
    fd->access_mode = amode;
    fd->flat_file   = NULL; /* flattend fileview in offset-length pairs */
    fd->io_buf      = NULL; /* collective buffer used by aggregators only */

    /* create and initialize info object */
    fd->hints = (ADIOI_Hints*) ADIOI_Calloc(1, sizeof(ADIOI_Hints));
    if (info == MPI_INFO_NULL)
        MPI_Info_create(&fd->info);
    else
        MPI_Info_dup(info, &fd->info);

    err = ADIO_File_SetInfo(fd, fd->info);
    if (err != NC_NOERR)
        return err;

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    { for (i=0; i<11; i++) fd->coll_write[i] = fd->coll_read[i] = 0; }
#endif

    ADIOI_Assert(fd->file_system != ADIO_FSTYPE_MPIIO);

    /* For Lustre, determining the I/O aggregators and constructing ranklist
     * requires file stripe count, which can only be obtained after file is
     * opened.
     */
    if (amode & MPI_MODE_CREATE)
        err = file_create(fd, amode);
    else
        err = file_open(fd);
    if (err != NC_NOERR) goto err_out;

    /* TODO: when hint no_indep_rw hint is set to true, only aggregators open
     * the file */
    fd->is_open = 1;

    /* construct cb_nodes rank list */
    if (fd->file_system == ADIO_LUSTRE) {
        ADIO_Lustre_set_cb_node_list(fd);
        MPI_Info_set(fd->info, "romio_filesystem_type", "LUSTRE:");
    }
    else if (fd->file_system == ADIO_UFS) {
        ADIO_GEN_set_cb_node_list(fd);
        MPI_Info_set(fd->info, "romio_filesystem_type", "UFS:");
    }

    /* set file striping hints */
    snprintf(value, sizeof(value), "%d", fd->hints->striping_unit);
    MPI_Info_set(fd->info, "striping_unit", value);

    snprintf(value, sizeof(value), "%d", fd->hints->striping_factor);
    MPI_Info_set(fd->info, "striping_factor", value);

    snprintf(value, sizeof(value), "%d", fd->hints->start_iodevice);
    MPI_Info_set(fd->info, "start_iodevice", value);

#ifdef HAVE_LUSTRE
    snprintf(value, sizeof(value), "%d", fd->hints->fs_hints.lustre.num_osts);
    MPI_Info_set(fd->info, "lustre_num_osts", value);

    snprintf(value, sizeof(value), "%d", fd->hints->fs_hints.lustre.overstriping_ratio);
    MPI_Info_set(fd->info, "lustre_overstriping_ratio", value);
#endif

    /* set file striping hints */
    snprintf(value, sizeof(value), "%d", fd->hints->cb_nodes);
    MPI_Info_set(fd->info, "cb_nodes", value);

    /* add hint "cb_node_list", list of aggregators' rank IDs */
    snprintf(value, sizeof(int_str), "%d", fd->hints->ranklist[0]);
    for (i=1; i<fd->hints->cb_nodes; i++) {
        snprintf(int_str, sizeof(int_str), " %d", fd->hints->ranklist[i]);
        if (strlen(value) + strlen(int_str) >= MPI_MAX_INFO_VAL-5) {
            strcat(value, " ...");
            break;
        }
        strcat(value, int_str);
    }
    MPI_Info_set(fd->info, "cb_node_list", value);

    /* collective buffer size must be at least file striping size */
    if (fd->hints->cb_buffer_size < fd->hints->striping_unit) {
        fd->hints->cb_buffer_size = fd->hints->striping_unit;
        snprintf(int_str, sizeof(int_str), " %d", fd->hints->cb_buffer_size);
        MPI_Info_set(fd->info, "cb_buffer_size", int_str);
    }

    /* collective buffer is used only by I/O aggregators only */
    if (fd->is_agg) {
        fd->io_buf = ADIOI_Calloc(1, fd->hints->cb_buffer_size);
        if (fd->io_buf == NULL)
            return NC_ENOMEM;
    }

err_out:
    MPI_Allreduce(&err, &min_err, 1, MPI_INT, MPI_MIN, comm);
    /* All NC errors are < 0 */
    if (min_err < 0) {
        if (err == 0) /* close file if opened successfully */
            close(fd->fd_sys);
        ADIOI_Free(fd->hints);
        if (fd->info != MPI_INFO_NULL)
            MPI_Info_free(&(fd->info));
        if (fd->io_buf != NULL)
            ADIOI_Free(fd->io_buf);
    }
    return err;
}

/*----< ADIO_File_close() >--------------------------------------------------*/
int ADIO_File_close(ADIO_File *fh)
{
    int err = NC_NOERR;

    err = close((*fh)->fd_sys);
    if (err != 0)
        err = ncmpii_error_posix2nc("close");

    if ((*fh)->hints->ranklist != NULL)
        ADIOI_Free((*fh)->hints->ranklist);
    if ((*fh)->hints->cb_config_list != NULL)
        ADIOI_Free((*fh)->hints->cb_config_list);
    if ((*fh)->hints != NULL)
        ADIOI_Free((*fh)->hints);
    if ((*fh)->info != MPI_INFO_NULL)
        MPI_Info_free(&((*fh)->info));
    if ((*fh)->io_buf != NULL)
        ADIOI_Free((*fh)->io_buf);
    ADIOI_Type_dispose(&(*fh)->filetype);

    if (ADIOI_Flattened_type_keyval != MPI_KEYVAL_INVALID) {
        MPI_Type_free_keyval(&ADIOI_Flattened_type_keyval);
        ADIOI_Flattened_type_keyval = MPI_KEYVAL_INVALID;
    }

#if defined(PNETCDF_PROFILING) && (PNETCDF_PROFILING == 1)
    int i, ntimers, rank;
    double tt[16], max_t[16];

    /* print two-phase I/O timing breakdown */
    MPI_Comm_rank((*fh)->comm, &rank);
    ntimers = 12;
    for (i=0; i<ntimers; i++) tt[i] = (*fh)->coll_write[i];
    MPI_Reduce(tt, max_t, ntimers, MPI_DOUBLE, MPI_MAX, 0, (*fh)->comm);
    if (rank == 0 && max_t[11] > 0)
        printf("%s: TWO-PHASE write init %5.2f pwrite %5.2f post %5.2f comm %5.2f collw %5.2f ntimes %d\n",
        __func__, max_t[1], max_t[2], max_t[4], max_t[3], max_t[0], (int)(max_t[11]));

    ntimers = 12;
    for (i=0; i<ntimers; i++) tt[i] = (*fh)->coll_write[i];
    MPI_Reduce(tt, max_t, ntimers, MPI_DOUBLE, MPI_MAX, 0, (*fh)->comm);
    if (rank == 0 && max_t[11] > 0)
        printf("%s: TWO-PHASE read  init %5.2f pread  %5.2f post %5.2f wait %5.2f collr %5.2f ntimes %d\n",
        __func__, max_t[1], max_t[2], max_t[4], max_t[3], max_t[0], (int)(max_t[11]));
#endif

    return err;
}
