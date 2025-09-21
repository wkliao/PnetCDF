/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <limits.h>
#include <stdarg.h> /* va_start(), va_end() */

#include <adio.h>

/* utility function to query a datatype for its combiner,
 * convenience wrapper around MPI_Type_get_envelope[_c] */

int PNCIO_Type_get_combiner(MPI_Datatype datatype, int *combiner)
{
    int ret;

    assert(datatype != MPI_DATATYPE_NULL);
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count ni, na, nc, nt;
    ret = MPI_Type_get_envelope_c(datatype, &ni, &na, &nc, &nt, combiner);
#else
    int ni, na, nt;
    ret = MPI_Type_get_envelope(datatype, &ni, &na, &nt, combiner);
#endif
    return ret;
}

/* utility function to determine whether a datatype is predefined:
 * a datatype is predefined if its combiner is MPI_COMBINER_NAMED
 * or MPI_COMBINER_F90_{INTEGER|REAL|COMPLEX} */

int PNCIO_Type_ispredef(MPI_Datatype datatype, int *flag)
{
    int ret, combiner;
    ret = PNCIO_Type_get_combiner(datatype, &combiner);
    switch (combiner) {
        case MPI_COMBINER_NAMED:
        case MPI_COMBINER_F90_INTEGER:
        case MPI_COMBINER_F90_REAL:
        case MPI_COMBINER_F90_COMPLEX:
            *flag = 1;
            break;
        default:
            *flag = 0;
            break;
    }
    return ret;
}

/* utility function for freeing user-defined datatypes,
 * MPI_DATATYPE_NULL and predefined datatypes are ignored,
 * datatype is set to MPI_DATATYPE_NULL upon return */

int PNCIO_Type_dispose(MPI_Datatype * datatype)
{
    int ret, flag;
    if (*datatype == MPI_DATATYPE_NULL)
        return MPI_SUCCESS;
    ret = PNCIO_Type_ispredef(*datatype, &flag);
    if (ret == MPI_SUCCESS && !flag)
        ret = MPI_Type_free(datatype);
    *datatype = MPI_DATATYPE_NULL;
    return ret;
}

#define CAST_INT32(func_name, count, bklen, disp, dType, newType) {          \
    int kk, iCount, *iBklen;                                                 \
    MPI_Aint *iDisp;                                                         \
    assert(count <= 2147483647); /* overflow 4-byte int */             \
    iCount = (int)count;                                                     \
    iBklen = (int*) NCI_Malloc(sizeof(int) * iCount);                      \
    for (kk=0; kk<iCount; kk++) {                                            \
        assert(bklen[kk] <= 2147483647); /* overflow 4-byte int */     \
        iBklen[kk] = (int)bklen[kk];                                         \
    }                                                                        \
    if (sizeof(MPI_Aint) != sizeof(MPI_Count)) {                             \
        iDisp = (MPI_Aint*) NCI_Malloc(sizeof(MPI_Aint) * iCount);         \
        for (kk=0; kk<iCount; kk++) {                                        \
            assert(disp[kk] <= 2147483647); /* overflow 4-byte int */  \
            iDisp[kk] = (MPI_Aint)disp[kk];                                  \
        }                                                                    \
    }                                                                        \
    else                                                                     \
        iDisp = (MPI_Aint*)disp;                                             \
    ret = func_name(iCount, iBklen, iDisp, dType, newType);                  \
    NCI_Free(iBklen);                                                      \
    if (sizeof(MPI_Aint) != sizeof(MPI_Count))                               \
        NCI_Free(iDisp);                                                   \
}

/* some systems do not have pread/pwrite, or requrie XOPEN_SOURCE set higher
 * than we would like.  see #1973 */
#if (HAVE_DECL_PWRITE == 0)

#include <sys/types.h>
#include <unistd.h>

ssize_t pread(int fd, void *buf, size_t count, off_t offset);
ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset);

ssize_t pread(int fd, void *buf, size_t count, off_t offset)
{
    off_t lseek_ret;
    off_t old_offset;
    ssize_t read_ret;

    old_offset = lseek(fd, 0, SEEK_CUR);
    lseek_ret = lseek(fd, offset, SEEK_SET);
    if (lseek_ret == -1)
        return lseek_ret;
    read_ret = read(fd, buf, count);
    if (read_ret < 0)
        return read_ret;
    /* man page says "file offset is not changed" */
    if ((lseek_ret = lseek(fd, old_offset, SEEK_SET)) < 0)
        return lseek_ret;

    return read_ret;
}

ssize_t pwrite(int fd, const void *buf, size_t count, off_t offset)
{
    off_t lseek_ret;
    off_t old_offset;
    ssize_t write_ret;

    old_offset = lseek(fd, 0, SEEK_CUR);
    lseek_ret = lseek(fd, offset, SEEK_SET);
    if (lseek_ret == -1)
        return lseek_ret;
    write_ret = write(fd, buf, count);
    if (write_ret < 0)
        return write_ret;
    /* man page says "file offset is not changed" */
    if ((lseek_ret = lseek(fd, old_offset, SEEK_SET)) < 0)
        return lseek_ret;

    return write_ret;
}
#endif

void PNCIO_Heap_merge(PNCIO_Access * others_req, MPI_Count * count,
                      MPI_Offset * srt_off, MPI_Count * srt_len, MPI_Count * start_pos,
                      int nprocs, int nprocs_recv, MPI_Count total_elements)
{
    typedef struct {
        MPI_Offset *off_list;
#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Offset *len_list;
#else
        int *len_list;
#endif
        MPI_Count nelem;
    } heap_struct;

    heap_struct *a, tmp;
    int i, j, heapsize, l, r, k, smallest;

    a = (heap_struct *) NCI_Malloc((nprocs_recv + 1) * sizeof(heap_struct));

    j = 0;
    for (i = 0; i < nprocs; i++)
        if (count[i]) {
            a[j].off_list = &(others_req[i].offsets[start_pos[i]]);
            a[j].len_list = &(others_req[i].lens[start_pos[i]]);
            a[j].nelem = count[i];
            j++;
        }

    /* build a heap out of the first element from each list, with
     * the smallest element of the heap at the root */

    heapsize = nprocs_recv;
    for (i = heapsize / 2 - 1; i >= 0; i--) {
        /* Heapify(a, i, heapsize); Algorithm from Cormen et al. pg. 143
         * modified for a heap with smallest element at root. I have
         * removed the recursion so that there are no function calls.
         * Function calls are too expensive. */
        k = i;
        for (;;) {
            l = 2 * (k + 1) - 1;
            r = 2 * (k + 1);

            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                tmp.off_list = a[k].off_list;
                tmp.len_list = a[k].len_list;
                tmp.nelem = a[k].nelem;

                a[k].off_list = a[smallest].off_list;
                a[k].len_list = a[smallest].len_list;
                a[k].nelem = a[smallest].nelem;

                a[smallest].off_list = tmp.off_list;
                a[smallest].len_list = tmp.len_list;
                a[smallest].nelem = tmp.nelem;

                k = smallest;
            } else
                break;
        }
    }

    for (i = 0; i < total_elements; i++) {
        /* extract smallest element from heap, i.e. the root */
        srt_off[i] = *(a[0].off_list);
        srt_len[i] = *(a[0].len_list);
        (a[0].nelem)--;

        if (!a[0].nelem) {
            a[0].off_list = a[heapsize - 1].off_list;
            a[0].len_list = a[heapsize - 1].len_list;
            a[0].nelem = a[heapsize - 1].nelem;
            heapsize--;
        } else {
            (a[0].off_list)++;
            (a[0].len_list)++;
        }

        /* Heapify(a, 0, heapsize); */
        k = 0;
        for (;;) {
            l = 2 * (k + 1) - 1;
            r = 2 * (k + 1);

            if ((l < heapsize) && (*(a[l].off_list) < *(a[k].off_list)))
                smallest = l;
            else
                smallest = k;

            if ((r < heapsize) && (*(a[r].off_list) < *(a[smallest].off_list)))
                smallest = r;

            if (smallest != k) {
                tmp.off_list = a[k].off_list;
                tmp.len_list = a[k].len_list;
                tmp.nelem = a[k].nelem;

                a[k].off_list = a[smallest].off_list;
                a[k].len_list = a[smallest].len_list;
                a[k].nelem = a[smallest].nelem;

                a[smallest].off_list = tmp.off_list;
                a[smallest].len_list = tmp.len_list;
                a[smallest].nelem = tmp.nelem;

                k = smallest;
            } else
                break;
        }
    }
    NCI_Free(a);
}

