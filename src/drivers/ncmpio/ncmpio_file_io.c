/*
 *  Copyright (C) 2018, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset() */

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< get_count() >--------------------------------------------------------*/
/* This subroutine is independent. On success, the number of bytes read/written
 * is returned (zero indicates nothing was read/written). Like POSIX read()/
 * write(), it is not an error if this number is smaller than the number of
 * bytes requested. On error, a negative value, an NC error code, is returned.
 */
static
MPI_Offset get_count(MPI_Status   *mpistatus,
                     MPI_Datatype  datatype)
{
    int mpireturn;

#ifdef HAVE_MPI_TYPE_SIZE_C
    MPI_Count type_size;
    /* MPI_Type_size_c is introduced in MPI 4.0 */
    MPI_Type_size_c(datatype, &type_size);
#elif defined(HAVE_MPI_TYPE_SIZE_X)
    MPI_Count type_size;
    /* MPI_Type_size_x is introduced in MPI 3.0 */
    MPI_Type_size_x(datatype, &type_size);
#else
    int type_size;
    MPI_Type_size(datatype, &type_size);
#endif

#ifdef HAVE_MPI_GET_COUNT_C
    MPI_Count count;
    mpireturn = MPI_Get_count_c(mpistatus, datatype, &count);
#else
    int count;
    mpireturn = MPI_Get_count(mpistatus, datatype, &count);
#endif

    if (mpireturn != MPI_SUCCESS || count == MPI_UNDEFINED)
        /* In case of partial read/write, MPI_Get_elements() is supposed to be
         * called to obtain the number of type map elements actually
         * read/written in order to calculate the true read/write amount. Below
         * skips this step and simply returns the partial read/write amount.
         * See an example usage of MPI_Get_count() in Example 5.12 from MPI
         * standard document.
         */
        return NC_EFILE;

    return (MPI_Offset)count * type_size;
}

/*----< ncmpio_file_read_at() >--------------------------------------------*/
/*
 * This function is independent.
 */
/* TODO: move check count against MAX_INT and call _c API */
MPI_Offset
ncmpio_file_read_at(NC           *ncp,
                    MPI_Offset    offset,
                    void         *buf,
                    MPI_Offset    count,
                    MPI_Datatype  buftype)
{
    int err=NC_NOERR, mpireturn;
    MPI_Status mpistatus;

    /* explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

#ifdef HAVE_MPI_LARGE_COUNT
        TRACE_IO(MPI_File_read_at_c, (fh, offset, buf, (MPI_Count)count,
                                      buftype, &mpistatus));
#else
        if (count > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW count="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,count);
#endif
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
        TRACE_IO(MPI_File_read_at, (fh, offset, buf, (int)count,
                                    buftype, &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
        }
    }
    else {
        err = PNCIO_File_read_at(ncp->adio_fh, offset, buf, count,
                                 buftype, &mpistatus);
    }

    /* update the number of bytes read since file open */
    if (err == NC_NOERR) {
        MPI_Offset amnt;
        amnt = get_count(&mpistatus, buftype);
        if (amnt >= 0) ncp->get_size += amnt;
        /* else: ignore if error, as this error is not fatal */
        return amnt;
    }
    else
        return err;
}

/*----< ncmpio_file_read_at_all() >-----------------------------------------*/
/*
 * This function is collective.
 */
MPI_Offset
ncmpio_file_read_at_all(NC           *ncp,
                        MPI_Offset    offset,
                        void         *buf,
                        MPI_Offset    count,
                        MPI_Datatype  buftype)
{
    int err=NC_NOERR, mpireturn;
    MPI_Status mpistatus;

    /* explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

#ifdef HAVE_MPI_LARGE_COUNT
        TRACE_IO(MPI_File_read_at_all_c, (fh, offset, buf, (MPI_Count)count,
                                          buftype, &mpistatus));
#else
        if (count > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW count="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,count);
#endif
            DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
            /* participate the collective call, but read nothing */
            count = 0;
        }
        TRACE_IO(MPI_File_read_at_all, (fh, offset, buf, (int)count,
                                        buftype, &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EREAD)
        }
    }
    else {
        err = PNCIO_File_read_at_all(ncp->adio_fh, offset, buf, count,
                                     buftype, &mpistatus);
    }

    /* update the number of bytes read since file open */
    if (err == NC_NOERR) {
        MPI_Offset amnt;
        amnt = get_count(&mpistatus, buftype);
        if (amnt >= 0) ncp->get_size += amnt;
        /* else: ignore if error, as this error is not fatal */
        return amnt;
    }
    else
        return err;
}

/*----< ncmpio_file_write_at() >--------------------------------------------*/
/*
 * This function is independent.
 */
MPI_Offset
ncmpio_file_write_at(NC           *ncp,
                     MPI_Offset    offset,
                     const void   *buf,
                     MPI_Offset    count,
                     MPI_Datatype  buftype)
{
    int err=NC_NOERR, mpireturn;
    MPI_Status mpistatus;

    /* explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

#ifdef HAVE_MPI_LARGE_COUNT
        TRACE_IO(MPI_File_write_at_c, (fh, offset, buf, (MPI_Count)count,
                                     buftype, &mpistatus));
#else
        if (count > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW count="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,count);
#endif
            DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
        TRACE_IO(MPI_File_write_at, (fh, offset, buf, (int)count,
                                     buftype, &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE)
        }
    }
    else {
        err = PNCIO_File_write_at(ncp->adio_fh, offset, buf, count,
                                  buftype, &mpistatus);
    }

    /* update the number of bytes written since file open */
    if (err == NC_NOERR) {
        MPI_Offset amnt;
        amnt = get_count(&mpistatus, buftype);
        if (amnt >= 0) ncp->put_size += amnt;
        /* else: ignore if error, as this error is not fatal */
        return amnt;
    }
    else
        return err;
}

/*----< ncmpio_file_write_at_all() >-----------------------------------------*/
/*
 * This function is collective.
 */
MPI_Offset
ncmpio_file_write_at_all(NC           *ncp,
                         MPI_Offset    offset,
                         const void   *buf,
                         MPI_Offset    count,
                         MPI_Datatype  buftype)
{
    int err=NC_NOERR, mpireturn;
    MPI_Status mpistatus;

    /* explicitly initialize mpistatus object to 0. For zero-length read/write,
     * MPI_Get_count may report incorrect result for some MPICH version,
     * due to the uninitialized MPI_Status object passed to MPI-IO calls.
     * Thus we initialize it above to work around. See MPICH ticket:
     * https://trac.mpich.org/projects/mpich/ticket/2332
     */
    memset(&mpistatus, 0, sizeof(MPI_Status));

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        MPI_File fh;

        fh = fIsSet(ncp->flags, NC_MODE_INDEP)
           ? ncp->independent_fh : ncp->collective_fh;

// printf("%s at %d: offset=%lld count=%lld\n",__func__,__LINE__,offset,count);
#ifdef HAVE_MPI_LARGE_COUNT
        TRACE_IO(MPI_File_write_at_all_c, (fh, offset, buf, (MPI_Count)count,
                                         buftype, &mpistatus));
#else
        if (count > NC_MAX_INT) {
#ifdef PNETCDF_DEBUG
            fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW count="OFFFMT"\n",
                    ncp->rank, __func__,__LINE__,count);
#endif
            DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
            /* participate the collective call, but write nothing */
            count = 0;
        }
        TRACE_IO(MPI_File_write_at_all, (fh, offset, buf, (int)count,
                                         buftype, &mpistatus));
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (err == NC_EFILE) DEBUG_ASSIGN_ERROR(err, NC_EWRITE)
        }
    }
    else {
// printf("%s at %d: offset=%lld count=%lld\n",__func__,__LINE__,offset,count);
        err = PNCIO_File_write_at_all(ncp->adio_fh, offset, buf, count,
                                      buftype, &mpistatus);
    }

    /* update the number of bytes written since file open */
    if (err == NC_NOERR) {
        MPI_Offset amnt;
        amnt = get_count(&mpistatus, buftype);
        if (amnt >= 0) ncp->put_size += amnt;
        /* else: ignore if error, as this error is not fatal */
        return amnt;
    }
    else
        return err;
}

/*----< ncmpio_getput_zero_req() >-------------------------------------------*/
/* This function is called when this process has zero-length I/O request and
 * must participate all the MPI collective calls involved in the collective
 * APIs and wait_all(), which include setting fileview, collective read/write,
 * another setting fileview.
 *
 * This function is collective.
 */
int
ncmpio_getput_zero_req(NC *ncp, int reqMode)
{
    int err, status=NC_NOERR;
    MPI_Offset rlen, wlen;

    /* When intra-node aggregation is enabled, non-aggregators do not access
     * the file.
     */
    if (ncp->num_aggrs_per_node > 0 && ncp->rank != ncp->my_aggr)
        return NC_NOERR;

    /* do nothing if this came from an independent API */
    if (fIsSet(reqMode, NC_REQ_INDEP)) return NC_NOERR;

    err = ncmpio_file_set_view(ncp, 0, MPI_BYTE, 0, NULL, NULL);
    if (status == NC_NOERR) status = err;

    if (fIsSet(reqMode, NC_REQ_RD)) {
        if (ncp->nprocs > 1)
            rlen = ncmpio_file_read_at_all(ncp, 0, NULL, 0, MPI_BYTE);
        else
            rlen = ncmpio_file_read_at(ncp, 0, NULL, 0, MPI_BYTE);
        if (status == NC_NOERR && rlen < 0) status = (int)rlen;
    }
    else { /* write request */
        if (ncp->nprocs > 1)
            wlen = ncmpio_file_write_at_all(ncp, 0, NULL, 0, MPI_BYTE);
        else
            wlen = ncmpio_file_write_at(ncp, 0, NULL, 0, MPI_BYTE);
        if (status == NC_NOERR && wlen < 0) status = (int)wlen;
    }

    /* No longer need to reset the file view, as the root's fileview includes
     * the whole file header.
     */

    return status;
}

/*----< ncmpio_read_write() >------------------------------------------------*/
int
ncmpio_read_write(NC           *ncp,
                  int           rw_flag,     /* NC_REQ_WR or NC_REQ_RD */
                  MPI_Offset    offset,
                  MPI_Offset    buf_count,
                  MPI_Datatype  buf_type,
                  void         *buf,
                  int           buftype_is_contig)
{
    char *mpi_name;
    int status=NC_NOERR, err=NC_NOERR, mpireturn, coll_indep;
    MPI_Offset req_size, rlen, wlen;
// printf("%s at %d offset=%lld buf_count=%lld\n",__func__,__LINE__,offset,buf_count);

#ifdef HAVE_MPI_TYPE_SIZE_C
    MPI_Count btype_size;
    /* MPI_Type_size_c is introduced in MPI 4.0 */
    mpireturn = MPI_Type_size_c(buf_type, &btype_size);
    mpi_name = "MPI_Type_size_c";
#elif defined(HAVE_MPI_TYPE_SIZE_X)
    MPI_Count btype_size;
    /* MPI_Type_size_x is introduced in MPI 3.0 */
    mpireturn = MPI_Type_size_x(buf_type, &btype_size);
    mpi_name = "MPI_Type_size_x";
#else
    int btype_size;
    mpireturn = MPI_Type_size(buf_type, &btype_size);
    mpi_name = "MPI_Type_size";
#endif
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        /* return the first encountered error if there is any */
        err = (err == NC_EFILE) ? NC_EREAD : err;
    }
    else if (btype_size == MPI_UNDEFINED) {
#ifdef PNETCDF_DEBUG
        fprintf(stderr,"%d: %s line %d: btype_size MPI_UNDEFINED buf_count="OFFFMT"\n",
                ncp->rank, __func__,__LINE__,buf_count);
#endif
        DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW)
    }

    coll_indep = NC_REQ_INDEP;
    if (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
        coll_indep = NC_REQ_COLL;

    if (err != NC_NOERR) {
        if (coll_indep == NC_REQ_COLL) {
            DEBUG_ASSIGN_ERROR(status, err)
            /* write nothing, but participate the collective call */
            buf_count = 0;
        }
        else
            DEBUG_RETURN_ERROR(err)
    }

    /* request size in bytes, may be > NC_MAX_INT */
    req_size = buf_count * btype_size;

    if (rw_flag == NC_REQ_RD) {
        void         *xbuf=buf;
        MPI_Datatype  xbuf_type=buf_type;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count xlen = (MPI_Count)buf_count;
#else
        int xlen = (int)buf_count;

        if (buf_count > NC_MAX_INT) {
            if (coll_indep == NC_REQ_COLL) {
#ifdef PNETCDF_DEBUG
                fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buf_count="OFFFMT"\n",
                        ncp->rank, __func__,__LINE__,buf_count);
#endif
                DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                /* write nothing, but participate the collective call */
                xlen = 0;
            }
            else
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
#endif

        if (xlen > 0 && !buftype_is_contig && req_size <= ncp->ibuf_size) {
            /* if read buffer is noncontiguous and size is < ncp->ibuf_size,
             * allocate a temporary buffer and use it to read, as some MPI,
             * e.g. Cray on KNL, can be significantly slow when read buffer is
             * noncontiguous.
             */
#ifdef HAVE_MPI_LARGE_COUNT
            xbuf_type = MPI_BYTE;
            xlen = (MPI_Count)req_size;
#else
            if (req_size > NC_MAX_INT) {
                mpireturn = MPI_Type_contiguous(xlen, buf_type, &xbuf_type);
                if (mpireturn != MPI_SUCCESS) {
                    err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_contiguous");
                    if (coll_indep == NC_REQ_COLL)
                        DEBUG_ASSIGN_ERROR(status, err)
                    else
                        DEBUG_RETURN_ERROR(err)
                }
                MPI_Type_commit(&xbuf_type);
                xlen = 1;
            }
            else {
                xbuf_type = MPI_BYTE;
                xlen = (int)req_size;
            }
#endif
            xbuf = NCI_Malloc((size_t)req_size);
        }

// printf("%s at %d: offset=%lld xlen=%lld\n",__func__,__LINE__,offset,xlen);
        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            rlen = ncmpio_file_read_at_all(ncp, offset, xbuf, xlen, xbuf_type);
        else
            rlen = ncmpio_file_read_at(ncp, offset, xbuf, xlen, xbuf_type);
        if (status == NC_NOERR && rlen < 0) status = (int)rlen;

        if (xbuf != buf) { /* unpack contiguous xbuf to noncontiguous buf */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count pos=0;
            mpireturn = MPI_Unpack_c(xbuf, xlen, &pos, buf, (MPI_Count)buf_count,
                                     buf_type, MPI_COMM_SELF);
            mpi_name = "MPI_Unpack_c";
#else
            int pos=0;
            mpireturn = MPI_Unpack(xbuf, xlen, &pos, buf, (int)buf_count,
                                   buf_type, MPI_COMM_SELF);
            mpi_name = "MPI_Unpack";
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                if (coll_indep == NC_REQ_COLL)
                    DEBUG_ASSIGN_ERROR(status, err)
                else
                    DEBUG_RETURN_ERROR(err)
            }
            NCI_Free(xbuf);
        }
        if (xbuf_type != buf_type && xbuf_type != MPI_BYTE)
            MPI_Type_free(&xbuf_type);
    } else { /* NC_REQ_WR */
        void         *xbuf=buf;
        MPI_Datatype  xbuf_type=buf_type;

#ifdef HAVE_MPI_LARGE_COUNT
        MPI_Count xlen = (MPI_Count)buf_count;
#else
        int xlen = (int)buf_count;
        if (buf_count > NC_MAX_INT) {
            if (coll_indep == NC_REQ_COLL) {
#ifdef PNETCDF_DEBUG
                fprintf(stderr,"%d: %s line %d:  NC_EINTOVERFLOW buf_count="OFFFMT"\n",
                        ncp->rank, __func__,__LINE__,buf_count);
#endif
                DEBUG_ASSIGN_ERROR(status, NC_EINTOVERFLOW)
                /* write nothing, but participate the collective call */
                xlen = 0;
            }
            else
                DEBUG_RETURN_ERROR(NC_EINTOVERFLOW)
        }
#endif

        if (xlen > 0 && !buftype_is_contig && req_size <= ncp->ibuf_size) {
            /* if write buffer is noncontiguous and size is < ncp->ibuf_size,
             * allocate a temporary buffer and use it to write, as some MPI,
             * e.g. Cray on KNL, can be significantly slow when write buffer is
             * noncontiguous.
             */
#ifdef HAVE_MPI_LARGE_COUNT
            MPI_Count pos=0;
            xbuf_type = MPI_BYTE;
            xlen = (MPI_Count)req_size;
            xbuf = NCI_Malloc(req_size);
            mpireturn = MPI_Pack_c(buf, (MPI_Count)buf_count, buf_type, xbuf,
                                   (MPI_Count)req_size, &pos, MPI_COMM_SELF);
            mpi_name = "MPI_Pack_c";
#else
            if (req_size > NC_MAX_INT) {
                /* skip packing write data into a temp buffer */
                xlen = (int)buf_count;
                xbuf_type = buf_type;
                mpireturn = MPI_SUCCESS;
            }
            else {
                int pos=0;
                xbuf_type = MPI_BYTE;
                xlen = (int)req_size;
                xbuf = NCI_Malloc(xlen);
                mpireturn = MPI_Pack(buf, (int)buf_count, buf_type, xbuf,
                                     xlen, &pos, MPI_COMM_SELF);
                mpi_name = "MPI_Pack";
            }
#endif
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
                if (coll_indep == NC_REQ_COLL)
                    DEBUG_ASSIGN_ERROR(status, err)
                else
                    DEBUG_RETURN_ERROR(err)
            }
        }

/*
int size; MPI_Type_size(xbuf_type, &size);
MPI_Aint extent; MPI_Type_extent(xbuf_type, &extent);
printf("%s at %d: offset=%lld xlen=%lld xbuf_type=%s size=%d extent=%ld\n",__func__,__LINE__,offset,xlen,(xbuf_type==MPI_BYTE)?"MPI_BYTE":"NOT MPI_BYTE",size,extent);
*/
        if (ncp->nprocs > 1 && coll_indep == NC_REQ_COLL)
            wlen = ncmpio_file_write_at_all(ncp, offset, xbuf, xlen, xbuf_type);
        else
            wlen = ncmpio_file_write_at(ncp, offset, xbuf, xlen, xbuf_type);
        if (status == NC_NOERR && wlen < 0) status = (int)wlen;

        if (xbuf != buf) NCI_Free(xbuf);
        if (xbuf_type != buf_type && xbuf_type != MPI_BYTE)
            MPI_Type_free(&xbuf_type);
    }

    /* fileview is never reused in PnetCDF */
/*
    MPI_Offset disp=0;
    ncmpio_file_set_view(ncp, &disp, MPI_BYTE, 0, NULL, NULL);
if (ncp->adio_fh != NULL) ncp->adio_fh->flat_file.count=0;
*/

    return status;
}

/*----< ncmpio_file_close() >------------------------------------------------*/
/*
 * This function is collective.
 */
int
ncmpio_file_close(NC *ncp)
{
    int err=NC_NOERR;

    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        int mpireturn;

        if (ncp->independent_fh != ncp->collective_fh &&
            ncp->independent_fh != MPI_FILE_NULL) {
            TRACE_IO(MPI_File_close, (&ncp->independent_fh));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }

        if (ncp->collective_fh != MPI_FILE_NULL) {
            TRACE_IO(MPI_File_close, (&ncp->collective_fh));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }
    }
    else {
        /* When intra-node aggregation is enabled, only aggregators have a
         * non-NULL ncp->adio_fh and non-aggregators has adio_fh == NULL.
         */
        if (ncp->adio_fh != NULL) {
            err = PNCIO_File_close(ncp->adio_fh);
            NCI_Free(ncp->adio_fh);
            ncp->adio_fh = NULL;
        }
    }

    return err;
}

/*----< ncmpio_file_delete() >-----------------------------------------------*/
/*
 * This function is collective.
 *
 * This subroutine is called only from ncmpi_abort. When the file is being
 * created and an error occurs, the program is still in define mode. In this
 * case, the file is deleted.
 */
int
ncmpio_file_delete(NC *ncp)
{
    int err=NC_NOERR;

    if (ncp->rank == 0) {
        if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
            char *mpi_name;
            int mpireturn;
            TRACE_IO(MPI_File_delete, ((char *)ncp->path, ncp->mpiinfo));
            if (mpireturn != MPI_SUCCESS)
                err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        }
        else
            err = PNCIO_File_delete(ncp->path);
    }

    if (ncp->nprocs > 1)
        MPI_Bcast(&err, 1, MPI_INT, 0, ncp->comm);

    return err;
}

/*----< ncmpio_file_sync() >-------------------------------------------------*/
/* This function must be called collectively, no matter if it is in collective
 * or independent data mode.
 */
int
ncmpio_file_sync(NC *ncp) {
    char *mpi_name;
    int mpireturn;

    if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
        if (ncp->adio_fh == NULL)
            return NC_NOERR;
        return PNCIO_File_sync(ncp->adio_fh);
    }

    /* the remaining of this subroutine are for when using MPI-IO */

    if (ncp->independent_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->independent_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }
    /* when nprocs == 1, ncp->collective_fh == ncp->independent_fh */
    if (ncp->nprocs == 1) return NC_NOERR;

    /* When intra-node aggregation is enabled, non-aggregator's
     * ncp->collective_fh is always MPI_FILE_NULL. When disabled,
     * ncp->collective_fh on all ranks is never MPI_FILE_NULL as collective
     * mode is default in PnetCDF.
     */
    if (ncp->collective_fh != MPI_FILE_NULL) {
        TRACE_IO(MPI_File_sync, (ncp->collective_fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }

    /* Barrier is not necessary ...
      TRACE_COMM(MPI_Barrier)(ncp->comm);
     */

    return NC_NOERR;
}

/*----< prepend_header() >---------------------------------------------------*/
/* Create a new data type by prepending the whole file header to old_type */
static
int prepend_header(const NC     *ncp,
                   MPI_Offset    disp,
                   MPI_Datatype  old_type,
                   MPI_Datatype *new_type)
{
    int err=NC_NOERR, mpireturn;
    MPI_Datatype ftypes[2];
#ifdef HAVE_MPI_LARGE_COUNT
    MPI_Count blocklens[2];
    MPI_Count disps[2];

    blocklens[0] = ncp->begin_var;
#else
    int blocklens[2];
    MPI_Aint disps[2];

    /* check if header size > 2^31 */
    if (ncp->begin_var > NC_MAX_INT) {
        DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW);
        goto err_out;
    }

    blocklens[0] = (int)ncp->begin_var;
#endif

    /* first block is the header extent */
        disps[0] = 0;
       ftypes[0] = MPI_BYTE;

    /* second block is old_type, the subarray request(s) to the variable */
    blocklens[1] = 1;
        disps[1] = disp;
       ftypes[1] = old_type;

#if !defined(HAVE_MPI_LARGE_COUNT) && (SIZEOF_MPI_AINT != SIZEOF_MPI_OFFSET)
    if (disp > NC_MAX_INT) {
        DEBUG_ASSIGN_ERROR(err, NC_EINTOVERFLOW);
        goto err_out;
    }
#endif

#ifdef HAVE_MPI_LARGE_COUNT
    mpireturn = MPI_Type_create_struct_c(2, blocklens, disps, ftypes,
                                         new_type);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_struct_c");
        goto err_out;
    }
#else
    mpireturn = MPI_Type_create_struct(2, blocklens, disps, ftypes,
                                       new_type);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_create_struct");
        goto err_out;
    }
#endif

err_out:
    if (err == NC_NOERR)
        MPI_Type_commit(new_type);
    else
        *new_type = MPI_BYTE;

    return err;
}

/*----< ncmpio_file_set_view() >---------------------------------------------*/
/* This subroutine is collective when using MPI-IO. When using internal I/O
 * driver, this subroutine can be call independently.
 *
 * This subroutine prepends an additional view of the entire file header for
 * root process, as root may have to update the number of records or attributes
 * in the file header while in data mode. In PnetCDF design, only root process
 * can read/write the file header. For non-root processes, no special treatment
 * is given.
 *
 * filetype: this argument passed into this subroutine does not include file
 *     header. This is why argument disp passed into this subroutine is always
 *     ncp->begin_var (i.e. file offset of data section).
 *
 * disp: this argument passed into this subroutine is always ncp->begin_var.
 *       Its value will be changed to the following when returned.
 *           For root, disp will be set to ncp->begin_var;
 *           For others, disp will be set to 0.
 *       It will be used in the successive MPI-IO calls to read/write variables
 *       in the data section (by all processes) or metadata in the header
 *       section (by root only).
 */
int
ncmpio_file_set_view(const NC     *ncp,
                     MPI_Offset   *disp,    /* IN/OUT */
                     MPI_Datatype  filetype,
                     MPI_Aint      npairs,
#ifdef HAVE_MPI_LARGE_COUNT
                     MPI_Count    *offsets,
                     MPI_Count    *lengths
#else
                     MPI_Offset   *offsets,
                     int          *lengths
#endif
)
{
    char *mpi_name;
    int err, mpireturn, status=NC_NOERR;
    MPI_File fh;

// printf("%s line %d: filetype = %s\n",__func__,__LINE__,(filetype == MPI_DATATYPE_NULL)?"NULL":"NOT NULL");

    if (ncp->fstype != PNCIO_FSTYPE_MPIIO) {
        /* Skip setting fileview for ranks whose adio_fh is NULL */
        if (ncp->adio_fh == NULL)
            return NC_NOERR;

        if (filetype == MPI_DATATYPE_NULL) {
            /* When PnetCDF's internal I/O driver is used and this is called
             * from intra_node_aggregation() which passes MPI_DATATYPE_NULL in
             * argument filetype to help this subroutine to tell this is called
             * from intra_node_aggregation(). In this case, this rank's
             * fileview has already been flattened into offset-length pairs
             * which can directly be used by the ADIO driver.
             *
             * Note when this subroutine is not called from
             * intra_node_aggregation(), i.e. intra-node aggregation is
             * disabled, this rank's fileview may not be flattened and both
             * offsets and lengths are NULL.
             */

            /* By this time, fd->flat_file should have need free by the
             * callback of MPI_Type_free(), i.e. previous round of PnetCDF I/O
             * operation.
             */
/*
            if (ncp->adio_fh->flat_file != NULL)
                NCI_Free(ncp->adio_fh->flat_file);
            ncp->adio_fh->flat_file = NULL;
*/

            /* Note: The passed-in offsets and lengths are not relative to any
             * MPI-IO fileview. They are flattened byte offsets and sizes.
             */
            *disp = (ncp->rank == 0) ? ncp->begin_var : 0;

            /* Pass the already flattened offsets and lengths to ADIO driver as
             * a flattened file type struct. This is avoid repeaated work of
             * constructing and flattening the same datatype.
             */
            return PNCIO_File_set_view(ncp->adio_fh, 0, filetype, npairs,
                                       offsets, lengths);
        }
        else { /* called from subroutines that are not using INA */
// printf("%s at %d disp=%lld filetype = %s\n",__func__,__LINE__,*disp, (filetype == MPI_BYTE)?"MPI_BYTE":"NOT MPI_BYTE");
/*
            if (ncp->adio_fh->flat_file != NULL) {
                NCI_Free(ncp->adio_fh->flat_file);
                ncp->adio_fh->flat_file = NULL;
            }
 */

            if (filetype == MPI_BYTE)
                /* make the whole file visible */
                return PNCIO_File_set_view(ncp->adio_fh, 0, MPI_BYTE, 0, NULL,
                                           NULL);

            if (ncp->rank == 0) {
                MPI_Datatype root_filetype;

                /* prepend the whole file header to filetype */
                err = prepend_header(ncp, *disp, filetype, &root_filetype);
                if (status == NC_NOERR) status = err;

                err = PNCIO_File_set_view(ncp->adio_fh, 0, root_filetype, 0,
                                          NULL, NULL);
                if (status == NC_NOERR) status = err;

                if (root_filetype != MPI_BYTE)
                    MPI_Type_free(&root_filetype);

                /* update the explicit disp to be used in MPI-IO call later */
                *disp = ncp->begin_var;
            }
            else {
                err = PNCIO_File_set_view(ncp->adio_fh, *disp, filetype, 0,
                                          NULL, NULL);
                if (status == NC_NOERR) status = err;

                /* the explicit disp is already set in fileview */
                *disp = 0;
            }
            return status;
        }
    }

    /* Now, ncp->fstype == PNCIO_FSTYPE_MPIIO */
    int to_free_filetype=0;
    MPI_Datatype root_filetype=MPI_BYTE;

    /* when ncp->nprocs == 1, ncp->collective_fh == ncp->independent_fh */
    fh = (ncp->nprocs > 1 && !fIsSet(ncp->flags, NC_MODE_INDEP))
       ? ncp->collective_fh : ncp->independent_fh;

    if (filetype == MPI_BYTE) {
        /* filetype is a contiguous space, make the whole file visible */
        TRACE_IO(MPI_File_set_view, (fh, 0, MPI_BYTE, MPI_BYTE, "native",
                                     MPI_INFO_NULL));
        return NC_NOERR;
    }

    if (filetype == MPI_DATATYPE_NULL) {
        /* This is called from an INA subroutine, or fillerup_aggregate(). We
         * construct (overwrite) filetype using offsets and lengths. Root's
         * offsets and lengths have already included the file header.
         */
// printf("%s at %d: disp=%lld npairs=%ld off=%lld len=%lld\n",__func__,__LINE__, *disp,npairs,offsets[0],lengths[0]);
#ifdef HAVE_MPI_LARGE_COUNT
        /* construct fileview */
        mpireturn = MPI_Type_create_hindexed_c(npairs, lengths, offsets,
                                               MPI_BYTE, &filetype);
#else
        assert(sizeof(*offsets) == sizeof(MPI_Aint));
        /* construct fileview */
        mpireturn = MPI_Type_create_hindexed(npairs, lengths,
                                             (MPI_Aint*)offsets,
                                             MPI_BYTE, &filetype);
#endif
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_create_hindexed");
            /* return the first encountered error if there is any */
            if (status == NC_NOERR) status = err;
        }
        else {
            mpireturn = MPI_Type_commit(&filetype);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn,"MPI_Type_commit");
                /* return the first encountered error if there is any */
                if (status == NC_NOERR) status = err;
            }
            else
                to_free_filetype = 1;
        }
    }
    else if (ncp->rank == 0) { /* NOT called from an INA subroutine */
        /* prepend the whole file header to filetype */
/*
printf("%s at %d: prepend_header disp=%lld\n",__func__,__LINE__, *disp);
int size; MPI_Type_size(filetype, &size);
MPI_Aint extent; MPI_Type_extent(filetype, &extent);
printf("%s at %d: disp=%lld filetype size=%d extent=%ld\n",__func__,__LINE__, *disp, size,extent);
*/
        err = prepend_header(ncp, *disp, filetype, &root_filetype);
        if (status == NC_NOERR) status = err;
        filetype = root_filetype;
        *disp = 0; /* root's fileview includes file header */
    }
/*
int size; MPI_Type_size(filetype, &size);
MPI_Aint extent; MPI_Type_extent(filetype, &extent);
printf("%s at %d: disp=%lld filetype size=%d extent=%ld\n",__func__,__LINE__, *disp, size,extent);
 */

    TRACE_IO(MPI_File_set_view, (fh, *disp, MPI_BYTE, filetype, "native",
                                 MPI_INFO_NULL));
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
        if (status == NC_NOERR) status = err;
    }

    if (root_filetype != MPI_BYTE)
        MPI_Type_free(&root_filetype);

    /* update the explicit disp to be used in MPI-IO call later */
    *disp = (ncp->rank == 0) ? ncp->begin_var : 0;

#if 0
    if (ncp->rank == 0) {
        MPI_Datatype root_filetype;

        /* prepend the whole file header to filetype */
        err = prepend_header(ncp, *disp, filetype, &root_filetype);
        if (status == NC_NOERR) status = err;

        TRACE_IO(MPI_File_set_view, (fh, 0, MPI_BYTE, root_filetype, "native",
                                     MPI_INFO_NULL));
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (status == NC_NOERR) status = err;
        }

        if (root_filetype != MPI_BYTE)
            MPI_Type_free(&root_filetype);

        /* update the explicit disp to be used in MPI-IO call later */
        *disp = ncp->begin_var;
    }
    else {
        TRACE_IO(MPI_File_set_view, (fh, *disp, MPI_BYTE, filetype, "native",
                                     MPI_INFO_NULL));
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
            if (status == NC_NOERR) status = err;
        }

        /* the explicit disp is already set in fileview */
        *disp = 0;
    }
#endif

    if (to_free_filetype)
        MPI_Type_free(&filetype);

    return status;
}

/*----< ncmpio_file_open() >-------------------------------------------------*/
int
ncmpio_file_open(NC         *ncp,
                 MPI_Comm    comm,
                 const char *path,
                 int         omode,
                 MPI_Info    info)
{
    int err=NC_NOERR;

    /* open file collectively */
    if (ncp->fstype == PNCIO_FSTYPE_MPIIO) {
        char *mpi_name;
        int mpireturn;
        MPI_File fh;

        TRACE_IO(MPI_File_open, (comm, path, omode, info, &fh));
        if (mpireturn != MPI_SUCCESS)
            return ncmpii_error_mpi2nc(mpireturn, mpi_name);

        /* Now the file has been successfully opened */
        ncp->collective_fh  = fh;
        ncp->independent_fh = (ncp->nprocs > 1) ? MPI_FILE_NULL : fh;

        /* get the I/O hints used/modified by MPI-IO */
        TRACE_IO(MPI_File_get_info, (fh, &ncp->mpiinfo));
        if (mpireturn != MPI_SUCCESS)
            err = ncmpii_error_mpi2nc(mpireturn, mpi_name);
    }
    else { /* ncp->fstype != PNCIO_FSTYPE_MPIIO */
        ncp->adio_fh = (PNCIO_File*) NCI_Calloc(1,sizeof(PNCIO_File));

        err = PNCIO_File_open(comm, path, omode, info, ncp->adio_fh);
        if (err != NC_NOERR) return err;

        /* Now the file has been successfully opened, obtain the I/O hints
         * used/modified by ADIO driver.
         */
        err = PNCIO_File_get_info(ncp->adio_fh, &ncp->mpiinfo);
    }

    return err;
}

