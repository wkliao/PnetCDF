dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2014, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id$ */

/*
 * This file implements the corresponding APIs defined in
 * src/dispatchers/var_getput.m4
 *
 * ncmpi_iget_varn_<type>() : dispatcher->iget_varn()
 * ncmpi_iput_varn_<type>() : dispatcher->iput_varn()
 * ncmpi_bput_varn_<type>() : dispatcher->bput_varn()
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <string.h> /* memcpy() */
#include <assert.h>

#include <mpi.h>

#include <pnc_debug.h>
#include <common.h>
#include "ncmpio_NC.h"

/*----< igetput_varn() >-----------------------------------------------------*/
/* The implementation of nonblocking varn APIs is to add "num" iget/iput
 * requests to the nonblocking request queue. These requests may be flattened
 * and sorted together with other nonblocking requests during the wait call.
 * All the "num" nonblocking requests posted by iget/iput/bput_varn() share
 * the same request ID.
 */
static int
igetput_varn(NC                *ncp,
             NC_var            *varp,
             int                num,
             MPI_Offset* const *starts, /* [num][varp->ndims], cannot be NULL */
             MPI_Offset* const *counts, /* [num][varp->ndims], may be NULL */
             void              *buf,
             MPI_Offset         bufcount,
             MPI_Datatype       buftype,  /* data type of the bufer */
             int               *reqid,    /* OUT: request ID */
             int                reqMode)
{
    int i, j, err=NC_NOERR, free_xbuf=0, isize, xsize, abuf_index=-1;
    int isContig=1, need_convert, need_swap, need_swap_back_buf=0, lead_off;
    int mpireturn, rem, new_nreqs;
    size_t memChunk;
    void *xbuf=NULL;
    char *xbufp;
    MPI_Offset nelems, *req_nelems=NULL, nbytes, *start_ptr;
    MPI_Datatype itype, xtype;
    NC_lead_req *lead_req;
    NC_req *req;

    /* if called from a bput API, check if buffer has been attached */
    if (fIsSet(reqMode, NC_REQ_NBB) && ncp->abuf == NULL) {
        DEBUG_ASSIGN_ERROR(err, NC_ENULLABUF)
        goto fn_exit;
    }

    /* validity of starts and counts has been checked at dispatcher layer */

    /* calculate nelems and the number of new non-lead and non-zero requests */
    req_nelems = (MPI_Offset*) NCI_Malloc(sizeof(MPI_Offset) * num);
    new_nreqs = 0;
    nelems    = 0;
    if (counts != NULL) {
        for (i=0; i<num; i++) {
            req_nelems[i] = 1;
            if (counts[i] != NULL) {
                for (j=0; j<varp->ndims; j++)
                    req_nelems[i] *= counts[i][j];
                if (req_nelems[i] == 0)
                    continue; /* ignore this 0-length request i */
                nelems += req_nelems[i];
                new_nreqs += (IS_RECVAR(varp)) ? (int)counts[i][0] : 1;
            }
            else { /* counts[i] == NULL, equivalent to all 1s */
                nelems++;
                new_nreqs++;
            }
        }
    } else { /* counts == NULL, equivalent to all 1s */
        for (i=0; i<num; i++)
            req_nelems[i] = 1;
        new_nreqs = num;
        nelems    = num;
    }

    /* xtype is the MPI element type in external representation, xsize is its
     * size in bytes. Similarly, itype and isize for internal representation.
     */
    xtype = ncmpii_nc2mpitype(varp->xtype);
    mpireturn = MPI_Type_size(xtype, &xsize);
    if (mpireturn != MPI_SUCCESS) {
        err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");
        goto fn_exit;
    }

    if (bufcount == NC_COUNT_IGNORE) {
        /* In this case, this subroutine is called from a high-level API.
         * buftype is one of the MPI predefined primitive data type. We set
         * itype to buftype. itype is the MPI element type in internal
         * representation. In addition, it means the user buf is contiguous.
         */
        itype = buftype;
        mpireturn = MPI_Type_size(itype, &isize); /* buffer element size */
        if (mpireturn != MPI_SUCCESS) {
            err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_size");
            goto fn_exit;
        }
    }
    else if (buftype == MPI_DATATYPE_NULL) {
        /* In this case, bufcount is ignored and the internal buffer data type
         * match the external variable data type. No data conversion will be
         * done. In addition, it means buf is contiguous. Hereinafter, buftype
         * is ignored.
         */
        itype = xtype;
        isize = xsize;
    }
    else { /* (bufcount > 0) */
        /* When bufcount > 0, this subroutine is called from a flexible API. If
         * buftype is noncontiguous, we pack buf into xbuf, a contiguous buffer.
         */
        MPI_Offset bnelems=0;

        /* itype (primitive MPI data type) from buftype
         * isize is the size of itype in bytes
         * bnelems is the number of itype elements in one buftype
         * isContig indicates whether user buffer, buf, is contiguous
         */
        err = ncmpii_dtype_decode(buftype, &itype, &isize, &bnelems,
                                  NULL, &isContig);
        if (err != NC_NOERR) goto fn_exit;

        /* size in bufcount * buftype must match with counts[] */
        if (bnelems * bufcount != nelems) {
            DEBUG_ASSIGN_ERROR(err, NC_EIOMISMATCH)
            goto fn_exit;
        }
    }

    /* nbytes is the amount of this varn request in bytes */
    nbytes = nelems * xsize;

    /* for nonblocking API, return now if request size is zero */
    if (nbytes == 0) goto fn_exit;

    memChunk = sizeof(MPI_Offset) * varp->ndims;

    /* check if type conversion and Endianness byte swap is needed */
    need_convert = ncmpii_need_convert(ncp->format, varp->xtype, itype);
    need_swap    = NEED_BYTE_SWAP(varp->xtype, itype);

    if (fIsSet(reqMode, NC_REQ_WR)) {
        /* check if in-place byte swap can be enabled */
        int can_swap_in_place = 1;
        if (need_swap) {
            if (fIsSet(ncp->flags, NC_MODE_SWAP_OFF))
                /* in-place byte swap is disabled by user through PnetCDF hint
                 * 'nc_in_place_swap'.
                 */
                can_swap_in_place = 0;
            else if (! fIsSet(ncp->flags, NC_MODE_SWAP_ON)) {
                /* auto mode, as user does not explicitly enable it */
                if (nbytes <= NC_BYTE_SWAP_BUFFER_SIZE)
                    /* If write amount is small, disable in-place swap.
                     * This is because the user buffer may be immutable. In
                     * this case, in-place swap will cause segmentation fault.
                     * Immutable buffers are usually small.  */
                    can_swap_in_place = 0;
            }
        }

        /* Because we are going to break the varn request into multiple vara
         * requests, we allocate a contiguous buffer, xbuf, if buftype is not
         * contiguous. So, we can do byte-swap and type-conversion on xbuf.
         */
        if (fIsSet(reqMode, NC_REQ_NBB)) {
            /* for bput call, check if the remaining buffer space is sufficient
             * to accommodate this varn request. If yes, allocate a space for
             * xbuf.
             */
            if (ncp->abuf->size_allocated - ncp->abuf->size_used < nbytes) {
                DEBUG_ASSIGN_ERROR(err, NC_EINSUFFBUF)
                goto fn_exit;
            }
            err = ncmpio_abuf_malloc(ncp, nbytes, &xbuf, &abuf_index);
            if (err != NC_NOERR) goto fn_exit;
        }
        else if (!need_convert && can_swap_in_place && isContig) {
            /* reuse buf and break it into multiple vara requests */
            xbuf = buf;
            if (need_swap) need_swap_back_buf = 1;
        }
        else { /* must allocate a buffer to convert/swap/pack */
            xbuf = NCI_Malloc((size_t)nbytes);
            free_xbuf = 1;
            if (xbuf == NULL) {
                DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
                goto fn_exit;
            }
        }

        /* when necessary, pack buf to xbuf and perform byte-swap and
         * type-conversion on xbuf, which will later be broken into num
         * sub-buffers (no additional malloc further), each to be added to the
         * nonblocking request queue.
         */
        err = ncmpio_pack_xbuf(ncp->format, varp, bufcount, buftype, isContig,
                               nelems, itype, isize, MPI_DATATYPE_NULL,
                               need_convert, need_swap, nbytes, buf, xbuf);
        if (err != NC_NOERR && err != NC_ERANGE) {
            if (fIsSet(reqMode, NC_REQ_NBB))
                ncmpio_abuf_dealloc(ncp, abuf_index);
            else if (free_xbuf)
                NCI_Free(xbuf);
            goto fn_exit;
        }

        /* allocate or expand the lead write request queue */
        if (ncp->numLeadPutReqs % NC_REQUEST_CHUNK == 0)
            ncp->put_lead_list = (NC_lead_req*) NCI_Realloc(ncp->put_lead_list,
                                 (ncp->numLeadPutReqs + NC_REQUEST_CHUNK) *
                                 sizeof(NC_lead_req));

        /* allocate or expand the non-lead write request queue */
        rem = ncp->numPutReqs % NC_REQUEST_CHUNK;
        if (rem) rem = NC_REQUEST_CHUNK - rem;

        if (ncp->put_list == NULL || new_nreqs > rem) {
            size_t req_alloc, nChunks;
            req_alloc = ncp->numPutReqs + new_nreqs;
            nChunks = req_alloc / NC_REQUEST_CHUNK;
            if (req_alloc % NC_REQUEST_CHUNK) nChunks++;
            req_alloc = nChunks * NC_REQUEST_CHUNK * sizeof(NC_req);
            ncp->put_list = (NC_req*) NCI_Realloc(ncp->put_list, req_alloc);
        }

#define SORT_LEAD_LIST_BASED_ON_VAR_BEGIN
#ifdef SORT_LEAD_LIST_BASED_ON_VAR_BEGIN
        /* add the new request to put_lead_list and keep put_lead_list sorted,
         * in an increasing order of variable begin offsets. The best is when
         * user makes varn API calls in an increasing order of variable begin
         * offsets, i.e. fixed-size variables first followed by record
         * variables and in an increasing order of variables IDs. If the
         * pending requests consist of multiple records, then keep the list
         * sorted based on the starting offsets of inidvidual records.
         */
        MPI_Offset req_off = varp->begin;
        if (IS_RECVAR(varp)) req_off += ncp->recsize * starts[0][0];

        for (i=ncp->numLeadPutReqs-1; i>=0; i--) {
            if (ncp->put_lead_list[i].varp->begin <= req_off)
                break;
            /* make space for new lead request */
            ncp->put_lead_list[i+1] = ncp->put_lead_list[i];
            ncp->put_lead_list[i+1].nonlead_off += new_nreqs;
        }
        lead_off = i + 1;
        lead_req = ncp->put_lead_list + lead_off;

        if (lead_off < ncp->numLeadPutReqs) {
            /* req is starting location to insert new non-lead requests */
            req = ncp->put_list + lead_req->nonlead_off;
            /* make space for new non-lead requests */
            for (i=ncp->numPutReqs-1; i>=lead_req->nonlead_off; i--) {
                ncp->put_list[i+new_nreqs] = ncp->put_list[i];
                ncp->put_list[i+new_nreqs].lead_off++;
            }
        }
        else {
            /* append new lead request at the end of ncp->put_lead_list */
            lead_req->nonlead_off = ncp->numPutReqs;
            /* append new non-lead requests at the end of ncp->put_list */
            req = ncp->put_list + ncp->numPutReqs;
        }
#else
        /* append new lead request at the end of ncp->put_lead_list */
        lead_req = ncp->put_lead_list + ncp->numLeadPutReqs;
        lead_req->nonlead_off = ncp->numPutReqs;
        /* append new non-lead requests at the end of ncp->put_list */
        req = ncp->put_list + ncp->numPutReqs;
        lead_off = ncp->numLeadPutReqs;
#endif

        /* the new request ID will be an even number (max of write ID + 2) */
        if (ncp->numLeadPutReqs == 0) {
            lead_req->id = 0;
            ncp->maxPutReqID = 0;
        } else {
            ncp->maxPutReqID += 2;
            lead_req->id = ncp->maxPutReqID;
        }

        ncp->numLeadPutReqs++;
        ncp->numPutReqs += new_nreqs;

        lead_req->flag     = 0;
        lead_req->bufcount = 0;
        lead_req->buftype  = MPI_DATATYPE_NULL;

        /* Only lead requests may free xbuf. For write, when xbuf == buf,
         * the user buffer, buf, may have been byte-swapped. In this case,
         * we need to swap it back after MPI-IO calls.
         */
        if (need_swap_back_buf) fSet(lead_req->flag, NC_REQ_BUF_BYTE_SWAP);
        if (free_xbuf)          fSet(lead_req->flag, NC_REQ_XBUF_TO_BE_FREED);

        /* when abuf_index >= 0 means called by bput_varn */
        lead_req->abuf_index = abuf_index; /* to mark space in abuf free */
    }
    else { /* read request */
        /* Type conversion and byte swap for read, if necessary, will be done
         * at the ncmpi_wait call */

        if (!need_convert && isContig) {
            /* reuse buf in later MPI file read */
            xbuf = buf;
        }
        else { /* must allocate a buffer for read/convert/swap/unpack */
            xbuf = NCI_Malloc((size_t)nbytes);
            if (xbuf == NULL) {
                DEBUG_ASSIGN_ERROR(err, NC_ENOMEM)
                goto fn_exit;
            }
            free_xbuf = 1;
        }

        /* allocate or expand the lead read request queue */
        if (ncp->numLeadGetReqs % NC_REQUEST_CHUNK == 0)
            ncp->get_lead_list = (NC_lead_req*) NCI_Realloc(ncp->get_lead_list,
                                 (ncp->numLeadGetReqs + NC_REQUEST_CHUNK) *
                                 sizeof(NC_lead_req));

        /* allocate or expand the non-lead read request queue */
        rem = ncp->numGetReqs % NC_REQUEST_CHUNK;
        if (rem) rem = NC_REQUEST_CHUNK - rem;

        if (ncp->get_list == NULL || new_nreqs > rem) {
            size_t req_alloc, nChunks;
            req_alloc = ncp->numGetReqs + new_nreqs;
            nChunks = req_alloc / NC_REQUEST_CHUNK;
            if (req_alloc % NC_REQUEST_CHUNK) nChunks++;
            req_alloc = nChunks * NC_REQUEST_CHUNK * sizeof(NC_req);
            ncp->get_list = (NC_req*) NCI_Realloc(ncp->get_list, req_alloc);
        }

#ifdef SORT_LEAD_LIST_BASED_ON_VAR_BEGIN
        /* add the new request to get_lead_list and keep get_lead_list sorted,
         * in an increasing order of variable begin offsets. The best is when
         * user makes varn API calls in an increasing order of variable begin
         * offsets, i.e. fixed-size variables first followed by record
         * variables and in an increasing order of variables IDs.
         */
        for (i=ncp->numLeadGetReqs-1; i>=0; i--) {
            if (ncp->get_lead_list[i].varp->begin <= varp->begin)
                break;
            /* make space for new lead request */
            ncp->get_lead_list[i+1] = ncp->get_lead_list[i];
            ncp->get_lead_list[i+1].nonlead_off += new_nreqs;
        }
        lead_off = i + 1;
        lead_req = ncp->get_lead_list + lead_off;

        if (lead_off < ncp->numLeadGetReqs) {
            /* req is starting location to insert new non-lead requests */
            req = ncp->get_list + lead_req->nonlead_off;
            /* make space for new non-lead requests */
            for (i=ncp->numGetReqs-1; i>=lead_req->nonlead_off; i--) {
                ncp->get_list[i+new_nreqs] = ncp->get_list[i];
                ncp->get_list[i+new_nreqs].lead_off++;
            }
        }
        else {
            /* append new lead request at the end of ncp->get_lead_list */
            lead_req->nonlead_off = ncp->numGetReqs;
            /* append new non-lead request at the end of ncp->get_list */
            req = ncp->get_list + ncp->numGetReqs;
        }
#else
        /* append new request at the end of ncp->get_lead_list */
        lead_req = ncp->get_lead_list + ncp->numLeadGetReqs;
        lead_req->nonlead_off = ncp->numGetReqs;
        /* append new non-lead requests at the end of ncp->get_list */
        req = ncp->get_list + ncp->numGetReqs;
        lead_off = ncp->numLeadGetReqs;
#endif

        /* the new request ID will be an odd number (max of read ID + 2) */
        if (ncp->numLeadGetReqs == 0) {
            lead_req->id = 1;
            ncp->maxGetReqID = 1;
        } else {
            ncp->maxGetReqID += 2;
            lead_req->id = ncp->maxGetReqID;
        }

        ncp->numLeadGetReqs++;
        ncp->numGetReqs += new_nreqs;

        lead_req->flag       = 0;
        lead_req->abuf_index = -1;

        /* Only lead requests may free xbuf. For read, only the lead requests
         * perform byte-swap, type-conversion, imap unpack, and buftype
         * unpacking from xbuf to buf.
         */
        if (need_convert) fSet(lead_req->flag, NC_REQ_BUF_TYPE_CONVERT);
        if (need_swap)    fSet(lead_req->flag, NC_REQ_BUF_BYTE_SWAP);
        if (free_xbuf)    fSet(lead_req->flag, NC_REQ_XBUF_TO_BE_FREED);

        if (isContig) {
            fSet(lead_req->flag, NC_REQ_BUF_TYPE_IS_CONTIG);
            lead_req->bufcount = 0;
            lead_req->buftype  = MPI_DATATYPE_NULL;
        } else {
            /* When read buftype is not contiguous, we duplicate buftype for
             * later used in the wait call to unpack xbuf using buftype to buf.
             */
            mpireturn = MPI_Type_dup(buftype, &lead_req->buftype);
            if (mpireturn != MPI_SUCCESS) {
                err = ncmpii_error_mpi2nc(mpireturn, "MPI_Type_dup");
                goto fn_exit;
            }
            lead_req->bufcount = bufcount;
        }
    }

    lead_req->varp        = varp;
    lead_req->itype       = itype;
    lead_req->xbuf        = xbuf;
    lead_req->buf         = buf;
    lead_req->nelems      = nelems;
    lead_req->imaptype    = MPI_DATATYPE_NULL;
    lead_req->status      = NULL;
    lead_req->max_rec     = -1;
    lead_req->nonlead_num = new_nreqs;

    /* varn APIs have no argument stride */
    fSet(lead_req->flag, NC_REQ_STRIDE_NULL);

    /* Lead request allocates a single array to store num start/count for all
     * non-lead requests, including individual record requests if
     * record variable.
     */
    lead_req->start = (MPI_Offset*) NCI_Malloc(memChunk * 2 * new_nreqs);

    /* break varn into multiple non-lead requests and buf/xbuf accordingly */
    start_ptr = lead_req->start;
    xbufp = (char*)xbuf;

    for (i=0; i<num; i++) {
        req->npairs = 0;

        if (req_nelems[i] == 0) continue; /* ignore this 0-length request i */

        req->nelems    = req_nelems[i];
        req->lead_off  = lead_off;
        req->xbuf      = xbufp;
        xbufp         += req_nelems[i] * xsize;

        /* copy starts[i] and counts[i] over to req */
        req->start = start_ptr;
        memcpy(start_ptr, starts[i], memChunk);
        start_ptr += varp->ndims; /* count[] */
        if (counts == NULL || counts[i] == NULL) {
            /* counts == NULL, equivalent to all 1s */
            for (j=0; j<varp->ndims; j++)
                 start_ptr[j] = 1; /* start_ptr is now counts[] */
        }
        else
            memcpy(start_ptr, counts[i], memChunk);
        start_ptr += varp->ndims;

        if (IS_RECVAR(varp)) {
            /* save the last record number accessed */
            MPI_Offset max_rec, num_rec;

            if (counts == NULL || counts[i] == NULL) num_rec = 1;
            else                                     num_rec = counts[i][0];

            /* calculate number of flattened offset-length pairs */
            req->npairs = 1;
            if (counts == NULL || counts[i] == NULL) {
                /* meaning var1 API, equivalent to all 1s */
                req->offset_start = 0;
                req->offset_end = varp->xsz;
            }
            else {
                for (j=1; j<varp->ndims-1; j++)
                    req->npairs *= counts[i][j];
                /* special treatment for when there is only one pair */
                if (req->npairs == 1) {
                    ncmpio_calc_off(ncp, varp, starts[i], counts[i],
                                    &req->offset_start);
                    req->offset_end = req->nelems * varp->xsz;
                }
            }

            max_rec = starts[i][0] + num_rec;
            lead_req->max_rec = MAX(lead_req->max_rec, max_rec);

            if (num_rec > 1) {
                /* If the number of requested records is more than 1, we split
                 * this request into multiple requests, one for each record.
                 * The number of records is only preserved in the lead request
                 * max_rec. All non-lead record-variable requests counts[i][0]
                 * are set to 1.
                 */
                NC_lead_req *lead_list;
                lead_list = (fIsSet(reqMode, NC_REQ_WR)) ? ncp->put_lead_list
                                                         : ncp->get_lead_list;


                req->nelems /= counts[i][0];
                if (req->npairs == 1)
                    req->offset_end = req->nelems * varp->xsz;

                /* append (counts[i][0]-1) number of requests to the queue */
                ncmpio_add_record_requests(lead_list, req, counts[i][0], NULL);
                start_ptr += (counts[i][0] - 1) * 2 * varp->ndims;
                req += counts[i][0];
            }
            else
                req++;
        }
        else {
            /* calculate number of flattened offset-length pairs */
            req->npairs = 1;
            if (counts == NULL || counts[i] == NULL) {
                /* meaning var1 API, equivalent to all 1s */
                req->offset_start = 0;
                req->offset_end = varp->xsz;
            }
            else {
                for (j=0; j<varp->ndims-1; j++)
                    req->npairs *= counts[i][j];
                /* special treatment for when there is only one pair */
                if (req->npairs == 1) {
                    ncmpio_calc_off(ncp, varp, starts[i], counts[i],
                                    &req->offset_start);
                    req->offset_end = req->nelems * varp->xsz;
                }
            }
            req++;
        }
    }

    if (reqid != NULL) *reqid = lead_req->id;

fn_exit:
    if (req_nelems != NULL) NCI_Free(req_nelems);

    return err;
}


include(`utils.m4')

dnl
define(`IsBput',    `ifelse(`$1',`bput', `1', `0')')dnl
define(`BufConst',  `ifelse(`$1',`get', , `const')')dnl
dnl
dnl VARN(iget/iput/bput)
dnl
define(`VARN',dnl
`dnl
/*----< ncmpio_$1_varn() >----------------------------------------------------*/
int
ncmpio_$1_varn(void               *ncdp,
               int                 varid,
               int                 num,
               MPI_Offset* const  *starts, /* cannot be NULL */
               MPI_Offset* const  *counts, /* may be NULL */
               BufConst(substr($1,1)) void  *buf,
               MPI_Offset          bufcount,
               MPI_Datatype        buftype,
               int                *reqid,
               int                 reqMode)
{
    NC *ncp=(NC*)ncdp;

    if (reqid != NULL) *reqid = NC_REQ_NULL;

    /* Note sanity check for ncdp and varid has been done in the dispatcher.
     * The same for zero-size request checking (return immediately)
     */

    if (fIsSet(reqMode, NC_REQ_ZERO)) return NC_NOERR;

    return igetput_varn(ncp, ncp->vars.value[varid], num, starts, counts,
                        (void*)buf, bufcount, buftype, reqid, reqMode);
}
')dnl
dnl

VARN(iput)
VARN(iget)
VARN(bput)

