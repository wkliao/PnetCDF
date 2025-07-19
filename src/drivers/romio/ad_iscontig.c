/*
 *  Copyright (C) 2025, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <adio.h>

void PNCIO_Datatype_iscontig(MPI_Datatype datatype, int *flag)
{
    int combiner;

    ADIOI_Type_get_combiner(datatype, &combiner);

    switch (combiner) {
        case MPI_COMBINER_NAMED:
            *flag = 1;
            break;
#ifdef MPIIMPL_HAVE_MPI_COMBINER_DUP
        case MPI_COMBINER_DUP:
#endif
        case MPI_COMBINER_CONTIGUOUS:
            {
                int *ints;
                MPI_Aint *adds;
                MPI_Count *cnts;
                MPI_Datatype *types;
                ADIOI_Assert(datatype != MPI_DATATYPE_NULL);
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Count nints, nadds, ncnts, ntypes;
                MPI_Type_get_envelope_c(datatype, &nints, &nadds, &ncnts, &ntypes, &combiner);
#else
                int nints, nadds, ncnts = 0, ntypes;
                MPI_Type_get_envelope(datatype, &nints, &nadds, &ntypes, &combiner);
#endif
                ints = (int *) ADIOI_Malloc((nints + 1) * sizeof(int));
                adds = (MPI_Aint *) ADIOI_Malloc((nadds + 1) * sizeof(MPI_Aint));
                cnts = (MPI_Count *) ADIOI_Malloc((ncnts + 1) * sizeof(MPI_Count));
                types = (MPI_Datatype *) ADIOI_Malloc((ntypes + 1) * sizeof(MPI_Datatype));
#ifdef HAVE_MPI_LARGE_COUNT
                MPI_Type_get_contents_c(datatype, nints, nadds, ncnts, ntypes, ints, adds, cnts,
                                        types);
#else
                MPI_Type_get_contents(datatype, nints, nadds, ntypes, ints, adds, types);
#endif
                PNCIO_Datatype_iscontig(types[0], flag);

                PNCIO_Type_dispose(types);
                ADIOI_Free(ints);
                ADIOI_Free(adds);
                ADIOI_Free(cnts);
                ADIOI_Free(types);
            }
            break;
        case MPI_COMBINER_F90_INTEGER:
        case MPI_COMBINER_F90_REAL:
        case MPI_COMBINER_F90_COMPLEX:
            *flag = 1;
            break;
        default:
            *flag = 0;
            break;
    }

    /* This function needs more work. It should check for contiguity
     * in other cases as well. */
}
