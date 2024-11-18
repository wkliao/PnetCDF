#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "pnc_lustre.h"

#define ERR { \
    if (err != NC_NOERR) \
        printf("%d: line %d: error (%s)\n", rank, __LINE__, ncmpi_strerrno(err)); \
}

int main(int argc, char *argv[])
{
    char *filename="dummy", buf[100];
    int i, err, rank, nprocs;
    int len, psizes[2]={0,0}, gsizes[2], sizes[2], starts[2];
    MPI_Datatype ftype;
    MPI_Info info;
    MPI_Status status;
    PNC_File fh;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    len = 10;

    err = MPI_Dims_create(nprocs, 2, psizes);
    ERR

    /* set global array sizes, local array sizes, starting offsets */
    starts[0] = len * (rank / psizes[1]);
    starts[1] = len * (rank % psizes[1]);
    gsizes[0] = len * psizes[0];
    gsizes[1] = len * psizes[1];
    sizes[0]  = sizes[1] = len;

    printf("rank %3d: gsizes=%4d %4d sizes=%4d %4d starts=%4d %4d\n", rank,
           gsizes[0],gsizes[1],sizes[0],sizes[1],starts[0],starts[1]);

    /* create file type: 2D subarray */
    err = MPI_Type_create_subarray(2, gsizes, sizes, starts, MPI_ORDER_C, MPI_BYTE, &ftype);
    ERR
    err = MPI_Type_commit(&ftype);
    ERR

    err = MPI_Info_create(&info); ERR

    err = MPI_Info_set(info, "romio_no_indep_rw", "true"); ERR

    if (rank == 0) {
        err = PNC_File_delete(filename); ERR
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int is_lustre = PNC_Check_Lustre(filename);
    printf("File %s is Lustre? %s\n",filename,(is_lustre)?"yes":"no");

    /* create a new file collectively */
    err = PNC_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE,
                        info, &fh); ERR

    /* set the file view */
    err = PNC_File_set_view(fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL); ERR
    err = MPI_Type_free(&ftype); ERR
    err = MPI_Info_free(&info); ERR

    for (i=0; i<100; i++) buf[i] = 'A'+i%52;

    err = PNC_File_seek(fh, 0, MPI_SEEK_SET); ERR

    err = PNC_File_write_all(fh, buf, 10, MPI_CHAR, &status); ERR

    err = PNC_File_seek(fh, 0, MPI_SEEK_CUR); ERR

    err = PNC_File_write(fh, buf, 10, MPI_CHAR, &status); ERR

    err = PNC_File_seek(fh, 0, MPI_SEEK_SET); ERR

    err = PNC_File_write_at_all(fh, 0, buf, 10, MPI_CHAR, &status); ERR

    err = PNC_File_seek(fh, 0, MPI_SEEK_SET); ERR

    err = PNC_File_read(fh, buf, 10, MPI_CHAR, &status); ERR

    err = PNC_File_seek(fh, 0, MPI_SEEK_SET); ERR

    err = PNC_File_read_all(fh, buf, 10, MPI_CHAR, &status); ERR

    err = PNC_File_seek(fh, 0, MPI_SEEK_SET); ERR

    err = PNC_File_read_at_all(fh, 0, buf, 10, MPI_CHAR, &status); ERR

    err = PNC_File_sync(fh); ERR

    MPI_Offset file_size;
    err = PNC_File_get_size(fh, &file_size); ERR
    printf("file size = %lld\n", file_size);

    MPI_Offset new_file_size = 10;
    err = PNC_File_set_size(fh, new_file_size); ERR
    err = PNC_File_get_size(fh, &file_size); ERR
    printf("After setting file size to %lld file size = %lld\n", new_file_size, file_size);

    err = PNC_File_close(&fh); ERR

    MPI_Finalize();
    return 0;
}

