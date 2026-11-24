#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* The editdist-mpi pair kernel exchanges two integers per row. Its following
   MPI_MAX reduction closes the timed region before the diagnostic gather. */
static double collective_seconds;
static int calls;

int MPI_Allgather(const void *send, int nsend, MPI_Datatype sendtype, void *recv,
                  int nrecv, MPI_Datatype recvtype, MPI_Comm comm) {
    const double start = PMPI_Wtime();
    int status = PMPI_Allgather(send, nsend, sendtype, recv, nrecv, recvtype, comm);
    if (nsend == 2 && nrecv == 2 && sendtype == MPI_INT && recvtype == MPI_INT) {
        collective_seconds += PMPI_Wtime() - start;
        ++calls;
    }
    return status;
}

int MPI_Reduce(const void *send, void *recv, int n, MPI_Datatype type, MPI_Op op,
               int root, MPI_Comm comm) {
    int status = PMPI_Reduce(send, recv, n, type, op, root, comm);
    if (comm != MPI_COMM_WORLD || n != 1 || type != MPI_DOUBLE || op != MPI_MAX || root != 0)
        return status;
    int rank, ranks;
    PMPI_Comm_rank(comm, &rank);
    PMPI_Comm_size(comm, &ranks);
    const double local[] = {*(const double *)send, collective_seconds, (double)calls};
    double *all = rank == root ? malloc((size_t)ranks * 3 * sizeof(double)) : NULL;
    if (rank == root && all == NULL) PMPI_Abort(comm, 2);
    PMPI_Gather(local, 3, MPI_DOUBLE, all, 3, MPI_DOUBLE, root, comm);
    if (rank == root) {
        int critical = 0;
        for (int r = 0; r < ranks; ++r) {
            if (all[3*r] > all[3*critical]) critical = r;
            if (all[3*r + 2] != calls) PMPI_Abort(comm, 3);
        }
        fprintf(stderr, "{\"critical_rank\":%d,\"calls\":%d,\"total\":%.12g,\"collective\":%.12g}\n",
                critical, calls, all[3*critical], all[3*critical+1]);
        free(all);
    }
    collective_seconds = 0;
    calls = 0;
    return status;
}
