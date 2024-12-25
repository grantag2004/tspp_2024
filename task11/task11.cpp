#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

std::vector<double> generate_vector(int N) 
{
    std::vector<double> b(N, 1.0); 
    for(int i = 0; i < N; ++i) 
        b[i] = static_cast<double>(i + 1); 
    
    return b;
}

std::vector<double> generate_matrix_block(int rows, int N, int local_start_row) 
{
    std::vector<double> A_local(rows * N, 1.0); 
    for(int i = 0; i < rows; ++i) 
        for(int j = 0; j < N; ++j) 
            A_local[i * N + j] = static_cast<double>((local_start_row + i) * N + j + 1);
    return A_local;
}

int main(int argc, char *argv[]) 
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N = 256; 

    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims); 


    if (size > 2) 
        assert(dims[0] > 1 && dims[1] > 1);

    int periods[2] = {0, 0}; 
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    int grid_rows = dims[0];
    int grid_cols = dims[1];
    int base_rows_per_grid_row = N / grid_rows;
    int remainder_rows = N % grid_rows;

    int rows_for_this_grid_row = base_rows_per_grid_row + (coords[0] < remainder_rows ? 1 : 0);
    int start_row = coords[0] * base_rows_per_grid_row + std::min(coords[0], remainder_rows);

    int rows_per_proc = rows_for_this_grid_row / grid_cols;
    int remainder_proc_rows = rows_for_this_grid_row % grid_cols;

    int local_rows = rows_per_proc + (coords[1] < remainder_proc_rows ? 1 : 0);
    int local_start_row = start_row + coords[1] * rows_per_proc + std::min(coords[1], remainder_proc_rows);

    std::vector<double> A_local = generate_matrix_block(local_rows, N, local_start_row);

    std::vector<double> b;
    if(rank == 0) 
        b = generate_vector(N);

    double *b_window = nullptr;
    MPI_Win win_b;
    if(rank == 0) 
    {
        b_window = new double[N];
        for(int i = 0; i < N; ++i) 
            b_window[i] = b[i];
    }

    MPI_Win_create(b_window, (rank == 0) ? sizeof(double) * N : 0, sizeof(double), MPI_INFO_NULL, cart_comm, &win_b);

    MPI_Win_fence(0, win_b);

    std::vector<double> b_local(N, 0.0);
    if(rank != 0) 
        MPI_Get(b_local.data(), N, MPI_DOUBLE, 0, 0, N, MPI_DOUBLE, win_b);
    else 
        b_local = b; 
    
    MPI_Win_fence(0, win_b);

    MPI_Win_free(&win_b);
    if(rank == 0) 
        delete[] b_window;

    MPI_Barrier(cart_comm);

    double start_time = MPI_Wtime();

    std::vector<double> c_local(local_rows, 0.0);
    for(int i = 0; i < local_rows; ++i)
        for(int j = 0; j < N; ++j) 
            c_local[i] += A_local[i * N + j] * b_local[j];
    
    MPI_Barrier(cart_comm);
    double end_time = MPI_Wtime();
    double local_elapsed = end_time - start_time;
    double max_elapsed;
    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);

    std::vector<int> recvcounts(size, 0);
    std::vector<int> displs(size, 0);

    for(int p = 0; p < size; ++p)
    {
        int p_coords[2];
        MPI_Cart_coords(cart_comm, p, 2, p_coords);
        int p_grid_row = p_coords[0];
        int p_grid_col = p_coords[1];
        int p_base_rows = base_rows_per_grid_row;
        int p_remainder = (p_grid_row < remainder_rows) ? 1 : 0;
        int p_rows_for_grid_row = p_base_rows + p_remainder;
        int p_rows_per_proc = p_rows_for_grid_row / grid_cols;
        int p_remainder_proc_rows = p_rows_for_grid_row % grid_cols;
        int p_local_rows = p_rows_per_proc + (p_grid_col < p_remainder_proc_rows ? 1 : 0);
        recvcounts[p] = p_local_rows;
    }

    displs[0] = 0;
    for(int p = 1; p < size; ++p)
        displs[p] = displs[p-1] + recvcounts[p-1];

    std::vector<double> c;
    if(rank == 0)
        c.resize(N, 0.0);

    MPI_Gatherv(c_local.data(), local_rows, MPI_DOUBLE, c.data(), recvcounts.data(), displs.data(), MPI_DOUBLE, 0, cart_comm);

    if(rank == 0)
        std::cout << "Максимальное время выполнения: " << max_elapsed << " секунд.\n";

    MPI_Comm_free(&cart_comm);

    MPI_Finalize();
    return 0;
}
