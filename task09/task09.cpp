#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <cmath>
#include <random>
#include <chrono>

void generateRandomMatrix(float* matrix, int size, int rank) 
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + rank;
    std::mt19937 generator(seed);  
    std::uniform_real_distribution<float> distribution(0.0, 10.0);  

    for (int i = 0; i < size; ++i) 
        matrix[i] = distribution(generator);
}

void sequentialMatrixMultiplication(const std::vector<float>& A, const std::vector<float>& B, std::vector<float>& C, int N)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                C[i * N + j] += A[i * N + k] * B[k * N + j];
}

bool compareMatrices(const std::vector<float>& C1, const std::vector<float>& C2, int size, float epsilon = 1e-5)
{
    for (int i = 0; i < size; ++i)
        if (std::fabs(C1[i] - C2[i]) > epsilon)
            return false;
    return true;
}

#define N 300

int main(int argc, char **argv) 
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int p = std::sqrt(size);

    if (p * p != size) 
    {
        if (rank == 0) 
            std::cerr << "Ошибка: число процессов должно быть квадратом.\n";
        MPI_Finalize();
        return 0;
    }

    int N_p = N / p;

    int block_size;

    if (rank == 0) 
    {
        if (argc < 2) 
        {
            std::cerr << "Ошибка: укажите размер блока в аргументах.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        block_size = std::atoi(argv[1]);
    }

    MPI_Bcast(&block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (N_p % block_size != 0) 
    {
        if (rank == 0) 
            std::cerr << "Ошибка: размер блока должен быть делителем N_p.\n";
        MPI_Finalize();
        return 0;
    }

    int num_blocks = N_p / block_size;

    std::vector<float> A(N_p * N_p);
    std::vector<float> B(N_p * N_p);
    std::vector<float> C(N_p * N_p, 0.0f);

    generateRandomMatrix(A.data(), N_p * N_p, rank);
    generateRandomMatrix(B.data(), N_p * N_p, rank);

    MPI_Comm comm, comm_row, comm_col;
    int dim[2] = {p, p}, period[2] = {0, 0}, reorder = 0;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

    int n_rank, coord[2];
    MPI_Comm_rank(comm, &n_rank);
    MPI_Cart_coords(comm, n_rank, 2, coord);

    int remain_dims[2] = {0, 1};
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    remain_dims[0] = 1; remain_dims[1] = 0;
    MPI_Cart_sub(comm, remain_dims, &comm_col);

    int row_rank, col_rank;
    MPI_Comm_rank(comm_row, &row_rank);
    MPI_Comm_rank(comm_col, &col_rank);

    std::vector<float> sequential_C;
    if (rank == 0) 
    {
        std::vector<float> full_A(size * N_p * N_p, 0.0f);
        std::vector<float> full_B(size * N_p * N_p, 0.0f);
        std::vector<float> ordered_A(N * N, 0.0f);
        std::vector<float> ordered_B(N * N, 0.0f);
        sequential_C.resize(N * N, 0.0f);

        MPI_Gather(A.data(), N_p * N_p, MPI_FLOAT, full_A.data(), N_p * N_p, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Gather(B.data(), N_p * N_p, MPI_FLOAT, full_B.data(), N_p * N_p, MPI_FLOAT, 0, MPI_COMM_WORLD);

        for (int proc = 0; proc < size; ++proc) 
        {
            int coords[2];
            MPI_Cart_coords(comm, proc, 2, coords);

            int row_start = coords[0] * N_p;
            int col_start = coords[1] * N_p;

            for (int i = 0; i < N_p; ++i)
                for (int j = 0; j < N_p; ++j) 
                {
                    ordered_A[(row_start + i) * N + (col_start + j)] = full_A[proc * N_p * N_p + i * N_p + j];
                    ordered_B[(row_start + i) * N + (col_start + j)] = full_B[proc * N_p * N_p + i * N_p + j];
                }
        }

        sequentialMatrixMultiplication(ordered_A, ordered_B, sequential_C, N);
    }
    else 
    {
        MPI_Gather(A.data(), N_p * N_p, MPI_FLOAT, nullptr, 0, MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Gather(B.data(), N_p * N_p, MPI_FLOAT, nullptr, 0, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }


    double start_time = MPI_Wtime();

    for (int k = 0; k < p; ++k) 
    {
        int rank_A = k, rank_B = k;

        for (int b = 0; b < num_blocks; ++b) 
        {
            std::vector<float> A_block(N_p * block_size);
            std::vector<float> B_block(N_p * block_size);

            if (row_rank == rank_A) 
                for (int i = 0; i < N_p; ++i) 
                    for (int j = 0; j < block_size; ++j) 
                        A_block[i * block_size + j] = A[i * N_p + j + b * block_size];

            if (col_rank == rank_B) 
                for (int i = 0; i < block_size; ++i) 
                    for (int j = 0; j < N_p; ++j) 
                        B_block[i * N_p + j] = B[(i + b * block_size) * N_p + j];

            MPI_Bcast(A_block.data(), N_p * block_size, MPI_FLOAT, rank_A, comm_row);
            MPI_Bcast(B_block.data(), block_size * N_p, MPI_FLOAT, rank_B, comm_col);

            for (int i = 0; i < N_p; ++i) 
                for (int j = 0; j < N_p; ++j) 
                    for (int z = 0; z < block_size; ++z) 
                        C[i * N_p + j] += A_block[i * block_size + z] * B_block[z * N_p + j];
        }
    }

    double end_time = MPI_Wtime();

    if (rank == 0) 
    {
        std::vector<float> full_C(size * N_p * N_p, 0.0f); 
        MPI_Gather(C.data(), N_p * N_p, MPI_FLOAT, full_C.data(), N_p * N_p, MPI_FLOAT, 0, MPI_COMM_WORLD);

        std::vector<float> ordered_C(N * N, 0.0f); 
        for (int proc = 0; proc < size; ++proc) 
        {
            int coords[2];
            MPI_Cart_coords(comm, proc, 2, coords);

            int row_start = coords[0] * N_p; 
            int col_start = coords[1] * N_p; 

            for (int i = 0; i < N_p; ++i)
                for (int j = 0; j < N_p; ++j)
                    ordered_C[(row_start + i) * N + (col_start + j)] = full_C[proc * N_p * N_p + i * N_p + j];
        }

        std::cout << "Время выполнения SUMMA: " << (end_time - start_time) << " секунд\n";
        bool correct = compareMatrices(ordered_C, sequential_C, N * N);
        if (correct)
            std::cout << "Результаты совпадают! Параллельное умножение выполнено корректно.\n";
        else
            std::cerr << "Ошибка: результаты не совпадают!\n";
    }
    else 
    {
        MPI_Gather(C.data(), N_p * N_p, MPI_FLOAT, nullptr, 0, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
