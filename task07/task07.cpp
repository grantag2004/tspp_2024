#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random> 

const int N = std::pow(2, 10);        
const int iterations = 1000;  

void initialization(std::vector<std::vector<double>>& table, int rank) 
{
    std::mt19937 gen(rank);
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (auto& row : table) 
        for (auto& val : row)
            val = dis(gen);  
}

void jacobi(const std::vector<std::vector<double>>& local_grid, std::vector<std::vector<double>>& changed_local_grid, int rows, int cols) 
{
    for (int i = 1; i < rows - 1; ++i) 
        for (int j = 1; j < cols - 1; ++j) 
            changed_local_grid[i][j] = 0.25 * (local_grid[i+1][j] + local_grid[i-1][j] + local_grid[i][j+1] + local_grid[i][j-1]);
}

double modul(const std::vector<std::vector<double>>& local_grid, const std::vector<std::vector<double>>& changed_local_grid, int rows, int cols) 
{
    double modul = 0.0;
    for (int i = 1; i < rows - 1; ++i) 
        for (int j = 1; j < cols - 1; ++j) 
            modul += std::pow(changed_local_grid[i][j] - local_grid[i][j], 2);
        
    return std::sqrt(modul);
}

int main(int argc, char** argv) 
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;

    int cols = N;
    int rows_count = N / size;

    std::vector<std::vector<double>> local_grid(rows_count + 2, std::vector<double>(cols, 0.0));
    std::vector<std::vector<double>> changed_local_grid(rows_count + 2, std::vector<double>(cols, 0.0));

    initialization(local_grid, rank);

    double time1 = MPI_Wtime();

    for (int iter = 0; iter < iterations; ++iter) 
    {
        if (rank > 0) 
        {
            MPI_Send(local_grid[1].data(), cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(local_grid[0].data(), cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
        }

        if (rank < size - 1) 
        {
            MPI_Send(local_grid[rows_count].data(), cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(local_grid[rows_count + 1].data(), cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
        }

        jacobi(local_grid, changed_local_grid, rows_count + 2, cols);

        if (iter == iterations - 1) 
        {
            double local_modul = modul(local_grid, changed_local_grid, rows_count + 2, cols);
            double global_modul;
            MPI_Reduce(&local_modul, &global_modul, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (rank == 0)
                std::cout << "Norm: " << std::sqrt(global_modul) << std::endl;
        }

        local_grid.swap(changed_local_grid);
    }

    double time2 = MPI_Wtime();

    if (rank == 0) 
        std::cout << "Execution time: " << time2 - time1 << " seconds" << std::endl;
    
    MPI_Finalize();
    
    return 0;
}
