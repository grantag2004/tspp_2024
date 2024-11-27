#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random> 

#define GRID_SIZE 1024
#define MAX_ITERATIONS 1000

void initialize_random(std::vector<std::vector<int>> &grid, int rows, int cols, double alive_probability, int rank) {
    std::mt19937 rng(std::random_device{}() + rank); 
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            grid[i][j] = (dist(rng) < alive_probability) ? 1 : 0;
}

int count_neighbors(const std::vector<std::vector<int> > &grid, int x, int y, int rows, int cols) 
{
    int count = 0;

    for (int i = -1; i <= 1; ++i) 
        for (int j = -1; j <= 1; ++j) 
        {
            if (i == 0 && j == 0) 
                continue; 
            int changed_x = (x + i + rows) % rows; 
            int changed_y = (y + j + cols) % cols; 
            count += grid[changed_x][changed_y];
        }

    return count;
}

void update_grid(std::vector<std::vector<int> > &grid, int rows, int cols) 
{
    std::vector<std::vector<int> > changed_grid = grid;

    for (int i = 0; i < rows; ++i) 
        for (int j = 0; j < cols; ++j) 
        {
            int neighbors = count_neighbors(grid, i, j, rows, cols);
            if (grid[i][j] == 1) {
                changed_grid[i][j] = (neighbors == 2 || neighbors == 3) ? 1 : 0;
            } else {
                changed_grid[i][j] = (neighbors == 3) ? 1 : 0;
            }
        }
    

    grid = changed_grid;
}

int count_alive(const std::vector<std::vector<int> > &grid) 
{
    int count = 0;

    for (int i = 0; i < grid.size(); ++i) 
        for (int j = 0; j < grid[i].size(); ++j)
            count += grid[i][j];

    return count;
}

int main(int argc, char **argv) 
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows = GRID_SIZE / size;
    int remaining_rows = GRID_SIZE % size;
    int local_rows = rows + (rank == size - 1 ? remaining_rows : 0);

    std::vector<std::vector<int>> local_grid(local_rows + 2, std::vector<int>(GRID_SIZE, 0));
    double alive_probability = 0.3;
    initialize_random(local_grid, local_rows, GRID_SIZE, alive_probability, rank);


    int iteration = 0;
    int global_alive = 0;
    bool stop = false;
    int K = 400; // Минимальное количество итераций игры

    double start_time = MPI_Wtime();

    while ((iteration < MAX_ITERATIONS) && (iteration < K || !stop)) 
    {
        std::vector<int> send_up(GRID_SIZE), send_down(GRID_SIZE);
        std::vector<int> recv_up(GRID_SIZE), recv_down(GRID_SIZE);

        MPI_Request send_requests[2], recv_requests[2];
        MPI_Isend(local_grid[1].data(), GRID_SIZE, MPI_INT, (rank - 1 + size) % size, 0, MPI_COMM_WORLD, &send_requests[0]);
        MPI_Isend(local_grid[rows].data(), GRID_SIZE, MPI_INT, (rank + 1) % size, 1, MPI_COMM_WORLD, &send_requests[1]);

        MPI_Irecv(recv_up.data(), GRID_SIZE, MPI_INT, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, &recv_requests[0]);
        MPI_Irecv(recv_down.data(), GRID_SIZE, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD, &recv_requests[1]);

        MPI_Waitall(2, recv_requests, MPI_STATUSES_IGNORE);

        local_grid[0] = recv_up;
        local_grid[rows + 1] = recv_down;

        int local_alive = count_alive(local_grid);

        update_grid(local_grid, rows + 2, GRID_SIZE);

        int new_local_alive = count_alive(local_grid);

        stop = (local_alive == new_local_alive);
        int stop_global = 0;
        MPI_Allreduce(&stop, &stop_global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        stop = stop_global;

        int local_total = count_alive(local_grid);
        MPI_Allreduce(&local_total, &global_alive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        iteration++;
    }

    double end_time = MPI_Wtime();

    if (rank == 0) 
    {
        std::cout << "Игра завершена после " << iteration << " итераций.\n";
        std::cout << "Общее число живых клеток: " << global_alive << "\n";
        std::cout << "Общее время выполнения: " << (end_time - start_time) << " секунд.\n";
    }

    MPI_Finalize();
    return 0;
}
