#include <mpi.h>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>

#define GRID_SIZE 64
#define ITERATIONS 10

void generate_random(std::vector<double>& matrix, int rank) 
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + rank;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (auto& element : matrix) 
        element = distribution(generator);
}


inline int idx(int x, int y, int z, int xx, int yy) 
{
    return z * xx * yy + y * xx + x;
}


void update_jacobi(const std::vector<double>& current, std::vector<double>& next, int xx, int yy, int zz)
{
    for (int z = 1; z < zz - 1; ++z) 
        for (int y = 1; y < yy - 1; ++y) 
            for (int x = 1; x < xx - 1; ++x) {
                int i = idx(x, y, z, xx, yy);
                double x_sum = current[i + 1] + current[i - 1];
                double y_sum = current[i + xx] + current[i - xx];
                double z_sum = current[i + xx*yy] + current[i - xx*yy];
                next[i] = (x_sum + y_sum + z_sum) / 6.0;
            }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[3] = {0, 0, 0};
    int periods[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims);   

    for (int i = 0; i < 3; ++i) 
        if (GRID_SIZE % dims[i] != 0) 
            if (rank == 0) 
                std::cerr << "[Warning] GRID_SIZE=" << GRID_SIZE
                          << " не делится нацело на dims[" << i 
                          << "]=" << dims[i] << std::endl;

    MPI_Comm cart_comm;
    int reorder = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &cart_comm);

    int cart_rank;
    int coords[3]; 
    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, 3, coords);

    int Nx = GRID_SIZE / dims[0];
    int Ny = GRID_SIZE / dims[1];
    int Nz = GRID_SIZE / dims[2];

    int xx = Nx + 2;
    int yy = Ny + 2;
    int zz = Nz + 2;

    std::vector<double> current_grid(xx * yy * zz, 0.0);
    std::vector<double> next_grid(xx * yy * zz, 0.0);

    generate_random(next_grid, rank);

    MPI_Datatype type_X; 
    MPI_Type_vector(yy * zz, 1, xx, MPI_DOUBLE, &type_X);
    MPI_Type_commit(&type_X);

    MPI_Datatype type_Y; 
    MPI_Type_vector(zz, xx, xx * yy, MPI_DOUBLE, &type_Y);
    MPI_Type_commit(&type_Y);

    MPI_Datatype type_Z; 
    MPI_Type_contiguous(xx * yy, MPI_DOUBLE, &type_Z);
    MPI_Type_commit(&type_Z);

    MPI_Datatype layer_types[3] = {type_X, type_Y, type_Z};

    auto get_neighbor = [&](int dim, int direction) 
	{
        int rank_source, rank_dest;
        MPI_Cart_shift(cart_comm, dim, 1, &rank_source, &rank_dest);
        if (direction == -1) return rank_source; 
        else if (direction == +1) return rank_dest;
    };

    //Определение смещений для отправки/приёма
    int offset_send_x_left = idx(1, 0, 0, xx, yy);  
    int offset_recv_x_left = idx(0, 0, 0, xx, yy);  

    int offset_send_x_right = idx(xx - 2, 0, 0, xx, yy);  
    int offset_recv_x_right = idx(xx - 1, 0, 0, xx, yy); 


    int offset_send_y_down = idx(0, 1, 0, xx, yy);  
    int offset_recv_y_down = idx(0, 0, 0, xx, yy); 

    int offset_send_y_up = idx(0, yy - 2, 0, xx, yy);  
    int offset_recv_y_up = idx(0, yy - 1, 0, xx, yy);  


    int offset_send_z_front = idx(0, 0, 1, xx, yy);  
    int offset_recv_z_front = idx(0, 0, 0, xx, yy);  

    int offset_send_z_back = idx(0, 0, zz - 2, xx, yy);  
    int offset_recv_z_back = idx(0, 0, zz - 1, xx, yy);  

    double t_start = MPI_Wtime();

    for (int iter = 0; iter < ITERATIONS; ++iter)
    {
        std::copy(next_grid.begin(), next_grid.end(), current_grid.begin());

        MPI_Request req[12];
        int rcount = 0;

        {
            int left_rank  = get_neighbor(0, -1);
            int right_rank = get_neighbor(0, +1);

            if (left_rank != MPI_PROC_NULL) 
			{
                MPI_Isend(&current_grid[offset_send_x_left], 1, layer_types[0],
                          left_rank,  101, cart_comm, &req[rcount++]);
                MPI_Irecv(&current_grid[offset_recv_x_left], 1, layer_types[0],
                          left_rank,  102, cart_comm, &req[rcount++]);
            }
            if (right_rank != MPI_PROC_NULL) 
			{
                MPI_Isend(&current_grid[offset_send_x_right], 1, layer_types[0],
                          right_rank, 201, cart_comm, &req[rcount++]);
                MPI_Irecv(&current_grid[offset_recv_x_right], 1, layer_types[0],
                          right_rank, 202, cart_comm, &req[rcount++]);
            }
        }

        {
            int down_rank = get_neighbor(1, -1);
            int up_rank   = get_neighbor(1, +1);

            if (down_rank != MPI_PROC_NULL) 
			{
                MPI_Isend(&current_grid[offset_send_y_down], 1, layer_types[1],
                          down_rank,  111, cart_comm, &req[rcount++]);
                MPI_Irecv(&current_grid[offset_recv_y_down], 1, layer_types[1],
                          down_rank,  112, cart_comm, &req[rcount++]);
            }
            if (up_rank != MPI_PROC_NULL) 
			{
                MPI_Isend(&current_grid[offset_send_y_up], 1, layer_types[1],
                          up_rank,    211, cart_comm, &req[rcount++]);
                MPI_Irecv(&current_grid[offset_recv_y_up], 1, layer_types[1],
                          up_rank,    212, cart_comm, &req[rcount++]);
            }
        }

        {
            int front_rank = get_neighbor(2, -1);
            int back_rank  = get_neighbor(2, +1);

            if (front_rank != MPI_PROC_NULL) 
			{
                MPI_Isend(&current_grid[offset_send_z_front], 1, layer_types[2],
                          front_rank,  121, cart_comm, &req[rcount++]);
                MPI_Irecv(&current_grid[offset_recv_z_front], 1, layer_types[2],
                          front_rank,  122, cart_comm, &req[rcount++]);
            }
            if (back_rank != MPI_PROC_NULL) 
			{
                MPI_Isend(&current_grid[offset_send_z_back], 1, layer_types[2],
                          back_rank,   221, cart_comm, &req[rcount++]);
                MPI_Irecv(&current_grid[offset_recv_z_back], 1, layer_types[2],
                          back_rank,   222, cart_comm, &req[rcount++]);
            }
        }

        MPI_Waitall(rcount, req, MPI_STATUS_IGNORE);

        update_jacobi(current_grid, next_grid, xx, yy, zz);
    }

    double t_end = MPI_Wtime();

    double local_norm = 0.0;
    for (size_t i = 0; i < current_grid.size(); ++i) 
        local_norm += std::fabs(current_grid[i] - next_grid[i]);

    local_norm /= (double)current_grid.size();


    double global_norm = 0.0;
    MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);
    global_norm /= size;

    if (cart_rank == 0) 
	{
        std::cout << "Общая усреднённая норма = " << global_norm << std::endl;
        std::cout << "Время выполнения: " << (t_end - t_start) << " секунд" << std::endl;
    }

    MPI_Type_free(&type_X);
    MPI_Type_free(&type_Y);
    MPI_Type_free(&type_Z);

    MPI_Finalize();
    return 0;
}
