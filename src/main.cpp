#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <zconf.h>
#include <memory>
#include <map>
#include <cassert>
#include <algorithm>
#include <mpi.h>

#include "zupply.hpp"

#define START_TIMING(v) \
auto v = MPI_Wtime()

#define RESTART_TIMING(v) \
v = MPI_Wtime() - v

#define STOP_TIMING(v) \
v = MPI_Wtime() - v

#define PAR_START_TIMING(v, comm) \
MPI_Barrier(comm); \
auto v = MPI_Wtime()

#define PAR_RESTART_TIMING(v, comm) \
MPI_Barrier(comm); \
v = MPI_Wtime() - v

#define PAR_STOP_TIMING(v, comm) \
MPI_Barrier(comm); \
v = MPI_Wtime() - v

inline long long position_to_cell(int msx, int msy, const std::pair<int, int> & position) {
    return position.first + msx * position.second;
}

inline long long position_to_cell(int msx, int msy, const int x, const int y) {
    return x + msx * y;
}

inline std::pair<int, int> cell_to_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msx);
}

struct Cell {
    int idx;
    int type;
};

#define WATER 1
#define ROCK  0

int main(int argc, char **argv) {
    int worldsize, rank;

    MPI_Init(&argc, &argv); 
    auto world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldsize);
    MPI_Comm_rank(world, &rank);

    int xprocs = std::atoi(argv[1]), yprocs = std::atoi(argv[2]);

    if(xprocs * yprocs != worldsize){
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int cell_per_process = std::atoi(argv[3]);

    int cell_per_row = std::sqrt(cell_per_process);
    int xcells = cell_per_row * xprocs, ycells = cell_per_row * yprocs; 
    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_position(xprocs, yprocs, rank);

    long total_cell = cell_per_process * worldsize;
    
    std::vector<Cell> my_cells; //local data
    
    


    
    return 0;
}
