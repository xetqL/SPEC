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
#include <numeric>
#include <random>

#include "../include/cell.hpp"
#include "../include/utils.hpp"
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

#define WATER_TYPE 1
#define ROCK_TYPE 0

int main(int argc, char **argv) {
    int worldsize, rank;

    MPI_Init(&argc, &argv);
    auto world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldsize);
    MPI_Comm_rank(world, &rank);

    // Generate data
    int xprocs = std::atoi(argv[1]), yprocs = std::atoi(argv[2]);

    if(xprocs * yprocs != worldsize){
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int cell_per_process = std::atoi(argv[3]);
    int cell_per_row = (int) std::sqrt(cell_per_process), cell_per_col = cell_per_row;

    std::vector<int> cols(cell_per_col);
    std::iota(cols.begin(), cols.end(), 0);
    std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
    gen.seed(0);
    std::shuffle(cols.begin(), cols.end(), gen);
    std::vector<int> water_cols(cols.begin(), cols.begin() + (int) (cols.size() / 3));
    std::sort(water_cols.begin(), water_cols.end());
    std::cout << water_cols.at(1) << std::endl;
    int xcells = cell_per_row * xprocs, ycells = cell_per_row * yprocs;
    int& msx = Cell::get_msx(); msx = cell_per_row * worldsize;
    int& msy = Cell::get_msy(); msy = cell_per_row * worldsize;
    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_position(xprocs, yprocs, rank);
    long total_cell = cell_per_process * worldsize;
    auto datatype   = Cell::register_datatype();
    std::vector<Cell> my_cells; //local data
    my_cells.reserve(cell_per_process);
    for(int i = 0; i < cell_per_row; ++i) {
        for(int j = 0; j < cell_per_col; ++j) {
            my_cells.emplace_back(
                    ((cell_per_row * worldsize) * i) + cell_per_row + j,
                    j + i * cell_per_row,
                    std::find(water_cols.begin(), water_cols.end(), j) == water_cols.end() ?
            ROCK_TYPE : WATER_TYPE, 0);
        }
    }

    datatype.free_datatypes();
    MPI_Finalize();
    return 0;
}
