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
#include <cnpy.h>
#include <complex>

#include "../include/cell.hpp"
#include "../include/utils.hpp"
#include "zupply.hpp"
#include "../include/zoltan_fn.hpp"

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

#define CHECKPOINT_TIMING(v, u) \
auto u = MPI_Wtime() - v


#define WATER_TYPE 1
#define ROCK_TYPE 0

int rank;

std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<float> ndist(90, 10);
std::uniform_real_distribution<float> udist(0, 1);


void dummy_erosion_computation(int msx, int msy,
                               std::vector<Cell>* _my_cells,
                               const std::vector<Cell>& remote_cells,
                               const std::vector<size_t>& data_pointers,
                               const std::tuple<int, int, int, int>& bbox) {


    std::vector<Cell>& my_cells = *(_my_cells);
    int x1,x2,y1,y2; std::tie(x1,x2,y1,y2) = bbox;
    int total_box = (x2-x1) * (y2-y1);
    const size_t my_cell_count = my_cells.size();
    std::vector<size_t> idx_neighbors(4);
    for(unsigned int i = 0; i < my_cell_count; ++i) {
        Cell& cell = my_cells[i];
        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
        if(cell.type) { // Do nothing if type == ROCK_TYPE
            std::this_thread::sleep_for(std::chrono::milliseconds( (int) std::max(ndist(gen), 0.0f)*0 ));
        } else {
            int _x,_y;
            std::tie(_x,_y) = cell_to_position((x2-x1), (y2-y1), lid);

            const Cell* neighbor_cell;

            if(lid+(x2-x1) < total_box) idx_neighbors[0] = (data_pointers[lid+(x2-x1)]);
            if(lid+1 < total_box)       idx_neighbors[1] = (data_pointers[lid+1]);
            if(lid-(x2-x1) >= 0)        idx_neighbors[2] = (data_pointers[lid-(x2-x1)]);
            if(lid-1 >= 0)              idx_neighbors[3] = (data_pointers[lid-1]);

            for(auto idx_neighbor : idx_neighbors){
                if(idx_neighbor >= my_cell_count) neighbor_cell = &remote_cells[idx_neighbor - my_cell_count];
                else neighbor_cell = &my_cells[idx_neighbor];
                if(neighbor_cell->type) {
                    auto p = udist(gen);
                    cell.type = p < cell.erosion_probability ? 1 : 0;
                    cell.weight = 1.0;
                    //break;
                }
            }
            std::fill(idx_neighbors.begin(), idx_neighbors.end(), 0);
        }
    }
}

void populate_data_pointers(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            const std::vector<Cell>& my_cells,
                            const std::vector<Cell>& remote_cells,
                            const std::tuple<int, int, int, int>& bbox) {
    std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1,x2,y1,y2) = bbox;
    int my_box = (x2-x1) * (y2-y1);
    data_pointers.resize(my_box);
    auto mine_size   = my_cells.size();
    auto remote_size = remote_cells.size();
    for (size_t i = 0; i < mine_size; ++i) {
        const Cell& cell = my_cells[i];
        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
        data_pointers[lid] = i;
    }
    for (size_t i = 0; i < remote_size; ++i) {
        const Cell& cell = remote_cells[i];
        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
        data_pointers[lid] = i+mine_size;
    }
}

template<class RealType, class Iter>
inline RealType skewness(Iter b, Iter e) {
    long N = std::distance(b, e);
    RealType mean = std::accumulate(b, e, (RealType) 0.0) / N;
    RealType diff1 = 0.0;
    std::for_each(b,e, [&diff1, mean](auto x){diff1 += std::pow(x - mean, 3.0);});
    diff1 /= N;
    RealType diff2 = 0.0;
    std::for_each(b,e, [&diff2, mean](auto x){diff2 += std::pow(x - mean, 2.0);});
    diff2 /= (N-1);
    return diff1 / std::pow(diff2, 3.0/2.0);
}

void print_bbox(std::tuple<int, int, int , int> bbox){
    std::cout << std::get<0>(bbox) << ", " << std::get<1>(bbox) << ", " << std::get<2>(bbox) << ", " << std::get<3>(bbox) << std::endl;
}

int main(int argc, char **argv) {
    int worldsize;

    MPI_Init(&argc, &argv);
    auto world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldsize);
    MPI_Comm_rank(world, &rank);
    std::string str_rank = "[RANK " + std::to_string(rank) + "] ";
    // Generate data
    int xprocs = std::atoi(argv[1]), yprocs = std::atoi(argv[2]);
    int cell_per_process = std::atoi(argv[3]);
    const int MAX_STEP = std::atoi(argv[4]);
    const int seed = std::atoi(argv[5]);
    zz::log::LoggerPtr perflogger, steplogger;

    zz::log::config_from_file("logger.cfg");
    perflogger = zz::log::get_logger("perf",  true);
    steplogger = zz::log::get_logger("steps", true);

    if(xprocs * yprocs != worldsize) {
        steplogger->fatal() << "Grid size does not match world size";
        MPI_Abort(world, 1);
        return EXIT_FAILURE;
    }

    if (xprocs < yprocs) {
        steplogger->fatal() << "Not a rectangular grid where X > Y";
        MPI_Abort(world, 1);
        return EXIT_FAILURE;
    }

    int cell_in_my_rows = (int) std::sqrt(cell_per_process), cell_in_my_cols = cell_in_my_rows;
    int xcells = cell_in_my_rows * xprocs, ycells = cell_in_my_rows * yprocs;
    std::vector<int> cols(ycells);
    std::iota(cols.begin(), cols.end(), 0);

    std::mt19937 gen(seed); ///predictable sequence of value

    std::shuffle(cols.begin(), cols.end(), gen);
    std::vector<int> water_cols(cols.begin(), cols.begin() + 1);
    std::sort(water_cols.begin(), water_cols.end());
    //std::for_each(water_cols.cbegin(), water_cols.cend(), [](auto v){ std::cout << v << std::endl; });

    int& msx = Cell::get_msx(); msx = xcells;
    int& msy = Cell::get_msy(); msy = ycells;
    int shape[2] = {msx,msy};
    if(!rank) cnpy::npz_save("out.npz", "shape", &shape[0], {2}, "w");

    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_global_position(xprocs, yprocs, rank);
    unsigned long total_cell = cell_per_process * worldsize;
    auto datatype   = Cell::register_datatype();

    if(!rank) {
        perflogger->info("CPU COUNT:")    << worldsize;
        perflogger->info("GRID PSIZE X:") << xprocs;
        perflogger->info("GRID PSIZE Y:") << yprocs;
    }

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);
    const int nb_bucket = 8;
    //gen.seed(1);
    if(!rank) steplogger->info() << cell_in_my_cols << " " << cell_in_my_rows;

    float minp, maxp;



    for(int j = 0; j < cell_in_my_cols; ++j) {

        if(j % 100 == 0){
            std::uniform_real_distribution<float> real_dist1(0.0f, 0.1f);
            minp = real_dist1(gen);
            std::uniform_real_distribution<float> real_dist2(0.2f, 0.5f);
            maxp = real_dist2(gen);
        }

        std::uniform_real_distribution<float> real_dist(minp, maxp);

        //min_bucket * 1/(1.5f*nb_bucket), max_bucket * 1/(1.5f*nb_bucket));

        for(int i = 0; i < cell_in_my_rows; ++i) {
            bool rock = std::find(water_cols.begin(), water_cols.end(), (i + (y_proc_idx * cell_in_my_cols))) == water_cols.end();
            int id = cell_in_my_rows * x_proc_idx + j + msx * (i + (y_proc_idx * cell_in_my_cols));
            my_cells.emplace_back(id, rock ? ROCK_TYPE : WATER_TYPE, rock ? 0.0f : 1.0f, rock ? real_dist(gen) : 0.0);
        }
    }

    const int my_cell_count = my_cells.size();
    std::vector<std::array<int,2>> all_types(total_cell);
    std::vector<std::array<int,2>> my_types(my_cell_count);
    for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};

    gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);

    if(!rank) {
        std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
        std::vector<int> types;
        std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
        assert((*(all_types.end() - 1))[0] == total_cell-1);
        cnpy::npz_save("out.npz", "step-"+std::to_string(0), &types[0], {total_cell}, "a");
    }

    int recv, sent;
    if(!rank) steplogger->info() << "End of map generation";
    /* Initial load balancing */
    auto zoltan_lb = zoltan_create_wrapper(true, world);
    zoltan_load_balance(&my_cells, zoltan_lb, true, true);

    /* lets make it fun now...*/
    int minx, maxx, miny, maxy;
    std::vector<size_t> data_pointers;
    if(!rank) all_types.resize(total_cell);
    PAR_START_TIMING(loop_time, world);
    for(unsigned int step = 0; step < MAX_STEP; ++step) {

        const int my_cell_count = my_cells.size();
        if(!rank) steplogger->info() << "Beginning step "<< step;
        PAR_START_TIMING(step_time, world);
        // Share my cells with my neighbors to propagate their state

        auto remote_cells = zoltan_exchange_data(zoltan_lb, my_cells, &recv, &sent, datatype.element_datatype, world, 1.0);

        auto bbox = get_bounding_box(my_cells, remote_cells);
        //print_bbox(bbox);
        populate_data_pointers(msx, msy, &data_pointers, my_cells, remote_cells, bbox);
        dummy_erosion_computation(msx, msy, &my_cells, remote_cells, data_pointers, bbox);

        CHECKPOINT_TIMING(step_time, my_step_time);
        PAR_STOP_TIMING(step_time, world);
        if(!rank) steplogger->info() << "Stop step "<< step;

        std::vector<double> timings(worldsize);
        MPI_Gather(&my_step_time, 1, MPI_DOUBLE, &timings.front(), 1, MPI_DOUBLE, 0, world);

        double max = *std::max_element(timings.cbegin(), timings.cend()),
               average = std::accumulate(timings.cbegin(), timings.cend(), 0.0) / worldsize,
               load_imbalance = (max / average - 1.0) * 100.0,
               skew = skewness<double>(timings.cbegin(), timings.cend());

        if(!rank) steplogger->info("stats: ") << "load imbalance: " << load_imbalance << " skewness: " << skew;
        std::vector<std::array<int,2>> my_types(my_cell_count);
        for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
        gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
        if(!rank){
            std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
            std::vector<int> types;
            std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
            cnpy::npz_save("out.npz", "step-"+std::to_string(step+1), &types[0], {total_cell}, "a");
        }
    }

    datatype.free_datatypes();
    MPI_Finalize();
    return 0;
}
