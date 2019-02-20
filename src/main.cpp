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
#include <complex>

#ifdef PRODUCE_OUTPUTS
    #include <cnpy.h>
#endif

#include "../include/cell.hpp"
#include "../include/utils.hpp"
#include "zupply.hpp"
#include "../include/zoltan_fn.hpp"
#include "../include/window.hpp"

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
#define EMPTY_TYPE -1

int rank;
zz::log::LoggerPtr perflogger, steplogger;

std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<float> ndist(9, 1);
std::normal_distribution<float> flopdist(145-52, 5); // https://arxiv.org/pdf/1703.08015.pdf
std::uniform_real_distribution<float> udist(0, 1);

std::vector<Cell> dummy_erosion_computation(int msx, int msy,
                               const std::vector<Cell>& _my_cells,
                               const std::vector<Cell>& remote_cells,
                               const std::vector<size_t>& data_pointers,
                               const std::tuple<int, int, int, int>& bbox) {


    std::vector<Cell> my_cells = _my_cells;
    const std::vector<Cell>& my_old_cells = _my_cells;
    int x1,x2,y1,y2; std::tie(x1,x2, y1,y2) = bbox;
    int total_box = (x2-x1) * (y2-y1);
    const size_t my_cell_count = my_cells.size();
    const size_t remote_count  = remote_cells.size();
    const size_t all_count     = my_cell_count + remote_count;
    std::vector<size_t> idx_neighbors(3, -1);
    for(unsigned int i = 0; i < all_count; ++i) {
        const Cell* cell;
        if(i < my_cell_count) cell = &my_old_cells[i];
        else cell = &remote_cells[i - my_cell_count];
        if(cell->type) { // Do nothing if type == ROCK_TYPE
            auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell->gid));
            if(cell->type) {
                std::fill(idx_neighbors.begin(), idx_neighbors.end(), -1);
                if(lid+(x2-x1)+1 < total_box) idx_neighbors[0] = (data_pointers[lid+(x2-x1)+1]);
                if(lid+1 < total_box)         idx_neighbors[1] = (data_pointers[lid+1]);
                if( ((lid-(x2-x1))+1) >= 0 )  idx_neighbors[2] = data_pointers[((lid-(x2-x1))+1)];
                for (auto idx_neighbor : idx_neighbors) {
                    if(idx_neighbor > my_cell_count) continue;
                    if(my_cells[idx_neighbor].type) continue;
                    auto erosion_proba = my_old_cells[idx_neighbor].erosion_probability;
                    auto p = udist(gen);
                    if (p < erosion_proba) {
                        my_cells[idx_neighbor].type   = 1;
                        my_cells[idx_neighbor].weight = 1;
                    }
                }
            }
        }
    }
    return my_cells;
}


std::vector<Cell> dummy_erosion_computation2(int msx, int msy,
                                            const std::vector<Cell>& _my_cells,
                                            const std::vector<Cell>& remote_cells,
                                            const std::vector<size_t>& data_pointers,
                                            const std::tuple<int, int, int, int>& bbox) {


    std::vector<Cell> my_cells = _my_cells;
    const std::vector<Cell>& my_old_cells = _my_cells;

    int x1,x2,y1,y2; std::tie(x1,x2, y1,y2) = bbox;
    int total_box = (x2-x1) * (y2-y1);
    const size_t my_cell_count = my_cells.size();
    const size_t remote_count  = remote_cells.size();
    const size_t all_count     = my_cell_count + remote_count;

    std::vector<size_t> idx_neighbors(8, -1);
    std::vector<float>  thetas(8, 0);

    for(unsigned int i = 0; i < all_count; ++i) {
        const Cell* cell;
        if(i < my_cell_count) cell = &my_old_cells[i];
        else cell = &remote_cells[i - my_cell_count];

        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell->gid));
        if(cell->type) { // Do nothing if type == ROCK_TYPE
            auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell->gid));
            if(cell->type) {
                std::fill(idx_neighbors.begin(), idx_neighbors.end(), -1);
                if(lid+1 < total_box) {
                    idx_neighbors[0] = (data_pointers[lid+1]);
                    thetas[0]        = 1.0f;
                }
                if((lid-(x2-x1))+1 >= 0) {
                    idx_neighbors[1] = (data_pointers[(lid-(x2-x1))+1]);
                    thetas[1]        = 1.0f/1.4142135f;
                }
                if(lid-(x2-x1) >= 0) {
                    idx_neighbors[2] = (data_pointers[lid-(x2-x1)]);
                    thetas[2]        = 0;
                }
                if(lid+(x2-x1) < total_box) {
                    idx_neighbors[6] = (data_pointers[lid+(x2-x1)]);
                    thetas[6]        = 0;
                }
                if(lid+(x2-x1)+1 < total_box) {
                    idx_neighbors[7] = (data_pointers[lid+(x2-x1)+1]);
                    thetas[7]        = 1.0f/1.4142135f;
                }
                for (int j = 0; j < 8; ++j) {
                    auto idx_neighbor = idx_neighbors[j];
                    auto theta = thetas[j];
                    if(idx_neighbor >= my_cell_count) continue;
                    if(my_cells[idx_neighbor].type) continue;
                    auto erosion_proba = my_old_cells[idx_neighbor].erosion_probability;
                    auto p = udist(gen);
                    if (p < (theta) * erosion_proba) {
                        my_cells[idx_neighbor].type   = 1;
                        my_cells[idx_neighbor].weight = 1;
                    }
                }

                /*DO NOT OPTIMIZE; SIMULATE COMPUTATION OF LBM FLUID WITH BGK D2Q9*/
                float flop = flopdist(gen);
                const int nb_flop = (int) flop;
                float res = 0.0;
                for(int i = 0; i < nb_flop; ++i) {
                    __asm__(""); res = flop * flop;
                }
                thetas[0] = res;
            }
        }
    }
    return my_cells;
}

/*
std::vector<Cell> invasion_percolation(int msx, int msy, float p,
                                    const std::vector<Cell>& _my_cells,
                                    const std::vector<Cell>& remote_cells,
                                    const std::vector<size_t>& data_pointers,
                                    const std::tuple<int, int, int, int>& bbox) {

    std::vector<Cell> my_cells = _my_cells;
    const size_t my_cell_count = my_cells.size();
    const std::vector<Cell>& my_old_cells = _my_cells;

    int x1,x2,y1,y2; std::tie(x1,x2, y1,y2) = bbox;
    int total_box = (x2-x1) * (y2-y1);
    const size_t remote_count  = remote_cells.size();


    const size_t all_count     = my_cell_count + remote_count;

    std::vector<size_t> idx_neighbors(4, -1);
    for(unsigned int i = 0; i < all_count; ++i) {
        const Cell* cell;

        if(i < my_cell_count) cell = &my_old_cells[i];
        else cell = &remote_cells[i - my_cell_count];

        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell->gid));
        if(cell->type) {
            int _x,_y; std::tie(_x,_y) = cell_to_position((x2-x1), (y2-y1), lid);
            const Cell* neighbor_cell;
            std::fill(idx_neighbors.begin(), idx_neighbors.end(), -1);

            if(lid+(x2-x1) < total_box) idx_neighbors[0] = (data_pointers[lid+(x2-x1)]);
            if(lid+1 < total_box)       idx_neighbors[1] = (data_pointers[lid+1]);
            if(lid-(x2-x1) >= 0)        idx_neighbors[2] = (data_pointers[lid-(x2-x1)]);
            if(lid-1 >= 0)              idx_neighbors[3] = (data_pointers[lid-1]);

            std::vector<std::pair<float, size_t>> neighbor_probability;
            for (unsigned long idx_neighbor : idx_neighbors) {
                if(idx_neighbor > (my_cell_count + remote_count)) continue;
                if(idx_neighbor >= my_cell_count) continue;
                if(my_old_cells[idx_neighbor].type) continue;

                auto gid = my_old_cells[idx_neighbor].gid;
                auto pos = cell_to_global_position(msx,msy,gid);
                float rel_pos = 1.0-((float) msx - (float) (pos.first+1)) / (float) (msx);
                //std::cout << udist(gen) * p + (1-p) * rel_pos << std::endl;

                neighbor_probability.emplace_back(udist(gen) * p + (1-p) * rel_pos, idx_neighbor);
            }
            std::sort(neighbor_probability.begin(),
                      neighbor_probability.end(), [](auto a, auto b){return a.first>b.first;});

            int idx_invasion_neighbor = -1;
            if(!neighbor_probability.empty()){
                idx_invasion_neighbor = neighbor_probability.front().second;
                my_cells[idx_invasion_neighbor].type   = 1;
                my_cells[idx_invasion_neighbor].weight = 1;
            }
        }
    }
    return my_cells;
}
*/

void populate_data_pointers(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            const std::vector<Cell>& my_cells,
                            const std::vector<Cell>& remote_cells,
                            const std::tuple<int, int, int, int>& bbox) {
    std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1,x2, y1,y2) = bbox;
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
    //if(!rank) std::cout << mine_size << std::endl;
    //if(!rank) std::for_each(data_pointers.cbegin(), data_pointers.cend(), [&](auto v){ std::cout << mine_size << ";" << v << ";" << my_cells[v].gid << std::endl; });
}
template<class RealType, class Iter>
inline RealType mean(Iter b, Iter e) {
    long N = std::distance(b, e);
    return std::accumulate(b, e, (RealType) 0.0) / N;
}

template<class RealType, class Iter>
inline RealType skewness(Iter b, Iter e) {
    long N = std::distance(b, e);
    RealType _mean = mean<RealType, Iter>(b,e);
    RealType diff1 = 0.0;
    std::for_each(b,e, [&diff1, mean=_mean](auto x){diff1 += std::pow(x - mean, 3.0);});
    diff1 /= N;
    RealType diff2 = 0.0;
    std::for_each(b,e, [&diff2, mean=_mean](auto x){diff2 += std::pow(x - mean, 2.0);});
    diff2 /= (N-1);
    return diff1 / std::pow(diff2, 3.0/2.0);
}

void print_bbox(std::tuple<int, int, int , int> bbox) {
    std::cout << std::get<0>(bbox) << ", " << std::get<1>(bbox) << ", " << std::get<2>(bbox) << ", " << std::get<3>(bbox) << std::endl;
}
std::vector<Cell> generate_lattice_CA_diffusion(int msx, int msy,
                      int x_proc_idx, int y_proc_idx, int cell_in_my_cols, int cell_in_my_rows, const std::vector<int>& water_cols){
    int cell_per_process = cell_in_my_cols*cell_in_my_rows;

    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

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

    return my_cells;
}

std::vector<Cell> generate_lattice_percolation_diffusion(int msx, int msy,
                                                int x_proc_idx, int y_proc_idx, int cell_in_my_cols, int cell_in_my_rows, const std::vector<int>& water_cols){
    const int smallest_coordinate = std::min(msx, msy);
    const int minr = (int) (smallest_coordinate*0.01f), maxr = (int) (smallest_coordinate*0.05f), avgr = maxr-minr;
    int circles = (int) ( 2*(msx * msy) / (M_PI*(avgr*avgr)));
    int nb_var = 4;
    int circle_data_space = circles*nb_var;
    std::vector<float> circles_var(circle_data_space);
    if(!rank){
        std::uniform_int_distribution<int> distx(0, msx), disty(0, msy), dist_r(minr, maxr);
        std::uniform_real_distribution<float> erosion_dist(0.01f, 0.8f);
        //std::gamma_distribution<float> erosion_dist(2, 0.01);
        for (int k = 0; k < circle_data_space; k+=4) {
            auto x = distx(gen); auto y = disty(gen);
            auto r = dist_r(gen);
            circles_var[k]   = ((float) x);
            circles_var[k+1] = ((float) y);
            circles_var[k+2] = ((float) r);
            circles_var[k+3] = (float)  erosion_dist(gen);
        }
    }
    MPI_Bcast(&circles_var.front(), circle_data_space, MPI_FLOAT, 0, MPI_COMM_WORLD);
    std::vector<std::tuple<int, int, int, float>> rock_fn;
    for (int l = 0; l < circle_data_space; l+=4) {
        rock_fn.emplace_back(circles_var[l], circles_var[l+1], circles_var[l+2], circles_var[l+3]);
    }

    int cell_per_process = cell_in_my_cols*cell_in_my_rows;
    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);

    for(int j = 0; j < cell_in_my_cols; ++j) {
        for(int i = 0; i < cell_in_my_rows; ++i) {
            int id = cell_in_my_rows * x_proc_idx + i + msx * (j + (y_proc_idx * cell_in_my_cols));
            bool is_rock = false;
            float eprob = 1.0;
            auto pos = cell_to_global_position(msx, msy, id);
            for ( auto circle : rock_fn) {
                int x2, y2, r; std::tie(x2, y2, r, eprob) = circle;
                is_rock = std::sqrt( std::pow(x2-pos.first,2) + std::pow(y2 - pos.second,2) ) < r;
                if(is_rock) break;
            }
            my_cells.emplace_back(id, is_rock ? ROCK_TYPE : WATER_TYPE, is_rock ? 0.0 : 1.0, is_rock ? eprob: 1.0);
        }
    }
    return my_cells;
}

void update_cell_weights(std::vector<Cell>* _my_cells, double relative_slope){
    std::vector<Cell>& my_cells = *(_my_cells);
    for(auto& cell : my_cells) {
        if(cell.type == ROCK_TYPE) { // rock
            cell.weight = (float) relative_slope;
        }
    }
}


#include "../include/io.hpp"

int main(int argc, char **argv) {
    int worldsize;

    MPI_Init(&argc, &argv);
    auto world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldsize);
    MPI_Comm_rank(world, &rank);
    std::string str_rank = "[RANK " + std::to_string(rank) + "] ";
    // Generate data
    int xprocs = std::atoi(argv[2]), yprocs = std::atoi(argv[1]);
    int cell_per_process = std::atoi(argv[3]);
    const int MAX_STEP = std::atoi(argv[4]);
    const int seed = std::atoi(argv[5]);

    zz::log::config_from_file("logger.cfg");
    perflogger = zz::log::get_logger("perf",  true);
    steplogger = zz::log::get_logger("steps", true);

    if(xprocs * yprocs != worldsize) {
        steplogger->fatal() << "Grid size does not match world size";
        MPI_Abort(world, 1);
        return EXIT_FAILURE;
    }

    int cell_in_my_rows = (int) std::sqrt(cell_per_process), cell_in_my_cols = cell_in_my_rows;
    int xcells = cell_in_my_rows * xprocs, ycells = cell_in_my_rows * yprocs;
    std::vector<int> cols(ycells);
    std::iota(cols.begin(), cols.end(), 0);

    std::mt19937 gen(seed); ///predictable sequence of value

    std::shuffle(cols.begin(), cols.end(), gen);
    cols = {0, 1};
    std::vector<int> water_cols(cols.begin(), cols.begin()+1 );
    std::sort(water_cols.begin(), water_cols.end());
    //std::for_each(water_cols.cbegin(), water_cols.cend(), [](auto v){ std::cout << v << std::endl; });

    int& msx = Cell::get_msx(); msx = xcells;
    int& msy = Cell::get_msy(); msy = ycells;

    int shape[2] = {msx,msy};
#ifdef PRODUCE_OUTPUTS
    //if(!rank) cnpy::npz_save("out.npz", "shape", &shape[0], {2}, "w");
    if(!rank) cnpy::npz_save("gids-out.npz", "shape", &shape[0], {2}, "w");
#endif
    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_global_position(xprocs, yprocs, rank);
    //std::cout << x_proc_idx << " ; " << y_proc_idx << std::endl;
    unsigned long total_cell = cell_per_process * worldsize;
    auto datatype   = Cell::register_datatype();

    if(!rank) {
        steplogger->info("CPU COUNT:")    << worldsize;
        steplogger->info("GRID PSIZE X:") << xprocs;
        steplogger->info("GRID PSIZE Y:") << yprocs;
        steplogger->info("GRID SIZE  X:") << msx;
        steplogger->info("GRID SIZE  Y:") << msy;
        steplogger->info("EACH SIZE  X:") << xcells;
        steplogger->info("EACH SIZE  Y:") << ycells;
    }

    if(!rank) steplogger->info() << cell_in_my_cols << " " << cell_in_my_rows;

    auto my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);

    const int my_cell_count = my_cells.size();
#ifdef PRODUCE_OUTPUTS
    std::vector<std::array<int,2>> all_types(total_cell);
    std::vector<std::array<int,2>> my_types(my_cell_count);
    for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
    gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
    if(!rank) {
        std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
        std::vector<int> types, water_gid;
        std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
        std::for_each(all_types.begin(), all_types.end(), [&water_gid](auto e){if(!e[1]) water_gid.push_back(e[0]);});
        assert((*(all_types.end() - 1))[0] == total_cell-1);
        //cnpy::npz_save("out.npz", "step-"+std::to_string(0), &types[0], {total_cell}, "a");
        cnpy::npz_save("gids-out.npz", "step-"+std::to_string(0), &water_gid[0], {water_gid.size()}, "a");
    }
    if(!rank) all_types.resize(total_cell);
#endif
    int recv, sent;
    if(!rank) steplogger->info() << "End of map generation";
    /* Initial load balancing */
    auto zoltan_lb = zoltan_create_wrapper(true, world);
    zoltan_load_balance(&my_cells, zoltan_lb, true, true);
    /* lets make it fun now...*/
    int minx, maxx, miny, maxy;
    std::vector<size_t> data_pointers;

    std::vector<double> lb_costs;
    SlidingWindow<double> window(20); //sliding window with max size = TODO: tune it?
    int ncall = 20, pcall=0;
    PAR_START_TIMING(loop_time, world);
    for(unsigned int step = 0; step < MAX_STEP; ++step) {

        const int my_cell_count = my_cells.size();
        if(!rank) steplogger->info() << "Beginning step "<< step;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// COMPUTATION START
        PAR_START_TIMING(step_time, world);
        auto remote_cells = zoltan_exchange_data(zoltan_lb,my_cells,&recv,&sent,datatype.element_datatype,world,1.0);
        auto bbox = get_bounding_box(my_cells, remote_cells);

        populate_data_pointers(msx, msy, &data_pointers, my_cells, remote_cells, bbox);
        my_cells = dummy_erosion_computation2(msx, msy, my_cells, remote_cells, data_pointers, bbox);

        CHECKPOINT_TIMING(step_time, my_step_time);
        window.add(my_step_time);
        PAR_STOP_TIMING(step_time, world);
        /// COMPUTATION STOP
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<double> timings(worldsize);
        std::vector<double> slopes (worldsize);
        std::vector<double> rel_slopes;
        std::vector<int>    PE_cells(worldsize);

        double my_time_slope = get_slope<double>(window.data_container), total_slope, relative_slope;

        MPI_Allgather(&my_step_time, 1, MPI_DOUBLE, &timings.front(), 1, MPI_DOUBLE, world);
        MPI_Allgather(&my_time_slope, 1, MPI_DOUBLE, &slopes.front(), 1, MPI_DOUBLE, world); //TODO: propagate information differently

        total_slope = std::max(std::accumulate(slopes.begin(), slopes.end(), 0.0), 0.0);
        std::for_each(slopes.cbegin(), slopes.cend(), [&](auto v){ rel_slopes.push_back(std::max(v/total_slope,0.0)); });
        double max = *std::max_element(timings.cbegin(), timings.cend()),
               average = std::accumulate(timings.cbegin(), timings.cend(), 0.0) / worldsize,
               load_imbalance = (max / average - 1.0) * 100.0,
               skew = skewness<double>(timings.cbegin(), timings.cend());

        if(!rank) {
            steplogger->info("stats: ") << "load imbalance: " << load_imbalance << " skewness: " << skew;
            perflogger->info("\"step\":") << step << ",\"load_imbalance\": " << load_imbalance << ",\"skewness\": " << skew << ",\"loads\": [" << timings << "],\"slopes\":["<<slopes<<"]";
        }

        int cell_cnt = my_cells.size();
#ifdef PRODUCE_OUTPUTS
        std::vector<std::array<int,2>> my_types(my_cell_count);
        for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
        gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
        if(!rank) {
            std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
            std::vector<int> types, rock_gid;
            std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
            std::for_each(all_types.begin(), all_types.end(), [&rock_gid](auto e){if(!e[1]) rock_gid.push_back(e[0]);});
            //cnpy::npz_save("out.npz", "step-"+std::to_string(step+1), &types[0], {total_cell}, "a");
            cnpy::npz_save("gids-out.npz", "step-"+std::to_string(step+1), &rock_gid[0], {rock_gid.size()}, "a");
        }
#endif
        if(pcall + ncall <= step && total_slope > 0 && skew > 0){
            if(!rank) steplogger->info("call LB");
            relative_slope = std::max(my_time_slope/total_slope, 0.0); //TODO: OVERLOADED PROCESSORS ARE THOSE WITH POSITIVE MYSLOPE-MEAN(ALL_WORKLOAD_SLOPE).
            PAR_START_TIMING(lb_cost, world);
            update_cell_weights(&my_cells, relative_slope);
            zoltan_load_balance(&my_cells, zoltan_lb, true, true);
            PAR_STOP_TIMING(lb_cost, world);

            window.data_container.clear();
            cell_cnt = my_cells.size();
        }
        MPI_Gather(&cell_cnt, 1, MPI_INT, &PE_cells.front(), 1, MPI_INT, 0, world);
        //if(!rank) perflogger->info("step:") << step << ",cells: (" << PE_cells << ")";
        if(!rank) steplogger->info() << "Stop step "<< step;
    }
    PAR_STOP_TIMING(loop_time, world);

    datatype.free_datatypes();
    MPI_Finalize();
    return 0;
}
