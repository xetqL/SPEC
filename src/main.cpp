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

const float bytes_per_flop = 2.77f;
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
    volatile float data_store = 0xAF;
    for(unsigned int i = 0; i < all_count; ++i) {
        const Cell* cell;
        if(i < my_cell_count) cell = &my_old_cells[i];
        else cell = &remote_cells[i - my_cell_count];

        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell->gid));
        if(cell->type) { // Do nothing if type == ROCK_TYPE
            auto __pos = cell_to_local_position(msx, msy, bbox, cell->gid);
            auto lid = position_to_cell(x2-x1, y2-y1, __pos);
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
                    if( erosion_proba < 1.0 ){
                        if (p < (theta) * erosion_proba) {
                            my_cells[idx_neighbor].type   = 1;
                            my_cells[idx_neighbor].weight = 1;
                        }
                    } else {
                        if (p < 1.0f - (msx - __pos.first) / (float) msx) {
                            my_cells[idx_neighbor].type   = 1;
                            my_cells[idx_neighbor].weight = 1;
                        }
                    }
                }

                /*DO NOT OPTIMIZE; SIMULATE COMPUTATION OF LBM FLUID WITH BGK D2Q9*/
                data_store = flopdist(gen);
                const int nb_flop = (int) data_store;
                const int nb_bytes = (int) (nb_flop * bytes_per_flop) + 1;
                float res = 0.0;
                for(int i = 0; i < nb_flop; ++i) {
                    __asm__(""); res = data_store * data_store;
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
    const int minr = (int) (smallest_coordinate*0.05f), maxr = (int) (smallest_coordinate*0.1f), avgr = maxr-minr;
    int circles = (int) ( 2*(msx * msy) / (M_PI*(avgr*avgr)));
    int nb_var = 4;
    int circle_data_space = circles*nb_var;
    std::vector<float> circles_var(circle_data_space);

    if(!rank) {
        std::uniform_int_distribution<int> distx(0, msx), disty(0, msy), dist_r(minr, maxr);
        std::uniform_real_distribution<float> erosion_dist(0.005f, 0.1f);
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
    int main_type;
    float main_weight;
    for(int j = 0; j < cell_in_my_cols; ++j) {

        for(int i = 0; i < cell_in_my_rows; ++i) {
            int id = cell_in_my_rows * x_proc_idx + i + msx * (j + (y_proc_idx * cell_in_my_cols));
            bool is_rock = false;
            float eprob = 1.0;
            auto pos = cell_to_global_position(msx, msy, id);

            if(x_proc_idx+i == 0){
                main_type = WATER_TYPE;
                main_weight = 1.0;
            } else {
                for ( auto circle : rock_fn) {
                    int x2, y2, r; std::tie(x2, y2, r, eprob) = circle;
                    is_rock = std::sqrt( std::pow(x2-pos.first,2) + std::pow(y2 - pos.second,2) ) < r;
                    if(is_rock) break;
                }
                main_type = ROCK_TYPE;
                main_weight = 0.0;
            }
            my_cells.emplace_back(id, is_rock ? ROCK_TYPE : main_type, is_rock ? 0.0 : main_weight, is_rock ? eprob : 1.0);
        }
    }
    return my_cells;
}

#ifdef PRODUCE_OUTPUTS
std::vector<Cell> generate_lattice_percolation_diffusion(int msx, int msy, int x_proc_idx, int y_proc_idx,
                                                         int cell_in_my_cols, int cell_in_my_rows,
                                                         const std::vector<int>& water_cols,
                                                         std::string rock_filename) {
    std::uniform_real_distribution<float> erosion_dist(0.005f, 0.1f);
    auto rock_shapes = cnpy::npy_load(std::move(rock_filename));
    auto rock_msx = rock_shapes.shape[1]; auto rock_msy = rock_shapes.shape[0];
    assert(cell_in_my_cols == rock_msx); assert(cell_in_my_rows == rock_msy);
    auto cell_per_process = cell_in_my_cols * cell_in_my_rows;
    std::vector<Cell> my_cells; my_cells.reserve(cell_per_process);
    auto lattice = rock_shapes.data<int>();
    int nb_labels = get_max(lattice, cell_per_process) + 1;
    std::vector<float> erosion_probabilities(nb_labels, 0);
    for(float& p : erosion_probabilities) p = erosion_dist(gen);
    int main_type;
    float main_weight;
    for(unsigned long k = 0; k < cell_per_process; ++k) {
        int i,j; std::tie(i,j) = cell_to_position(cell_in_my_cols, cell_in_my_rows, k);
        int id = cell_in_my_rows * x_proc_idx + i + msx * (j + (y_proc_idx * cell_in_my_cols));
        if(x_proc_idx+i == 0){
            main_type = WATER_TYPE;
            main_weight = 1.0;
        } else {
            main_type = ROCK_TYPE;
            main_weight = 0.0;
        }
        auto val = lattice[k];
        my_cells.emplace_back(id, val > 0 ? ROCK_TYPE : main_type, val > 0 ? 0.0 : main_weight, val > 0 ? erosion_probabilities[val] : 1.0);
    }
    std::sort(my_cells.begin(), my_cells.end(), [](auto a, auto b){return a.gid < b.gid;});
    return my_cells;
}
#endif

float compute_load(const std::vector<Cell>& _my_cells) {
    float load = 0;
    for (const auto& cell : _my_cells) load += cell.weight;
    return load;
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
    if(!rank) cnpy::npz_save("gids-out.npz", "shape", &shape[0], {2}, "w");
#endif

    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_global_position(xprocs, yprocs, rank);
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
    std::vector<Cell> my_cells;
#ifndef PRODUCE_OUTPUTS
    my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);
#else
    if(argc == 7){
        auto lattice_fname = argv[6];
        my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols, lattice_fname);
    } else{
        my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);
    }
#endif

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
    std::vector<double> lb_costs;
    auto zoltan_lb = zoltan_create_wrapper(true, world);

    PAR_START_TIMING(current_lb_cost, world);
    zoltan_load_balance(&my_cells, zoltan_lb, true, true);
    PAR_STOP_TIMING(current_lb_cost, world);
    lb_costs.push_back(current_lb_cost);

    /* lets make it fun now...*/
    int minx, maxx, miny, maxy;
    std::vector<size_t> data_pointers;

#ifdef LB_METHOD
    auto avg_lb_cost = mean<double>(lb_costs.begin(), lb_costs.end());
    auto tau = 25; // ???
#endif

    SlidingWindow<double> window(20); //sliding window with max size = TODO: tune it?
    int ncall = 20, pcall=0;
#if LB_METHOD==1
    std::cout << "SALE PUTE"<< std::endl;
    ncall = 100;
#endif

    float average_load_at_lb;
    double skew = 0, relative_slope, my_time_slope, total_slope;
    std::vector<double> timings(worldsize);
    PAR_START_TIMING(loop_time, world);
    for(unsigned int step = 0; step < MAX_STEP; ++step) {
        const int my_cell_count = my_cells.size(); //TODO: impl. iterative mean formula..

#ifdef LB_METHOD
        PAR_START_TIMING(step_time, world);
#if LB_METHOD==1   // load balance every 100 iterations
        if((pcall + ncall) == step) {
            if(!rank) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            zoltan_load_balance(&my_cells, zoltan_lb, true, true);
            PAR_STOP_TIMING(current_lb_cost, world);
            lb_costs.push_back(current_lb_cost);
            ncall = 100;
            if(!rank) steplogger->info("next LB call at: ") << (step+ncall);
            pcall = step;
        }
#elif LB_METHOD == 2 // http://sc16.supercomputing.org/sc-archive/tech_poster/poster_files/post247s2-file3.pdf
        double my_time_slope = get_slope<double>(window.data_container), total_slope, relative_slope;
        if (ncall + pcall == step) {
            if(!rank) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            zoltan_load_balance(&my_cells, zoltan_lb, true, true);
            PAR_STOP_TIMING(current_lb_cost, world);
            lb_costs.push_back(current_lb_cost);
            avg_lb_cost = mean<double>(lb_costs.begin(), lb_costs.end());
            std::vector<double> slopes(worldsize);
            MPI_Allgather(&my_time_slope, 1, MPI_DOUBLE, &slopes.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
            auto max_slope  = *std::max_element(slopes.begin(), slopes.end());
            if(max_slope > 0) ncall = std::sqrt((2 * avg_lb_cost) / max_slope);
            else ncall = MAX_STEP;
            window.data_container.clear();
            if(!rank) steplogger->info("next LB call at: ") << (step+ncall);
            pcall = step;
        }
#elif LB_METHOD == 3 // Unloading Model
        double my_slope_before_lb;
        double my_last_step_time;

        if(pcall + ncall <= step) {
            my_slope_before_lb = my_time_slope;
            std::vector<double> slopes(worldsize);
            my_time_slope = get_slope<double>(window.data_container);

            MPI_Allgather(&my_time_slope, 1, MPI_DOUBLE, &slopes.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
            auto avg_slope = mean<double>(slopes.begin(), slopes.end());
            int overloaded = std::count_if(slopes.begin(), slopes.end(), [&](auto s){ return s > avg_slope; });
            if(overloaded < worldsize / 2) {
                if(!rank) steplogger->info("call LB at: ") << step;
                my_last_step_time  = window.data_container.back();
                total_slope = std::max(std::accumulate(slopes.begin(), slopes.end(), 0.0), 0.0);
                relative_slope = std::max(my_time_slope / total_slope, 0.0); // TODO: OVERLOADED PROCESSORS ARE THOSE WITH POSITIVE MYSLOPE-MEAN(ALL_WORKLOAD_SLOPE).
                update_cell_weights(&my_cells, relative_slope);
                PAR_START_TIMING(current_lb_cost, world);
                zoltan_load_balance(&my_cells, zoltan_lb, true, true);
                PAR_STOP_TIMING(current_lb_cost, world);
                lb_costs.push_back(current_lb_cost);
                window.data_container.clear();
            } else {
                if(!rank) steplogger->info("call LB at: ") << step;
                std::vector<double> slopes(worldsize);
                MPI_Allgather(&my_time_slope, 1, MPI_DOUBLE, &slopes.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
                auto max_slope  = *std::max_element(slopes.begin(), slopes.end());
                if(max_slope > 0){
                    PAR_START_TIMING(current_lb_cost, world);
                    zoltan_load_balance(&my_cells, zoltan_lb, true, true);
                    PAR_STOP_TIMING(current_lb_cost, world);
                    avg_lb_cost = mean<double>(lb_costs.begin(), lb_costs.end());
                    ncall = (int) std::sqrt((2 * avg_lb_cost) / max_slope);
                } else ncall = 100;

                window.data_container.clear();
                if(!rank) steplogger->info("next LB call at: ") << (step+ncall);
                pcall = step;
            }
        }
#endif
        PAR_STOP_TIMING(step_time, world);

#endif

        if(!rank) steplogger->info() << "Beginning step "<< step;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// COMPUTATION START
        PAR_RESTART_TIMING(step_time, world);
        auto remote_cells = zoltan_exchange_data(zoltan_lb,my_cells,&recv,&sent,datatype.element_datatype,world,1.0);
        auto bbox = get_bounding_box(my_cells, remote_cells);
        populate_data_pointers(msx, msy, &data_pointers, my_cells, remote_cells, bbox);
        my_cells = dummy_erosion_computation2(msx, msy, my_cells, remote_cells, data_pointers, bbox);
        CHECKPOINT_TIMING(step_time, my_step_time);
        window.add(my_step_time);

#if LB_METHOD==3 // compute ncall
        if(pcall + ncall <= step) {
            std::vector<double> slopes(worldsize), times_and_slopes(worldsize * 2);
            auto my_slope_after_lb   = my_last_step_time - my_step_time;
            auto my_slope_prediction = my_slope_before_lb - my_slope_after_lb;
            double time_and_slope[2] ={my_step_time, my_slope_prediction};
            MPI_Allgather(time_and_slope, 2, MPI_DOUBLE, &times_and_slopes.front(), 2, MPI_DOUBLE, world); // TODO: propagate information differently
            for (unsigned int i = 0; i < worldsize; ++i) {
                timings[i] = times_and_slopes[2 * i];
                slopes[i]  = times_and_slopes[2 * i + 1];
            }
            auto mu = mean<double>(timings.begin(), timings.end());
            // auto min_timing = *std::min_element(timings.begin(), timings.end());
            // auto max_slope  = *std::max_element(slopes.begin(), slopes.end());
            ncall = 1; for (int j = 0; j < worldsize; ++j) ncall = std::max(ncall, (int) ((mu-timings[j]) / slopes[j]));
            total_slope = std::max(std::accumulate(slopes.begin(), slopes.end(), 0.0), 0.0);
            skew = skewness<double>(timings.cbegin(), timings.cend());
            pcall = step;
            if(!rank) steplogger->info("next LB call at: ") << (step+ncall);

        }
#endif
        PAR_STOP_TIMING(step_time, world);
        /// COMPUTATION STOP
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<int> PE_cells(worldsize);

        MPI_Allgather(&my_step_time, 1, MPI_DOUBLE, &timings.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
        std::vector<double> slopes(worldsize);
        my_time_slope = get_slope<double>(window.data_container);
        MPI_Allgather(&my_time_slope, 1, MPI_DOUBLE, &slopes.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
        double max = *std::max_element(timings.cbegin(), timings.cend()),
               average = std::accumulate(timings.cbegin(), timings.cend(), 0.0) / worldsize,
               load_imbalance = (max / average - 1.0) * 100.0;

        if(!rank) {
            steplogger->info("stats: ") << "load imbalance: " << load_imbalance << " skewness: " << skew;
            perflogger->info("\"step\":") << step << ",\"time\": " << step_time << ",\"lb_cost\": " << mean<double>(lb_costs.begin(), lb_costs.end()) << ",\"load_imbalance\": " << load_imbalance << ",\"skewness\": " << skew << ",\"loads\": [" << timings << "],\"slopes\":["<<slopes<<"]";
        }


#ifdef PRODUCE_OUTPUTS
        int cell_cnt = my_cells.size();
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
        if(!rank) steplogger->info() << "Stop step "<< step;
    }
    PAR_STOP_TIMING(loop_time, world);
    if(!rank) perflogger->info("\"total_time\":") << loop_time;
    datatype.free_datatypes();
    MPI_Finalize();
    return 0;
}
