//
// Created by xetql on 5/3/19.
//
//
// Created by xetql on 2/25/19.
//

#ifndef SPEC_MAIN_HPP
#define SPEC_MAIN_HPP

#include "Utils.hpp"

#include <vector>
#include <random>
#include <algorithm>
#include <thread>
#include "Cell.hpp"

//#define WATER_TYPE 1
//#define ROCK_TYPE 0
//#define EMPTY_TYPE -1

static std::random_device rd;
static std::mt19937 gen(rd());
static std::normal_distribution<float> ndist(9, 1);
static std::uniform_real_distribution<float> udist(0, 1);
static volatile double res = 0.0;
static volatile double one = 1.0;

/// SIMULATION
void consume_cpu_flops(double flop_to_consume, int i);
void consume_cpu_flops(float flop_to_consume, int i);
/*
/// GENERATION
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
}*/

std::vector<Cell> generate_lattice_single_type( int msx, int msy,
                                                int x_proc_idx, int y_proc_idx,
                                                int cell_in_my_cols, int cell_in_my_rows,
                                                int type, float weight, float erosion_probability);

void generate_lattice_rocks( const int rocks_per_stripe, int msx, int msy,
                             std::vector<Cell>* _cells,
                             float erosion_probability,
                             int begin_stripe, int end_stripe, float radius);
int checkpoint(int h, int k, int x, int y, int a, int b);



void generate_lattice_ellipse(int msx, int msy,
                              std::vector<Cell>* _cells,
                              float erosion_probability,
                              int begin_stripe, int end_stripe);
/*
std::vector<Cell> generate_lattice_percolation_diffusion(int msx, int msy,
                                                         int x_proc_idx, int y_proc_idx, int cell_in_my_cols, int cell_in_my_rows, const std::vector<int>& water_cols) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int smallest_coordinate = std::min(msx, msy);
    const int minr = (int) (smallest_coordinate*0.05f), maxr = (int) (smallest_coordinate*0.1f), avgr = maxr-minr;
    int circles = (int) (0.5f * (msx * msy) / (M_PI * (avgr * avgr)));
    int nb_var = 4;
    int circle_data_space = circles*nb_var;
    std::vector<float> circles_var(circle_data_space);

    if(!rank) {
        std::uniform_int_distribution<int> distx(0, msx), disty(0, msy), dist_r(minr, maxr);
        std::uniform_real_distribution<float> erosion_dist(0.01f, 0.1f);
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

            if(pos.first < (msx / 4)){
                main_type = WATER_TYPE;
                main_weight = WATER_WEIGHT;
            } else {
                main_type = ROCK_TYPE;
                main_weight = ROCK_WEIGHT;
            }
            for ( auto circle : rock_fn) {
                int x2, y2, r; std::tie(x2, y2, r, eprob) = circle;
                is_rock = std::sqrt( std::pow(x2-pos.first,2) + std::pow(y2 - pos.second,2) ) < r;
                if(is_rock) break;
            }
            my_cells.emplace_back(id, is_rock ? ROCK_TYPE : main_type, is_rock ? ROCK_WEIGHT : main_weight, is_rock ? eprob : 1.0);
        }
    }
    return my_cells;
}*/

#ifdef PRODUCE_OUTPUTS
/*
#include <cnpy.h>
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
            main_weight = WATER_WEIGHT;
        } else {
            main_type = ROCK_TYPE;
            main_weight = ROCK_WEIGHT;
        }
        auto val = lattice[k];
        my_cells.emplace_back(id, val > 0 ? ROCK_TYPE : main_type, val > 0 ? ROCK_WEIGHT : main_weight, val > 0 ? erosion_probabilities[val] : 1.0);
    }
    std::sort(my_cells.begin(), my_cells.end(), [](auto a, auto b){return a.gid < b.gid;});
    return my_cells;
}
 */
#endif

/// COMPUTATION
/*
std::vector<Cell> dummy_erosion_computation(int msx, int msy,
                                            const std::vector<Cell>& _my_cells,
                                            const std::vector<Cell>& remote_cells,
                                            const std::vector<size_t>& data_pointers,
                                            const std::tuple<int, int, int, int>& bbox) {


    std::vector<Cell> my_cells = _my_cells;
    const std::vector<Cell> &my_old_cells = _my_cells;
    int x1, x2, y1, y2;
    std::tie(x1, x2, y1, y2) = bbox;
    int total_box = (x2 - x1) * (y2 - y1);
    const size_t my_cell_count = my_cells.size();
    const size_t remote_count  = remote_cells.size();
    const size_t all_count = my_cell_count + remote_count;
    std::vector<size_t> idx_neighbors(3, -1);
    for (unsigned int i = 0; i < all_count; ++i) {
        const Cell *cell;
        if (i < my_cell_count) cell = &my_old_cells[i];
        else cell = &remote_cells[i - my_cell_count];
        if (cell->type) { // Do nothing if type == ROCK_TYPE
            auto lid = position_to_cell(x2 - x1, y2 - y1, cell_to_local_position(msx, msy, bbox, cell->gid));
            if (cell->type) {
                std::fill(idx_neighbors.begin(), idx_neighbors.end(), -1);
                if (lid + (x2 - x1) + 1 < total_box) idx_neighbors[0] = (data_pointers[lid + (x2 - x1) + 1]);
                if (lid + 1 < total_box) idx_neighbors[1] = (data_pointers[lid + 1]);
                if (((lid - (x2 - x1)) + 1) >= 0) idx_neighbors[2] = data_pointers[((lid - (x2 - x1)) + 1)];
                for (auto idx_neighbor : idx_neighbors) {
                    if (idx_neighbor > my_cell_count) continue;
                    if (my_cells[idx_neighbor].type) continue;
                    auto erosion_proba = my_old_cells[idx_neighbor].erosion_probability;
                    auto p = udist(gen);
                    if (p < erosion_proba) {
                        my_cells[idx_neighbor].type = 1;
                        my_cells[idx_neighbor].weight = 1;
                    }
                }
            }
        }
    }
    return my_cells;
}*/

std::pair<unsigned long, std::vector<unsigned long>> create_water_ptr_vector(const std::vector<Cell>& cells);

const std::vector<const Cell*> create_water_ptr(const std::vector<Cell>& cells);

/*
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

        if(cell->type) { // Do nothing if type == ROCK_TYPE
            auto __pos = cell_to_local_position(msx, msy, bbox, cell->gid);
            auto lid = position_to_cell(x2-x1, y2-y1, __pos);
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
                //auto x = cell_to_global_position(msx, msy, idx_neighbor);
                auto p = udist(gen);
                if( erosion_proba < 1.0 ) {
                    if (p < (theta) * erosion_proba) {
                        my_cells[idx_neighbor].type   = 1;
                        my_cells[idx_neighbor].weight = 1;
                    }
                } else {
                    if (p < (msx - __pos.first) / (float) msx) {
                        my_cells[idx_neighbor].type   = 1;
                        my_cells[idx_neighbor].weight = 1;
                    }
                }
            }

            //DO NOT OPTIMIZE; SIMULATE COMPUTATION OF LBM FLUID WITH BGK D2Q9
            consume_cpu_flops(flops);


        }
    }
    return my_cells;
}
*/
/*std::tuple<std::vector<Cell>, std::vector<unsigned long>, double> dummy_erosion_computation3(int step,
                                                                                             int msx, int msy,
                                                                                             const std::vector<Cell>& my_old_cells,
                                                                                             const std::vector<unsigned long>& my_water_ptr,
                                                                                             const std::vector<Cell>& remote_cells,
                                                                                             const std::vector<unsigned long>& remote_water_ptr,
                                                                                             const std::vector<size_t>& data_pointers,
                                                                                             const std::tuple<int, int, int, int>& bbox);*/

std::vector<Cell>
dummy_erosion_computation3(int step,
                           int msx, int msy,
                           const std::vector<Cell>& my_old_cells,
                           const std::vector<unsigned long>& my_water_ptr,
                           const std::vector<Cell>& remote_cells,
                           const std::vector<unsigned long>& remote_water_ptr,
                           const size_t *data_pointers,
                           const std::tuple<int, int, int, int>& bbox,
                           std::vector<unsigned long>* new_water_cells, double* total_weight);

void compute_fluid(const std::vector<Cell>& my_old_cells);

void compute_fluid(float total_cells);

void compute_fluid_time(float total_cells);

#endif //SPEC_PROBLEM_HPP
