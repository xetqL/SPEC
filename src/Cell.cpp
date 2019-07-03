//
// Created by xetql on 5/3/19.
//
#include "Utils.hpp"
#include "Cell.hpp"

void add_remote_data_to_arr(int msx, int msy,
                            std::vector<size_t>* data_pointers,
                            long mine_size,
                            const std::vector<Cell>& remote_cells,
                            const std::tuple<int, int, int, int>& bbox){
    //std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1, x2, y1, y2) = bbox;
    auto remote_size = remote_cells.size();
    for (size_t i = 0; i < remote_size; ++i) {
        const Cell& cell = remote_cells[i];
        auto lid = position_to_cell(x2 - x1 + 1, y2 - y1 + 1, cell_to_local_position(msx, msy, bbox, cell.gid));
        assert(lid < data_pointers->size());
        data_pointers->operator[](lid) = i+mine_size;
    }
}

void add_remote_data_to_arr(int msx, int msy,
                            size_t* data_pointers,
                            long mine_size,
                            const std::vector<Cell>& remote_cells,
                            const std::tuple<int, int, int, int>& bbox){
    //std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1, x2, y1, y2) = bbox;
    auto remote_size = remote_cells.size();
    for (size_t i = 0; i < remote_size; ++i) {
        const Cell& cell = remote_cells[i];
        auto lid = position_to_cell(x2 - x1 + 1, y2 - y1 + 1, cell_to_local_position(msx, msy, bbox, cell.gid));
        data_pointers[lid] = i+mine_size;
    }
}

void populate_data_pointers(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            const std::vector<Cell>& my_cells,
                            const std::vector<Cell>& remote_cells,
                            const std::tuple<int, int, int, int>& bbox,
                            bool create) {
    std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1, x2, y1, y2) = bbox;
    int my_box = (x2-x1) * (y2-y1);
    auto mine_size   = my_cells.size();
    if(create) {
        data_pointers.resize(my_box);
        std::fill(data_pointers.begin(), data_pointers.end(), msx * msy + 1);
        for (size_t i = 0; i < mine_size; ++i) {
            const Cell& cell = my_cells[i];
            auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
            data_pointers[lid] = i;
        }
    }

    add_remote_data_to_arr(msx, msy, _data_pointers->data(), mine_size, remote_cells, bbox);
}

void init_populate_data_pointers(int msx, int msy,
                                 std::vector<size_t>* _data_pointers,
                                 const std::vector<Cell>& my_cells,
                                 const std::tuple<int, int, int, int>& bbox) {
    std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1, x2, y1, y2) = bbox;
    int my_box = (x2-x1) * (y2-y1);
    auto mine_size   = my_cells.size();

    data_pointers.resize(my_box + 4 * (msx < msy ? msx : msy));
    std::fill(data_pointers.begin(), data_pointers.end(), msx * msy + 1);
    for (size_t i = 0; i < mine_size; ++i) {
        const Cell& cell = my_cells[i];
        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
        data_pointers[lid] = i;
    }
}


void populate_data_pointers(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            const std::vector<Cell>& cells,
                            int displ,
                            const std::tuple<int, int, int, int>& bbox,
                            bool create) {
    std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1, x2, y1, y2) = bbox;
    int my_box = (x2-x1) * (y2-y1);
    if(create) {
        data_pointers.clear();
        data_pointers.resize(my_box);
        std::fill(data_pointers.begin(), data_pointers.end(), msx * msy + 1);
    }
    auto mine_size = cells.size();
    for (size_t i = 0; i < mine_size; ++i) {
        const Cell& cell = cells[i];
        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
        data_pointers[lid] = i+displ;
    }
}

float compute_estimated_workload(const std::vector<Cell>& _my_cells, int type) {
    float load = 0;
    for (const auto& cell : _my_cells) load += cell.weight * (type == cell.type);
    return load;
}

float compute_estimated_workload(const std::vector<Cell>& _my_cells) {
    float load = 0;
    for (const auto& cell : _my_cells) load += cell.weight * (cell.type);
    return load;
}

long compute_effective_workload(const std::vector<Cell>& _my_cells, int type) {
    return std::count_if(_my_cells.begin(), _my_cells.end(), [&](auto c){return c.type == type;});
}

/**
 * Divide the speed at which I am loading among all the work units that can load
 * @param _my_cells
 * @param slope
 */

void update_cell_weights(std::vector<Cell>* _my_cells, double slope, int type){
    std::vector<Cell>& my_cells = *(_my_cells);
    int nb_type = 0;
    slope = std::max(slope, 0.0); // wtf is a negative slope
    for(auto& cell : my_cells) if(cell.type == type) nb_type++;
    for(auto& cell : my_cells) if(cell.type == type) cell.weight = cell.weight - (float) slope * 1.0f / nb_type;
}

void reset_cell_weights(std::vector<Cell>* _my_cells){
    std::vector<Cell>& my_cells = *(_my_cells);
    for(auto& cell : my_cells) cell.weight = cell.type ? WATER_WEIGHT : ROCK_WEIGHT;
}
