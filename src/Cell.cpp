//
// Created by xetql on 5/3/19.
//
#include "Utils.hpp"
#include "Cell.hpp"

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
        data_pointers.clear();
        data_pointers.resize(my_box);
        std::fill(data_pointers.begin(), data_pointers.end(), msx * msy + 1);
        for (size_t i = 0; i < mine_size; ++i) {
            const Cell& cell = my_cells[i];
            auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
            data_pointers[lid] = i;
        }
    }

    auto remote_size = remote_cells.size();

    for (size_t i = 0; i < remote_size; ++i) {
        const Cell& cell = remote_cells[i];
        auto lid = position_to_cell(x2-x1, y2-y1, cell_to_local_position(msx, msy, bbox, cell.gid));
        data_pointers[lid] = i+mine_size;
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
    int x,y;
    for (size_t i = 0; i < mine_size; ++i) {
        const Cell& cell = cells[i];
        cell_to_local_position(msx, msy, bbox, cell.gid, &x, &y);
        auto lid = position_to_cell(x2-x1, y2-y1, x, y);
        data_pointers[lid] = i+displ;
    }
}
// O(n) => very long
template<class IntegerType>
IntegerType compute_estimated_workload(const std::vector<Cell>& _my_cells) {
    IntegerType load = 0;
    for (const auto& cell : _my_cells)
        load += (IntegerType) cell.weight;
    return load;
}
template<>
int64_t compute_estimated_workload(const std::vector<Cell>& _my_cells) {
    int64_t load = 0;
    for (const auto& cell : _my_cells)
        load += (int64_t) cell.weight;
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
