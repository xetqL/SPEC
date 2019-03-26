//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_CELL_HPP
#define SPEC_CELL_HPP

#include <ostream>
#include "communication.hpp"
#include "utils.hpp"

#define ROCK_WEIGHT 0.0
#define WATER_WEIGHT 1.0

struct CellStatistics {
    double average_load;
};

struct Cell {
    int gid, type; //type = -1:empty, 0:rock, 1:water
    float weight, erosion_probability;
    double average_load = 2000;

    Cell() : gid(0), type(0), weight(ROCK_WEIGHT), erosion_probability(0) {};
    Cell(int gid, int type, float weight, float erosion_probability)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability) {};
    Cell(int gid, int type, float weight, float erosion_probability, double average_load)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability), average_load(average_load) {};

    template<class NumericType>
    std::array<NumericType, 2> get_position_as_array() const {
        auto position_as_pair = cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
        std::array<NumericType, 2> array = {(NumericType) position_as_pair.first, (NumericType) position_as_pair.second};
        return array;
    }

    std::pair<int, int> get_position_as_pair() const {
        return cell_to_global_position(Cell::get_msx(), Cell::get_msy(), gid);
    }

    static CommunicationDatatype register_datatype() {

        MPI_Datatype cell_datatype, gid_type_datatype;

        MPI_Aint intex, lb, floatex;

        const int number_of_int_elements = 2;
        const int number_of_float_elements = 2;
        const int number_of_double_elements = 1;

        int blockcount_element[3];

        blockcount_element[0] = number_of_int_elements; // gid, lid, exit, waiting_time
        blockcount_element[1] = number_of_float_elements; // position <x,y>
        blockcount_element[2] = number_of_double_elements; // position <x,y>
        //int
        MPI_Type_contiguous(number_of_int_elements, MPI_INT, &gid_type_datatype); // position
        MPI_Type_commit(&gid_type_datatype);

        MPI_Datatype blocktypes[3];
        blocktypes[0] = MPI_INT;
        blocktypes[1] = MPI_FLOAT;
        blocktypes[2] = MPI_DOUBLE;

        MPI_Type_get_extent(MPI_INT, &lb, &intex);
        MPI_Type_get_extent(MPI_FLOAT, &lb, &floatex);

        MPI_Aint offset[3];
        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = 2*intex;
        offset[2] = 2*intex + 2*floatex;

        MPI_Type_create_struct(3, blockcount_element, offset, blocktypes, &cell_datatype);
        MPI_Type_commit(&cell_datatype);

        return CommunicationDatatype(cell_datatype, gid_type_datatype);
    }

    static int set_msx(int _msx){
        static int msx = _msx;
    }
    static int set_msy(int _msy){
        static int msy = _msy;
    }
    static int& get_msx(){
        static int msx;
        return msx;
    }
    static int& get_msy(){
        static int msy;
        return msy;
    }

    friend std::ostream &operator<<(std::ostream &os, const Cell &cell) {
        os << "gid: " << cell.gid << " type: " << cell.type << " weight: " << cell.weight << " erosion_probability: "
           << cell.erosion_probability;
        return os;
    }

};

void populate_data_pointers(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            const std::vector<Cell>& my_cells,
                            const std::vector<Cell>& remote_cells,
                            const std::tuple<int, int, int, int>& bbox,
                            bool create = false) {
    std::vector<size_t>& data_pointers = *(_data_pointers);
    int x1, x2, y1, y2; std::tie(x1, x2, y1, y2) = bbox;
    int my_box = (x2-x1) * (y2-y1);

    if(create){
        data_pointers.clear();
        data_pointers.resize(my_box);
        std::fill(data_pointers.begin(), data_pointers.end(), msx * msy + 1);
    }

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

void populate_data_pointers(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            const std::vector<Cell>& cells,
                            int displ,
                            const std::tuple<int, int, int, int>& bbox,
                            bool create = false) {
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

float compute_estimated_workload(const std::vector<Cell>& _my_cells) {
    float load = 0;
    for (const auto& cell : _my_cells) load += cell.weight;
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
template<class Predicate>
void update_cell_weights(std::vector<Cell>* _my_cells, double slope, int type, Predicate f){
    std::vector<Cell>& my_cells = *(_my_cells);
    int nb_type = 0;
    slope = std::max(slope, 0.0); // wtf is a negative slope
    for(auto& cell : my_cells) if(cell.type == type) nb_type++;
    for(auto& cell : my_cells) if(cell.type == type) {
        cell.weight = f(cell.weight, (float) slope);
    }

}
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

#endif //SPEC_CELL_HPP
