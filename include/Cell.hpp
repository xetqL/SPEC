//
// Created by xetql on 2/11/19.
//

#ifndef CELL_HPP
#define CELL_HPP

#include <algorithm>
#include <ostream>
#include <vector>
#include "Communication.hpp"
#include "Utils.hpp"

#define ROCK_WEIGHT 0.0
#define WATER_WEIGHT 1.0

struct Cell {
    static const int WATER_TYPE = 1;
    static const int ROCK_TYPE  = 0;

    int gid, type; //type = -1:empty, 0:rock, 1:water
    float weight, erosion_probability;
    mutable float fakeInnerData = 1.0f;
    double slope;

    Cell() : gid(0), type(0), weight(ROCK_WEIGHT), erosion_probability(0), slope(1.0) {};
    Cell(int gid, int type, float weight, float erosion_probability)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability), slope(1.0-type) {};
    Cell(int gid, int type, float weight, float erosion_probability, double slope)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability), slope(slope) {};

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
        const int number_of_float_elements = 3;
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
        offset[1] = number_of_int_elements*intex;
        offset[2] = number_of_int_elements*intex + number_of_float_elements*floatex;

        MPI_Type_create_struct(3, blockcount_element, offset, blocktypes, &cell_datatype);
        MPI_Type_commit(&cell_datatype);

        return {cell_datatype, gid_type_datatype};
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
                            bool create = false);

void populate_data_pointers(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            const std::vector<Cell>& cells,
                            int displ,
                            const std::tuple<int, int, int, int>& bbox,
                            bool create = false);
void init_populate_data_pointers(int msx, int msy,
                                 std::vector<size_t>* _data_pointers,
                                 const std::vector<Cell>& my_cells,
                                 const std::tuple<int, int, int, int>& bbox);

void add_remote_data_to_arr(int msx, int msy,
                            std::vector<size_t>* _data_pointers,
                            int mine_size,
                            const std::vector<Cell>& remote_cells,
                            const std::tuple<int, int, int, int>& bbox);
float compute_estimated_workload(const std::vector<Cell>& _my_cells);

long compute_effective_workload(const std::vector<Cell>& _my_cells, int type) ;

/**
 * Divide the speed at which I am loading among all the work units that can load
 * @param _my_cells
 * @param slope
 */
template<class MapFunc>
void update_cell_weights(std::vector<Cell>* _my_cells, double slope, int type, MapFunc f) {
    std::vector<Cell>& my_cells = *(_my_cells);
    int nb_type = 0;
    slope = std::max(slope, 0.0); // wtf is a negative slope
    for(auto& cell : my_cells) if(cell.type == type) nb_type++;
    for(auto& cell : my_cells) if(cell.type == type) {
        cell.weight = f(cell.weight, (float) slope);
    }

}

template<class Applier>
void update_cells(std::vector<Cell>* _my_cells, double slope, Applier f) {
    std::vector<Cell>& my_cells = *(_my_cells);
    slope = std::max(slope, 0.0); // wtf is a negative slope
    for(auto& cell : my_cells) f(cell, slope);
}

void update_cell_weights(std::vector<Cell>* _my_cells, double slope, int type);

void reset_cell_weights(std::vector<Cell>* _my_cells);

#endif