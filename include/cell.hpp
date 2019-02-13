//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_CELL_HPP
#define SPEC_CELL_HPP

#include <ostream>
#include "communication.hpp"
#include "utils.hpp"

struct Cell {
    int gid, type;
    float weight, erosion_probability;

    Cell() : gid(0), type(0), weight(0), erosion_probability(0) {};
    Cell(int gid, int type, float weight, float erosion_probability)
            : gid(gid), type(type), weight(weight), erosion_probability(erosion_probability) {};

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

        MPI_Datatype element_datatype,
                     gid_type_datatype,
                     weight_prob_datatype,
                     oldtype_element[2];
        MPI_Aint offset[2], intex, pos_offset;

        const int number_of_int_elements = 2;
        const int number_of_float_elements = 2;

        int blockcount_element[2];

        //int
        MPI_Type_contiguous(number_of_int_elements, MPI_INT, &gid_type_datatype); // position
        MPI_Type_commit(&gid_type_datatype);
        //float
        MPI_Type_contiguous(number_of_float_elements, MPI_FLOAT, &weight_prob_datatype); // position
        MPI_Type_commit(&weight_prob_datatype);

        blockcount_element[0] = 1; // gid, lid, exit, waiting_time
        blockcount_element[1] = 1; // position <x,y>

        oldtype_element[0] = gid_type_datatype;
        oldtype_element[1] = weight_prob_datatype;

        MPI_Type_extent(gid_type_datatype, &pos_offset);

        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = pos_offset;

        MPI_Type_struct(2, blockcount_element, offset, oldtype_element, &element_datatype);

        MPI_Type_commit(&element_datatype);

        return CommunicationDatatype(element_datatype, gid_type_datatype);
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

#endif //SPEC_CELL_HPP
