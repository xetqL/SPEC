//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_CELL_HPP
#define SPEC_CELL_HPP

#include "communication.hpp"

struct Cell {
    int gid, lid, type;
    float erosion_probability;
    Cell() : gid(0), lid(0), type(0), erosion_probability(0) {};
    Cell(int gid, int lid, int type, float erosion_probability)
            : gid(gid), lid(lid), type(type), erosion_probability(erosion_probability) {};

    static CommunicationDatatype register_datatype() {

        MPI_Datatype element_datatype,
                     position_datatype,
                     oldtype_element[2];
        MPI_Aint offset[2], intex, pos_offset;

        const int number_of_int_elements = 3;
        const int number_of_proba_elements = 1;

        int blockcount_element[2];

        // register particle element type
        MPI_Type_contiguous(number_of_int_elements, MPI_INT, &position_datatype); // position
        MPI_Type_commit(&position_datatype);

        blockcount_element[0] = 1; // gid, lid, exit, waiting_time
        blockcount_element[1] = 1; // position <x,y>

        oldtype_element[0] = position_datatype;
        oldtype_element[1] = MPI_FLOAT;

        MPI_Type_extent(MPI_INT, &intex);
        MPI_Type_extent(position_datatype, &pos_offset);

        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = pos_offset;

        MPI_Type_struct(2, blockcount_element, offset, oldtype_element, &element_datatype);

        MPI_Type_commit(&element_datatype);

        return CommunicationDatatype(position_datatype, element_datatype);
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

};

#endif //SPEC_CELL_HPP
