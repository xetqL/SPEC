//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_COMM_HPP
#define SPEC_COMM_HPP

#include <mpi.h>

struct CommunicationDatatype {

    MPI_Datatype element_datatype, vec_datatype;

    CommunicationDatatype(const MPI_Datatype &el, const MPI_Datatype &vec) : element_datatype(el), vec_datatype(vec){}

    void free_datatypes() {
        MPI_Type_free(&vec_datatype);
        MPI_Type_free(&element_datatype);
    }
};


#endif //SPEC_COMM_HPP
