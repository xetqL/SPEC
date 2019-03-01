//
// Created by xetql on 3/1/19.
//

#ifndef SPEC_BAND_PARTITIONER_HPP
#define SPEC_BAND_PARTITIONER_HPP

#include <mpi.h>
#include <vector>
#include <algorithm>
#include "communication.hpp"
#include "cell.hpp"


class StripePartitioner {
    int sizeX, sizeY;
    int worldsize, myrank;
    MPI_Comm world;
    MPI_Datatype datatype;
public:
    StripePartitioner(int &X, int &Y, const MPI_Datatype &datatype, const MPI_Comm& world) : sizeX(X), sizeY(Y), world(world) {
        MPI_Comm_rank(world, &myrank);
        MPI_Comm_size(world, &worldsize);
    }
    void load_balance(std::vector<Cell>* _data, double &my_slope) {
        std::vector<Cell>& data = *_data;
        auto mesh   = gather_elements_on(data, 0, datatype, world);
        auto slopes = gather_elements_on({my_slope}, 0, datatype, world);
        std::sort(mesh.begin(), mesh.end(), [](auto a,auto b) { return a.gid < b.gid;} );
        auto stripes_load = get_stripes_load(data);

    }
private:
    const std::vector<double> get_stripes_load(const std::vector<Cell>& data) {
        std::vector<double> stripesload(sizeY, 0.0);
        for(const auto& cell : data) {
            long stripe = cell_to_position(sizeX, sizeY, cell.gid).second;
            stripesload[stripe] += cell.average_load;
        }
        return stripesload;
    }

    std::vector<std::vector<int>> partition_and_map(const std::vector<int>& desired_workload) {

    }
};

#endif //SPEC_BAND_PARTITIONER_HPP
