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
#include "../src/zupply.hpp"
#include <random>


class StripeLoadBalancer {
    typedef unsigned int PEIndex;
    int sizeX, sizeY;
    int worldsize, myrank;
    MPI_Comm world;
    MPI_Datatype datatype;
    std::vector<PEIndex> partition;
    zz::log::LoggerPtr loggerPtr = nullptr;
public:
    StripeLoadBalancer(int &X, int &Y, const MPI_Datatype &datatype, const MPI_Comm& world) : sizeX(X), sizeY(Y), datatype(datatype), world(world) {
        MPI_Comm_rank(world, &myrank);
        MPI_Comm_size(world, &worldsize);
        partition.resize(2*worldsize);
    }

    void setLoggerPtr(const zz::log::LoggerPtr &loggerPtr) {
        StripeLoadBalancer::loggerPtr = loggerPtr;
    }

    /**
     * Procedure used to re-balance the mesh onto all the processes.
     * The mesh is spatially partitionned into stripes, each stripes having the same workload.
     * Alpha is as an extra parameter to tell then partitioner that (alpha*100)% of the total load of the process
     * should be "scattered".
     *
     * Implementation is as follows:
     * 1. Data are gather onto process 0
     * 2. P0 computes the load of each row in the mesh
     * 3. P0 attributes the same load minus (alpha*100)% (greedy algorithm) to each process, producing "stripes".
     * 4. P0 broadcasts the partitioning among processes. Done.
     * /!\ It is a cyclic process, p0 will always get the top part of the mesh etc.
     * @param _data
     * @param alpha
     */
    void load_balance(std::vector<Cell>* _data, double alpha = 0.0) {
        if(worldsize == 1) return;
        std::vector<Cell>& data = *_data;

        // Gather
        auto mesh   = gather_elements_on(data,    0, datatype,   world);
        auto alphas = gather_elements_on({alpha}, 0, MPI_DOUBLE, world);

        // Partition
        execute_partitioning_algorithm(mesh, alphas);

        // Mapping
        execute_mapping_algorithm(&mesh);

        // Affectation
        _data->assign(mesh.begin(), mesh.end());

    }

    std::vector<Cell> share_frontier_with_neighbors(const std::vector<Cell> &data,
                                                    int* nb_elements_recv,
                                                    int* nb_elements_sent,
                                                    double cell_size = 1.0) {

        std::vector<Cell> remote_data;
        if(worldsize == 1) return remote_data;
        std::vector<std::vector<Cell>> data_to_migrate(2);

        auto neighbors = get_my_neighbors();
        auto my_top_frontier    = partition[2 * myrank];
        auto my_bottom_frontier = partition[2 * myrank + 1];

        for(const auto& cell : data) {
            long row = cell_to_position(sizeX, sizeY, cell.gid).second;
            if(row == my_top_frontier && neighbors.first > -1) {
                data_to_migrate[0].push_back(cell);
            } else if (row == my_bottom_frontier && neighbors.second > -1) {
                data_to_migrate[1].push_back(cell);
            }
        }

        MPI_Request reqs[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

        if(neighbors.first > -1) {
            MPI_Isend(&(data_to_migrate[0].front()), data_to_migrate[0].size(), datatype, neighbors.first, 547, world,
                      &reqs[0]);
        }
        if(neighbors.second > -1) {
            MPI_Isend(&(data_to_migrate[1].front()), data_to_migrate[1].size(), datatype, neighbors.second, 547, world,
                      &reqs[1]);
        }

        std::vector<Cell> buffer(sizeX);
        if(neighbors.first > -1) {
            MPI_Recv(&buffer.front(), sizeX, datatype, neighbors.first, 547, world, MPI_STATUS_IGNORE);
            remote_data.insert(remote_data.end(), std::make_move_iterator(buffer.begin()), std::make_move_iterator(buffer.end()));
        }
        if(neighbors.second > -1) {
            MPI_Recv(&buffer.front(), sizeX, datatype, neighbors.second,547, world, MPI_STATUS_IGNORE);
            remote_data.insert(remote_data.end(), std::make_move_iterator(buffer.begin()), std::make_move_iterator(buffer.end()));
        }

        MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

        return remote_data;

    }

    std::pair<unsigned int, unsigned int> get_domain(PEIndex rank) {
        auto top_frontier    = partition[2 * rank];
        auto bottom_frontier = partition[2 * rank + 1];

        return std::make_pair(top_frontier, bottom_frontier);
    };
private:
    /**
     * Partition data into stripe of equal workload (minus alpha*100%).
     * Broadcast the resulting partition to the world.
     * @param mesh
     * @param alphas
     */
    void execute_partitioning_algorithm(const std::vector<Cell>& mesh, std::vector<double> alphas){
        if(!myrank) {
            auto rows_load = get_rows_load(mesh);
            auto total_workload = std::accumulate(rows_load.cbegin(), rows_load.cend(), 0.0);
            auto average_workload = total_workload / (double) worldsize;

            // Count the number of processes that load up.
            auto N = std::count_if(alphas.cbegin(), alphas.cend(), [](auto a){return a > 0.0;});

            if(N > worldsize / 2) return;

            double alpha = *std::max_element(alphas.begin(), alphas.end());

            // Compute all the workload for each PEs
            decltype(alphas) desired_workloads, effective_workloads(worldsize);

            std::transform(alphas.cbegin(), alphas.cend(), std::back_inserter(desired_workloads),
                           [&](auto a) { return a > 0.0 ? average_workload * (1.0 - alpha) : (1.0 + (alpha*N)/(worldsize-N)) * average_workload; } );

            std::for_each(desired_workloads.begin(), desired_workloads.end(), [](auto v){std::cout << v << std::endl;});

            // Greedy algorithm to compute the stripe for each process
            unsigned int begin_stripe=0, end_stripe=0;
            for(PEIndex p = 0; p < worldsize; ++p)
            {
                double current_process_workload = 0.0;

                // While we have not ~reached~ the desired workload or we reached the end of the mesh
                while( ((current_process_workload < desired_workloads[p]) || p == worldsize-1) && end_stripe < sizeY) {
                    current_process_workload += rows_load[end_stripe];
                    end_stripe++;
                }

                assert(begin_stripe != end_stripe);
                assert(current_process_workload > 0.0);

                effective_workloads[p] = current_process_workload;
                partition[2 * p]   = begin_stripe;
                partition[2 * p + 1] = end_stripe - 1;
                std::cout << p << "  " << begin_stripe << " -> " << (end_stripe-1) << " = " << current_process_workload << std::endl;
                begin_stripe = end_stripe; // begin at next stripe
                end_stripe = begin_stripe;
            }
            double desired   = std::accumulate(desired_workloads.begin(),   desired_workloads.end(), 0.0);
            double effective = std::accumulate(effective_workloads.begin(), effective_workloads.end(), 0.0);
            assert(almost_equal(desired, effective, 2));
        }

        // Broadcast the new partition
        MPI_Bcast(&partition.front(), 2*worldsize, MPI_INT, 0, world);
    }

    void execute_mapping_algorithm(std::vector<Cell>* _mesh){
        std::vector<Cell>& mesh = *_mesh;
        if(!myrank) {
            size_t data_id = 0;
            PEIndex PE;
            std::vector<std::vector<Cell>> data_to_migrate(worldsize);
            while (data_id < mesh.size()) {
                resolve_process_membership(mesh[data_id].gid, &PE);
                if (PE != myrank) {
                    std::iter_swap(mesh.begin() + data_id, mesh.end() - 1);
                    data_to_migrate.at(PE).push_back(*(mesh.end()- 1));
                    mesh.pop_back();
                } else data_id++; //if the element must stay with me then check the next one
            }
            PE = 0;
            std::vector<MPI_Request> reqs(worldsize, MPI_REQUEST_NULL);
            for(auto& buf_to_send : data_to_migrate){
                if(PE != myrank){
                    auto send_size = buf_to_send.size();
                    if (send_size) {
                        MPI_Isend(&buf_to_send.front(), send_size, datatype, PE, 300, world, &reqs[PE]);
                    }
                }
                PE++;
            }
            MPI_Waitall(worldsize, &reqs.front(), MPI_STATUSES_IGNORE);
        } else {
            MPI_Status status;
            MPI_Probe(0, 300, world, &status);
            int size;
            MPI_Get_count(&status, datatype, &size);
            mesh.resize(size);
            MPI_Recv(&mesh.front(), size, datatype, 0, 300, world, MPI_STATUS_IGNORE);
        }
    }

    void resolve_process_membership(unsigned int gid, PEIndex* p){
        auto row = cell_to_position(sizeX, sizeY, gid).second;
        for((*p) = 0; (*p) < worldsize; ++(*p)) {
            auto b = partition[2 * (*p)];
            auto e = partition[2 * (*p) + 1];
            if(b <= row && row <= e) return;
        }
        std::cout << gid <<" " <<row << std::endl;
        if(loggerPtr != nullptr) {
            loggerPtr->error("error in resolving process membership.");
        }
        MPI_Abort(world,  MPI_ERR_OTHER);
    }


    /**
     * Compute the average load of a given row. It is used to determine the load of a stripe.
     * @param data the mesh
     * @return a vector containing the load of each row
     */
    const std::vector<double> get_rows_load(const std::vector<Cell>& data) {
        std::vector<double> rowsload(sizeY, 0.0);
        for(const auto& cell : data) {
            long row = cell_to_position(sizeX, sizeY, cell.gid).second;
            rowsload[row] += cell.type;
            /*if(cell.average_load > 0 && cell.type)
                rowsload[row] += cell.average_load;
            else
                rowsload[row] += cell.weight;
                */
        }
        return rowsload;
    }

    /**
     * Get the process id of neighboring processes (i.e., the processes that owns my begin_stripe-1 or my end_stripe+1)
     * Only one of the two returned process id can be negative (i.e., there is no neighbor)
     * @return a pair of process ids
     */
    std::pair<int,int> get_my_neighbors(){
        auto my_top_frontier    = partition[2 * myrank];
        auto my_bottom_frontier = partition[2 * myrank + 1];
        int bottom_neighbor = -1, top_neighbor = -1;
        for(int p = 0; p < worldsize; ++p){
            auto his_top_frontier    = partition[2 * p];
            auto his_bottom_frontier = partition[2 * p + 1];
            if(his_bottom_frontier + 1  == my_top_frontier) {
                top_neighbor = p;
            } else if(my_bottom_frontier + 1  == his_top_frontier) {
                bottom_neighbor = p;
            }
            if(top_neighbor > -1 && bottom_neighbor > -1) break;
        }
        assert(top_neighbor != bottom_neighbor);
        assert(bottom_neighbor != -1 || top_neighbor != -1);
        return std::make_pair(top_neighbor, bottom_neighbor);
    };

};

#endif //SPEC_BAND_PARTITIONER_HPP
