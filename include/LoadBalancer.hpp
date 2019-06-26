//
// Created by xetql on 5/3/19.
//

#ifndef SPEC_LOAD_BALANCER_HPP
#define SPEC_LOAD_BALANCER_HPP

#include <vector>
#include <thread>
#include <algorithm>
#include <mpi.h>
#include "Time.hpp"
#include "Utils.hpp"
#include "LoadBalancingApproach.hpp"
#include "StdApproach.hpp"
#include <zupply.hpp>
template<typename Data>
class LoadBalancer {
protected:
    unsigned int pcall = 0;
    std::vector<double> lb_costs;
    MPI_Comm world;
    MPI_Datatype datatype;
    int rank, worldsize;
    std::unique_ptr<const LoadBalancingApproach> approach;
    const zz::log::LoggerPtr logger;

public:

    LoadBalancer(MPI_Comm world, MPI_Datatype datatype, LoadBalancingApproach* approach) :
        world(world),
        datatype(datatype),
        approach(approach),
        logger(zz::log::get_logger("LB", true)){
        MPI_Comm_rank(world, &rank);
        MPI_Comm_size(world, &worldsize);
    };

    std::tuple<int, int, int, int> activate_load_balance(int msx, int msy, int step,
                                                         std::vector<Data>* _data, std::vector<size_t> *_data_pointers){
        if(rank == 0) logger->info("Calling load balancer at step {} with {}", step, this->approach->to_string());

        START_TIMING(current_lb_cost);
        load_balance(_data);
        auto bbox = get_bounding_box(msx, msy, *(_data));
        init_populate_data_pointers(msx, msy, _data_pointers, *(_data), bbox);
        STOP_TIMING(current_lb_cost);

        MPI_Allreduce(MPI_IN_PLACE, &current_lb_cost, 1, MPI_DOUBLE, MPI_MAX, world);

        pcall = (unsigned int) step;
        lb_costs.push_back(current_lb_cost);

        return bbox;
    }

    virtual std::vector<Data> propagate(const std::vector<Data> &data,
                                        int* nb_elements_recv, int* nb_elements_sent,
                                        double cell_size) = 0;

    double get_average_cost() {
        return stats::mean<double>(lb_costs.begin(), lb_costs.end());
    }

    unsigned int get_last_call() {
        return pcall;
    }

    double estimate_best_ncall(double weight, double slope) {
        return sqrt(2.0 * weight * get_average_cost() / slope);
    }

    void set_approach(LoadBalancingApproach* approach) {
        this->approach.reset(approach);
    }

private:
    virtual void load_balance(std::vector<Data>* _data) = 0;
};
#endif //SPEC_LOAD_BALANCER_HPP
