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

template<typename Data>
class LoadBalancer {
protected:
    unsigned int pcall = 0;
    std::vector<double> lb_costs;
    MPI_Comm world;
    MPI_Datatype datatype;
public:
    LoadBalancer(MPI_Comm world, MPI_Datatype datatype) : world(world), datatype(datatype){};

    void activate_load_balance(int step, std::vector<Data>* _data, double alpha) {
        PAR_START_TIMING(current_lb_cost, world);
        load_balance(_data, alpha);
        PAR_STOP_TIMING(current_lb_cost, world);
        MPI_Allreduce(&current_lb_cost, &current_lb_cost, 1, MPI_DOUBLE, MPI_MAX, world);
        lb_costs.push_back(current_lb_cost);
        pcall = (unsigned int) step;
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

private:
    virtual void load_balance(std::vector<Data>* _data, double alpha) = 0;
};
#endif //SPEC_LOAD_BALANCER_HPP
