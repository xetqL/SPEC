//
// Created by xetql on 5/3/19.
//

#ifndef SPEC_SIMULATED_LBM_HPP
#define SPEC_SIMULATED_LBM_HPP

#include <memory>
#include <mpi.h>

#include "Cell.hpp"
#include "LoadBalancer.hpp"
#include "SimulationParams.hpp"
#include "zupply.hpp"

class SimulatedLBM {
    using LBAlgorithm = LoadBalancer<Cell>;
    const int FOREMAN = 0;

    SimulationParams params;

    MPI_Comm comm;
    std::unique_ptr<LBAlgorithm> load_balancer;

    zz::log::LoggerPtr perflogger, steplogger, proctime;
public:
    SimulatedLBM(SimulationParams params, const MPI_Comm comm, LBAlgorithm *load_balancer);

    virtual void run(float alpha);

    void set_loggers(zz::log::LoggerPtr perflogger, zz::log::LoggerPtr steplogger, zz::log::LoggerPtr proctime) {
        this->perflogger = perflogger;
        this->steplogger = steplogger;
        this->proctime   = proctime;
    }

};

#endif //SPEC_SIMULATED_LBM_HPP
