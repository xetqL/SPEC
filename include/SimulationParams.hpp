//
// Created by xetql on 5/3/19.
//

#ifndef SPEC_SIMULATIONPARAMS_HPP
#define SPEC_SIMULATIONPARAMS_HPP

#include <string>
#include <cmath>
struct SimulationParams {
    unsigned int xprocs              = 0;
    unsigned int yprocs              = 0;
    unsigned int N                   = 0;
    unsigned int MAX_STEP            = 0;
    unsigned int cell_per_process    = 0;
    unsigned int interval = MAX_STEP;
    unsigned int xcells, ycells;
    int seed = 0;
    bool load_lattice = false;
    bool verbose      = false;
    float alpha       = 0.0;
    float zthreshold  = 0.0;
    float strong_erosion_probability = 0.0;
    std::string filename, outputfname;


    SimulationParams() {}

    SimulationParams(const unsigned int xprocs, const unsigned int yprocs, const unsigned int N,
                     const unsigned int MAX_STEP, const unsigned int cell_per_process,
                     const unsigned int interval, const int seed,
                     const bool load_lattice, const bool verbose,
                     const std::string &filename, const std::string& outputfname) :
            xprocs(xprocs), yprocs(yprocs), N(N), MAX_STEP(MAX_STEP), cell_per_process(cell_per_process), interval(interval), seed(seed),
            load_lattice(load_lattice), verbose(verbose), filename(filename), outputfname(outputfname) {
        xcells = xprocs * (int) std::sqrt(cell_per_process);
        ycells = yprocs * (int) std::sqrt(cell_per_process);
    }
};

#endif //SPEC_SIMULATIONPARAMS_HPP
