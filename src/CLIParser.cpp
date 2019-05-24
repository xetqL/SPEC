//
// Created by xetql on 5/6/19.
//

#include "CLIParser.hpp"

std::pair<SimulationParams, bool>  parse_cli(int argc, char** argv) {
    unsigned int xprocs, yprocs;
    unsigned int cell_per_process;
    unsigned int MAX_STEP, interval = 100;
    int seed;
    unsigned int N;
    bool load_lattice, verbose;
    float alpha;
    std::string lattice_fname, outputfname;

    zz::cfg::ArgParser parser;
    parser.add_opt_version('V', "version", "0.1");
    parser.add_opt_help('h', "help"); // use -h or --help

    parser.add_opt_value('x', "xprocs", xprocs, (unsigned int) 0, "set the number of PE in x dimension", "INT").require(); // require this option
    parser.add_opt_value('y', "yprocs", yprocs, (unsigned int) 0, "set the number of PE in y dimension", "INT").require(); // require this option
    parser.add_opt_value('c', "cellPerCPU", cell_per_process, (unsigned int) 0, "set the number of cell per CPU must be a perfect square", "INT").require(); // require this option
    parser.add_opt_value('N', "overloading", N, (unsigned int) 0, "set the number of overloading CPU", "INT").require(); // require this option
    parser.add_opt_value('t', "steps", MAX_STEP, (unsigned int) 0,"set the number of PE in y dimension", "INT").require(); // require this option
    parser.add_opt_value('s', "seed" , seed, 0, "set the random seed", "INT");
    parser.add_opt_flag('v', "verbose", "verbosity", &verbose);
    parser.add_opt_flag('l', "loadLattice", "load an external lattice instead of random generation", &load_lattice);
    parser.add_opt_value('f', "filename" , lattice_fname, std::string(""), "load this lattice file", "STRING");
    parser.add_opt_value('a', "alpha" , alpha, 0.1f, "set the alpha-value (unloading model)", "FLOAT");

#ifdef CYCLIC_LOAD_BALANCING
    parser.add_opt_value('I', "interval" , interval, (unsigned int) 0, "load balancing interval (CYCLIC LB)", "INT");
#endif

#ifdef PRODUCE_OUTPUTS
    parser.add_opt_value('o', "output", outputfname, std::string("gids-out.npz"),
                         "produce the resulting lattice in NPZ archive", "STRING");
#endif
    parser.parse(argc, argv);
    if (parser.count_error() > 0) {
        std::cout << parser.get_error() << std::endl;
        std::cout << parser.get_help() << std::endl;
        return std::make_pair<SimulationParams>({}, true);
    } else {
        SimulationParams p = {xprocs, yprocs, N, MAX_STEP, cell_per_process, interval, seed, load_lattice, verbose, lattice_fname, outputfname};
        p.alpha = alpha;
        return {p, false};
   
    }
}
