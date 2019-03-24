#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <zconf.h>
#include <memory>
#include <map>
#include <cassert>
#include <algorithm>
#include <mpi.h>
#include <numeric>
#include <random>
#include <complex>

#ifdef PRODUCE_OUTPUTS
    #include <cnpy.h>
#endif

#include "../include/cell.hpp"
#include "../include/utils.hpp"

#ifdef WITH_ZOLTAN
#include "../include/zoltan_fn.hpp"
#endif

#include "../include/window.hpp"
#include "../include/io.hpp"
#include "../include/CPULoadDatabase.hpp"

#include "../include/time.hpp"
#include "../include/band_partitioner.hpp"

#include "zupply.hpp"
#include "main.hpp"

zz::log::LoggerPtr perflogger, steplogger, proctime;

#define FOREMAN 0

int main(int argc, char **argv) {

    int worldsize;
    int rank;

    MPI_Init(&argc, &argv);
    auto world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldsize);
    MPI_Comm_rank(world, &rank);
    const bool i_am_foreman = rank == FOREMAN;

    auto datatype   = Cell::register_datatype();

    std::string str_rank = "[RANK " + std::to_string(rank) + "] ";
    // Generate data
    //int xprocs = std::atoi(argv[2]), yprocs = std::atoi(argv[1]);
    unsigned int xprocs, yprocs;
    unsigned int cell_per_process;
    unsigned int MAX_STEP;
    int seed;
    unsigned int N;
    bool load_lattice, verbose;
    float alpha;
    std::string lattice_fname, outputfname;
    ///-----------------------------------------------------------------------------------------------------------------
    /// Parser

    zz::log::config_from_file("logger.cfg");
    perflogger = zz::log::get_logger("perf",  true);
    steplogger = zz::log::get_logger("steps", true);
    proctime   = zz::log::get_logger("proctime", true);

    zz::cfg::ArgParser parser;
    parser.add_opt_version('V', "version", "0.1");
    parser.add_opt_help('h', "help"); // use -h or --help

    parser.add_opt_value('x', "xprocs", xprocs, (unsigned int) 0, "set the number of PE in x dimension", "INT").require(); // require this option
    parser.add_opt_value('y', "yprocs", yprocs, (unsigned int) 0, "set the number of PE in y dimension", "INT").require(); // require this option
    parser.add_opt_value('c', "cellPerCPU", cell_per_process, (unsigned int) 0, "set the number of cell per CPU must be a perfect square", "INT").require(); // require this option
    parser.add_opt_value('N', "overloading", N, (unsigned int) 0, "set the number of overloading CPU", "INT").require(); // require this option
    parser.add_opt_value('t', "steps", MAX_STEP, (unsigned int) 0, "set the number of PE in y dimension", "INT").require(); // require this option
    parser.add_opt_value('s', "seed" , seed, 0, "set the random seed", "INT");
    parser.add_opt_flag('v', "verbose", "verbosity", &verbose);
    parser.add_opt_flag('l', "loadLattice", "load an external lattice instead of random generation", &load_lattice);
    parser.add_opt_value('f', "filename" , lattice_fname, std::string(""), "load this lattice file", "STRING");
    parser.add_opt_value('a', "alpha" , alpha, 0.1f, "set the alpha-value (unloading model)", "FLOAT");

#ifdef PRODUCE_OUTPUTS
    parser.add_opt_value('o', "output", outputfname, std::string("gids-out.npz"),
                         "produce the resulting lattice in NPZ archive", "STRING");
#endif
    parser.parse(argc, argv);

    if (parser.count_error() > 0) {
        if(i_am_foreman) {
            std::cout << parser.get_error() << std::endl;
            std::cout << parser.get_help() << std::endl;
        }
        MPI_Abort(world, MPI_ERR_UNKNOWN);
        return EXIT_FAILURE;
    }

    if(i_am_foreman) {
        steplogger->info("xprocs: \t") << xprocs;
        steplogger->info("yprocs: ") << yprocs;
        steplogger->info("cellPerCPU: ") << cell_per_process;
        steplogger->info("N: ") << N;
        steplogger->info("steps: ") << MAX_STEP;
        steplogger->info("seed: ") << seed;
        steplogger->info("alpha: ") << alpha;
    }

    if(verbose) {
        steplogger->info() << N;
    }

    /// finish parsing cli arguments
    ///-----------------------------------------------------------------------------------------------------------------

    std::mt19937 gen(seed); // common seed
    std::uniform_int_distribution<>  proc_dist(0, worldsize-1);

    if(xprocs * yprocs != worldsize) {
        steplogger->fatal() << "Grid size does not match world size! " << (xprocs*yprocs) << " " << worldsize;
        MPI_Abort(world, MPI_ERR_UNKNOWN);
        return EXIT_FAILURE;
    }

    int cell_in_my_rows = (int) std::sqrt(cell_per_process), cell_in_my_cols = cell_in_my_rows;
    int xcells = cell_in_my_rows * xprocs, ycells = cell_in_my_rows * yprocs;
    std::vector<int> cols(ycells);
    std::iota(cols.begin(), cols.end(), 0);

    std::shuffle(cols.begin(), cols.end(), gen);
    cols = {0, 1};
    std::vector<int> water_cols(cols.begin(), cols.begin()+1 );
    std::sort(water_cols.begin(), water_cols.end());
    //std::for_each(water_cols.cbegin(), water_cols.cend(), [](auto v){ std::cout << v << std::endl; });

    int& msx = Cell::get_msx(); msx = xcells;
    int& msy = Cell::get_msy(); msy = ycells;
    int inner_type = WATER_TYPE;
    int shape[2] = {msx,msy};
    std::vector<int> loading_procs;
    //const int loading_proc = proc_dist(gen) % worldsize; //one randomly chosen load proc
    std::uniform_int_distribution<>  lproc_dist((int)N/2, worldsize-1-(int)N/2);

    int median_lproc = lproc_dist(gen);
    if(N >= 1) loading_procs.push_back(median_lproc);
    for(int i = 1; i <= (int) N/2; ++i) {
        if(loading_procs.size() < N) {
            loading_procs.push_back(median_lproc-i);
        } else break;
        if(loading_procs.size() < N) {
            loading_procs.push_back(median_lproc+i);
        } else break;
    }
    std::for_each(loading_procs.begin(), loading_procs.end(), [](auto v){std::cout << v << std::endl;});

    const bool i_am_loading_proc = std::find(loading_procs.begin(), loading_procs.end(), rank) != loading_procs.end();

#ifdef PRODUCE_OUTPUTS
    if(i_am_foreman) cnpy::npz_save("gids-out.npz", "shape", &shape[0],   {2}, "w");
    if(i_am_foreman) cnpy::npz_save("gids-out.npz", "type",  &inner_type, {1}, "a");
#endif

    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_global_position(xprocs, yprocs, rank);
    unsigned long total_cell = cell_per_process * worldsize;

    if(i_am_foreman) {
        steplogger->info("CPU COUNT:")    << worldsize;
        steplogger->info("GRID PSIZE X:") << xprocs;
        steplogger->info("GRID PSIZE Y:") << yprocs;
        steplogger->info("GRID SIZE  X:") << msx;
        steplogger->info("GRID SIZE  Y:") << msy;
        steplogger->info("EACH SIZE  X:") << xcells;
        steplogger->info("EACH SIZE  Y:") << ycells;
    }

    if(i_am_foreman) steplogger->info() << cell_in_my_cols << " " << cell_in_my_rows;
    std::vector<Cell> my_cells;
#ifndef PRODUCE_OUTPUTS
    //my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);
    my_cells = generate_lattice_single_type(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, WATER_TYPE, 1.0, 0.0);
#else
    if(load_lattice) {
        my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols, lattice_fname);
    } else {
        my_cells = generate_lattice_single_type(msx, msy, x_proc_idx, y_proc_idx, cell_in_my_cols,cell_in_my_rows, WATER_TYPE, 1.0, 0.0);
    }
#endif
    const int my_cell_count = my_cells.size();

    int recv, sent;
    // if(i_am_foreman) steplogger->info("End of map generation");
    /* Initial load balancing */
    std::vector<double> lb_costs;
    // auto zoltan_lb = zoltan_create_wrapper(true, world);
    StripeLoadBalancer stripe_lb(FOREMAN, msx, msy, datatype.element_datatype, world);
    stripe_lb.setLoggerPtr(steplogger);
    PAR_START_TIMING(current_lb_cost, world);
    stripe_lb.load_balance(&my_cells, 0.0);
    PAR_STOP_TIMING(current_lb_cost, world);

    MPI_Allreduce(&current_lb_cost, &current_lb_cost, 1, MPI_DOUBLE, MPI_MAX, world);
    lb_costs.push_back(current_lb_cost);
    auto avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());

    auto my_domain = stripe_lb.get_domain(rank);

    generate_lattice_rocks(1, msx, msy, &my_cells, i_am_loading_proc ? 0.5f : 0.01f, my_domain.first, my_domain.second);

#ifdef PRODUCE_OUTPUTS
    std::vector<std::array<int,2>> all_types(total_cell);
    std::vector<std::array<int,2>> my_types(my_cell_count);
    for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
    gather_elements_on(my_types, FOREMAN, &all_types, datatype.minimal_datatype, world);
    if(i_am_foreman) {
        std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
        std::vector<int> types, water_gid;
        std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
        std::for_each(all_types.begin(), all_types.end(), [&water_gid, &inner_type](auto e){if(e[1] == inner_type) water_gid.push_back(e[0]);});
        assert((*(all_types.end() - 1))[0] == total_cell-1);
        cnpy::npz_save("gids-out.npz", "step-"+std::to_string(0), &water_gid[0], {water_gid.size()}, "a");
    }
    if(i_am_foreman) all_types.resize(total_cell);
#endif

    if(i_am_foreman) perflogger->info("LB_time: ") << current_lb_cost;

    std::vector<unsigned long> my_water_ptr;
    int n;
    std::vector<size_t> data_pointers, remote_data_pointers;
    std::tie(n, my_water_ptr) = create_water_ptr_vector(my_cells);
    std::tuple<int, int, int, int> bbox; // = add_to_bbox(msx, msy, get_bounding_box(my_cells), -10, 10, -10, 10);
    //populate_data_pointers(msx, msy, &data_pointers, my_cells, 0, bbox, true);
    /* lets make it fun now...*/
    CPULoadDatabase gossip_workload_db(worldsize,  2, 9999, world),
                    gossip_waterslope_db(worldsize,  2, 8888, world);

    unsigned int ncall = 10, pcall=0;
#ifndef LB_METHOD
    ncall = 0;
#endif
#if LB_METHOD==1
    ncall = 25;
#endif

    double skew = 0, degradation_since_last_lb = 0.0;

    std::vector<double> timings(worldsize), all_degradations, water;

    SlidingWindow<double> window_step_time(15);  // sliding window with max size = TODO: tune it?
    SlidingWindow<double> window_water(ncall);   // sliding window with max size = TODO: tune it?
    SlidingWindow<double> window_my_time(100);   // sliding window with max size = TODO: tune it?

    water.push_back(my_water_ptr.size()); window_water.add(my_water_ptr.size());

    //double time_since_start;
    PAR_START_TIMING(loop_time, world);
    for(unsigned int step = 0; step < MAX_STEP; ++step) {
        if(i_am_foreman) steplogger->info() << "Beginning step "<< step;

        PAR_START_TIMING(step_time, world);
        bool lb_condition = false;
#ifdef LB_METHOD
#if LB_METHOD==1   // load balance every 100 iterations
        lb_condition = (pcall + ncall) == step;
        if(lb_condition) {
            if(i_am_foreman) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            stripe_lb.load_balance(&my_cells, 0.0);
            PAR_STOP_TIMING(current_lb_cost, world);
            if(i_am_foreman) perflogger->info("LB_time: ") << current_lb_cost;
            lb_costs.push_back(current_lb_cost);
            ncall = 100;
            if(i_am_foreman) steplogger->info("next LB call at: ") << (step+ncall);
            pcall = step;
        }
#elif LB_METHOD == 2
        // http://sc16.supercomputing.org/sc-archive/tech_poster/poster_files/post247s2-file3.pdf +
        if(i_am_foreman) steplogger->info("degradation method 2: ") << (degradation_since_last_lb*(step-pcall))/2.0 << " avg_lb_cost " << avg_lb_cost;
        lb_condition = pcall + ncall <= step;// || (degradation_since_last_lb*(step-pcall))/2.0 > avg_lb_cost;
        if(lb_condition) {
            auto total_slope = get_slope<double>(window_step_time.data_container);
            if(i_am_foreman) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            stripe_lb.load_balance(&my_cells, 0.0);
            PAR_STOP_TIMING(current_lb_cost, world);
            MPI_Allreduce(&current_lb_cost, &current_lb_cost, 1, MPI_DOUBLE, MPI_MAX, world);
            lb_costs.push_back(current_lb_cost);
            if(i_am_foreman) perflogger->info("LB_time: ") << current_lb_cost;
            avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());
            if(total_slope > 0) {
                ncall = (unsigned int) std::floor(std::sqrt((2.0 * avg_lb_cost) / total_slope));
                ncall = std::max((unsigned int) 1, ncall);
            } else ncall = MAX_STEP;

            gossip_workload_db.reset();
            water.clear();
            degradation_since_last_lb = 0.0;
            window_my_time.data_container.clear();
            window_step_time.data_container.clear();
            window_water.data_container.clear();
            pcall = step;

            if(i_am_foreman) steplogger->info("next LB call at: ") << (step+ncall);

            MPI_Bcast(&ncall, 1, MPI_INT, 0, world);
        }
#elif LB_METHOD == 3
        // http://sc16.supercomputing.org/sc-archive/tech_poster/poster_files/post247s2-file3.pdf +
        // http://ics2018.ict.ac.cn/essay/ics18-final62.pdf
        auto total_slope = get_slope<double>(window_step_time.data_container);
        if(i_am_foreman) steplogger->info("degradation method 3: ") << degradation_since_last_lb << " avg_lb_cost " << avg_lb_cost << " total slope: " << total_slope;

        lb_condition = pcall + ncall <= step || degradation_since_last_lb > avg_lb_cost;
        if(lb_condition) {

            if(i_am_foreman) steplogger->info("call LB at: ") << step;
            std::cout << lb_condition << std::endl;
            PAR_START_TIMING(current_lb_cost, world);
            stripe_lb.load_balance(&my_cells, 0.0);
            PAR_STOP_TIMING(current_lb_cost, world);
            MPI_Allreduce(&current_lb_cost, &current_lb_cost, 1, MPI_DOUBLE, MPI_MAX, world);
            lb_costs.push_back(current_lb_cost);
            perflogger->info("LB_time: ") << current_lb_cost;
            avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());

            if(total_slope > 0) {
                if(i_am_foreman) steplogger->info("DEBUG")<< (2.0 * avg_lb_cost) << "/" << total_slope;
                ncall = (unsigned int) std::floor(std::sqrt((2.0 * avg_lb_cost) / total_slope));
                ncall = std::max((unsigned int) 1, ncall);
            } else ncall = MAX_STEP;

            gossip_workload_db.reset();
            water.clear();
            degradation_since_last_lb = 0.0;
            window_my_time.data_container.clear();
            window_step_time.data_container.clear();
            window_water.data_container.clear();
            pcall = step;

            if(i_am_foreman) steplogger->info("next LB call at: ") << (step+ncall);

            MPI_Bcast(&ncall, 1, MPI_INT, 0, world);
        }
#elif LB_METHOD == 4 // Unloading Model
        constexpr double alpha = 1.0/4.0;
        if(i_am_foreman) steplogger->info("degradation: ") << degradation_since_last_lb << " avg_lb_cost " << avg_lb_cost;
        lb_condition = pcall + ncall <= step || degradation_since_last_lb > avg_lb_cost;
        if( lb_condition ) {
            double my_slope           = (std::floor(get_slope<double>(window_my_time.data_container)*1000.0)) / 1000.0;
            double step_slope         = (std::floor(get_slope<double>(window_step_time.data_container)*1000.0)) / 1000.0;
            auto total_last_step_time = *(window_step_time.end()-1);
            auto my_last_step_time    = *(window_my_time.end()-1);

            PAR_START_TIMING(current_lb_cost, world);
            // if(i_am_overloaded && i_am_increasing) {
            // std::cout << rank << " !!!!!!!!!!!!!!!!!!! " << step_slope << " " << my_slope << std::endl;
            int parts_num[1] = {rank}, weight_per_obj[1] = {0};
            float part_size[1] = {(1.0f - (float) my_slope)};
            // Zoltan_LB_Set_Part_Sizes(zoltan_lb, 1, 1, parts_num, weight_per_obj, part_size);
            update_cell_weights(&my_cells, part_size[0], WATER_TYPE, [](auto a, auto b){return a * 1.0/b;});
            /*} else {
                int parts_num[1] = {rank}, weight_per_obj[1] = {0};
                float part_size[1] = {1.0};
                Zoltan_LB_Set_Part_Sizes(zoltan_lb, 1, 1, parts_num, weight_per_obj, part_size);
            }*/
            // auto data = zoltan_load_balance(&my_cells, zoltan_lb, true, true);
            // std::cout << rank << " " << data << std::endl;
            PAR_STOP_TIMING(current_lb_cost, world);
            if(i_am_foreman) perflogger->info("LB_time: ") << current_lb_cost;
            lb_costs.push_back(current_lb_cost);
            avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());
            //std::cout << rank << " number of cells: " << compute_effective_workload(my_cells, WATER_TYPE) << " effective load " << compute_estimated_workload(my_cells) << std::endl ;
            degradation_since_last_lb = 0.0;
            window_my_time.data_container.clear();
            window_step_time.data_container.clear();
            pcall = step;
            ncall = 10;
        }
#elif LB_METHOD == 5
        if(i_am_foreman) steplogger->info("degradation method 4: ") << degradation_since_last_lb << " avg_lb_cost " << avg_lb_cost;
        lb_condition = pcall + ncall <= step || degradation_since_last_lb > avg_lb_cost*1.1;
        if(lb_condition) {
            bool overloading = gossip_waterslope_db.zscore(rank) > 3.0;
            if(overloading) std::cout << "I WILL BE UNLOADED" << std::endl;

            PAR_START_TIMING(current_lb_cost, world);
            stripe_lb.load_balance(&my_cells, overloading ? alpha : 0.0);
            PAR_STOP_TIMING(current_lb_cost, world);
            MPI_Allreduce(&current_lb_cost, &current_lb_cost, 1, MPI_DOUBLE, MPI_MAX, world);
            lb_costs.push_back(current_lb_cost);
            avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());
            perflogger->info("LB_time: ") << current_lb_cost;

            gossip_workload_db.reset();
            water.clear();
            degradation_since_last_lb = 0.0;
            window_my_time.data_container.clear();
            window_step_time.data_container.clear();
            window_water.data_container.clear();
            pcall = step;
            ncall = MAX_STEP;
        }
#endif
        if(lb_condition) {
            std::tie(n, my_water_ptr) = create_water_ptr_vector(my_cells);
            water.push_back(my_water_ptr.size());
            window_water.add(my_water_ptr.size());
        }
#endif
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// COMPUTATION START

        PAR_START_TIMING(comp_time, world);

        auto remote_cells = stripe_lb.share_frontier_with_neighbors(my_cells, &recv, &sent);//zoltan_exchange_data(zoltan_lb,my_cells,&recv,&sent,datatype.element_datatype,world,1.0);
        decltype(my_water_ptr) remote_water_ptr;

        std::tie(std::ignore, remote_water_ptr) = create_water_ptr_vector(remote_cells);
        decltype(my_water_ptr) new_water_ptr;

        if(lb_condition || step == 0) bbox = get_bounding_box(my_cells, remote_cells);
        populate_data_pointers(msx, msy, &data_pointers, my_cells, remote_cells, bbox, lb_condition || step == 0);

        std::tie(my_cells, new_water_ptr) = dummy_erosion_computation3(msx, msy, my_cells, my_water_ptr, remote_cells, remote_water_ptr, data_pointers, bbox);
        my_water_ptr.insert(my_water_ptr.end(), std::make_move_iterator(new_water_ptr.begin()), std::make_move_iterator(new_water_ptr.end()));
        n += 4 * new_water_ptr.size(); // adapt the number of cell to compute

        water.push_back(n);
        window_water.add(n);

        CHECKPOINT_TIMING(comp_time, my_comp_time);
        PAR_STOP_TIMING(comp_time, world);
        PAR_STOP_TIMING(step_time, world);
        CHECKPOINT_TIMING(loop_time, time_since_start);
        PAR_STOP_TIMING(loop_time, world);

        if(i_am_foreman) steplogger->info("time for step ") << step << " = " << step_time<< " time for comp. = "<< comp_time << " total: " << time_since_start;
        //if(i_am_loading_proc) perflogger->info(str_rank.c_str())<< " is loading with a current load of " << n;

        MPI_Allreduce(&comp_time, &comp_time, 1, MPI_DOUBLE, MPI_MAX, world); // i should not need that!

        window_step_time.add(comp_time);  // monitor evolution of load in time with a window
        window_my_time.add(my_comp_time); // monitor evolution of my load in time with a window

        if(worldsize > 2)
        {
            if(step > 0)
            {
                gossip_waterslope_db.finish_gossip_step();
                gossip_workload_db.finish_gossip_step();
            }
            gossip_waterslope_db.gossip_update(rank, get_slope<double>(water.begin(), water.end()));
            gossip_waterslope_db.gossip_propagate();
            gossip_workload_db.gossip_update(rank, my_comp_time);
            gossip_workload_db.gossip_propagate();
        }

        if(pcall + 1 < step) {
            degradation_since_last_lb +=
                    (stats::median<double>(window_step_time.end()-3, window_step_time.end())
                            - stats::mean<double>(window_step_time.begin(), window_step_time.end()));
            //std::for_each(window_step_time.newest(), window_step_time.window_step_time.newest()-2(), [](auto v){ std::cout << v << std::endl; });
        }
        /// COMPUTATION STOP
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<double> exch_timings(worldsize);
        std::vector<double> slopes(worldsize);
        std::vector<int> tloads(worldsize);

        MPI_Gather(&my_comp_time, 1, MPI_DOUBLE, &timings.front(), 1, MPI_DOUBLE, FOREMAN, world);
        if(i_am_foreman) {
            double max = *std::max_element(timings.cbegin(), timings.cend()),
                   average = std::accumulate(timings.cbegin(), timings.cend(), 0.0) / worldsize,
                   load_imbalance = (max / average - 1.0) * 100.0;
            perflogger->info("\"step\":") << step << ",\"LI\": " << load_imbalance;
            proctime->info("\"step\":") << step << ",\"proctime\": " << timings;
        }
#ifdef PRODUCE_OUTPUTS
        if(step % 10 == 0){
            unsigned long cell_cnt = my_cells.size();
            std::vector<std::array<int,2>> my_types(cell_cnt);
            for (unsigned int i = 0; i < cell_cnt; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
            gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
            if(i_am_foreman) {
                std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
                std::vector<int> types, rock_gid;
                std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
                std::for_each(all_types.begin(), all_types.end(), [&rock_gid, &inner_type](auto e){if(e[1] == inner_type) rock_gid.push_back(e[0]);});
                cnpy::npz_save("gids-out.npz", "step-"+std::to_string(step+1), &rock_gid[0], {rock_gid.size()}, "a");
            }
        }
#endif
        if(i_am_foreman) steplogger->info() << "Stop step "<< step;
        PAR_RESTART_TIMING(loop_time, world);
    }
    PAR_STOP_TIMING(loop_time, world);
    if(i_am_foreman) perflogger->info("\"total_time\":") << loop_time;
    // if(i_am_foreman) steplogger->info("\"total_time\":") << loop_time;
    datatype.free_datatypes();
    MPI_Finalize();
    return 0;
}
