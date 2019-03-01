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
#include "zupply.hpp"
#include "../include/zoltan_fn.hpp"
#include "../include/window.hpp"
#include "../include/io.hpp"
#include "../include/CPULoadDatabase.hpp"
#include "main.hpp"
#include "../include/time.hpp"


zz::log::LoggerPtr perflogger, steplogger;

int main(int argc, char **argv) {
    int worldsize;
    int rank;

    MPI_Init(&argc, &argv);
    auto world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldsize);
    MPI_Comm_rank(world, &rank);
    std::string str_rank = "[RANK " + std::to_string(rank) + "] ";
    // Generate data
    int xprocs = std::atoi(argv[2]), yprocs = std::atoi(argv[1]);
    int cell_per_process = std::atoi(argv[3]);
    const int MAX_STEP = std::atoi(argv[4]);
    const int seed = std::atoi(argv[5]);

    zz::log::config_from_file("logger.cfg");
    perflogger = zz::log::get_logger("perf",  true);
    steplogger = zz::log::get_logger("steps", true);

    if(xprocs * yprocs != worldsize) {
        steplogger->fatal() << "Grid size does not match world size";
        MPI_Abort(world, 1);
        return EXIT_FAILURE;
    }

    int cell_in_my_rows = (int) std::sqrt(cell_per_process), cell_in_my_cols = cell_in_my_rows;
    int xcells = cell_in_my_rows * xprocs, ycells = cell_in_my_rows * yprocs;
    std::vector<int> cols(ycells);
    std::iota(cols.begin(), cols.end(), 0);

    std::mt19937 gen(seed); ///predictable sequence of value

    std::shuffle(cols.begin(), cols.end(), gen);
    cols = {0, 1};
    std::vector<int> water_cols(cols.begin(), cols.begin()+1 );
    std::sort(water_cols.begin(), water_cols.end());
    //std::for_each(water_cols.cbegin(), water_cols.cend(), [](auto v){ std::cout << v << std::endl; });

    int& msx = Cell::get_msx(); msx = xcells;
    int& msy = Cell::get_msy(); msy = ycells;

    int shape[2] = {msx,msy};

#ifdef PRODUCE_OUTPUTS
    if(!rank) cnpy::npz_save("gids-out.npz", "shape", &shape[0], {2}, "w");
#endif

    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_global_position(xprocs, yprocs, rank);
    unsigned long total_cell = cell_per_process * worldsize;
    auto datatype   = Cell::register_datatype();
    if(!rank) {
        steplogger->info("CPU COUNT:")    << worldsize;
        steplogger->info("GRID PSIZE X:") << xprocs;
        steplogger->info("GRID PSIZE Y:") << yprocs;
        steplogger->info("GRID SIZE  X:") << msx;
        steplogger->info("GRID SIZE  Y:") << msy;
        steplogger->info("EACH SIZE  X:") << xcells;
        steplogger->info("EACH SIZE  Y:") << ycells;
    }
    if(!rank) steplogger->info() << cell_in_my_cols << " " << cell_in_my_rows;
    std::vector<Cell> my_cells;
#ifndef PRODUCE_OUTPUTS
    my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);
#else
    if(argc == 7) {
        auto lattice_fname = argv[6];
        my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols, lattice_fname);
    } else {
        my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);
    }
#endif

    const int my_cell_count = my_cells.size();

#ifdef PRODUCE_OUTPUTS
    std::vector<std::array<int,2>> all_types(total_cell);
    std::vector<std::array<int,2>> my_types(my_cell_count);
    for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
    gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
    if(!rank) {
        std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
        std::vector<int> types, water_gid;
        std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
        std::for_each(all_types.begin(), all_types.end(), [&water_gid](auto e){if(!e[1]) water_gid.push_back(e[0]);});
        assert((*(all_types.end() - 1))[0] == total_cell-1);
        cnpy::npz_save("gids-out.npz", "step-"+std::to_string(0), &water_gid[0], {water_gid.size()}, "a");
    }
    if(!rank) all_types.resize(total_cell);
#endif

    int recv, sent;
    if(!rank) steplogger->info() << "End of map generation";

    /* Initial load balancing */

    std::vector<double> lb_costs;
    auto zoltan_lb = zoltan_create_wrapper(true, world);

    PAR_START_TIMING(current_lb_cost, world);
    zoltan_load_balance(&my_cells, zoltan_lb, true, true);
    PAR_STOP_TIMING(current_lb_cost, world);
    if(!rank) perflogger->info("LB_time: ") << current_lb_cost;

    auto my_water_ptr = create_water_ptr_vector(my_cells);
    //lb_costs.push_back(current_lb_cost/2.0);

    /* lets make it fun now...*/
    int minx, maxx, miny, maxy;
    std::vector<size_t> data_pointers;

#ifdef LB_METHOD
    double avg_lb_cost = 100.0; //stats::mean<double>(lb_costs.begin(), lb_costs.end());
    auto tau = 25; // ???
#endif

    SlidingWindow<double> window_step_time(100); //sliding window with max size = TODO: tune it?
    SlidingWindow<double> window_my_time(100); //sliding window with max size = TODO: tune it?

    int ncall = 10, pcall=0;
#if LB_METHOD==1
    ncall = 25;
#endif

    float average_load_at_lb;
    double skew = 0, relative_slope, my_time_slope, total_slope, degradation=0.0,
           degradation_since_last_lb = 0.0;

    std::vector<double> timings(worldsize);
    std::vector<double> all_degradations;
    CPULoadDatabase gossip_workload_db(world);

    PAR_START_TIMING(loop_time, world);
    for(unsigned int step = 0; step < MAX_STEP; ++step) {
        if(!rank) steplogger->info() << "Beginning step "<< step;
        PAR_START_TIMING(step_time, world);
#ifdef LB_METHOD
        bool lb_condition = false;
#if LB_METHOD==1   // load balance every 100 iterations
        lb_condition = (pcall + ncall) == step;
        if(lb_condition) {
            if(!rank) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            zoltan_load_balance(&my_cells, zoltan_lb, true, true);
            PAR_STOP_TIMING(current_lb_cost, world);
            if(!rank) perflogger->info("LB_time: ") << current_lb_cost;

            lb_costs.push_back(current_lb_cost);
            ncall = 25;
            if(!rank) steplogger->info("next LB call at: ") << (step+ncall);
            pcall = step;
        }
#elif LB_METHOD == 2 // http://sc16.supercomputing.org/sc-archive/tech_poster/poster_files/post247s2-file3.pdf
        if(!rank) steplogger->info("degradation: ") << (degradation_since_last_lb*(step-pcall))/2.0 << " avg_lb_cost " << avg_lb_cost;
        lb_condition = (pcall + ncall) == step;
        if(lb_condition) {
            double total_slope = get_slope<double>(window_step_time.data_container);
            if(!rank) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            zoltan_load_balance(&my_cells, zoltan_lb, true, true);
            PAR_STOP_TIMING(current_lb_cost, world);
            if(!rank) perflogger->info("LB_time: ") << current_lb_cost;

            lb_costs.push_back(current_lb_cost);
            avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());
            if(total_slope > 0) ncall = std::sqrt((2 * avg_lb_cost) / total_slope);
            else ncall = MAX_STEP;
            degradation_since_last_lb = 0.0;
            window_my_time.data_container.clear();
            window_step_time.data_container.clear();
            if(!rank) steplogger->info("next LB call at: ") << (step+ncall);
            pcall = step;
        }
#elif LB_METHOD == 3 // Unloading Model
        constexpr double alpha = 1.0/4.0;
        if(!rank) steplogger->info("degradation: ") << degradation_since_last_lb << " avg_lb_cost " << avg_lb_cost;
        lb_condition = pcall + ncall <= step || (degradation_since_last_lb*(step-pcall))/2.0 > avg_lb_cost;
        if( lb_condition ) {
            double my_slope =  (std::floor(get_slope<double>(window_my_time.data_container)*1000.0)) / 1000.0;
            double step_slope =  (std::floor(get_slope<double>(window_step_time.data_container)*1000.0)) / 1000.0;
            auto total_last_step_time = *(window_step_time.end()-1);
            auto my_last_step_time    = *(window_my_time.end()-1);

            PAR_START_TIMING(current_lb_cost, world);
            //if(i_am_overloaded && i_am_increasing) {
                //std::cout << rank << " !!!!!!!!!!!!!!!!!!! " << step_slope << " " << my_slope << std::endl;
            int parts_num[1] = {rank}, weight_per_obj[1] = {0};
            float part_size[1] = {(1.0f - (float) my_slope)};
            Zoltan_LB_Set_Part_Sizes(zoltan_lb, 1, 1, parts_num, weight_per_obj, part_size);
            update_cell_weights(&my_cells, part_size[0], WATER_TYPE, [](auto a, auto b){return a * 1.0/b;});
            /*} else {
                int parts_num[1] = {rank}, weight_per_obj[1] = {0};
                float part_size[1] = {1.0};
                Zoltan_LB_Set_Part_Sizes(zoltan_lb, 1, 1, parts_num, weight_per_obj, part_size);
            }*/
            auto data = zoltan_load_balance(&my_cells, zoltan_lb, true, true);
            std::cout << rank << " " << data << std::endl;
            PAR_STOP_TIMING(current_lb_cost, world);
            if(!rank) perflogger->info("LB_time: ") << current_lb_cost;
            lb_costs.push_back(current_lb_cost);
            avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());
            //std::cout << rank << " number of cells: " << compute_effective_workload(my_cells, WATER_TYPE) << " effective load " << compute_estimated_workload(my_cells) << std::endl ;
            degradation_since_last_lb = 0.0;
            window_my_time.data_container.clear();
            window_step_time.data_container.clear();
            pcall = step;
            ncall = 10;
        }
#endif
        if(lb_condition) {
            my_water_ptr = create_water_ptr_vector(my_cells);
        }
#endif
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// COMPUTATION START

        PAR_START_TIMING(comp_time, world);
        auto remote_cells = zoltan_exchange_data(zoltan_lb,my_cells,&recv,&sent,datatype.element_datatype,world,1.0);
        CHECKPOINT_TIMING(comp_time, my_exchange_time);
        auto remote_water_ptr = create_water_ptr_vector(remote_cells);
        auto bbox = get_bounding_box(my_cells, remote_cells);
        populate_data_pointers(msx, msy, &data_pointers, my_cells, remote_cells, bbox);
        //my_cells = dummy_erosion_computation2(msx, msy, my_cells,  remote_cells,  data_pointers, bbox);
        decltype(my_water_ptr) new_water_ptr;
        std::tie(my_cells, new_water_ptr) = dummy_erosion_computation3(msx, msy, my_cells, my_water_ptr, remote_cells, remote_water_ptr, data_pointers, bbox);
        my_water_ptr.insert(my_water_ptr.end(), std::make_move_iterator(new_water_ptr.begin()), std::make_move_iterator(new_water_ptr.end()));
        CHECKPOINT_TIMING(comp_time, my_comp_time);
        PAR_STOP_TIMING(comp_time, world);
        PAR_STOP_TIMING(step_time, world);

        window_step_time.add(comp_time); // monitor evolution of load in time with a window
        window_my_time.add(my_comp_time);    // monitor evolution of my load in time with a window

        if(step > 0)
            gossip_workload_db.finish_gossip_step();
        gossip_workload_db.gossip_update(my_comp_time);
        gossip_workload_db.gossip_propagate();

        if(window_step_time.size() > 2)
            degradation += (window_step_time.mean() - window_step_time.median(window_step_time.end() - 3, window_step_time.end() - 1)); // median among [cts-2, cts]

        if(pcall+1 < step)
            degradation_since_last_lb += *(window_step_time.end() - 1) - *(window_step_time.end() - 2) ;

        /// COMPUTATION STOP
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<double> exch_timings(worldsize);
        MPI_Allgather(&my_exchange_time, 1, MPI_DOUBLE, &exch_timings.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
        MPI_Allgather(&my_comp_time,     1, MPI_DOUBLE, &timings.front(),      1, MPI_DOUBLE, world); // TODO: propagate information differently
        std::vector<double> slopes(worldsize);
        std::vector<int> tloads(worldsize);
        my_time_slope = (std::floor(get_slope<double>(window_my_time.data_container)*1000.0)) / 1000.0;
        MPI_Allgather(&my_time_slope, 1, MPI_DOUBLE, &slopes.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
        int eload = compute_effective_workload(my_cells, 1);
        MPI_Allgather(&eload, 1, MPI_INT, &tloads.front(), 1, MPI_INT, world); // TODO: propagate information differently

        double max = *std::max_element(timings.cbegin(), timings.cend()),
               average = std::accumulate(timings.cbegin(), timings.cend(), 0.0) / worldsize,
               load_imbalance = (max / average - 1.0) * 100.0;

        if(!rank) {
            steplogger->info("stats: ") << "load imbalance: " << load_imbalance << " skewness: " << skew;
            perflogger->info("\"step\":") << step
            << ",\"step_time\": " << step_time
            << ",\"total_time\": " << std::accumulate(timings.begin(), timings.end(), 0.0)
            << ",\"lb_cost\": " << stats::mean<double>(lb_costs.begin(), lb_costs.end())
            << ",\"load_imbalance\": " << load_imbalance
            << ",\"skewness\": " << stats::skewness<double>(timings.begin(), timings.end())
            << ",\"loads\": [" << timings
            << ",\"exch\": ["  << exch_timings
            << "],\"tloads\":["<<tloads
            << "],\"slopes\":["<<slopes<<"]";
        }


#ifdef PRODUCE_OUTPUTS
        unsigned long cell_cnt = my_cells.size();
        std::vector<std::array<int,2>> my_types(cell_cnt);
        for (int i = 0; i < cell_cnt; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
        gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
        if(!rank) {
            std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
            std::vector<int> types, rock_gid;
            std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
            std::for_each(all_types.begin(), all_types.end(), [&rock_gid](auto e){if(!e[1]) rock_gid.push_back(e[0]);});
            cnpy::npz_save("gids-out.npz", "step-"+std::to_string(step+1), &rock_gid[0], {rock_gid.size()}, "a");
        }
#endif
        if(!rank) steplogger->info() << "Stop step "<< step;
    }
    PAR_STOP_TIMING(loop_time, world);
    if(!rank) perflogger->info("\"total_time\":") << loop_time;
    if(!rank) steplogger->info("\"total_time\":") << loop_time;
    datatype.free_datatypes();
    MPI_Finalize();
    return 0;
}
