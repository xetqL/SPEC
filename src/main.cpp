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

#include "../include/zoltan_fn.hpp"
#include "../include/window.hpp"
#include "../include/io.hpp"
#include "../include/CPULoadDatabase.hpp"

#include "../include/time.hpp"
#include "../include/band_partitioner.hpp"

#include "zupply.hpp"
#include "main.hpp"

zz::log::LoggerPtr perflogger, steplogger;

int main(int argc, char **argv) {
    int worldsize;
    int rank;

    MPI_Init(&argc, &argv);
    auto world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldsize);
    MPI_Comm_rank(world, &rank);


    auto datatype   = Cell::register_datatype();

    std::string str_rank = "[RANK " + std::to_string(rank) + "] ";
    // Generate data
    int xprocs = std::atoi(argv[2]), yprocs = std::atoi(argv[1]);
    int cell_per_process = std::atoi(argv[3]);
    const int MAX_STEP = std::atoi(argv[4]);
    const int seed = std::atoi(argv[5]);
    std::mt19937 gen(seed); // common seed
    std::uniform_int_distribution<>  proc_dist(0, worldsize-1);


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

    std::shuffle(cols.begin(), cols.end(), gen);
    cols = {0, 1};
    std::vector<int> water_cols(cols.begin(), cols.begin()+1 );
    std::sort(water_cols.begin(), water_cols.end());
    //std::for_each(water_cols.cbegin(), water_cols.cend(), [](auto v){ std::cout << v << std::endl; });

    int& msx = Cell::get_msx(); msx = xcells;
    int& msy = Cell::get_msy(); msy = ycells;
    int inner_type = WATER_TYPE;
    int shape[2] = {msx,msy};

#ifdef PRODUCE_OUTPUTS
    if(!rank) cnpy::npz_save("gids-out.npz", "shape", &shape[0],   {2}, "w");
    if(!rank) cnpy::npz_save("gids-out.npz", "type",  &inner_type, {1}, "a");
#endif

    int x_proc_idx, y_proc_idx;
    std::tie(x_proc_idx, y_proc_idx) = cell_to_global_position(xprocs, yprocs, rank);
    unsigned long total_cell = cell_per_process * worldsize;

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
    //my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);
    my_cells = generate_lattice_single_type(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, WATER_TYPE, 1.0, 0.0);
#else
    if(argc == 7) {
        auto lattice_fname = argv[6];
        my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols, lattice_fname);
    } else {
        //my_cells = generate_lattice_percolation_diffusion(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, water_cols);
        my_cells = generate_lattice_single_type(msx, msy, x_proc_idx, y_proc_idx,cell_in_my_cols,cell_in_my_rows, inner_type, 1.0, 0.0);
    }
#endif
    const int my_cell_count = my_cells.size();


    int recv, sent;
    if(!rank) steplogger->info("End of map generation");

    int loading_proc = 1;//proc_dist(gen);
    bool i_am_loading_proc = rank == loading_proc;
    /* Initial load balancing */
    std::vector<double> lb_costs;
    // auto zoltan_lb = zoltan_create_wrapper(true, world);
    StripeLoadBalancer stripe_lb(i_am_loading_proc, msx, msy, datatype.element_datatype, world);
    stripe_lb.setLoggerPtr(steplogger);
    PAR_START_TIMING(current_lb_cost, world);
    stripe_lb.load_balance(&my_cells, 0.0);
    PAR_STOP_TIMING(current_lb_cost, world);
    lb_costs.push_back(current_lb_cost);
    auto my_domain = stripe_lb.get_domain(rank);

    generate_lattice_rocks(4, msx, msy, &my_cells, i_am_loading_proc ? 0.5f : 0.001f, my_domain.first, my_domain.second);

    //stripe_lb.load_balance(&my_cells, i_am_loading_proc ? 0 : 0.0);

#ifdef PRODUCE_OUTPUTS
    std::vector<std::array<int,2>> all_types(total_cell);
    std::vector<std::array<int,2>> my_types(my_cell_count);
    for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
    gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
    if(i_am_loading_proc) {
        std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
        std::vector<int> types, water_gid;
        std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
        std::for_each(all_types.begin(), all_types.end(), [&water_gid, &inner_type](auto e){if(e[1] == inner_type) water_gid.push_back(e[0]);});
        assert((*(all_types.end() - 1))[0] == total_cell-1);
        cnpy::npz_save("gids-out.npz", "step-"+std::to_string(0), &water_gid[0], {water_gid.size()}, "a");
    }
    if(i_am_loading_proc) all_types.resize(total_cell);
#endif

    if(i_am_loading_proc) perflogger->info("LB_time: ") << current_lb_cost;

    auto my_water_ptr = create_water_ptr_vector(my_cells);
    std::vector<size_t> data_pointers, remote_data_pointers;
    my_water_ptr = create_water_ptr_vector(my_cells);
    std::tuple<int, int, int, int> bbox; // = add_to_bbox(msx, msy, get_bounding_box(my_cells), -10, 10, -10, 10);
    //populate_data_pointers(msx, msy, &data_pointers, my_cells, 0, bbox, true);
    /* lets make it fun now...*/
    int minx, maxx, miny, maxy;

#ifdef LB_METHOD
    auto avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());
#endif

    int ncall = 10, pcall=0;
#if LB_METHOD==1
    ncall = 100;
#endif

    float average_load_at_lb;
    double skew = 0, relative_slope, my_time_slope, total_slope, degradation=0.0, degradation_since_last_lb = 0.0;

    std::vector<double> timings(worldsize), all_degradations, water;

    SlidingWindow<double> window_step_time(100); // sliding window with max size = TODO: tune it?
    SlidingWindow<double> window_water(ncall);   // sliding window with max size = TODO: tune it?
    SlidingWindow<double> window_my_time(100);   // sliding window with max size = TODO: tune it?

    water.push_back(my_water_ptr.size()); window_water.add(my_water_ptr.size());

    CPULoadDatabase gossip_workload_db(world);
    CPULoadDatabase gossip_waterslope_db(world);
    //double time_since_start;
    PAR_START_TIMING(loop_time, world);
    for(unsigned int step = 0; step < MAX_STEP; ++step) {
        if(i_am_loading_proc) steplogger->info() << "Beginning step "<< step;

        PAR_START_TIMING(step_time, world);
        bool lb_condition = false;
#ifdef LB_METHOD

#if LB_METHOD==1   // load balance every 100 iterations
        lb_condition = (pcall + ncall) == step;
        if(lb_condition) {
            if(!rank) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            stripe_lb.load_balance(&my_cells, 0.0);
            PAR_STOP_TIMING(current_lb_cost, world);
            if(!rank) perflogger->info("LB_time: ") << current_lb_cost;
            lb_costs.push_back(current_lb_cost);
            ncall = 100;
            if(!rank) steplogger->info("next LB call at: ") << (step+ncall);
            pcall = step;
        }
#elif LB_METHOD == 2 // http://sc16.supercomputing.org/sc-archive/tech_poster/poster_files/post247s2-file3.pdf
        //if(!rank) steplogger->info("degradation method 2: ") << (degradation_since_last_lb*(step-pcall))/2.0 << " avg_lb_cost " << avg_lb_cost;
        lb_condition = pcall + ncall <= step;// || (degradation_since_last_lb*(step-pcall))/2.0 > avg_lb_cost;
        if(lb_condition) {
            auto total_slope = get_slope<double>(window_step_time.data_container);
            if(!rank) steplogger->info("call LB at: ") << step;
            PAR_START_TIMING(current_lb_cost, world);
            stripe_lb.load_balance(&my_cells, 0.0);
            PAR_STOP_TIMING(current_lb_cost, world);
            lb_costs.push_back(current_lb_cost);
            if(!rank) perflogger->info("LB_time: ") << current_lb_cost;
            avg_lb_cost = stats::mean<double>(lb_costs.begin(), lb_costs.end());
            if(total_slope > 0) {
                ncall = (int) std::floor(std::sqrt((2.0 * avg_lb_cost) / total_slope));
                std::cout << ncall << std::endl;
            } else
                ncall = MAX_STEP;
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
#elif LB_METHOD==4
        if(i_am_loading_proc) steplogger->info("degradation method 4: ") << ((degradation_since_last_lb*(step-pcall))/2.0) << " avg_lb_cost " << avg_lb_cost;
        lb_condition = pcall + ncall <= step || ((degradation_since_last_lb*(step-pcall))/2.0 > avg_lb_cost && gossip_waterslope_db.has_converged(7));
        if(lb_condition) {
            bool overloading = gossip_waterslope_db.zscore(rank) > 3.0;
            PAR_START_TIMING(current_lb_cost, world);
            stripe_lb.load_balance(&my_cells, overloading ? 0.05 : 0.0);
            PAR_STOP_TIMING(current_lb_cost, world);
            lb_costs.push_back(current_lb_cost);
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
            my_water_ptr = create_water_ptr_vector(my_cells);
            water.push_back(my_water_ptr.size());
            window_water.add(my_water_ptr.size());
        }
#endif
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /// COMPUTATION START

        PAR_START_TIMING(comp_time, world);
        if(i_am_loading_proc) steplogger->info("start communication");
        auto remote_cells = stripe_lb.share_frontier_with_neighbors(my_cells, &recv, &sent);//zoltan_exchange_data(zoltan_lb,my_cells,&recv,&sent,datatype.element_datatype,world,1.0);
        if(i_am_loading_proc) steplogger->info("create water ptr");
        auto remote_water_ptr = create_water_ptr_vector(remote_cells);
        decltype(my_water_ptr) new_water_ptr;
        MPI_Barrier(world);
        if(i_am_loading_proc) steplogger->info("get_bounding_box");
        if(lb_condition || step == 0) bbox = get_bounding_box(my_cells, remote_cells);
        MPI_Barrier(world);
        if(i_am_loading_proc) steplogger->info("populate_data_pointers");
        populate_data_pointers(msx, msy, &data_pointers, my_cells, remote_cells, bbox, lb_condition || step == 0);
        MPI_Barrier(world);
        if(i_am_loading_proc) steplogger->info("dummy_erosion_computation3");
        std::tie(my_cells, new_water_ptr) = dummy_erosion_computation3(msx, msy, my_cells, my_water_ptr, remote_cells, remote_water_ptr, data_pointers, bbox);
        my_water_ptr.insert(my_water_ptr.end(), std::make_move_iterator(new_water_ptr.begin()), std::make_move_iterator(new_water_ptr.end()));

        CHECKPOINT_TIMING(comp_time, my_comp_time);
        PAR_STOP_TIMING(comp_time, world);
        PAR_STOP_TIMING(step_time, world);

        CHECKPOINT_TIMING(loop_time, time_since_start);
        if(i_am_loading_proc) steplogger->info("time until step ") << step << " = " << time_since_start;
        MPI_Barrier(world);
        if(i_am_loading_proc) steplogger->info("start pushing latest data");
        water.push_back(my_water_ptr.size());
        window_water.add(my_water_ptr.size());

        window_step_time.add(comp_time); // monitor evolution of load in time with a window
        window_my_time.add(my_comp_time);    // monitor evolution of my load in time with a window
        MPI_Barrier(world);
        if(i_am_loading_proc) steplogger->info("start gossip step");
        if(step > 0) {
            gossip_workload_db.finish_gossip_step();
            gossip_waterslope_db.finish_gossip_step();
        }

        gossip_waterslope_db.gossip_update(get_slope<double>(water.begin(), water.end()));
        gossip_waterslope_db.gossip_propagate();

        gossip_workload_db.gossip_update(my_comp_time);
        gossip_workload_db.gossip_propagate();
        MPI_Barrier(world);
        if(i_am_loading_proc) steplogger->info("start degradation computation");

        if(window_step_time.size() > 2)
            degradation += (window_step_time.mean() - window_step_time.median(window_step_time.end() - 3, window_step_time.end() - 1)); // median among [cts-2, cts]

        if(pcall + 1 < step)
            degradation_since_last_lb += *(window_step_time.end() - 1) - *(window_step_time.end() - 2) ;
        MPI_Barrier(world);
        if(i_am_loading_proc) steplogger->info("End of data processing");
        /// COMPUTATION STOP
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef PRODUCE_OUTPUTS
        std::vector<double> exch_timings(worldsize);
        //MPI_Allgather(&cpt, 1, MPI_DOUBLE, &exch_timings.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
        MPI_Allgather(&my_comp_time, 1, MPI_DOUBLE, &timings.front(),      1, MPI_DOUBLE, world); // TODO: propagate information differently
        std::vector<double> slopes(worldsize);
        std::vector<int> tloads(worldsize);
        auto my_water_slope = get_slope<double>(water.begin(), water.end());
        MPI_Allgather(&my_water_slope, 1, MPI_DOUBLE, &slopes.front(), 1, MPI_DOUBLE, world); // TODO: propagate information differently
        int eload = compute_effective_workload(my_cells, 1);
        MPI_Allgather(&eload, 1, MPI_INT, &tloads.front(), 1, MPI_INT, world); // TODO: propagate information differently

        double max = *std::max_element(timings.cbegin(), timings.cend()),
               average = std::accumulate(timings.cbegin(), timings.cend(), 0.0) / worldsize,
               load_imbalance = (max / average - 1.0) * 100.0;

        if(i_am_loading_proc) {
            steplogger->info("stats: ") << "load imbalance: " << load_imbalance << " skewness: " << skew;
            steplogger->info("Water Size: ") << water;
            perflogger->info("\"step\":") << step
            << ",\"step_time\": " << step_time
            << ",\"total_time\": " << std::accumulate(timings.begin(), timings.end(), 0.0)
            << ",\"lb_cost\": "  << stats::mean<double>(lb_costs.begin(), lb_costs.end())
            << ",\"load_imbalance\": " << load_imbalance
            << ",\"skewness\": " << stats::skewness<double>(timings.begin(), timings.end())
            << ",\"loads\": ["   << timings
            //<< ",\"exch\": ["    << exch_timings
            << "],\"tloads\":["  << tloads
            << "],\"slopes\":["  << gossip_waterslope_db.get_all_data()<<"]";
        }

        unsigned long cell_cnt = my_cells.size();
        std::vector<std::array<int,2>> my_types(cell_cnt);
        for (int i = 0; i < cell_cnt; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};
        gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
        if(i_am_loading_proc) {
            std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
            std::vector<int> types, rock_gid;
            std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
            std::for_each(all_types.begin(), all_types.end(), [&rock_gid, &inner_type](auto e){if(e[1] == inner_type) rock_gid.push_back(e[0]);});
            cnpy::npz_save("gids-out.npz", "step-"+std::to_string(step+1), &rock_gid[0], {rock_gid.size()}, "a");
        }
#endif
        if(i_am_loading_proc) steplogger->info() << "Stop step "<< step;
    }
    PAR_STOP_TIMING(loop_time, world);
    if(i_am_loading_proc) perflogger->info("\"total_time\":") << loop_time;
    if(i_am_loading_proc) steplogger->info("\"total_time\":") << loop_time;
    datatype.free_datatypes();
    MPI_Finalize();
    return 0;
}
