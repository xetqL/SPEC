//
// Created by xetql on 5/3/19.
//

#ifndef SPEC_SIMULATED_LBM_HPP
#define SPEC_SIMULATED_LBM_HPP

#include <memory>
#include <mpi.h>
#include <random>

#include "Cell.hpp"
#include "LoadBalancer.hpp"
#include "SimulationParams.hpp"
#include "zupply.hpp"
#include <GossipDatabase.hpp>
#include <Window.hpp>
#include "Simflow.hpp"
#include "SimulatedLBM.hpp"
#include "Utils.hpp"
#include "io.hpp"

#include "ULBA.hpp"

#ifdef PRODUCE_OUTPUTS
#include <cnpy.h>
#endif

class SimulatedLBM {
    using LBAlgorithm = LoadBalancer<Cell>;
    const int FOREMAN = 0;

    SimulationParams params;

    MPI_Comm comm;
    std::unique_ptr<LBAlgorithm> load_balancer;

    zz::log::LoggerPtr perflogger, steplogger, proctime;
public:
    SimulatedLBM(SimulationParams params, MPI_Comm comm, SimulatedLBM::LBAlgorithm *load_balancer) :
            params(std::move(params)), comm(comm), load_balancer(load_balancer) {}

    template<class Approach> void run(float alpha) {
        auto world = this->comm;
        int worldsize;
        int rank;
        MPI_Comm_size(world, &worldsize);
        MPI_Comm_rank(world, &rank);

        const unsigned int xprocs = params.xprocs,
                yprocs = params.yprocs,
                cell_per_process = params.cell_per_process,
                MAX_STEP = params.MAX_STEP;

        const int seed = params.seed;
        unsigned int N = params.N;
        bool load_lattice = params.load_lattice, verbose = params.verbose;

        SlidingWindow<double> window_step_time(15);  // sliding window with max size = TODO: tune it?
        SlidingWindow<double> window_my_time(100);   // sliding window with max size = TODO: tune it?

        GossipDatabase<double> gossip_workload_db(worldsize,    2, 9999, world),
                gossip_waterslope_db(worldsize,  2, 8888, world);

        const bool i_am_foreman = rank == FOREMAN;
        auto datatype   = Cell::register_datatype();

        std::string str_rank = "[RANK " + std::to_string(rank) + "] ";
        std::mt19937 gen(seed); // common seed
        std::uniform_int_distribution<>  proc_dist(0, worldsize-1);

        int cell_in_my_rows = (int) std::sqrt(cell_per_process), cell_in_my_cols = cell_in_my_rows;
        int xcells = cell_in_my_rows * xprocs, ycells = cell_in_my_rows * yprocs;

        auto total_cell = xcells * ycells;
        std::vector<int> cols(ycells);
        std::iota(cols.begin(), cols.end(), 0);

        std::shuffle(cols.begin(), cols.end(), gen);
        cols = {0, 1};
        std::vector<int> water_cols(cols.begin(), cols.begin()+1 );

        std::sort(water_cols.begin(), water_cols.end());
        //std::for_each(water_cols.cbegin(), water_cols.cend(), [](auto v){ std::cout << v << std::endl; });

        int& msx = Cell::get_msx(); msx = xcells;
        int& msy = Cell::get_msy(); msy = ycells;

        std::vector<int> loading_procs;
        //const int loading_proc = proc_dist(gen) % worldsize; //one randomly chosen load proc
        std::uniform_int_distribution<>  lproc_dist((int)N/2, worldsize-1-(int)N/2);

        const int center = worldsize-1;//lproc_dist(gen);
        loading_procs.push_back(center);

        for(int i = 1; i < N; i++) {
            loading_procs.push_back(center - i);
            loading_procs.push_back(center + i);
        }

        //std::for_each(loading_procs.begin(), loading_procs.end(), [](auto v){std::cout << v << std::endl;});

        const bool i_am_loading_proc = std::find(loading_procs.begin(), loading_procs.end(), rank) != loading_procs.end();

        int x_proc_idx, y_proc_idx;
        std::tie(x_proc_idx, y_proc_idx) = cell_to_global_position(xprocs, yprocs, rank);

        if(i_am_foreman) steplogger->info() << cell_in_my_cols << " " << cell_in_my_rows;
        std::vector<Cell> my_cells;
        my_cells = generate_lattice_single_type(msx, msy, x_proc_idx, y_proc_idx, cell_in_my_cols, cell_in_my_rows, Cell::WATER_TYPE, 1.0, 0.0);

        int bottom, top;

        bottom = y_proc_idx * cell_in_my_rows;
        top    = (y_proc_idx+1) * cell_in_my_rows;
        generate_lattice_rocks(1, msx, msy, &my_cells, i_am_loading_proc ? 0.4f : 0.00f, bottom, top, (top-bottom) / 3);

        int recv, sent;
        std::vector<double> lb_costs;

        this->load_balancer->activate_load_balance(0, &my_cells);

        if(std::is_same<Approach , ULBA>::value)
            this->load_balancer->set_approach(new ULBA(world, &gossip_waterslope_db, 3.0, alpha));

#ifdef PRODUCE_OUTPUTS
        int inner_type = Cell::ROCK_TYPE;
        int shape[2] = {msx, msy};

        if(i_am_foreman) cnpy::npz_save("gids-out.npz", "shape", &shape[0],   {2}, "w");
        if(i_am_foreman) cnpy::npz_save("gids-out.npz", "type",  &inner_type, {1}, "a");

        auto my_cell_count = my_cells.size();

        std::vector<std::array<int,2>> all_types(worldsize*cell_per_process), my_types(my_cell_count);

        for (int i = 0; i < my_cell_count; ++i) my_types[i] = {my_cells[i].gid, my_cells[i].type};

        gather_elements_on(my_types, FOREMAN, &all_types, datatype.minimal_datatype, world);

        if(i_am_foreman) {

            std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});

            std::vector<int> types, water_gid;

            //std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});

            std::for_each(all_types.begin(), all_types.end(), [&water_gid, &inner_type](auto e){if(e[1] == inner_type) water_gid.push_back(e[0]);});

            //assert((*(all_types.end() - 1))[0] == total_cell-1);
            cnpy::npz_save("gids-out.npz", "step-"+std::to_string(0), &water_gid[0], {water_gid.size()}, "a");

        }
        MPI_Barrier(world);
        if(i_am_foreman) all_types.resize(total_cell);

#endif

        std::vector<unsigned long> my_water_ptr;
        int64_t n;
        std::vector<size_t> data_pointers, remote_data_pointers;
        std::tie(n, my_water_ptr) = create_water_ptr_vector(my_cells);
        std::tuple<int, int, int, int> bbox; // = add_to_bbox(msx, msy, get_bounding_box(my_cells), -10, 10, -10, 10);
        //populate_data_pointers(msx, msy, &data_pointers, my_cells, 0, bbox, true);

        /* lets make it fun now...*/
        double degradation_since_last_lb = 0.0;
        double perfect_time_value = 0.0;
        std::vector<double> timings(worldsize), all_degradations, water, compTimes, stepTimes, deltaWorks, loadImbalance;
        unsigned int ncall = MAX_STEP;
        water.push_back(my_water_ptr.size());
        double my_gossip_time = 0.0;

        n = compute_estimated_workload<int64_t>(my_cells);

        if(!rank) steplogger->info("Starting simulation: ") << params << " with " << load_balancer->get_approach_name();

        PAR_START_TIMING(loop_time, world);
        for(unsigned int step = 0; step < MAX_STEP; ++step) {
            if(i_am_foreman) steplogger->info() << "Beginning step "<< step;
            PAR_START_TIMING(step_time, world);

#ifdef AUTONOMIC_LOAD_BALANCING
            bool lb_condition = this->load_balancer->get_last_call() + ncall <= step || degradation_since_last_lb > this->load_balancer->get_average_cost();
#elif  CYCLIC_LOAD_BALANCING
            bool lb_condition = this->load_balancer->get_last_call() + params.interval <= step;
#else   // NO_LOAD_BALANCING
        bool lb_condition = false;
#endif
            if(lb_condition) {

                // bool overloading = gossip_waterslope_db.zscore(rank) > 3.0;
                this->load_balancer->activate_load_balance(step, &my_cells);
#ifdef AUTONOMIC_LOAD_BALANCING
                double median;
                if(std::distance(window_step_time.begin(), window_step_time.end() - 3) < 0)
                    median  = stats::median<double>(window_step_time.begin(), window_step_time.end());
                else median = stats::median<double>(window_step_time.end() - 3, window_step_time.end());
                auto mean   = stats::mean<double>(window_step_time.begin(), window_step_time.end());
                unsigned int pcall = this->load_balancer->get_last_call();
                ncall = (unsigned int) this->load_balancer->estimate_best_ncall(((step - pcall) - 1), median-mean);
#elif  CYCLIC_LOAD_BALANCING
                ncall = params.interval;
#endif
                gossip_workload_db.reset();
                water.clear();
                degradation_since_last_lb = 0.0;
                window_my_time.data_container.clear();
                window_step_time.data_container.clear();
                std::tie(n, my_water_ptr) = create_water_ptr_vector(my_cells);
                water.push_back(n);
                deltaWorks.clear();
                n = compute_estimated_workload<int64_t>(my_cells);
            }

            //PAR_STOP_TIMING(step_time, world);
            //STOP_TIMING(loop_time);
            double add_weight = 0;

            auto remote_cells = this->load_balancer->propagate(my_cells, &recv, &sent, 1.0);

            PAR_START_TIMING(comp_time, world);


            decltype(my_water_ptr) remote_water_ptr;

            //it is actually cheap because remote_cells is small
            std::tie(std::ignore, remote_water_ptr) = create_water_ptr_vector(remote_cells);

            decltype(my_water_ptr) new_water_ptr;

            STOP_TIMING(comp_time);
            if(lb_condition || step == 0) bbox = get_bounding_box(my_cells, remote_cells);
            populate_data_pointers(msx, msy, &data_pointers, my_cells, remote_cells, bbox, lb_condition || step == 0);
            RESTART_TIMING(comp_time);
            // compute_estimated_workload(my_cells);

            //RESTART_TIMING(loop_time);
            //PAR_RESTART_TIMING(step_time, world);

            my_cells = dummy_erosion_computation3(step, msx, msy, my_cells, my_water_ptr, remote_cells, remote_water_ptr, data_pointers.data(), bbox, &new_water_ptr, &add_weight);
            compute_fluid_time(n);

            my_water_ptr.insert(my_water_ptr.end(), std::make_move_iterator(new_water_ptr.begin()), std::make_move_iterator(new_water_ptr.end()));
            n += (int64_t) add_weight; // adapt the number of cell to compute
            water.push_back(n);

            CHECKPOINT_TIMING(comp_time, my_comp_time);

            PAR_STOP_TIMING(comp_time, world);
            PAR_STOP_TIMING(step_time, world);
            CHECKPOINT_TIMING(loop_time, time_since_start);
            STOP_TIMING(loop_time);
            MPI_Allreduce(&my_comp_time, &comp_time, 1, MPI_DOUBLE, MPI_MAX, world); // i should not need that!

            // this is the reference time for further steps
            if(this->load_balancer->get_last_call() == step) perfect_time_value = comp_time;
            double currDegradation = std::max((comp_time - perfect_time_value), 0.0);

            deltaWorks.push_back(currDegradation);
            compTimes.push_back(comp_time);
            stepTimes.push_back(step_time);
            window_step_time.add(comp_time);  // monitor evolution of computing time with a window
            window_my_time.add(my_comp_time); // monitor evolution of my workload    with a window

            RESTART_TIMING(loop_time);
            if(std::is_same<Approach , ULBA>::value) {
                RESTART_TIMING(my_gossip_time);
                gossip_waterslope_db.execute(rank, get_slope<double>(water.begin(), water.end()));
                gossip_workload_db.execute(  rank,   my_comp_time);
                STOP_TIMING(my_gossip_time);
            }

            //update_cells(&my_cells, gossip_waterslope_db.get(rank), [](Cell &c, auto v){c.slope = v;});

            if(this->load_balancer->get_last_call() + 1 < step) {
                degradation_since_last_lb +=
                        stats::median<double>(window_step_time.end()-3, window_step_time.end()) - perfect_time_value;
            }
            STOP_TIMING(loop_time);

            /// COMPUTATION STOP
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            std::vector<double> exch_timings(worldsize), slopes(worldsize);
            std::vector<int> tloads(worldsize);
            std::vector<int64_t> all_the_workloads(worldsize);

            int64_t my_workload = n;

            MPI_Gather(&my_workload,  1, MPI_INT64_T,  &all_the_workloads.front(), 1, MPI_INT64_T,  FOREMAN, world);
            MPI_Gather(&my_comp_time, 1, MPI_DOUBLE,   &timings.front(),           1, MPI_DOUBLE,   FOREMAN, world);

            if(i_am_foreman) {
                auto total_slope = get_slope<double>(window_step_time.data_container);
                steplogger->info("degradation since last LB ") << degradation_since_last_lb << ", avg_lb_cost " << this->load_balancer->get_average_cost() << ", total slope: " << total_slope;
                steplogger->info("time for step ") << step << " = " << step_time << " computation_time = "<< comp_time
                                                   << " total = " << time_since_start
                                                   << " dW: "<< stats::mean<double>(deltaWorks.begin(), deltaWorks.end());

                double max = *std::max_element(timings.cbegin(), timings.cend()),
                        average = std::accumulate(timings.cbegin(), timings.cend(), 0.0) / worldsize,
                        load_imbalance = (max / average - 1.0) * 100.0;

                loadImbalance.push_back(load_imbalance);

                perflogger->info("\"step\":") << step << ",\"LI\": " << load_imbalance;
                perflogger->info("\"step\":") << step << ",\"workloads\": " << all_the_workloads;
                proctime->info("\"step\":") << step << ",\"proctime\": " << timings;
            }

#ifdef PRODUCE_OUTPUTS
            unsigned long cell_cnt = my_cells.size();
            std::vector<std::array<int,2>> my_types(cell_cnt);

            if(step % 5 == 0) {
                //std::cout << cell_cnt << std::endl;
                for (unsigned int i = 0; i < cell_cnt; ++i) {
                    my_types[i] = {my_cells[i].gid, my_cells[i].type};
                    //this->steplogger->info() << my_cells[i].type;
                }
                gather_elements_on(my_types, 0, &all_types, datatype.minimal_datatype, world);
                if(i_am_foreman) {
                    std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
                    std::vector<int> types, rock_gid;
                    std::for_each(all_types.begin(), all_types.end(), [&types](auto e){types.push_back(e[1]);});
                    std::for_each(all_types.begin(), all_types.end(), [&rock_gid, &inner_type](auto e){if(e[1] == inner_type) rock_gid.push_back(e[0]);});
                    cnpy::npz_save("gids-out.npz", "step-"+std::to_string(step+1), &rock_gid[0], {rock_gid.size()}, "a");
                }

                for (unsigned int i = 0; i < cell_cnt; ++i) my_types[i] = {my_cells[i].gid, rank};

                gather_elements_on(my_types, FOREMAN, &all_types, datatype.minimal_datatype, world);

                if(i_am_foreman) {
                    std::sort(all_types.begin(), all_types.end(), [](auto a, auto b){return a[0] < b[0];});
                    std::vector<int> part_idx;
                    std::for_each(all_types.begin(), all_types.end(), [&part_idx](auto e){part_idx.push_back(e[1]);});
                    cnpy::npz_save("gids-out.npz", "partition-step-"+std::to_string(step+1), &part_idx[0], {part_idx.size()}, "a");
                }
            }
#endif
            RESTART_TIMING(loop_time);
        }
        PAR_STOP_TIMING(loop_time, world);
        double total_gossip_time;
        double total_loop_time;
        MPI_Allreduce(&my_gossip_time, &total_gossip_time, 1, MPI_DOUBLE, MPI_MAX, world);
        MPI_Allreduce(&loop_time,      &total_loop_time, 1,   MPI_DOUBLE, MPI_MAX, world);
        if(i_am_foreman) perflogger->info("total_time for ") << load_balancer->to_string() << " with " << load_balancer->get_approach_name() << " :"        << total_loop_time;
        if(i_am_foreman) perflogger->info("\"step times\":")        << stepTimes;
        if(i_am_foreman) perflogger->info("\"comp times\":")        << compTimes;
        if(i_am_foreman) perflogger->info("\"load imbalance\":")    << loadImbalance;
        if(i_am_foreman) perflogger->info("\"total_gossip_time\":") << total_gossip_time;
        datatype.free_datatypes();
    }

    void set_loggers(zz::log::LoggerPtr perflogger, zz::log::LoggerPtr steplogger, zz::log::LoggerPtr proctime) {
        this->perflogger = perflogger;
        this->steplogger = steplogger;
        this->proctime   = proctime;
    }

};

#endif //SPEC_SIMULATED_LBM_HPP
