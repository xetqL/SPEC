//
// Created by xetql on 5/3/19.
//

#ifndef SPEC_ZOLTANLOADBALANCER_HPP
#define SPEC_ZOLTANLOADBALANCER_HPP

#include <functional>
#include "LoadBalancer.hpp"
#include "ZoltanUtils.hpp"
#include "Utils.hpp"

template<class Data>
class ZoltanLoadBalancer : public LoadBalancer<Data> {
    using ZoltanCreateFunc = std::function<Zoltan_Struct* (bool, MPI_Comm)>;
    using ZoltanLBFunc     = std::function<void (std::vector<Data>*, Zoltan_Struct*, bool)>;
    Zoltan_Struct* zoltan_lb;
    const ZoltanLBFunc lb_func;
    const GossipDatabase<unsigned long>* workdb;
    MPI_Comm comm;
public:
    ZoltanLoadBalancer(MPI_Comm world, MPI_Datatype datatype,
                       const GossipDatabase<unsigned long>* workdb,
                       const ZoltanCreateFunc &zoltan_create_wrapper,
                       const ZoltanLBFunc &lb_func) :
        LoadBalancer<Data>(world, datatype, new StdApproach(world)),
        lb_func(lb_func),
        workdb(workdb), comm(world) {
        zoltan_lb = zoltan_create_wrapper(true, world);
    }

    ZoltanLoadBalancer(MPI_Comm world, MPI_Datatype datatype, const GossipDatabase<unsigned long>* workdb,
                       LoadBalancingApproach* approach,
                       const ZoltanCreateFunc& zoltan_create_wrapper,
                       const ZoltanLBFunc& lb_func) :
                       LoadBalancer<Data>(world, datatype, approach),
                       lb_func(lb_func),
                       workdb(workdb), comm(world) {
        zoltan_lb = zoltan_create_wrapper(true, world);
    }

    std::vector<Data> propagate(const std::vector<Data> &data,
                                int *nb_elements_recv, int *nb_elements_sent,
                                double cell_size) override;

private:
    void load_balance(std::vector<Data> *_data) override;

    std::vector<int> neighbors, cell_per_neighbors;
    std::vector<std::vector<unsigned long>> neighboring_cells;

    std::tuple<std::vector<int>, std::vector<std::vector<unsigned long>>, std::vector<int>> compute_neighborhood(const std::vector<Data>* data) {
        int wsize;
        MPI_Comm_size(comm, &wsize);
        int caller_rank;
        MPI_Comm_rank(comm, &caller_rank);

        std::vector<Data> buffer;

        long nb_reqs = 0;

        if (wsize == 1) return std::make_tuple(std::vector<int>(), std::vector<std::vector<unsigned long>>(), std::vector<int>());

        std::vector<std::vector<unsigned long>> data_to_migrate(wsize);
        size_t data_id = 0;
        std::vector<int> PEs(wsize, -1), parts(wsize, -1);
        int num_found, num_known = 0;
        std::vector<int> export_gids, export_lids, export_procs;

        // so much memory could be allocated here... potentially PE * n * DIM * 44 bytes => so linear in N
        // as DIM << PE <<<< n
        while (data_id < data->size()) {
            auto pos_in_double = data->at(data_id).template get_position_as_array<double>();
            int numparts;
            Zoltan_LB_Box_PP_Assign(zoltan_lb,
                                    pos_in_double.at(0) - 1,
                                    pos_in_double.at(1) - 1,
                                    pos_in_double.size() == 3 ? pos_in_double.at(2) - 1 : 0.0,
                                    pos_in_double.at(0) + 1,
                                    pos_in_double.at(1) + 1,
                                    pos_in_double.size() == 3 ? pos_in_double.at(2) + 1 : 0.0,
                                    &PEs.front(), &num_found, parts.data(), &numparts);

            for (int PE_idx = 0; PE_idx < num_found; PE_idx++) {
                int PE = PEs[PE_idx];
                if (PE >= 0 && PE != caller_rank) {
                    export_gids.push_back(data->at(data_id).gid);
                    export_lids.push_back(data_id);
                    export_procs.push_back(PE);

                    //get the value and copy it into the "to migrate" vector
                    data_to_migrate.at(PE).push_back((data_id));
                    num_known++;

                }
            }
            data_id++; //if the element must stay with me then check the next one
        }

        std::vector<int> num_import_from_procs(wsize);
        std::vector<int> import_from_procs;

        {
            nb_reqs = std::count_if(data_to_migrate.cbegin(), data_to_migrate.cend(), [](auto buf){return !buf.empty();});
            std::vector<MPI_Request> send_reqs(nb_reqs), rcv_reqs(wsize);
            std::vector<MPI_Status> statuses(wsize);

            int nb_neighbor = 0;
            for(int comm_pe = 0; comm_pe < wsize; ++comm_pe) {
                MPI_Irecv(&num_import_from_procs[comm_pe], 1, MPI_INT, comm_pe, 666, comm, &rcv_reqs[comm_pe]);
                int send_size = data_to_migrate.at(comm_pe).size();
                if (send_size) {
                    MPI_Isend(&send_size, 1, MPI_INT, comm_pe, 666, comm, &send_reqs[nb_neighbor]);
                    nb_neighbor++;
                }
            }

            MPI_Waitall(send_reqs.size(), &send_reqs.front(), MPI_STATUSES_IGNORE);
            MPI_Barrier(comm);
            for(int comm_pe = 0; comm_pe < wsize; ++comm_pe) {
                int flag; MPI_Status status;
                MPI_Test(&rcv_reqs[comm_pe], &flag, &status);
                if(!flag) MPI_Cancel(&rcv_reqs[comm_pe]);
            }

            import_from_procs.clear();
            for(int i = 0; i < wsize; ++i){
                if(num_import_from_procs[i] > 0) import_from_procs.push_back(i);
            }
        }

        return std::make_tuple(import_from_procs, data_to_migrate, num_import_from_procs);
    }
};

//TODO: alpha should be a struct holding the algorithm and the data?
template<class Data> void ZoltanLoadBalancer<Data>::load_balance(std::vector<Data> *_data) {


    lb_func(_data, zoltan_lb, true);

    std::tie(neighbors, neighboring_cells, cell_per_neighbors) = compute_neighborhood(_data);
}

template<class Data> std::vector<Data> ZoltanLoadBalancer<Data>::propagate(const std::vector<Data> &data,
                            int *nb_elements_recv, int *nb_elements_sent, double cell_size) {
    std::vector<Data> r;
    std::vector<std::vector<unsigned long>> i;

    std::tie(r,i) = zoltan_exchange_data<Data>(data, neighbors, cell_per_neighbors, neighboring_cells, this->datatype, this->world);

    assert(std::equal(i.begin(), i.end(), neighboring_cells.begin(), neighboring_cells.end()));

    return r;
    //return zoltan_exchange_data<Data>(zoltan_lb, data, &r, &s, this->datatype, this->world);
}


#endif //SPEC_ZOLTANLOADBALANCER_HPP
