//
// Created by xetql on 5/3/19.
//

#ifndef SPEC_ZOLTANLOADBALANCER_HPP
#define SPEC_ZOLTANLOADBALANCER_HPP

#include <functional>
#include "LoadBalancer.hpp"
#include "ZoltanUtils.hpp"

template<class Data>
class ZoltanLoadBalancer : public LoadBalancer<Data> {
    using ZoltanCreateFunc = std::function<Zoltan_Struct* (bool, MPI_Comm)>;
    using ZoltanLBFunc     = std::function<void (std::vector<Data>*, Zoltan_Struct*, bool)>;

    Zoltan_Struct* zoltan_lb;
    ZoltanLBFunc lb_func;

public:
    ZoltanLoadBalancer(MPI_Comm world, MPI_Datatype datatype,
                       ZoltanCreateFunc zoltan_create_wrapper,
                       ZoltanLBFunc lb_func) :
            LoadBalancer<Data>(world, datatype, new StdApproach(world)), lb_func(lb_func) {
        zoltan_lb = zoltan_create_wrapper(true, world);
    }

    ZoltanLoadBalancer(MPI_Comm world, MPI_Datatype datatype, LoadBalancingApproach* approach,
                       ZoltanCreateFunc& zoltan_create_wrapper,
                       ZoltanLBFunc& lb_func) :
    LoadBalancer<Data>(world, datatype, approach), lb_func(lb_func) {
        zoltan_lb = zoltan_create_wrapper(true, world);
    }

    std::vector<Data> propagate(const std::vector<Data> &data,
                                int *nb_elements_recv, int *nb_elements_sent,
                                double cell_size) override;
private:
    void load_balance(std::vector<Data> *_data) override;
};

//TODO: alpha should be a struct holding the algorithm and the data?
template<class Data> void ZoltanLoadBalancer<Data>::load_balance(std::vector<Data> *_data) {
    float share, alpha;
    std::tie(share, alpha) = this->approach->compute_share(this->rank);
    int parts_num[2] = {this->rank, this->rank}, weight_per_obj[2] = {0, 1};
    float part_size[2] = {share, 1};
    //Zoltan_LB_Set_Part_Sizes(zoltan_lb, 1, 2, parts_num, weight_per_obj, part_size);
    lb_func(_data, zoltan_lb, true);
}

template<class Data> std::vector<Data> ZoltanLoadBalancer<Data>::propagate(const std::vector<Data> &data,
                            int *nb_elements_recv, int *nb_elements_sent, double cell_size) {
    return zoltan_exchange_data<Data>(zoltan_lb, data, nb_elements_recv, nb_elements_sent, this->datatype, this->world, cell_size);
}


#endif //SPEC_ZOLTANLOADBALANCER_HPP
