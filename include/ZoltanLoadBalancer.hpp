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
    const ZoltanLBFunc lb_func;
    const GossipDatabase<unsigned long>* workdb;
public:
    ZoltanLoadBalancer(MPI_Comm world, MPI_Datatype datatype,
                       const GossipDatabase<unsigned long>* workdb,
                       const ZoltanCreateFunc &zoltan_create_wrapper,
                       const ZoltanLBFunc &lb_func) :
        LoadBalancer<Data>(world, datatype, new StdApproach(world)),
        lb_func(lb_func),
        workdb(workdb) {
        zoltan_lb = zoltan_create_wrapper(true, world);
    }

    ZoltanLoadBalancer(MPI_Comm world, MPI_Datatype datatype, const GossipDatabase<unsigned long>* workdb,
                       LoadBalancingApproach* approach,
                       const ZoltanCreateFunc& zoltan_create_wrapper,
                       const ZoltanLBFunc& lb_func) :
                       LoadBalancer<Data>(world, datatype, approach),
                       lb_func(lb_func),
                       workdb(workdb) {
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
    lb_func(_data, zoltan_lb, true);
}

template<class Data> std::vector<Data> ZoltanLoadBalancer<Data>::propagate(const std::vector<Data> &data,
                            int *nb_elements_recv, int *nb_elements_sent, double cell_size) {
    return zoltan_exchange_data<Data>(zoltan_lb, data, nb_elements_recv, nb_elements_sent, this->datatype, this->world, cell_size);
}


#endif //SPEC_ZOLTANLOADBALANCER_HPP
