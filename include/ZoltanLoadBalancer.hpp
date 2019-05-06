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
    //using ZoltanCreateFunc = Zoltan_Struct* (bool, MPI_Comm);
    using ZoltanLBFunc     = std::function<void (std::vector<Data>*, Zoltan_Struct*, bool)>;
    //using ZoltanLBFunc     = void (std::vector<Data>*, Zoltan_Struct*, bool);
    Zoltan_Struct* zoltan_lb;
    //ZoltanLBFunc   lb_func;
    void (*lb_func)(std::vector<Data>*, Zoltan_Struct*, bool);
public:
    ZoltanLoadBalancer(MPI_Comm world, MPI_Datatype datatype, Zoltan_Struct* (*zoltan_create_wrapper)(bool, MPI_Comm),
                       void (*lb_func)(std::vector<Data>*, Zoltan_Struct*, bool)) :
    LoadBalancer<Data>(world, datatype), lb_func(lb_func) {
        zoltan_lb = zoltan_create_wrapper(true, world);
    }

    std::vector<Data> propagate(const std::vector<Data> &data,
                                int *nb_elements_recv, int *nb_elements_sent,
                                double cell_size) override;
private:
    void load_balance(std::vector<Data> *_data, double alpha) override;
};

template<class Data> void ZoltanLoadBalancer<Data>::load_balance(std::vector<Data> *_data, double alpha) {
    lb_func(_data, zoltan_lb, true);
}

template<class Data> std::vector<Data> ZoltanLoadBalancer<Data>::propagate(const std::vector<Data> &data,
                            int *nb_elements_recv, int *nb_elements_sent, double cell_size) {
    return zoltan_exchange_data<Data>(zoltan_lb, data, nb_elements_recv, nb_elements_sent, this->datatype, this->world, cell_size);
}


#endif //SPEC_ZOLTANLOADBALANCER_HPP
