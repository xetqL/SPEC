//
// Created by xetql on 5/8/19.
//

#ifndef SPEC_LOADBALANCINGAPPROACH_HPP
#define SPEC_LOADBALANCINGAPPROACH_HPP

#include <mpi.h>
#include <ostream>
class LoadBalancingApproach {
protected:
    MPI_Comm world;
public:
    using WorkloadShare  = float;
    using WorkloadWeight = float;

    explicit LoadBalancingApproach(MPI_Comm world);

    virtual std::pair<WorkloadShare, WorkloadWeight> compute_share(int rank) const = 0;

    virtual std::string to_string() const {
        return "Undefined approach.";
    }

    //virtual void adapt_weights(T* data) = 0;
};

#endif //SPEC_LOADBALANCINGAPPROACH_HPP
