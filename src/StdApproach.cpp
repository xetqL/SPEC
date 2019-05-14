//
// Created by xetql on 5/8/19.
//

#include "StdApproach.hpp"

StdApproach::StdApproach(MPI_Comm world) : LoadBalancingApproach(world) {}

std::pair<LoadBalancingApproach::WorkloadShare, LoadBalancingApproach::WorkloadWeight>
StdApproach::compute_share(int rank) const {
    return {1, 0};
}

std::string StdApproach::to_string() const {
    return "Standard";
}


