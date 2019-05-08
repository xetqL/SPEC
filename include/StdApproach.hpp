//
// Created by xetql on 5/8/19.
//

#ifndef SPEC_STDAPPROACH_HPP
#define SPEC_STDAPPROACH_HPP

#include "LoadBalancingApproach.hpp"
#include "GossipDatabase.hpp"

class StdApproach : public LoadBalancingApproach {
public:
    explicit StdApproach(MPI_Comm world) ;

    std::pair<WorkloadShare, WorkloadWeight> compute_share(int rank) const override;

    std::string to_string() const override;
};

#endif //SPEC_STDAPPROACH_HPP
