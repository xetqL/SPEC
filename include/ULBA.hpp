//
// Created by xetql on 5/8/19.
//

#ifndef SPEC_ULBA_HPP
#define SPEC_ULBA_HPP

#include "LoadBalancingApproach.hpp"
#include "GossipDatabase.hpp"

class ULBA : public LoadBalancingApproach {
    GossipDatabase<double> * const wirdb;
    const double threshold, alpha;
    int P;
public:

    ULBA(MPI_Comm world, GossipDatabase<double> *wirdb, double threshold, double alpha);
    std::pair<WorkloadShare, WorkloadWeight> compute_share(int rank) const override;

    std::string to_string() const override;
};

#endif //SPEC_ULBA_HPP
