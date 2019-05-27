//
// Created by xetql on 5/8/19.
//

#include "ULBA.hpp"
#include <iostream>
ULBA::ULBA(MPI_Comm world, GossipDatabase<double> *wirdb, double threshold, double alpha) :
        LoadBalancingApproach(world),
        wirdb(wirdb),
        threshold(threshold),
        alpha(alpha) {
    MPI_Comm_size(world, &P);
}

std::pair<LoadBalancingApproach::WorkloadShare, LoadBalancingApproach::WorkloadWeight> ULBA::compute_share(int rank) const {
    int N = 0, overloading = wirdb->zscore(rank) > threshold ? 1 : 0;
    auto data = wirdb->get_all_data();
    for(auto pute: data){
       std::cout << rank << " " << pute << std::endl;
    } 
    // Count the number of overloading processes
    MPI_Allreduce(&overloading, &N, 1, MPI_INT, MPI_SUM, this->world);

    double share = 1;
    
    if(overloading) {
        std::cout << "outlier..." << std::endl;
        share = (1 - this->alpha) * share;
        return {share, this->alpha};
    } else {
        std::cout << N << " outliers" << std::endl;
        share = (1 + this->alpha * ((double) N / (double) (P-N))) * share;
        return {share, 0.0};
    }

}

std::string ULBA::to_string() const {
    return "ULBA";
}
