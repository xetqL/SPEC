//
// Created by xetql on 6/27/19.
//

#ifndef SPEC_WEIGHTUPDATER_HPP
#define SPEC_WEIGHTUPDATER_HPP

#include <vector>
#include <cstdlib>
#include "ULBA.hpp"
#include "StdApproach.hpp"
#include "Cell.hpp"

template<class Data>
class WeightUpdater {
public:
    virtual void update_weight(
            std::vector<Data>* data, std::vector<unsigned long> potential_targets,
            const LoadBalancingApproach* approach,
            float W, float my_load) = 0;
};

template<class Data>
class RandomWeightUpdater : public WeightUpdater<Data> {
    const int target_cnt;
    template<class BidiIter >
    BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random) {
        size_t left = std::distance(begin, end);
        while (num_random--) {
            BidiIter r = begin;
            std::advance(r, std::rand() % left);
            std::swap(*begin, *r);
            ++begin;
            --left;
        }
        return begin;
    }
public:
    RandomWeightUpdater(int target_cnt) : target_cnt(target_cnt) {}

    void update_weight(
                    std::vector<Data>* data, std::vector<unsigned long> potential_targets,
                    const LoadBalancingApproach* approach,
                    float W, float my_load) override {
        float share, alpha;
        std::tie(share, alpha) = approach->compute_share(get_rank());

        //now I want alpha*W of it, i.e., remove (1-alpha)*W.
        //1. Compute the difference between my workload and the average
        if(alpha > 0) { //I should over-estimate my workload to get less cells
            share *= W; //use the estimated mean
            auto diff_with_share = my_load - share;

            float weight_amount = diff_with_share / target_cnt;
            //1. select randomly among potential_targets
            auto random_targets = random_unique(potential_targets.begin(), potential_targets.end(), target_cnt);
            //2. update the weight of data[i] with a part of value
            for (auto target = random_targets; random_targets != (random_targets+target_cnt); random_targets++) {
                data->at(*target).weight += (weight_amount);
            }
        }

    }
};

template<class Data>
class TypeOnlyWeightUpdater : public WeightUpdater<Data> {
    const int target_cnt;
    const std::function<bool(Data)> predicate;
public:
    TypeOnlyWeightUpdater(int target_cnt, std::function<bool(Data)>&& predicate) : target_cnt(target_cnt), predicate(predicate) {}

    void update_weight(
            std::vector<Data>* _data, std::vector<unsigned long> potential_targets,
            const LoadBalancingApproach* approach,
            float W, float my_load) override {
        float share, alpha;
        auto rank = get_rank();
        std::tie(share, alpha) = approach->compute_share(rank);
        std::vector<Data>& data = *_data;
        //now I want alpha*W of it, i.e., remove (1-alpha)*W.
        //1. Compute the difference between my workload and the average
        float weight_amount = 0.0f;
        std::vector<unsigned long> targets;
        std::copy_if(potential_targets.begin(), potential_targets.end(), std::back_inserter(targets), [this, &data](auto tar){return this->predicate(data[tar]);});
        if(alpha > 0) { //I should over-estimate my workload to get less cells
            share *= W; //use the estimated mean
            auto diff_with_share = my_load - share;
            weight_amount = diff_with_share / targets.size();
            std::cout << "Desired share for "<< get_rank() << " " << share << std::endl;
        }
        for(const auto id : targets) {
            Data& cell  = _data->at(id);
            cell.weight = weight_amount;
        }
    }
};
#endif //SPEC_WEIGHTUPDATER_HPP
