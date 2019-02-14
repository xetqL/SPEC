//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_UTILS_HPP
#define SPEC_UTILS_HPP

#include <utility>
#include "cell.hpp"

inline long long position_to_cell(int msx, int msy, const std::pair<int, int> & position) {
    return position.first + msx * position.second;
}

inline long long position_to_cell(int msx, int msy, const int x, const int y) {
    return x + msx * y;
}

inline std::pair<int, int> cell_to_global_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msx);
}

inline std::pair<int, int> cell_to_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msx);
}

inline std::pair<int, int> cell_to_local_position(int msx, int msy, std::tuple<int,int,int,int> bounding_box, long long position){
    int minx, maxx, miny, maxy; std::tie(minx, maxx, miny, maxy) = bounding_box;
    int gidx =  position % msx, gidy = (int) position / msx;
    return std::make_pair(gidx - minx,  gidy - miny);
}

template<class A>
std::tuple<int, int, int, int> get_bounding_box(
        const std::vector<A>& my_data,
        const std::vector<A>& remote_data) {
    int x, y, minx= std::numeric_limits<int>::max(), miny = std::numeric_limits<int>::max(), maxx=-1, maxy=-1;

    //create boundaries from vehicles
    for(const auto& v : my_data) {
        std::tie(x, y) = v.get_position_as_pair();
        minx = std::min(x, minx); miny = std::min(y, miny);
        maxx = std::max(x, maxx); maxy = std::max(y, maxy);
    }
    for(const auto& v : remote_data) {
        std::tie(x, y) = v.get_position_as_pair();
        minx = std::min(x, minx); miny = std::min(y, miny);
        maxx = std::max(x, maxx); maxy = std::max(y, maxy);
    }
    maxx++;maxy++;
    assert(minx >= 0);
    assert(miny >= 0);
    assert((maxx-minx) * (maxy-miny) >= (my_data.size() + remote_data.size()));
    return std::make_tuple(minx, maxx, miny, maxy);
}

#endif //SPEC_UTILS_HPP
