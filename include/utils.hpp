//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_UTILS_HPP
#define SPEC_UTILS_HPP

#include <utility>

inline long long position_to_cell(int msx, int msy, const std::pair<int, int> & position) {
    return position.first + msx * position.second;
}

inline long long position_to_cell(int msx, int msy, const int x, const int y) {
    return x + msx * y;
}

inline std::pair<int, int> cell_to_position(int msx, int msy, long long position){
    return std::make_pair(position % msx, (int) position / msx);
}
#endif //SPEC_UTILS_HPP
