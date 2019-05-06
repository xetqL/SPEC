//
// Created by xetql on 5/3/19.
//

//
// Created by xetql on 2/11/19.
//

#include "Utils.hpp"

 long long position_to_cell(int msx, int msy, const std::pair<int, int> & position) {
    return position.first + msx * position.second;
}

 long long position_to_cell(int msx, int msy, const int x, const int y) {
    return x + msx * y;
}

 std::pair<int, int> cell_to_global_position(int msx, int msy, long long position) {
    return std::make_pair(position % msx, (int) position / msx);
}

 std::pair<int, int> cell_to_local_position(int msx, int msy, std::tuple<int,int,int,int> bounding_box, long long position) {
    int minx, maxx, miny, maxy; std::tie(minx, maxx, miny, maxy) = bounding_box;
    int gidx =  position % msx, gidy = (int) position / msx;
    return std::make_pair(gidx - minx,  gidy - miny);
}