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

 std::pair<int, int> cell_to_local_position(int msx, int msy, const std::tuple<int,int,int,int>& bounding_box, long long position) {
    int minx = std::get<0>(bounding_box), maxx = std::get<1>(bounding_box), miny = std::get<2>(bounding_box), maxy = std::get<3>(bounding_box);
    int gidx =  position % msx, gidy = (int) position / msx;
    return std::make_pair(gidx - minx,  gidy - miny);
}
void cell_to_local_position(int msx, int msy, const std::tuple<int,int,int,int>& bounding_box, long long position, int* pX, int* pY ) {
    int minx = std::get<0>(bounding_box), maxx = std::get<1>(bounding_box), miny = std::get<2>(bounding_box), maxy = std::get<3>(bounding_box);
    int gidx =  position % msx, gidy = (int) position / msx;
    *pX = gidx - minx;
    *pY = gidy - miny;
}

void to_global_position(int msx, int msy, long long position, int *x, int *y) {
*x = position % msx;
*y = (int) position / msx;
}