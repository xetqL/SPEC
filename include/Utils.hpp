//
// Created by xetql on 2/11/19.
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <cmath>
#include <cassert>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <map>
#include <vector>
#include <mpi.h>

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    std::size_t total_size = 0;
    for (const auto& sub : v)
        total_size += sub.size(); // I wish there was a transform_accumulate
    std::vector<T> result;
    result.reserve(total_size);
    for (const auto& sub : v)
        result.insert(result.end(), sub.begin(), sub.end());
    return result;
}

inline int get_rank(){
    int r;
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    return r;
}

long long position_to_cell(int msx, int msy, const std::pair<int, int> & position) ;

long long position_to_cell(int msx, int msy, const int x, const int y);

std::pair<int, int> cell_to_global_position(int msx, int msy, long long position);
void cell_to_global_position(int msx, int msy, long long position, int *px, int *py);


std::pair<int, int> cell_to_local_position(int msx, int msy, std::tuple<int,int,int,int> bounding_box, long long position);

template<typename K, typename V>
std::vector<std::pair<K,V>> mapToVector(const std::unordered_map<K,V> &map) {
    return std::vector<std::pair<K,V>>(map.begin(), map.end());
}

template<typename K, typename V>
std::vector<std::pair<K,V>> mapToVector(const std::map<K,V> &map) {
    return std::vector<std::pair<K,V>>(map.begin(), map.end());
}

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
           // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}

template<class IntegerType=unsigned int>
inline std::pair<IntegerType, IntegerType> cell_to_position(IntegerType msx, IntegerType msy, long long position) {
    return std::make_pair(position % msx, (IntegerType) position / msx);
}

template<class A>
std::tuple<int, int, int, int> get_bounding_box(
        const std::vector<A>& my_data,
        const std::vector<A>& remote_data) {
    int x, y, minx= std::numeric_limits<int>::max(), miny = std::numeric_limits<int>::max(), maxx=-1, maxy=-1;

    // create boundaries from vehicles
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
    assert((unsigned int)  (maxx-minx) * (maxy-miny) >= (my_data.size() + remote_data.size()));
    return std::make_tuple(minx, maxx, miny, maxy);
}


template<class A>
std::tuple<int, int, int, int> get_bounding_box(const std::vector<A>& my_data) {
    int x, y, minx= std::numeric_limits<int>::max(), miny = std::numeric_limits<int>::max(), maxx=-1, maxy=-1;

    // create boundaries from vehicles
    for(const auto& v : my_data) {
        std::tie(x, y) = v.get_position_as_pair();
        minx = std::min(x, minx); miny = std::min(y, miny);
        maxx = std::max(x, maxx); maxy = std::max(y, maxy);
    }

    assert(minx >= 0);
    assert(miny >= 0);
    assert((unsigned int)  (maxx-minx+1) * (maxy-miny+1) >= (my_data.size()));
    //return std::make_tuple(std::max(0, minx), std::min(msx, maxx+2), std::max(0, miny-1), std::min(msy, maxy+1));
    return std::make_tuple(minx, maxx, miny, maxy);
}

template<class A>
std::tuple<int, int, int, int> update_bounding_box(const std::vector<A>& my_data, const std::tuple<int, int, int, int>& bbox) {
    int x1,y1,x2,y2,x,y;
    std::tie(x1, x2, y1, y2) = bbox;

    // create boundaries from vehicles
    for(const auto& v : my_data) {
        v.get_position(&x, &y);
        x1 = std::min(x, x1); y1 = std::min(y, y1);
        x2 = std::max(x, x2); y2 = std::max(y, y2);
    }

    //return std::make_tuple(std::max(0, minx), std::min(msx, maxx+2), std::max(0, miny-1), std::min(msy, maxy+1));
    return std::make_tuple(x1, x2, y1, y2);
}
/*
template<class A>
std::tuple<int, int, int, int> get_bounding_box(const std::vector<A>& my_data) {
    int x, y, minx= std::numeric_limits<int>::max(), miny = std::numeric_limits<int>::max(), maxx=-1, maxy=-1;
    for(const auto& v : my_data) { //create boundaries from vehicles
        std::tie(x, y) = v.get_position_as_pair();
        minx = std::min(x, minx); miny = std::min(y, miny);
        maxx = std::max(x, maxx); maxy = std::max(y, maxy);
    }
    maxx++; maxy++;
    assert(minx >= 0); assert(miny >= 0);
    assert((maxx-minx) * (maxy-miny) >= (my_data.size()));
    return std::make_tuple(minx, maxx, miny, maxy);
}
*/
template<typename Realtype, typename ContainerA>
Realtype get_slope(const ContainerA& y) {

    std::vector<Realtype> x(y.size());
    std::iota(x.begin(),x.end(), 0);

    int n = x.size();

    Realtype avgX = std::accumulate(x.begin(), x.end(), 0.0) / n;
    Realtype avgY = std::accumulate(y.begin(), y.end(), 0.0) / n;

    Realtype numerator = 0.0;
    Realtype denominator = 0.0;

    for(int i=0; i < n; ++i) {
        numerator   += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    return denominator == 0 ? (Realtype) 0 : numerator / denominator;
}

template<typename Realtype, typename Iter>
Realtype get_slope(const Iter& beginy, const Iter& endy) {

    if(std::distance(beginy, endy) < 1) return (Realtype) 0.0;

    const auto n = std::distance(beginy, endy);

    std::vector<Realtype> x(n);
    std::iota(x.begin(), x.end(), 0);

    Realtype avgX = std::accumulate(x.begin(), x.end(), 0.0) / (Realtype) n;
    Realtype avgY = std::accumulate(beginy, endy, 0.0) / (Realtype) n;
    Realtype numerator = 0.0, denominator = 0.0;

    for(int i=0; i < n; ++i) {
        numerator   += (x[i] - avgX) * (*(beginy+i) - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    return denominator == 0 ? (Realtype) 0 : numerator / denominator;
}



namespace stats
{
template<class RealType, class Iter>
inline RealType mean(Iter b, Iter e) {
    const long N = std::distance(b, e);
    if(N > 0)
        return std::accumulate(b, e, (RealType) 0.0) / N;
    else
        return (RealType) 0.0;
}

template<class RealType, class Iter>
inline RealType skewness(Iter b, Iter e) {
    long N = std::distance(b, e);
    RealType _mean = mean<RealType, Iter>(b,e);
    RealType diff1 = 0.0;
    std::for_each(b,e, [&diff1, mean=_mean](auto x){diff1 += std::pow(x - mean, 3.0);});
    diff1 /= N;
    RealType diff2 = 0.0;
    std::for_each(b,e, [&diff2, mean=_mean](auto x){diff2 += std::pow(x - mean, 2.0);});
    diff2 /= (N-1);
    return diff1 / std::pow(diff2, 3.0/2.0);
}

template<class Realtype, class Iter>
Realtype median(Iter begin, Iter  end) {
    std::vector<typename Iter::value_type> tmp(begin, end);
    std::sort(tmp.begin(), tmp.end());
    if(std::distance(begin, end) == 0)
        return (Realtype) 0.0;
    else
        return (Realtype) tmp[(int) (tmp.size() / 2)];
}


} //end of namespace: stats

#endif //SPEC_UTILS_HPP
