//
// Created by xetql on 2/13/19.
//

#ifndef SPEC_IO_HPP
#define SPEC_IO_HPP


#include <ostream>
#include <vector>
#include <deque>
#include <iostream>

template<class A>
std::ostream &operator<<(std::ostream &os, const std::vector<A> &data) {
    const long sz = data.size();
    os << "" ;
    for (long i = 0; i < sz-1; ++i) {
        os << data.at(i)<<", ";
    }
    os << data.at(sz-1);
    return os;
}

template<class A, class B>
std::ostream &operator<<(std::ostream &os, const std::pair<A,B> &data) {

    os << "("<<data.first<<","<<data.second<<")";
    return os;
}

template<class A>
std::ostream &operator<<(std::ostream &os, const std::deque<A> &data) {
    const long sz = data.size();
    os << "" ;
    for (long i = 0; i < sz-1; ++i) {
        os << data.at(i)<<", ";
    }
    os << data.at(sz-1);
    return os;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
static std::ostream &operator<<(std::ostream &os, const std::tuple<int, int, int , int>& bbox) {
    os << std::get<0>(bbox) << ", " << std::get<1>(bbox) << ", " << std::get<2>(bbox) << ", " << std::get<3>(bbox);
    return os;
}

static void print_bbox(std::tuple<int, int, int , int> bbox) {
    std::cout << std::get<0>(bbox) << ", " << std::get<1>(bbox) << ", " << std::get<2>(bbox) << ", " << std::get<3>(bbox) << std::endl;
}

static void print_bbox(int rank, std::tuple<int, int, int , int> bbox) {
    std::cout << rank << " " << std::get<0>(bbox) << ", " << std::get<1>(bbox) << ", " << std::get<2>(bbox) << ", " << std::get<3>(bbox) << std::endl;
}
#pragma GCC diagnostic pop
#endif //SPEC_IO_HPP
