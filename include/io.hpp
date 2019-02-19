//
// Created by xetql on 2/13/19.
//

#ifndef SPEC_IO_HPP
#define SPEC_IO_HPP


#include <ostream>
#include <vector>

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

#endif //SPEC_IO_HPP
