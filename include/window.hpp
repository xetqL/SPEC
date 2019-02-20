//
// Created by xetql on 1/23/19.
//

#ifndef CA_ROAD_WINDOW_HPP
#define CA_ROAD_WINDOW_HPP

#include <deque>

template<class T>
struct SlidingWindow {
    std::deque<T> data_container;
    size_t window_max_size;

    SlidingWindow(size_t window_max_size, T init_value) : data_container(window_max_size, init_value), window_max_size(window_max_size) {};
    SlidingWindow(size_t window_max_size) : window_max_size(window_max_size) {};

    inline void add(const T &data) {
        if (data_container.size() < window_max_size)
            data_container.push_back(data); // add while not full
        else {                              // when full
            data_container.pop_front();     // delete oldest data
            data_container.push_back(data); // push new data
        }
    }

    typename std::deque<T>::iterator begin(){
        return data_container.begin();
    }

    typename std::deque<T>::iterator end(){
        return data_container.end();
    }

    unsigned long size(){
        return data_container.size();
    }

};

#endif //CA_ROAD_WINDOW_HPP
