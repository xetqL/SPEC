#!/bin/bash
rm bin/*
git pull origin master
rm CMakeCache.txt
~/cmake-3.11.0-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DBINARY_NAME=SPEC_NO_LB . && make
rm CMakeCache.txt
~/cmake-3.11.0-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DLB_METHOD=1 -DBINARY_NAME=SPEC_1 . && make
