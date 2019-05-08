#!/bin/bash
#make clean
git pull origin with-zoltan

#rm CMakeCache.txt
#~/cmake-3.11.0-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DBINARY_NAME=SPEC_NO_LB . && make
#rm CMakeCache.txt
~/cmake-3.11.0-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DWITH_ZOLTAN=1 -DBINARY_NAME=SPEC_1 -DLOAB_BALANCING_CALL=Autonomic . && make
#rm CMakeCache.txt
#~/cmake-3.11.0-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DLB_METHOD=3 -DBINARY_NAME=SPEC_3 . && make
#rm CMakeCache.txt
#~/cmake-3.11.0-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DLB_METHOD=5 -DBINARY_NAME=SPEC_5 . && make
#rm CMakeCache.txt
#~/cmake-3.11.0-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release -DLB_METHOD=2 -DBINARY_NAME=SPEC_2 . && make
#rm CMakeCache.txt
#~/cmake-3.11.0-Linux-x86_64/bin/cmake -DPRODUCE_OUTPUTS=1 -DCMAKE_BUILD_TYPE=Release -DLB_METHOD=2 -DBINARY_NAME=SPEC_2_WITH_OUTPUTS . && make
#rm CMakeCache.txt
#~/cmake-3.11.0-Linux-x86_64/bin/cmake -DPRODUCE_OUTPUTS=1 -DCMAKE_BUILD_TYPE=Release -DLB_METHOD=3 -DBINARY_NAME=SPEC_3_WITH_OUTPUTS . && make
#rm CMakeCache.txt
#~/cmake-3.11.0-Linux-x86_64/bin/cmake -DPRODUCE_OUTPUTS=1 -DCMAKE_BUILD_TYPE=Release -DLB_METHOD=5 -DBINARY_NAME=SPEC_5_WITH_OUTPUTS . && make
