//
// Created by xetql on 2/25/19.
//

#ifndef SPEC_TIME_HPP
#define SPEC_TIME_HPP
#include <mpi.h>

#define START_TIMING(v) \
auto v = MPI_Wtime()

#define RESTART_TIMING(v) \
v = MPI_Wtime() - v

#define STOP_TIMING(v) \
v = MPI_Wtime() - v

#define PAR_START_TIMING(v, comm) \
MPI_Barrier(comm); \
auto v = MPI_Wtime()

#define PAR_RESTART_TIMING(v, comm) \
MPI_Barrier(comm); \
v = MPI_Wtime() - v

#define PAR_STOP_TIMING(v, comm) \
MPI_Barrier(comm); \
v = MPI_Wtime() - v

#define CHECKPOINT_TIMING(v, u) \
auto u = MPI_Wtime() - v

#endif //SPEC_TIME_HPP
