// This library takes care of the multithreaded torque calculation using OMP.

#ifndef _PARALLEL_H
#define _PARALLEL_H


#include "parameters.h"
#include "utils.h"

void calculateNeighborInteraction(double *dAlign, int *neighborCount, int *max_neighbor_count, Vector *positions, double *prevVelocities);

void initThreads(int _thread_num);
void killThreads();

#endif