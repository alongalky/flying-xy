#include "stdafx.h"
#include "collision.h"
#include "parameters.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include <omp.h>

typedef struct tagthreadInfo {
	pairIterator *it;
	double *dAlign;
	int *neighborCount;
	Vector *positions;
	double *prevVelocities;
} threadInfo;

// constants
int thread_num;
double **dAlign_partial;
int **neighborCount_partial;
pairIterator *iterators;

// Concurrency data structure
threadInfo *inf_array;

threadInfo initThreadInfo(pairIterator *it, double *dAlign, int *neighborCount, Vector *positions, double *prevVelocities){
	threadInfo inf;
	inf.it = it;
	inf.dAlign = dAlign;
	inf.neighborCount = neighborCount;
	inf.positions = positions;
	inf.prevVelocities = prevVelocities;
	return inf;
}
void initThreads(int _thread_num){
	time_printf("Initializing %d threads...\n", _thread_num);
	int count;
	thread_num = _thread_num;
	dAlign_partial = (double**)malloc(thread_num * sizeof(double*));
	neighborCount_partial = (int**)malloc(thread_num * sizeof(int*));
	iterators = (pairIterator*)malloc(thread_num * sizeof(pairIterator));

	inf_array = (threadInfo*)malloc(thread_num * sizeof(threadInfo));
	for (count = 0; count < thread_num; count++){
		dAlign_partial[count] = (double*)malloc(NPART*sizeof(double));
		neighborCount_partial[count] = (int*)malloc(NPART*sizeof(int));
		int first = (NPART * count) / thread_num;
		int next = (NPART * (count + 1)) / thread_num;
		iterators[count] = initPairIterator(first, next);
		inf_array[count] = initThreadInfo(&(iterators[count]), dAlign_partial[count], neighborCount_partial[count], NULL, NULL);
	}
}
void killThreads(){
	time_printf("Killing all threads!\n");
	int count;
	free(inf_array);
	// free all allocated memory
	for (count = 0; count < thread_num; count++){
		free(dAlign_partial[count]);
		free(neighborCount_partial[count]);
	}
	free(dAlign_partial);
	free(iterators);
}
void calculateNeighborInteraction(double *dAlign, int *neighborCount, int *max_neighborcount, Vector *positions, double *prevVelocities){
	int i, count;
	//time_printf("Entering parallel region with %d threads\n", thread_num);
#pragma omp parallel num_threads(thread_num)
	{
		int thread_id = omp_get_thread_num();
		int i, j;
		threadInfo *inf = &(inf_array[thread_id]);
		memset(inf->dAlign, 0, sizeof(double)*NPART);
		memset(inf->neighborCount, 0, sizeof(int)*NPART);
		inf->positions = positions;
		inf->prevVelocities = prevVelocities;
		resetPairIterator(inf->it);
		while (getNextNeighbors(&i, &j, inf->positions, inf->it) != -1){
			//time_printf("thread with first = %d found pair %d %d\n", inf->it->first, i, j);
			double temp = sin(inf->prevVelocities[j] - inf->prevVelocities[i]);
			inf->dAlign[i] += temp;
			inf->dAlign[j] -= temp;
			(inf->neighborCount[i])++;
			(inf->neighborCount[j])++;
			//time_printf("neighborcount[%d] = %d;\n", i, inf->neighborCount[i]);
		}
	}
	for (i = 0; i < NPART; i++){
		//time_printf("neighborcount[%d] = %d; max = %d\n", i, neighborCount[i], *max_neighborcount);
		for (count = 0; count < thread_num; count++){
			//time_printf("partial dAlign is %f\n", dAlign_partial[count][i]);
			dAlign[i] += dAlign_partial[count][i];
			neighborCount[i] += neighborCount_partial[count][i];
		}
		//time_printf("dAlign[%d] = %f\n", i, dAlign[i]);
	}
	*max_neighborcount = 0;
	for (i = 0; i < NPART; i++){
		if (neighborCount[i] > *max_neighborcount)
			*max_neighborcount = neighborCount[i];
	}
}

