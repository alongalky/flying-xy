#include "stdafx.h"
#include "collision.h"
#include "parameters.h"
#include <stdlib.h>
#include <omp.h>

// ********* optimization parameter *********
int switch_axis_freq = 2000;

// ********* global variables *********
int *xOrder, *yOrder;
int current, next_neighbor;
unsigned long long perf_good_checks, perf_total_checks, perf_good_checks_iter, perf_total_checks_iter;
double last_x_efficiency;
double last_y_efficiency;
int sorting_axis;	// x = 0, y = 1
// uses iterNum for periodic debugging
extern int iterNum;

// ********* internal function headers *********
void inline swap(int *A, int i, int j);
void xSort(int *A, Vector *positions);
void ySort(int *A, Vector *positions);
inline int areNeighbors(int i_part, int j_part, Vector *positions);
void updatePerformanceStats();
// debugging
int debugVerifyOrder(Vector *positions);
int debugVerifyPairCount(Vector *positions);

// ********* internal function definitions *********
void inline swap(int *A, int i, int j){
	int temp = A[i];
	A[i] = A[j];
	A[j] = temp;
}
// insertion sort on X
void xSort(int *A, Vector *positions){
	int i,j;
	double x;
	int perf_counter = 0;
	for (i = 1; i <= NPART - 1; i++){
		x = positions[A[i]].x;
		j = i;
		while (j > 0 && positions[A[j - 1]].x > x){
			swap(A, j, j - 1);
			j--;
		}
	}
}
// insertion sort on Y
void ySort(int *A, Vector *positions){
	int i, j;
	double x;
	for (i = 1; i <= NPART - 1; i++){
		x = positions[A[i]].y;
		j = i;
		while (j > 0 && positions[A[j - 1]].y > x){
			swap(A, j, j - 1);
			j--;
		}
	}
	// test
	for (i = 1; i < NPART - 1; i++){
		if (positions[A[i - 1]].y > positions[A[i]].y){
			panic("disorder");
		}
	}
}
void updatePerformanceStats(){
	double efficiency = perf_total_checks_iter == 0 ? 0 : (double)perf_good_checks_iter / (double)perf_total_checks_iter;
	if (sorting_axis == 0){
		last_x_efficiency = efficiency;
	}
	else{
		last_y_efficiency = efficiency;
	}
}
// verifies exactly whether two particles are neighbors, using simulation index
inline int areNeighbors(int i_part, int j_part, Vector *positions){
	return distSq(positions[i_part], positions[j_part]) < RSTAR_SQ;
}

// ********* interface function definitions *********
void initCollision(Vector *positions){
	xOrder = (int *)malloc(sizeof(int) * NPART);
	yOrder = (int *)malloc(sizeof(int) * NPART);
	int i;
	for (i = 0; i < NPART; i++){
		xOrder[i] = i;
		yOrder[i] = i;
	}
	updatePositions(positions);
	perf_good_checks = 0;
	perf_total_checks = 0;
	last_x_efficiency = 1.0;
	last_y_efficiency = 1.0;
	sorting_axis = 0;
}
void killCollision(){
	time_printf("Sweep efficiency: %llu / %llu, %4.2f %%\n", perf_good_checks, perf_total_checks, (float)perf_good_checks*100 / perf_total_checks);
	free(xOrder);
	free(yOrder);
}
void updatePositions(Vector *positions){
	xSort(xOrder, positions);
	ySort(yOrder, positions);
	
	if (iterNum % switch_axis_freq == 1){
		if (!debugVerifyOrder(positions)){
			panic("wrong order");
		}
		updatePerformanceStats();
		sorting_axis = 1 - sorting_axis;	// switch between y and x, just for one iteration
	}
	if (iterNum % switch_axis_freq == 2){
		updatePerformanceStats();
		if (last_x_efficiency > last_y_efficiency){
			sorting_axis = 0;
		}
		else{
			sorting_axis = 1;
		}
		time_printf("iter%d x_eff = %4.2f %% ; y_eff = %4.2f %%. Choosing %s\n", iterNum, last_x_efficiency * 100, last_y_efficiency * 100, sorting_axis == 0 ? "x" : "y");
	}
	// update global counters
	perf_good_checks += perf_good_checks_iter;
	perf_total_checks += perf_total_checks_iter;

	current = 0;
	next_neighbor = 1;
	perf_good_checks_iter = 0;
	perf_total_checks_iter = 0;
}
pairIterator initPairIterator(int first, int last){
	pairIterator ret;
	ret.first			= first;
	ret.last			= last;
	ret.current			= first;
	ret.next_neighbor	= first + 1;
	ret.perf_good_checks = 0;
	ret.perf_total_checks = 0;
	return ret;
}
void resetPairIterator(pairIterator *it){
	it->current			= it->first;
	it->next_neighbor	= it->first + 1;
	it->perf_good_checks = 0;
	it->perf_total_checks = 0;
}
// iterates over possibly neighboring pairs
int getNextNeighbors(int *i_part, int *j_part, Vector *positions, pairIterator *it){
	int current = it->current;
	int last = it->last;
	if (sorting_axis == 0){	
		// sort on x
		if (it->current >= last){
			// finished scanning. Update performance statistics
#pragma omp atomic
			perf_good_checks_iter += it->perf_good_checks;
#pragma omp atomic
			perf_total_checks_iter += it->perf_total_checks;
			return -1;
		}
		double pos = positions[xOrder[it->current]].x + RSTAR;
		// checked all pairs, now check over the boundary
		if (it->next_neighbor == NPART){
			it->next_neighbor = 0;
		}
		// interactions without crossing the wall
		while (it->current < it->next_neighbor && positions[xOrder[it->next_neighbor]].x < pos){
			// intersect on x axis, let's check them:
			(it->perf_total_checks)++;
			if (areNeighbors(xOrder[it->current], xOrder[it->next_neighbor], positions)){
				(it->perf_good_checks)++;
				*i_part = xOrder[it->current];
				*j_part = xOrder[it->next_neighbor];
				it->next_neighbor = (it->next_neighbor + 1) % NPART;
				return 1;
			}
			it->next_neighbor = (it->next_neighbor + 1) % NPART;
		}
		// check interactions across the border
		if (pos > BOX_X){
			pos = pos - BOX_X;
			while (it->current > it->next_neighbor && positions[xOrder[it->next_neighbor]].x < pos){
				// intersect on x axis, let's check them:
				(it->perf_total_checks)++;
				if (areNeighbors(xOrder[it->current], xOrder[it->next_neighbor], positions)){
					(it->perf_good_checks)++;
					*i_part = xOrder[it->current];
					*j_part = xOrder[it->next_neighbor];
					it->next_neighbor = (it->next_neighbor + 1) % NPART;
					return 1;
				}
				it->next_neighbor = (it->next_neighbor + 1) % NPART;
			}
		}
		(it->current)++;
		it->next_neighbor = (it->current + 1) % NPART;
		return getNextNeighbors(i_part, j_part, positions, it);
	}
	else	{
		// sort on y
		if (it->current >= last){
			// finished scanning. Update performance statistics
#pragma omp atomic
			perf_good_checks_iter += it->perf_good_checks;
#pragma omp atomic
			perf_total_checks_iter += it->perf_total_checks;
			return -1;
		}
		double pos = positions[yOrder[it->current]].y + RSTAR;
		// checked all pairs, now check over the boundary
		if (it->next_neighbor == NPART){
			it->next_neighbor = 0;
		}
		// interactions without crossing the wall
		while (it->current < it->next_neighbor && positions[yOrder[it->next_neighbor]].y < pos){
			// intersect on y axis, let's check them:
			(it->perf_total_checks)++;
			if (areNeighbors(yOrder[it->current], yOrder[it->next_neighbor], positions)){
				(it->perf_good_checks)++;
				*i_part = yOrder[it->current];
				*j_part = yOrder[it->next_neighbor];
				it->next_neighbor = (it->next_neighbor + 1) % NPART;
				return 1;
			}
			it->next_neighbor = (it->next_neighbor + 1) % NPART;
		}
		// check interactions across the border
		if (pos > BOX_X){
			pos = pos - BOX_X;
			while (it->current > it->next_neighbor && positions[yOrder[it->next_neighbor]].y < pos){
				// intersect on y axis, let's check them:
				(it->perf_total_checks)++;
				if (areNeighbors(yOrder[it->current], yOrder[it->next_neighbor], positions)){
					(it->perf_good_checks)++;
					*i_part = yOrder[it->current];
					*j_part = yOrder[it->next_neighbor];
					it->next_neighbor = (it->next_neighbor + 1) % NPART;
					return 1;
				}
				it->next_neighbor = (it->next_neighbor + 1) % NPART;
			}
		}
		(it->current)++;
		it->next_neighbor = (it->current + 1) % NPART;
		return getNextNeighbors(i_part, j_part, positions, it);
	}
}

// ********* debugging functions *********
int debugVerifyOrder(Vector *positions){
	int i;
	for (i = 1; i < NPART; i++){
		if (positions[xOrder[i - 1]].x > positions[xOrder[i]].x){
			return 0;
		}

		if (positions[yOrder[i - 1]].y > positions[yOrder[i]].y){
			return 0;
		}
	}
	return 1;
}
int debugVerifyPairCount(Vector *positions){
	double RSTAR_SQ = SQR(RSTAR);
	int i, j, count = 0;
	for (i = 0; i < NPART - 1; i++){
		for (j = i + 1; j < NPART; j++){
			if (distSq(positions[i], positions[j]) < RSTAR_SQ){
				//time_printf("iter%d Found real pair %d %d\n", iterNum, i, j);
				count++;
			}
		}
	}
	return count;
}
