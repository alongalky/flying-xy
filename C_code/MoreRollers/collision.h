// This library takes care of the iteration over neighbor particle pairs (using Sweep and Prune algorithm)

#ifndef __COLLISION_H
#define __COLLISION_H

#include "utils.h"

typedef  struct tagpairIterator{
	int first; // inclusive
	int last;  // exclusive
	int current;
	int next_neighbor;
	// To assess performance (calculate hit rate):
	int perf_good_checks;
	int perf_total_checks;
} pairIterator;

// module constructor and destructor
void initCollision(Vector *positions);
void killCollision();

// interface
void updatePositions(Vector *positions);
int getNextNeighbors(int *i, int *j, Vector *positions, pairIterator *it);	// iterates over neighboring pairs
pairIterator initPairIterator(int first, int last);
void resetPairIterator(pairIterator *it);


// debugging
// ---------------

// runs a naive pair count to see how many pairs there are in total (slow, only for debugging)
int debugVerifyPairCount(Vector *positions);

#endif