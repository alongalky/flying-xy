// horrible windows stuff:
#include "stdafx.h"

#include "random.h"
#include "utils.h"
#include "parameters.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

// GSL stuff. comment it out if you fail to link gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
gsl_rng * r;

#ifndef PI
#define PI (3.14159265358979323846)
#endif

void initRandom(){
	unsigned int seed = (unsigned int)(time(NULL) + (int)NOISE*200);
	time_printf("seed = %d\n", seed);
	//srand(seed);
	// select random number generator 
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, seed);
}
double nextRand(){
	return gsl_rng_uniform(r);
//	return (double)rand() / (double)RAND_MAX;
}
double nextNormalRand(double sigma, double mu){
	return gsl_ran_gaussian(r, sigma);
//	double rand1 = nextRand();
//	if (rand1 < 1e-100L) rand1 = 1e-100;
//	return (sigma*sqrt(-2 * log(rand1))*cos(2 * PI*nextRand())) + mu;
}

void killRandom(){
	gsl_rng_free(r);
}
