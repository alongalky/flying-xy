// this file handles random number generation. 
// If the gsl library is giving you trouble for some reason you can go back to rand(), which always works but is less random.

#ifndef _RANDOM_H
#define _RANDOM_H

void initRandom();
double nextRand();
double nextNormalRand(double sigma, double mu);
void killRandom();

#endif