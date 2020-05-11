// horrible windows stuff:
#include "stdafx.h"



// This file declares the parameters and gives them default values.
// These should be overwritten by the input file data.
#include "parameters.h"

// constants
char * formatNames[] = { "csv", "bin" };
double RSTAR_SQ = 1.0;

// file names
char LOG_FILE[256] = "log$.dat";
char MSD_FILE[256] = "msd.dat";
char VCOR_FILE[256] = "vcor.dat";
char INPUT_FILE[256] = "input";
char OUTPUT_DIR[256] = "output";

// output format
int OUTPUT_FORMAT = FORMAT_CSV;


// physical parameters
int NPART = 100;		// number of particles
int NITER = 1000;	// number of iterations

double NOISE = 1.0;		// noise strength
double BOX_X = 10.0;		// box X length
double BOX_Y = 10.0;		// box Y length
double V0 = 1.0;		// particle speed
double originalDT = 0.01;		// default sampling time
double ALIGN = 1.0;		// alignment strength
double DEFLECT = 1.0;		// deflect strength
double RSTAR = 0.1;		// interparticle interaction radius

// AUTOCODE START
int THREADS = 30;
int MAX_RUNTIME_MINUTES = 2000000;
double HYSTER_NOISE = -1.0;
int HYSTER = FALSE;
char INITIAL_DATA_FILE[256] = "";
double SIMTIME = 10.0;
double TIME_BETWEEN_SAMPLES = 1.0;
int NOISE_MODE = ANGULAR;
int NORMALIZE_NEIGHBORS = TRUE;
int PERIODIC_Y = TRUE;
int PERIODIC_X = TRUE;
double OBSTACLE_RADIUS = 1.0;
double OBSTACLEDEFLECTION = 1.0;
int OBSTACLENUM = 5;
double SOFT_WALLS_STRENGTH = 1.0;
int SOFT_WALLS = FALSE;
double RDEFLECTION = 0.1;
// AUTOCODE END

