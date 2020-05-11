// Roller.cpp : Defines the entry point for the console application.
//

// horrible windows stuff:
#include "stdafx.h"
#ifdef _USING_WINDOWS_
#include <conio.h>
#include <direct.h>

#define MKDIR(x) _mkdir(x)
#else
#include <sys/stat.h> // mkdir
#define MKDIR(x) mkdir(x,0777)
#endif

// my code
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "parameters.h"
#include "utils.h"
#include "log.h"
#include "random.h"
#include "profiler.h"
#include "collision.h"
#include "parallel.h"
// global variables
Vector *positions;
Vector *obstaclePositions;
double *velocities;	// only angular
double *prevVelocities;	// previous time step, for torque calculations
double elapsedTime;
double dt;				// adaptative time steps: this guys changes all the time
double noise;				// adaptative noise, for hysterisis calculations
int iterations_in_memory;
int iteration_size;
int iterNum;
char outputDir_LongTime[256];
int resume_previous_simulation = FALSE;
char prev_simulation_lastfile[256] = "";
int filesWrittenLongTime;

// global observables and statistical measures
Vector **positionsHistoryLongTime;
double **velocitiesHistoryLongTime;
double *timeHistoryLongTime;

// profiling timers
timer_type noiseCounter;
timer_type alignmentCounter1;
timer_type alignmentCounter2;
timer_type torqueCounterMem;
timer_type torqueNeighborCount;


// function headers
void initialize();					// initializes random variables, log files, and puts particles in random positions and velocities 
void simulate();					// runs the entire simulation
void cleanup();						// cleans up
void MoveParticles();		// moves all the particles along their direction of motion, for one dt
void applyTorques();
void printSummary();
double wallDecayFunction(double r);
void updateHistory(Vector **positionsHistory, double **velocitiesHistory, double *timeHistory, int bufferIterNum);
void clearHistory();
void makeRandomState(Vector *positions, double *velocities);
void makeRandomObstacles(Vector *obstaclePositions);
void initArrays();

// main program
int main(int argc, char* argv[])
{
	time_printf("Welcome to Alon\'s Roller simulation!\n");
	parse_commandLineArgs(argc, argv);
	initialize();
	simulate();
	cleanup();

	printSummary();

	// exit cleanly
	time_printf("bye!\n");
	//_getch();
	return 0;
}

void initialize(){
	char filename[256];

	time_printf("initalizing\n");

	// read parameter file
	readParameterFile(INPUT_FILE);

	if (PERIODIC_Y && SOFT_WALLS){
		time_printf("Error: doesn\'t make sense to have both periodic Y boundary and soft walls.");
			exit(1);
	}
	// set up output directories
	strcpy(outputDir_LongTime, OUTPUT_DIR);	 strcat(outputDir_LongTime, separator);	strcat(outputDir_LongTime, "longtime"); strcat(outputDir_LongTime, separator);
	MKDIR(OUTPUT_DIR);
	MKDIR(outputDir_LongTime);

	// calculate number of iterations which fit in memory
	iteration_size = 3 * sizeof(double) * NPART;	// memory space a single configuration occupies, in bytes
	iteration_size += 1 * sizeof(double);			// plus <elapsed_time> in each frame
	iterations_in_memory = (MAX_MEMORY * 1024) / iteration_size;
	time_printf("Iterations in each output file: %d\n", iterations_in_memory);

	// initialize random number generator
	initRandom();
	// init logging
	initLog();

	// allocate position and velocity vectors
	initArrays();

	// initalize particle positions:
	filesWrittenLongTime = 0;
	if (resume_previous_simulation == TRUE){
		filesWrittenLongTime	= find_last_logfile(outputDir_LongTime, filename)	+ 1;
		time_printf("debug: filesWrittenLongTime = %d\n", filesWrittenLongTime);
		int goodLastFrame = -1;
		if (filesWrittenLongTime != 0){
			goodLastFrame = readLastFrame(positions, velocities, obstaclePositions, &elapsedTime, filename);
			time_printf("filename = %s\n", filename);
		}
		if (goodLastFrame == -1){
			time_printf("Warning: couldn\'t find previous simulation data, starting from random configuration or initial data file (if given)\n");
			filesWrittenLongTime = 0;
		}
		
	}
	if (filesWrittenLongTime == 0) {				// not resuming simulation
		if (STRCMPI(INITIAL_DATA_FILE,"") != 0){	// load initial configuration from file
			if (readInitConfFile(INITIAL_DATA_FILE, positions, velocities, obstaclePositions) == -1){
				cleanup();
				panic("No luck with input file");
			}
			time_printf("Successfully read data from configuration file!\n");
			//time_printf("positions1 = %f %f\n", positions[0].x, positions[0].y);
		}
		else{
			makeRandomState(positions, velocities);
		}
		makeRandomObstacles(obstaclePositions);
		// write obstacles to files
		logObstacles(obstaclePositions, outputDir_LongTime);	// write obstacle locations to file

		elapsedTime = 0;
		filesWrittenLongTime = 0;
	}
	// initalize time
	dt = originalDT;
	iterNum = 0;

	// init profiling
	noiseCounter = newTimer();
	torqueNeighborCount = newTimer();
	alignmentCounter1 = newTimer();
	alignmentCounter2 = newTimer();
	torqueCounterMem = newTimer();

	// initialize collision detector
	initCollision(positions);

	// init hysterysis
	noise = NOISE;
	
	// init Neighbor counting threads :)
	if (THREADS > 0){
		initThreads(THREADS);
		time_printf("number of CPUs is: %d, number of threads = %d\n", omp_get_num_procs(), THREADS);
	}
	
}
void moveParticles(){
	int i;
	Vector newPos;
	for (i = 0; i < NPART; i++){
		// initalize position
		newPos = positions[i];
		newPos.x += dt*V0*cos(velocities[i]);
		newPos.y += dt*V0*sin(velocities[i]);

		// apply boundary conditions on X
		if (PERIODIC_X){
			// apply periodic boundary conditions on x axis
			while (newPos.x < 0) newPos.x += BOX_X;
			while (newPos.x > BOX_X) newPos.x -= BOX_X;
		}
		else{
			// apply (hard) reflecting boundary conditions on X axis
			if (newPos.x < 0){
				newPos.x = 0;								// return particle to the box
				velocities[i] = XMIRROR(velocities[i]);		// and flip its velocity
			}
			if (newPos.x > BOX_X){
				newPos.x = BOX_X;							// return particle to the box
				velocities[i] = XMIRROR(velocities[i]);		// and flip its velocity
			}

		}
		// apply boundary conditions on Y
		if (PERIODIC_Y){
			// apply periodic boundary conditions on Y axis
			while (newPos.y < 0) newPos.y += BOX_Y;
			while (newPos.y > BOX_Y) newPos.y -= BOX_Y;
		}
		else{
			// apply (hard) reflecting boundary conditions on Y axis
			if (newPos.y < 0){
				newPos.y = 0;							// return particle to the box
				velocities[i] = YMIRROR(velocities[i]);		// and flip its velocity
			}
			if (newPos.y > BOX_Y){
				newPos.y = BOX_Y;							// return particle to the box
				velocities[i] = YMIRROR(velocities[i]);		// and flip its velocity
			}
		}

		if (newPos.x > BOX_X || newPos.x < 0 || newPos.y > BOX_Y || newPos.y < 0){
			time_printf("Error! particle %d escaped the box!\n", i);
			time_printf("position (%f,%f)\n", newPos.x, newPos.y);
			exit(1);
		}

		// update array
		positions[i] = newPos;
	}
	updatePositions(positions); // for collision detection
}
void applyTorques(){
	int i, j, a=0;
	int max_neighborcount = 0;
	double dNoise;
	//double dDeflect_i, dDeflect_j, dObstacle;
	double temp;
	double noiseFactor;

	startTimer(&torqueCounterMem);

	double *dAlign = (double*)(malloc(sizeof(double) * NPART));
	memset(dAlign, 0, sizeof(double)*NPART);
	// for Vicksek behavior and for determining the value of "dt"
	int *neighborCount = (int*)(malloc(sizeof(int) * NPART));
	memset(neighborCount, 0, sizeof(int)*NPART);
	// keep old velocities values
	memcpy(prevVelocities, velocities, sizeof(double)*NPART);
	stopTimer(&torqueCounterMem);

	// Loop over neighbors, calculate dt, and calculate alignment
	// TODO: deflection should probably also go here
	startTimer(&alignmentCounter1);
	calculateNeighborInteraction(dAlign, neighborCount, &max_neighborcount, positions, prevVelocities);
	stopTimer(&alignmentCounter1);

	startTimer(&alignmentCounter2);
	// calculate dt by making sure (max_alignment_torque * dt) < original_dt
	if (NORMALIZE_NEIGHBORS == FALSE && max_neighborcount > 0){	// not relevant for Vicsek
		double maxAlign = 0.0;
		for (i = 0; i < NPART; i++){
			maxAlign = MAX(maxAlign, abs(dAlign[i]));
		}
		maxAlign *= ALIGN;
		
		if (maxAlign > originalDT){
			dt = originalDT / maxAlign;
		}
		else{
			dt = originalDT;
		}
		if (dt * noise > 1.0){
			time_printf("Warning: dt * D = %f; neighbor count = %d\n", dt * NOISE, max_neighborcount);
			dt = originalDT / noise;
		}
	}
	stopTimer(&alignmentCounter2);
	noiseFactor = sqrt(2 * noise * dt);
	for (i = 0; i < NPART; i++){
		// add noise term
		startTimer(&noiseCounter);
		if (NOISE_MODE == ANGULAR){
			dNoise = noiseFactor * nextNormalRand(1, 0);
		}
		else if (NOISE_MODE == VECTORIAL){
			panic("Vectorial noise mode not supported");
			dNoise = noiseFactor *(nextNormalRand(1, 0)*cos(prevVelocities[i]) + nextNormalRand(1, 0)*sin(prevVelocities[i]));
		}
		velocities[i] = putInUnitCircle(velocities[i] + dNoise);
		/*
		// add soft wall term
		double dSoftWall;
		if (SOFT_WALLS){	
			// deflection from top wall:
			dSoftWall = dt* SOFT_WALLS_STRENGTH * (-1) * cos(prevVelocities[i]) * wallDecayFunction(BOX_Y - positions[i].y);
			// bottom wall:
			dSoftWall += dt*SOFT_WALLS_STRENGTH  * cos(prevVelocities[i]) *wallDecayFunction(positions[i].y);
			// update velocity:
			velocities[i] = putInUnitCircle(velocities[i] + dSoftWall);
		}
		*/

		// add alignment term
		if (NORMALIZE_NEIGHBORS && neighborCount[i] > 0){
			velocities[i] = putInUnitCircle(velocities[i] + ((ALIGN*dt*dAlign[i]) / (double)neighborCount[i]));
		}
		else{
			velocities[i] = putInUnitCircle(velocities[i] + (ALIGN*dt*dAlign[i]));
		}
		// debug
		if (neighborCount < 0){
			panic("Error! negative neighbor count.\n");
		}
		// add obstacle deflection term
		/*
		double dSquared
		if (OBSTACLEDEFLECTION > 0){
			dObstacle = 0;
			for (a = 0; a < OBSTACLENUM; a++){
				dSquared = distSq(positions[i], obstaclePositions[a]);
				if (dSquared < SQR(OBSTACLE_RADIUS)){
					Vector r_ia = vecSub(obstaclePositions[a], positions[i]);
					double r_ia_angle = vecDirection(r_ia);
					// TODO: ask Denis about units etc.
					dObstacle += OBSTACLEDEFLECTION * dt * sin(r_ia_angle - prevVelocities[i]);
				}
			}
			velocities[i] = putInUnitCircle(velocities[i] + dObstacle);
		}
		*/
	}
	stopTimer(&noiseCounter);
	startTimer(&torqueCounterMem);
	free(neighborCount);
	free(dAlign);
	stopTimer(&torqueCounterMem);

}
void simulate(){
	int progress_update_frequency = 500;// MAX((int)(NITER / 30), 1);
	char ch;
	// variable that control updating the user about the simulation process
	clock_t max_time_between_updates = 1000 * 10;	// maximum time between two user updates in miliseconds
	clock_t startTime = getTime();
	time_t timeout_startTime = time(NULL);
	int minute_diff;
	int speed_assessement_every_x_iterations = 100;
	double average_iterTime = -1.0;
	timer_type lastXiterTimer = newTimer();
	// profiling clocks
	timer_type applyTorqueTimer = newTimer();
	timer_type moveParticlesTimer = newTimer();
	timer_type updateHistoryTimer = newTimer();
	timer_type logHistoryTimer = newTimer();
	timer_type totalTime = newTimer();

	// control long and short time logging frequency
	double lastLogTime_LongTime = elapsedTime;
	double nextLogTime_LongTime = ((int)(elapsedTime / TIME_BETWEEN_SAMPLES) + 1) * TIME_BETWEEN_SAMPLES;
	// TODO: make this configurable
	// log 0.5 full unit time every 10 unit times, to get a local picture of what's going on.

	int bufferIterNumLongTime = 0;


	startTimer(&totalTime);
	startTimer(&lastXiterTimer);
	for (; iterNum < NITER || elapsedTime < SIMTIME; iterNum++){
#ifdef _USING_WINDOWS_
		if (kbhit()){
			if ((ch =_getch()) == 'S'){
				ungetc(ch, stdin);
				time_printf("Interrupted by user after %d iterations, exiting gracefully\n", iterNum);
				break;
			}
		}
#endif
		minute_diff = (int)(difftime(time(NULL), timeout_startTime) / 60);
		if (minute_diff > MAX_RUNTIME_MINUTES){
			time_printf("Timed out after %d minutes and %d iterations, exiting gracefully\n", MAX_RUNTIME_MINUTES, iterNum);
			break;
		}
		// this is the actual simulation
		startTimer(&applyTorqueTimer);
		applyTorques();		// apply torques and update dt
		stopTimer(&applyTorqueTimer);

		startTimer(&moveParticlesTimer);
		moveParticles();
		stopTimer(&moveParticlesTimer);

		if (HYSTER){
			noise = noiseFunction(elapsedTime);
		}
		elapsedTime += dt;

		// Sample frame into "Long Time" every <time_between_logs_LongTime>
		if (elapsedTime > nextLogTime_LongTime){
			startTimer(&updateHistoryTimer);
			updateHistory(positionsHistoryLongTime, velocitiesHistoryLongTime, timeHistoryLongTime, bufferIterNumLongTime);
			stopTimer(&updateHistoryTimer);
			bufferIterNumLongTime++;
			nextLogTime_LongTime = ((int)(elapsedTime / TIME_BETWEEN_SAMPLES) + 1) * TIME_BETWEEN_SAMPLES;
			lastLogTime_LongTime = elapsedTime;
		}
		if (bufferIterNumLongTime == iterations_in_memory){		// time to write a LongTime file
			startTimer(&updateHistoryTimer);
			logHistory(positionsHistoryLongTime, velocitiesHistoryLongTime, timeHistoryLongTime, bufferIterNumLongTime, filesWrittenLongTime, outputDir_LongTime);
			bufferIterNumLongTime = 0;
			filesWrittenLongTime++;
			stopTimer(&updateHistoryTimer);
		}
		// every X frames measure how quickly we're going, to estimate ETA
		if (iterNum % speed_assessement_every_x_iterations == 0){
			stopTimer(&lastXiterTimer);
			average_iterTime = getTimerValue(lastXiterTimer) / speed_assessement_every_x_iterations;
			lastXiterTimer = newTimer();
		}
		if (iterNum > 0 && iterNum % progress_update_frequency == 0 || getTime() - startTime > max_time_between_updates){
			int iterations_left = MAX((NITER - iterNum), (int)((SIMTIME - elapsedTime) / dt));
			int secs_left = (int)(iterations_left*average_iterTime / 1000.0);
			int min_left = (int)(secs_left / 60.0);
			int hours_left = (int)(min_left / 60.0);
			time_printf("%4.1f%%: %10d / %10d iterations. Time = %5.1f, time left: %3d:%02d:%02d, %d iteration time: %f sec, Noise = %5.3f, dt = %5.4f\n", 
				(iterNum*100.0f) / (iterations_left + iterNum), iterNum, iterations_left + iterNum, elapsedTime,
				hours_left, min_left % 60, secs_left % 60, speed_assessement_every_x_iterations, average_iterTime*speed_assessement_every_x_iterations / 1000.0, noise, dt);
			startTime = getTime();
		}

	}

	startTimer(&updateHistoryTimer);
	logHistory(positionsHistoryLongTime, velocitiesHistoryLongTime, timeHistoryLongTime, bufferIterNumLongTime, filesWrittenLongTime, outputDir_LongTime);
	stopTimer(&updateHistoryTimer);

	stopTimer(&totalTime);
	// print profiling information:
	double total_time = getTimerValue(totalTime);
	double total_torquetime = getTimerValue(applyTorqueTimer);
	time_printf("Total simulation time:    %10.3f sec\n", total_time / 1000.0);
	time_printf("Total applyTorque time:   %10.3f sec  %6.3f %%\n", total_torquetime / 1000.0, total_torquetime * 100 / total_time);
	time_printf("     alignment1:          %10.3f sec  %6.3f %%\n", getTimerValue(alignmentCounter1) / 1000.0, getTimerValue(alignmentCounter1) * 100 / total_time);
	time_printf("     alignment2:          %10.3f sec  %6.3f %%\n", getTimerValue(alignmentCounter2) / 1000.0, getTimerValue(alignmentCounter2) * 100 / total_time);
	time_printf("     memory operations:   %10.3f sec  %6.3f %%\n", getTimerValue(torqueCounterMem) / 1000.0, getTimerValue(torqueCounterMem) * 100 / total_time);
	time_printf("     noise:               %10.3f sec  %6.3f %%\n", getTimerValue(noiseCounter) / 1000.0, getTimerValue(noiseCounter) * 100 / total_time);
	time_printf("Total moveParticles time: %10.3f sec  %6.3f %%\n", getTimerValue(moveParticlesTimer) / 1000.0, getTimerValue(moveParticlesTimer) * 100 / total_time);
	time_printf("Total updateHistory time: %10.3f sec  %6.3f %%\n", getTimerValue(updateHistoryTimer) / 1000.0, getTimerValue(updateHistoryTimer) * 100 / total_time);
	time_printf("Total logHistory time:    %10.3f sec  %6.3f %%\n", getTimerValue(logHistoryTimer) / 1000.0, getTimerValue(logHistoryTimer) * 100 / total_time);
}
void cleanup(){
	int i;
	// free 1d array memory
	free(positions);
	free(velocities);
	free(prevVelocities);
	free(obstaclePositions);

	// free 2d array memory
	for (i = 0; i < iterations_in_memory; i++){
		free(positionsHistoryLongTime[i]);
		free(velocitiesHistoryLongTime[i]);
	}
	free(positionsHistoryLongTime);
	free(velocitiesHistoryLongTime);
	free(timeHistoryLongTime);

	killRandom();
	killLog();
	killCollision();
	if (THREADS > 0)
		killThreads();
}
void printSummary(){
	logSummary();	// print logging time
}
double wallDecayFunction(double r){
	if (r < 0){
		printf("Error! negative distance\n");
		exit(1);
	}
	double NormalizedDistanceFromWall = (r / BOX_Y);
	double wall_spatial_extent = 0.1;
	double closest_distance = 0.001;
	// TODO: change to something generic
	if (NormalizedDistanceFromWall < closest_distance){
		time_printf("Warning: particle very close to wall in iteration %d\n", iterNum);
		return pow(closest_distance / wall_spatial_extent, -4);
	}
	else
		return pow(NormalizedDistanceFromWall / wall_spatial_extent, -4);
}
void updateHistory(Vector **positionsHistory, double **velocitiesHistory, double *timeHistory, int bufferIterNum){
	// update history before moving anything
	memcpy(positionsHistory[bufferIterNum], positions, sizeof(Vector)*NPART);
	memcpy(velocitiesHistory[bufferIterNum], velocities, sizeof(double)*NPART);
	timeHistory[bufferIterNum] = elapsedTime;
}
void clearHistory(){
	/*
	int i, j;
	for (i = 0; i < iterations_in_memory; i++){
		for (j = 0; j < NPART; j++){
			positionsHistory[i][j].x = 0;
			positionsHistory[i][j].y = 0;
			velocitiesHistory[i][j] = 0;
		}
	}
	*/
}
void makeRandomState(Vector *positions, double *velocities){
	int i;
	for (i = 0; i < NPART; i++){
		// initalize position
		positions[i].x = nextRand() * BOX_X;
		//time_printf("debug: x = %f\n",positions[i].x);
		positions[i].y = nextRand() * BOX_Y;
		//positions[i].y = BOX_Y * (nextRand() * 0.8 + 0.1);
		// initialize velocity
		velocities[i] = nextRand() * 2 * PI;
	}
}

void makeRandomObstacles(Vector *obstaclePositions){
	int i;
	// initialize obstacles positions
	for (i = 0; i < OBSTACLENUM; i++){
		// initalize position
		obstaclePositions[i].x = nextRand() * BOX_X;
		obstaclePositions[i].y = nextRand() * BOX_Y;
	}
}

void initArrays(){
	int i;
	positions = (Vector*)malloc(NPART * sizeof(Vector));
	velocities = (double*)malloc(NPART * sizeof(double));
	prevVelocities = (double*)malloc(NPART * sizeof(double));
	obstaclePositions = (Vector*)malloc(OBSTACLENUM * sizeof(Vector));
	// allocate history arrays:
	timeHistoryLongTime = (double*)malloc(iterations_in_memory * sizeof(double));
	positionsHistoryLongTime = (Vector**)malloc(iterations_in_memory * sizeof(Vector*));
	velocitiesHistoryLongTime = (double**)malloc(iterations_in_memory * sizeof(double*));
	for (i = 0; i < iterations_in_memory; i++){
		positionsHistoryLongTime[i] = (Vector*)malloc(NPART * sizeof(Vector));
		velocitiesHistoryLongTime[i] = (double*)malloc(NPART * sizeof(double));
	}

}
