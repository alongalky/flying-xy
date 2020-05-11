// These functions will handle reading and writing data files, to store the simulation data at the end of the run.
// The files are binary, comprised only of "double" numbers, and the format is: <time> <x_1> <y_1> <theta_1>  ... <x_N> <y_N> <theta_N>

#ifndef _LOG_H
#define _LOG_H


// logging functions
void logHistory(Vector **positionsHistory, double **velocitiesHistory, double *elapsed_times, int iterations_to_write, int filesWritten, char *output_dir);
void killLog();
void initLog();
void logSummary();
void logObstacles(Vector *obstaclePositions, char *output_dir);
int readLastFrame(Vector *positions, double *velocities, Vector *obstaclePositions, double *elapsedTime, char *filename);
int readInitConfFile(char *filename, Vector *positions, double *velocities, Vector *obstaclePositions);

#endif