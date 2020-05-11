// horrible windows stuff:
#include "stdafx.h"

#ifdef _USING_WINDOWS_
#include <string.h>
#else
#include <strings.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "log.h"
#include "parameters.h"
#include "profiler.h"


// external variables
extern double elapsedTime;
extern Vector **positionsHistory; 
extern double **velocitiesHistory;
extern double *dtHistory;
extern int iterations_in_memory;
extern int iteration_size;
extern int iterNum;


// profiling timers
timer_type openTimer;
timer_type closeTimer;
timer_type writeTimer;
timer_type bufTimer;
timer_type totalTimer;

void initLog(){
	timer_type openTimer = newTimer();
	timer_type closeTimer = newTimer();
	timer_type writeTimer = newTimer();
	timer_type bufTimer = newTimer();
	timer_type totalTimer = newTimer();
}
void killLog(){
}
void logObstacles(Vector *obstaclePositions, char *output_dir){
	FILE *obstacleFile;
	int i;
	char name[256] = "obstacles.dat";
	char full_name[256];
	double *buf;

	int counter;

	if (OUTPUT_FORMAT == FORMAT_CSV){
		printf("error! csv no longer supported");
		return;
	}

	buf = (double*)malloc(sizeof(double)*2*OBSTACLENUM);
	if (buf == NULL){
		printf("error! cannot allocate memory\n");
		exit(1);
	}
	if (((sizeof(double) * 2 * OBSTACLENUM) / 1000) > MAX_MEMORY){
		printf("Warning: obstacles are taking up a lot of space, writing may be slow\n");
	}
	counter = 0;
	for (i = 0; i < OBSTACLENUM;i++){
		// write two values: x, y
		buf[counter++] = obstaclePositions[i].x;
		buf[counter++] = obstaclePositions[i].y;
	}
	if (counter != (2 * OBSTACLENUM)){
		printf("output error: counter = %d, wanted it to be: %d\n!", counter, (2 * OBSTACLENUM));
		exit(1);
	}
	// construct file name
	strcpy(full_name, output_dir);
	strcat(full_name, separator);
	strcat(full_name, name);
	// open file
	obstacleFile = fopen(full_name, "wb");	// this "wb" is important, otherwise a bunch of random spaces are inserted and the files are unreadable
	if (obstacleFile == NULL){
		printf("Error opening output file. Probably output directory doesn\'t exist.\n");
		exit(1);
	}
	// write to file
	fwrite(buf, sizeof(double), counter, obstacleFile);
	// close file
	fclose(obstacleFile);
	// free buffer
	free(buf);
}

void logHistory(Vector **positionsHistory, double **velocitiesHistory, double *elapsed_times, int iterations_to_write, int filesWritten, char *output_dir){
	if (iterations_to_write <= 0){
		time_printf("Warning: trying to log empty file\n");
		return;
	}
	FILE *timeSliceFile;
	int i, j;
	char name[256];
	char full_name[256];
	double *buf;
	int counter;
	if (OUTPUT_FORMAT == FORMAT_CSV){
		printf("error! csv no longer supported");
		return;
	}
	buf = (double*)malloc(iterations_in_memory * iteration_size);
	if (buf == NULL){
		printf("error! cannot allocate memory\n");
		exit(1);
	}

	startTimer(&totalTimer);
	if (iterations_to_write != iterations_in_memory && iterNum != NITER){
		time_printf("Warning: is this the last iteration?\n");
	}
	counter = 0;
	for (j = 0; j < iterations_to_write; j++){
		startTimer(&bufTimer);
		buf[counter++] = elapsed_times[j];
		for (i = 0; i < NPART; i++){
			// write three values: x, y, theta
			buf[counter++] = positionsHistory[j][i].x;
			buf[counter++] = positionsHistory[j][i].y;
			buf[counter++] = velocitiesHistory[j][i];
		}
		stopTimer(&bufTimer);
		if (counter != (3 * NPART)*(j + 1) + (j + 1)){
			printf("output error: counter = %d, wanted it to be: %d\n!", counter, (3 * NPART)*(j + 1));
			exit(1);
		}
	}
	// construct file name
	replaceDollarWithNumber(name, LOG_FILE, filesWritten);
	strcpy(full_name, output_dir);
	strcat(full_name, separator);
	strcat(full_name, name);
	// open file
	startTimer(&openTimer);
	timeSliceFile = fopen(full_name, "wb");	// this "wb" is important, otherwise a bunch of random spaces are inserted and the files are unreadable
	stopTimer(&openTimer);
	if (timeSliceFile == NULL){
		printf("Error opening output file. Probably output directory doesn\'t exist.\n");
		exit(1);
	}
	// write to file
	startTimer(&writeTimer);
	fwrite(buf, sizeof(double), counter, timeSliceFile);
	stopTimer(&writeTimer);
	// close file
	startTimer(&closeTimer);
	fclose(timeSliceFile);
	stopTimer(&closeTimer);
	// free buffer
	free(buf);
	stopTimer(&totalTimer);
}

void logSummary(){
	int totaltime = getTimerValue(totalTimer);
	time_printf("Time filling buffers: %5.3f%%\n", getTimerValue(bufTimer)*100.0 / totaltime);
	time_printf("Time writing to files: %5.3f%%\n", getTimerValue(writeTimer)*100.0 / totaltime);
	time_printf("Time opening files: %5.3f%%\n", getTimerValue(openTimer)*100.0 / totaltime);
	time_printf("Time closing files: %5.3f%%\n", getTimerValue(closeTimer)*100.0 / totaltime);
}
int readLastFrame(Vector *positions, double *velocities, Vector *obstaclePositions, double *elapsedTime, char *filename){
	if (OUTPUT_FORMAT == FORMAT_CSV){
		time_printf("error! csv no longer supported");
		return -1;
	}
	char full_name[256];
	strcpy(full_name, OUTPUT_DIR);
	strcat(full_name, separator);
	strcat(full_name, "longtime");
	strcat(full_name, separator);
	strcat(full_name, filename);
	double *buf;
	FILE *lastFrameFile = fopen(full_name, "rb");
	int fsize = fseek(lastFrameFile, 0L, SEEK_END);
	fsize = ftell(lastFrameFile);
	if (fsize % iteration_size != 0){
		panic("Bad first frame");
	}
	if (fsize == 0){	// empty file
		fclose(lastFrameFile);
		return -1;
	}
	int iterations_in_file = fsize / iteration_size;
	fseek(lastFrameFile, iteration_size * (iterations_in_file - 1), SEEK_SET);
	buf = (double*)malloc(iterations_in_memory * iteration_size);
	if (buf == NULL){
		panic("error! cannot allocate memory\n");
	}
	fread(buf, 1, iteration_size, lastFrameFile);
	fclose(lastFrameFile);
	int i=0, j=0,counter = 1;
	*elapsedTime = buf[0];
	while(counter < (NPART*3 +1)){
		positions[i].x = buf[counter++];
		positions[i].y = buf[counter++];
		velocities[i] = buf[counter++];
		if (positions[i].x < 0 || positions[i].x > BOX_X ||
			positions[i].y < 0 || positions[i].y > BOX_Y ||
			velocities[i] < 0 || velocities[i] > 2 * PI){
			panic("Bad frame");
		}
		i++;
		if (i > NPART){
			panic("Bad frame bug");
		}
	}
	free(buf);
	// initialize obstacles positions
	strcpy(full_name, OUTPUT_DIR);
#ifdef _USING_WINDOWS_
	strcat(full_name, "\\longtime\\obstacles.dat");
#else
	strcat(full_name, "/longtime/obstacles.dat");
#endif
	FILE *obstacleFile = fopen(full_name, "rb");
	if (obstacleFile == NULL){
		panic("cannot find obstacle file");
	}
	buf = (double*)malloc(sizeof(double) * 2 * OBSTACLENUM);
	fread(buf, sizeof(double), 2 * OBSTACLENUM, obstacleFile);
	fclose(obstacleFile);
	counter = 0;
	while (counter < 2 * OBSTACLENUM){
		// initalize position
		obstaclePositions[i].x = buf[counter++];
		obstaclePositions[i].y = buf[counter++];
		if (obstaclePositions[i].x < 0 || obstaclePositions[i].x > BOX_X ||
			obstaclePositions[i].y < 0 || obstaclePositions[i].y > BOX_Y){
			panic("Bad obstacle file");
		}
	}

	free(buf);
	return 1;
}
char *getExtenstion(char *filename){
	char *p = filename;
	while (*p && *p != '.') p++;
	p++;
	return p;
}
char *getNextFloat (char *line, double *ans){
	if (!line) return NULL;
	char *p = line;
	char *start;
	char num[50];
	while (*p && !isdigit(*p)) p++;
	if (!*p) return NULL;
	start = p;
	while (*p && isdigit(*p)) p++;
	if (*p != '.'){
		memcpy(num, start, sizeof(char) * (p - start));
		*ans = atof(num);
		return p;
	}
	p++;
	while (*p && isdigit(*p)) p++;
	memcpy(num, start, sizeof(char) * (p - start));
	*ans = atof(num);
	return p;
}
int readInitConfFile(char *filename, Vector *positions, double *velocities, Vector *obstaclePositions){
	char *ext = getExtenstion(filename);
	if (STRCMPI(ext, "csv") == 0 || STRCMPI(ext, "txt") == 0){
		FILE *f = fopen(filename, "r");
		char line[256];
		char *p;
		int counter = 0;
		double x, y, theta;
		if (f == NULL){
			time_printf("Error: cannot find initial configuration file\n");
			return -1;
		}
		time_printf("Reading parameters...\n");
		while (fgets(line, sizeof(line), f)) {
			if (line[0] == '/')	// skip comments
				continue;
			p = line;
			//sscanf(line, "%f%*[ \n\t],%*[ \n\t]%f%*[ \n\t],%*[ \n\t]%f%*[ \n\t]", &x, &y, &theta);
			//sscanf(line, "%f,%f,%f", &x, &y, &theta);
			p = getNextFloat(p, &x);
			p = getNextFloat(p, &y);
			p = getNextFloat(p, &theta);
			if (!p) continue;
			theta = putInUnitCircle(theta * PI / 180.0);
			//time_printf("parsed values: %f %f %f\n", x, y, theta);
			positions[counter].x = x;
			positions[counter].y = y;
			velocities[counter] = theta;
			counter++;
		}
		// particle number doesn't match input
		if (counter != NPART)
			return -1;
		return 1;
	}
	else { // binary file
		FILE *lastFrameFile = fopen(filename, "rb");
		if (lastFrameFile == NULL){
			time_printf("First frame file does not exist: %s\n", filename);
			panic("Bad first frame");
		}
		int fsize = fseek(lastFrameFile, 0L, SEEK_END);
		fsize = ftell(lastFrameFile);
		if (fsize % iteration_size != 0){
			panic("Bad first frame");
		}
		if (fsize == 0){	// empty file
			fclose(lastFrameFile);
			return -1;
		}
		int iterations_in_file = fsize / iteration_size;
		fseek(lastFrameFile, iteration_size * (iterations_in_file - 1), SEEK_SET);
		double *buf = (double*)malloc(iterations_in_memory * iteration_size);
		if (buf == NULL){
			panic("error! cannot allocate memory\n");
		}
		fread(buf, 1, iteration_size, lastFrameFile);
		fclose(lastFrameFile);
		int i = 0, j = 0, counter = 1;
		//*elapsedTime = buf[0];
		while (counter < (NPART * 3 + 1)){
			positions[i].x = buf[counter++];
			positions[i].y = buf[counter++];
			velocities[i] = buf[counter++];
			if (positions[i].x < 0 || positions[i].x > BOX_X ||
				positions[i].y < 0 || positions[i].y > BOX_Y ||
				velocities[i] < 0 || velocities[i] > 2 * PI){
				panic("Bad frame");
			}
			i++;
			if (i > NPART){
				panic("Bad frame bug");
			}
		}
		free(buf);
		return 1;
	}
	
}
