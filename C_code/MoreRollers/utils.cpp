// horrible windows stuff:
#include "stdafx.h"
#ifdef _USING_WINDOWS_
#include <io.h>
#include <direct.h>
#include <string.h>
#else
#include <strings.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "random.h"
#include "utils.h"
#include "parameters.h"
#include <ctype.h>
#include <stdarg.h>
#include <time.h>

extern int resume_previous_simulation;
extern char prev_simulation_lastfile[256];
extern double NOISE;	// original noise, for hysterisis

// Hysterisis!
double noiseFunction(double time){
	// triangle function
	if (HYSTER_NOISE < 0 || time > SIMTIME){
		time_printf("SIMTIME = %f, HYSTER_NOISE = %f, time = %f\n", SIMTIME, HYSTER_NOISE, time);
		panic("Bad hysterisis setup");
	}
	double half_sim_time = SIMTIME / 2;
	if (time < half_sim_time){
		return NOISE + (time / half_sim_time)*(HYSTER_NOISE - NOISE);
	}
	else{
		return HYSTER_NOISE + ((time-half_sim_time) / half_sim_time)*(NOISE - HYSTER_NOISE);
	}
}
void panic(char *msg){
	time_printf(msg);
	exit(1);
}
double XMIRROR(double angle){
	double ans;
	//if (angle < PI / 4)
		ans = PI - angle;
	while (ans < 0) ans += 2 * PI;
	return ans;
}
double YMIRROR(double angle){
	double ans;
	ans = 2*PI - angle;
	return ans;
}
double putInUnitCircle(double angle){
	double ret = angle;
	while (ret < 0) ret += 2 * PI;
	while (ret > 2 * PI) ret -= 2 * PI;
	return ret;
}
Vector vecAdd(Vector a, Vector b){
	Vector res;
	res.x = a.x + b.x;
	res.y = a.y + b.y;
	return res;
}
Vector vecSub(Vector a, Vector b){
	Vector res;
	res.x = b.x - a.x;
	res.y = b.y - a.y;
	return res;
}
double vecDot(Vector a, Vector b){
	return (a.x*b.x) + (a.y*b.y);
}
double vecSq(Vector a){
	return vecDot(a, a);
}
double distSq(Vector position_a, Vector position_b){
	// calculate the distance between two particles. Periodic boundary condition in X
	double dist_x, dist_y, delta;
	if (PERIODIC_X){
		delta = abs(position_a.x - position_b.x);
		dist_x = MIN(delta, BOX_X - delta);
	}
	else{
		dist_x = position_a.x - position_b.x;
	}
	if (PERIODIC_Y){
		delta = abs(position_a.y - position_b.y);
		dist_y = MIN(delta, BOX_Y - delta);
	}
	else{
		dist_y = position_a.y - position_b.y;
	}
	return SQR(dist_x) + SQR(dist_y);
}
double vecDirection(Vector a){
	return putInUnitCircle(atan2(a.y, a.x));
}

// internal parsing functions
void cleanWhitespace(char *str){
	char *start = str;
	int i = 0;
	while (*start == ' ' || *start == '\n' || *start == '\t') start++;
	while (*start){
		str[i] = *start;
		i++; start++;
	}
	str[i] = 0;
	int len = strlen(str) - 1;
	while (str[len] == ' ' || str[len] == '\n' || str[len] == '\t'){
		str[len] = NULL;
		len--;
	}
}
void cleanComment(char *str){
	char *p = str;
	while (*p){
		if (*p == '/' && *(p + 1) == '/'){
			*p = NULL;
		}
		else p++;
	}
}
void getFieldName(char *output, char *line){
	char *pos;
	pos = strchr(line, '=');
	*output = NULL;
	if (pos != NULL){
		memcpy(output, line, pos - line);
		output[pos - line] = NULL;
		cleanWhitespace(output);
	}
}
void getFieldValue(char *output, char *line){
	char *pos;
	*output = NULL;
	pos = strchr(line, '=');
	if (pos != NULL){
		pos++;
		memcpy(output, pos, strlen(pos));
		output[strlen(pos) - 1] = NULL;
		cleanComment(output);
		cleanWhitespace(output);
	}
}
void readParameterFile(char *str){
	FILE *fInput = fopen(str, "r");
	char line[256];
	char field_name[256];
	char field_val[256];
	double val;
	int val_int;
	int padding = 15;

	if (fInput == NULL)	{
		time_printf("Error opening input file! Cannot find %s\n", str);
		exit(1);
	}

	time_printf("Reading parameters...\n");
	while (fgets(line, sizeof(line), fInput)) {
		if (line[0] == '/')	// skip comments
			continue;

		getFieldName(field_name, line);
		getFieldValue(field_val, line);
		//time_printf("field = %s ; value = %s\n", field_name, field_val);
		if (STRCMPI(field_name, "LOG_FILE") == 0){
			time_printf("log file is %s\n", field_val);
			strcpy(LOG_FILE, field_val);
		}
		else if (STRCMPI(field_name, "OUTPUT_DIR") == 0){
			time_printf("%-30s %s\n", "OUTPUT_DIR file is ", field_val);
			strcpy(OUTPUT_DIR, field_val);
		}
		else if (STRCMPI(field_name, "NPART") == 0){
			time_printf("%-30s %d\n","NPART is ", atoi(field_val));
			val_int = atoi(field_val);
			NPART = val_int;
		}
		else if (STRCMPI(field_name, "NITER") == 0){
			time_printf("%-30s %d\n", "NITER is ", atoi(field_val));
			val_int = atoi(field_val);
			NITER = val_int;
		}
		else if (STRCMPI(field_name, "NOISE") == 0){
			val = atof(field_val);
			time_printf("%-30s %5.3f\n", "NOISE is ", val);
			NOISE = val;
		}
		else if (STRCMPI(field_name, "BOX_X") == 0){
			val = atof(field_val);
			time_printf("%-30s %5.3f\n", "BOX_X is ", val);
			BOX_X = val;
		}
		else if (STRCMPI(field_name, "BOX_Y") == 0){
			val = atof(field_val);
			time_printf("%-30s %5.3f\n", "BOX_Y is ", val);
			BOX_Y = val;
		}
		else if (STRCMPI(field_name, "V0") == 0){
			val = atof(field_val);
			time_printf("%-30s %5.3f\n", "V0 is ", val);
			V0 = val;
		}
		else if (STRCMPI(field_name, "DT") == 0){
			val = atof(field_val);
			time_printf("%-30s %5.3f\n", "DT is ", val);
			originalDT = val;
		}
		else if (STRCMPI(field_name, "ALIGN") == 0){
			val = atof(field_val);
			time_printf("%-30s %5.3f\n", "ALIGN is ", val);
			ALIGN = val;
		}
		else if (STRCMPI(field_name, "DEFLECT") == 0){
			val = atof(field_val);
			if (val != 0){
				panic("Deflection isn\'t tested yet, please use zero.");
			}
			time_printf("%-30s %5.3f\n", "DEFLECT is ", val);
			DEFLECT = val;
		}
		else if (STRCMPI(field_name, "RSTAR") == 0){
			val = atof(field_val);
			time_printf("%-30s %5.3f\n", "RSTAR is ", val);
			RSTAR = val;
			RSTAR_SQ = SQR(val);
		}
		else if (STRCMPI(field_name, "OUTPUT_FORMAT") == 0){
			if (STRCMPI(field_val, formatNames[FORMAT_CSV]) == 0)
				val_int = FORMAT_CSV;
			else if (STRCMPI(field_val, formatNames[FORMAT_BIN]) == 0)
				val_int = FORMAT_BIN;
			else{
				time_printf("OUTPUT_FORMAT is invalid, defaulting to bin\n");
				val_int = FORMAT_BIN;
			}
			time_printf("%-30s %s\n", "OUTPUT_FORMAT file is ", formatNames[val_int]);
			OUTPUT_FORMAT = val_int;
		}
		// AUTOCODE START
		else if (STRCMPI(field_name, "THREADS") == 0){
			time_printf("%-30s %d\n","THREADS is ", atoi(field_val));
		    val_int = atoi(field_val);
		    THREADS = val_int;
		}
		else if (STRCMPI(field_name, "MAX_RUNTIME_MINUTES") == 0){
			time_printf("%-30s %d\n", "MAX_RUNTIME_MINUTES is ", atoi(field_val));
			val_int = atoi(field_val);
		    MAX_RUNTIME_MINUTES = val_int;
		}
		else if (STRCMPI(field_name, "HYSTER_NOISE") == 0){
		    val = atof(field_val);
			time_printf("%-30s %5.3f\n", "HYSTER_NOISE is ", val);
		    HYSTER_NOISE = val;
		}
		else if (STRCMPI(field_name, "HYSTER") == 0){
		    if (STRCMPI(field_val, "TRUE") == 0)
		        val = 1;
		    else
		        val = 0;
			time_printf("%-30s %s\n", "HYSTER is ", val ? "true" : "false");
		    HYSTER = val;
		}
		else if (STRCMPI(field_name, "INITIAL_DATA_FILE") == 0){
			time_printf("%-30s %s\n", "INITIAL_DATA_FILE file is ", field_val);
		    strcpy(INITIAL_DATA_FILE, field_val);
		}
		else if (STRCMPI(field_name, "SIMTIME") == 0){
		    val = atof(field_val);
			time_printf("%-30s %5.3f\n", "SIMTIME is ", val);
		    SIMTIME = val;
		}
		else if (STRCMPI(field_name, "TIME_BETWEEN_SAMPLES") == 0){
		    val = atof(field_val);
			time_printf("%-30s %5.3f\n", "TIME_BETWEEN_SAMPLES is ", val);
		    TIME_BETWEEN_SAMPLES = val;
		}
		else if (STRCMPI(field_name, "NOISE_MODE") == 0){
			if (STRCMPI(field_val, "VECTORIAL") == 0){
				panic("Are you sure you want vectorial noise? it\'s really confusing");
				val = VECTORIAL;
			}
			else if (STRCMPI(field_val, "ANGULAR") == 0){
				val = ANGULAR;
			}
			else
				continue;
			time_printf("%-30s %s\n", "NOISE_MODE is ", field_val);
		    NOISE_MODE = val;
		}
		else if (STRCMPI(field_name, "NORMALIZE_NEIGHBORS") == 0){
		    if (STRCMPI(field_val, "TRUE") == 0)
		        val = 1;
		    else
		        val = 0;
			time_printf("%-30s %s     (%s Model)\n", "NORMALIZE_NEIGHBORS is ", val ? "true" : "false", val ? "Vicsek" : "Flying XY");
		    NORMALIZE_NEIGHBORS = val;
		}
		else if (STRCMPI(field_name, "PERIODIC_Y") == 0){
		    if (STRCMPI(field_val, "TRUE") == 0)
		        val = 1;
		    else
		        val = 0;
			time_printf("%-30s %s\n", "PERIODIC_Y is ", val ? "true" : "false");
		    PERIODIC_Y = val;
		}
		else if (STRCMPI(field_name, "PERIODIC_X") == 0){
		    if (STRCMPI(field_val, "TRUE") == 0)
		        val = 1;
		    else
		        val = 0;
			time_printf("%-30s %s\n", "PERIODIC_X is ", val ? "true" : "false");
		    PERIODIC_X = val;
		}
		else if (STRCMPI(field_name, "OBSTACLE_RADIUS") == 0){
		    val = atof(field_val);
		    time_printf("OBSTACLE_RADIUS is %f\n", val);
		    OBSTACLE_RADIUS = val;
		}
		else if (STRCMPI(field_name, "OBSTACLEDEFLECTION") == 0){
		    val = atof(field_val);
			time_printf("%-30s %5.3f\n", "OBSTACLEDEFLECTION is ", val);
		    OBSTACLEDEFLECTION = val;
		}
		else if (STRCMPI(field_name, "OBSTACLENUM") == 0){
			time_printf("%-30s %d\n","OBSTACLENUM is ", atoi(field_val));
		    val_int = atoi(field_val);
		    OBSTACLENUM = val_int;
		}
		else if (STRCMPI(field_name, "SOFT_WALLS_STRENGTH") == 0){
		    val = atof(field_val);
			time_printf("%-30s %5.3f\n", "SOFT_WALLS_STRENGTH is ", val);
			SOFT_WALLS_STRENGTH = val;
		}
		else if (STRCMPI(field_name, "SOFT_WALLS") == 0){
		    if (STRCMPI(field_val, "TRUE") == 0)
		        val = 1;
		    else
		        val = 0;
			time_printf("%-30s %s\n", "SOFT_WALLS is ", val ? "true" : "false");
		    SOFT_WALLS = val;
		}
		else if (STRCMPI(field_name, "RDEFLECTION") == 0){
		    val = atof(field_val);
			time_printf("%-30s %5.3f\n", "RDEFLECTION is ", val);
		    RDEFLECTION = val;
		}
		// AUTOCODE END
	}

	fclose(fInput);
}
// replaces the first $ in the string with num. if no dollar is found, does nothing
void replaceDollarWithNumber(char *dest, char *source, int num){
	char *p = dest;
	char num_str[256];
	while (*(source) != '$' && source){
		*(dest++) = *(source++);
	}
	if (!source){	// no dollar sign
		time_printf("warning: no dollar found\n");
		*dest = 0;
		return;
	}
	sprintf(num_str, "%09d", num);
	strcpy(dest, num_str);
	strcat(dest, source+1);
}
int extract_int(char *str){
	char *p = str;
	char num[256];
	int i = 0;
	while (!isdigit(*p)) p++;
	while (isdigit(*p)){
		num[i++] = *p++;
	}
	num[i] = '\0';
	return atoi(num);
}
int endwith(char *str, char *suffix){
	int len = strlen(str);
	int suf_len = strlen(suffix);
	if (len < suf_len){
		return 0;
	}
	return (strcmp(str+len-suf_len, suffix) == 0);
}
int startswith(char *str, char *prefix){
	int len = strlen(str);
	int pref_len = strlen(prefix);
	int i;
	if (len < pref_len){
		return 0;
	}
	for (i = 0; i < pref_len; i++){
		if (str[i] != prefix[i])
			return 0;
	}
	return 1;
}
int find_last_logfile(char *dir, char *output){
#ifdef _USING_WINDOWS_
	struct _finddata_t file_info;
	char path[256];
	strcpy(path, dir);
	strcat(path, separator);
	strcat(path, "log*.dat");

	char last_file[256];
	int last_file_num = -1;
	int file_num;

	intptr_t handle = 0;
	memset(&file_info, 0, sizeof(file_info));
	handle = _findfirst(path, &file_info);

	if (handle != -1)
	{
		do
		{
			file_num = extract_int(file_info.name);
			if (file_num > last_file_num){
				last_file_num = file_num;
				strcpy(last_file, file_info.name);
			}
		} while (_findnext(handle, &file_info) != -1);
	}
	else{
		file_num = -1;
	}
	_findclose(handle);
	if (output != NULL && file_num != -1)
		strcpy(output, last_file);
	return last_file_num;
#else //linux

#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
	int last_file_num = -1;
	int file_num;
	char last_file[256];
	int found_any = 0;

	struct dirent *dp;
	DIR *dfd;

	if ((dfd = opendir(dir)) == NULL)
	{
		fprintf(stderr, "Can't open %s\n", dir);
		return 0;
	}

	char filename_qfd[256];
	char new_name_qfd[256];
	while ((dp = readdir(dfd)) != NULL)
	{
		struct stat stbuf;
		sprintf(filename_qfd, "%s/%s", dir, dp->d_name);
		//printf("debug: %s\n", filename_qfd);
		if (stat(filename_qfd, &stbuf) == -1)
		{
			printf("Unable to stat file: %s\n", filename_qfd);
			continue;
		}

		if ((stbuf.st_mode & S_IFMT) == S_IFDIR)
		{
			continue;
			// Skip directories
		}
		else
		{
			if (endwith(dp->d_name, ".dat") && startswith(dp->d_name, "log")){
				time_printf("debug: found logfile: %s\n",dp->d_name);
				file_num = extract_int(dp->d_name);
				if (file_num > last_file_num){
					last_file_num = file_num;
					strcpy(last_file, dp->d_name);
				}
				found_any = 1;
			}
		}
	}
	if (found_any == 0)
		return -1;
	if (output != NULL && file_num != -1)
		strcpy(output, last_file);
	return last_file_num;
#endif
}
void parse_commandLineArgs(int argc, char* argv[]){
	int i;
	for (i = 1; i < argc; i++){
		if (argv[i][0] == '-'){
			if (STRCMPI(argv[i], "-r") == 0){
				resume_previous_simulation = TRUE;
				time_printf("Trying to resume previous simulation\n");
			}
			else if (STRCMPI(argv[i], "--resume") == 0){
				resume_previous_simulation = TRUE;
				time_printf("Trying to resume previous simulation\n");
			}
		}
		else {
			strcpy(INPUT_FILE, (char*)argv[1]);
			time_printf("Input file: %s\n", INPUT_FILE);
		}
	}
}

void time_printf(const char *fmt, ...){
	time_t timer;
	char buffer[255];
	struct tm* tm_info;
	time(&timer);
	tm_info = localtime(&timer);
	strftime(buffer, 25, "%Y/%m/%d %H:%M:%S> ", tm_info);
	// finished making time stamp
	strcat(buffer, fmt);
	va_list args;
	va_start(args, fmt);
	vprintf(buffer, args);
	va_end(args);
	fflush(stdout);
}
// returns time
clock_t getTime(){
#ifdef _USING_WINDOWS_
	return clock() / (CLOCKS_PER_SEC / 1000);
#else
	struct timespec ts;
	unsigned theTick = 0U;
	clock_gettime(CLOCK_REALTIME, &ts);
	theTick = ts.tv_nsec / 1000000;
	theTick += ts.tv_sec * 1000;
	return theTick;
#endif
}
