#ifndef __UTILS_H
#define __UTILS_H

#include <time.h>

#ifdef _USING_WINDOWS_
#include <string.h>
#define STRCMPI(a,b) _strcmpi(a, b)
#else
#include <strings.h>
#define STRCMPI(a,b) strcasecmp(a, b)
#endif



// macros
#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

// types
typedef struct tagVector {
	double x, y;
} Vector;
// constants
#define PI (3.14159265358979323846)

// Simulation functions
double noiseFunction(double time);		// In the case of a hysteresis calculation

// Error handling
void panic(char *msg);

// angle functions
double XMIRROR(double angle);
double YMIRROR(double angle);
double putInUnitCircle(double angle1);

// vector functions
double distSq(Vector position_a, Vector position_b);
Vector vecAdd(Vector a, Vector b);
Vector vecSub(Vector a, Vector b);
double vecDot(Vector a, Vector b);
double vecSq(Vector a);
double vecDirection(Vector a);

// parsing functions
void readParameterFile(char *str);
void replaceDollarWithNumber(char *dest, char *source, int num);
void parse_commandLineArgs(int argc, char* argv[]);
int find_last_logfile(char *dir, char *output);

// time functions
void time_printf(const char *fmt, ...);
clock_t getTime();



#endif