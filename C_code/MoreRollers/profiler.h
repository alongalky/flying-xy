#ifndef __PROFILER_H
#define __PROFILER_H

#include "utils.h"

#ifdef _USING_WINDOWS_	// windows section

#include <windows.h>

typedef struct _timer_type{
	LARGE_INTEGER prev_time;
	LARGE_INTEGER total_time;
} timer_type;

#else // using linux
#include <time.h>
typedef struct _timer_type{
	double prev_time;
	double total_time;
} timer_type;

#endif

timer_type newTimer();
void startTimer(timer_type *timer);
void stopTimer(timer_type *timer);
double getTimerValue(timer_type timer);


#endif