// horrible windows stuff:
#include "stdafx.h"

#include "profiler.h"
#include "utils.h"

#ifdef _USING_WINDOWS_
#include <windows.h>
#else
#include <time.h>
#endif
double PCFreq = 0.0;

timer_type newTimer(){
	timer_type timer;

#ifdef _USING_WINDOWS_
	LARGE_INTEGER freq;
	if (!QueryPerformanceFrequency(&freq))
		panic("QueryPerformanceFrequency failed!\n");
	timer.total_time.QuadPart = 0;
	PCFreq = (double)(freq.QuadPart) / 1000.0;
	QueryPerformanceCounter(&(timer.prev_time));
#else // linux
	struct timespec now;
	clock_gettime(CLOCK_MONOTONIC, &now);
	double dNow = now.tv_sec * 1000 + now.tv_nsec / 1000000.0;
	timer.prev_time = dNow;
	timer.total_time = 0;
#endif
	return timer;
}

void startTimer(timer_type *timer){
#ifdef _USING_WINDOWS_
	QueryPerformanceCounter(&(timer->prev_time));
#else
	struct timespec now;
	clock_gettime(CLOCK_MONOTONIC, &now);
	timer->prev_time = now.tv_sec*1000 + now.tv_nsec / 1000000.0;
#endif
}
void stopTimer(timer_type *timer){
#ifdef _USING_WINDOWS_
	LARGE_INTEGER now;
	QueryPerformanceCounter(&now);
	__int64 diff = now.QuadPart - timer->prev_time.QuadPart;
	timer->total_time.QuadPart = timer->total_time.QuadPart + diff;
#else
	struct timespec now;
	clock_gettime(CLOCK_MONOTONIC, &now);
	double dNow = now.tv_sec*1000 + now.tv_nsec / 1000000.0;
	double diff = dNow - timer->prev_time;
	timer->total_time = timer->total_time + diff;
#endif
}
double getTimerValue(timer_type timer){
#ifdef _USING_WINDOWS_
	return (double)(timer.total_time.QuadPart) / PCFreq;	// in miliseconds
#else
	return timer.total_time;	// in miliseconds
#endif
}

