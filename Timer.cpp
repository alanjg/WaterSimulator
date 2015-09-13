#include "stdafx.h"
#include "Timer.h"

#ifdef _WIN32
#pragma comment(lib, "winmm.lib")
Timer::Timer()
{
	last_time = 0;
	if (QueryPerformanceFrequency((LARGE_INTEGER*)&perf_cnt))
	{
		perf_flag = true;
		time_count = DWORD(perf_cnt); //perf_cnt counts per second
		QueryPerformanceCounter((LARGE_INTEGER*)&last_time);
		time_scale = 1.0 / perf_cnt;
		QPC = true;
	}
	else
	{
		perf_flag = false;
		last_time = timeGetTime();
		time_scale = 0.001;
		time_count = 33;
	}
	Reset();
}

float Timer::GetElapsedTime()
{
	if (perf_flag)
		QueryPerformanceCounter((LARGE_INTEGER*)&cur_time);
	else
		cur_time = timeGetTime();

	float time_elapsed = float((cur_time - last_time)*time_scale);
	//last_time=cur_time;
	return time_elapsed;
}

void Timer::Reset()
{
	GetElapsedTime();
	last_time = cur_time;
}


#else
//***********************************unix specific*********************************

Timer::Timer()
{
	gettimeofday(&cur_time, 0);
}

void Timer::Reset()
{
	gettimeofday(&cur_time, 0);
}

float Timer::GetElapsedTime()
{
	float dif;
	timeval newtime;
	gettimeofday(&newtime, 0);
	dif = (newtime.tv_sec - cur_time.tv_sec);
	dif += (newtime.tv_usec - cur_time.tv_usec) / 1000000.0;
	return dif;
}

#endif
