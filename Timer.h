#pragma once

#ifdef _WIN32
#include <windows.h>

class Timer
{
public:
	Timer();

	//in seconds
	float GetElapsedTime();
	void Reset();
private:
	LONGLONG cur_time;

	DWORD time_count;
	LONGLONG perf_cnt;
	bool perf_flag;
	LONGLONG last_time;
	double time_scale;

	bool QPC;
};
#else
//*****************************unix stuff****************************
#include <sys/time.h>


class Timer
{
public:
	Timer();

	void Reset();
	float GetElapsedTime();

private:
	timeval cur_time;

};

#endif //unix
