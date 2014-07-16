/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file Timer.h
 * @brief Wall clock timing.
 */

#ifndef __UTIL_OSTIMER_H__
#define __UTIL_OSTIMER_H__

// Created by bdl 5. april 2002
// The purpose of this file is to make a timer function that is as 
// precise as posible on any given platform
#if (_MSC_VER >= 1200)
#pragma warning (disable: 4244)
#endif
#ifdef WIN32
#include <windows.h>
static LARGE_INTEGER	largeInteger;
#else
#include <sys/time.h>
#endif

namespace Util{
	class Timer {
#ifdef WIN32
		double freq;
		double start_count;
#else
		timeval start_time;
#endif
	public:
		void start() {
#ifdef WIN32
			QueryPerformanceFrequency(&largeInteger);
			freq = largeInteger.QuadPart;
			QueryPerformanceCounter(&largeInteger);
			start_count = largeInteger.QuadPart;
#else
			gettimeofday(&start_time,0);
#endif
		}

		float get_secs() {
#ifdef WIN32
			QueryPerformanceCounter(&largeInteger);
			double now_count = largeInteger.QuadPart;
			return (float)((now_count-start_count)/freq);
#else
			timeval now;
			gettimeofday(&now,0);
			return (now.tv_sec-start_time.tv_sec) + 
				(now.tv_usec-start_time.tv_usec)/1.0e6;
#endif
		}
	};
}

#endif
