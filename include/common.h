#ifndef _COMMON_H_
#define _COMMON_H_

#include <iostream>
#include <fstream>
#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include<vector>
#include<algorithm>
#include<queue>
#include <malloc.h>
#include <chrono>

typedef chrono::high_resolution_clock Clock;



#define TDEF(x_) chrono::high_resolution_clock::time_point x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = Clock::now();
#define TEND(x_) x_##_t1 = Clock::now();
#define TPRINT(x_, str) printf("%-20s \t%.6f\t sec\n", str, chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()/1e6);
