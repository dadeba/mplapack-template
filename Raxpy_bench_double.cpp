/*
 * Copyright (c) 2008-2022
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */
#include <cstdint>
#include <cstring>
#include <random>
#include <iostream>

#include <time.h>
#include "mpblas.hpp"

#define MFLOPS 1e-6
#define NANOSECOND 1e-9

int main(int argc, char *argv[]) {
    int64_t n;
    int64_t incx = 1, incy = 1, STEP = 97, N0 = 1, LOOP = 3, TOTALSTEPS = 3000;
    double alpha;
    double elapsedtime;
    long elapsedtime_l;
    long t1, t2;
    int i, p;
    int check_flag;
    struct timespec ts;

    if (argc != 1) {
        for (i = 1; i < argc; i++) {
            if (strcmp("-N", argv[i]) == 0) {
                N0 = atoi(argv[++i]);
            } else if (strcmp("-STEP", argv[i]) == 0) {
                STEP = atoi(argv[++i]);
            } else if (strcmp("-NOCHECK", argv[i]) == 0) {
                check_flag = 0;
            } else if (strcmp("-LOOP", argv[i]) == 0) {
                LOOP = atoi(argv[++i]);
            } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
                TOTALSTEPS = atoi(argv[++i]);
            }
        }
    }

    n = N0;
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::uniform_real_distribution<> urdist(-1.0, 1.0);

    for (p = 0; p < TOTALSTEPS; p++) {
        double *x = new double[n];
        double *y = new double[n];
        for (i = 0; i < n; i++) {
            x[i] = urdist(engine);
            y[i] = urdist(engine);
        }
        alpha = urdist(engine);
        elapsedtime_l = 0;
        for (int j = 0; j < LOOP; j++) {
            clock_gettime(CLOCK_REALTIME, &ts);
            t1 = ts.tv_nsec;
            mpblas::Raxpy<double>(n, alpha, x, incx, y, incy);
            clock_gettime(CLOCK_REALTIME, &ts);
            t2 = ts.tv_nsec;
            elapsedtime_l = elapsedtime_l + t2 - t1;
        }
        elapsedtime = (double)elapsedtime_l * NANOSECOND / (double)LOOP;
        printf("         n       MFLOPS\n");
        printf("%10d   %10.3f\n", (int)n, (2.0 * (double)n) / elapsedtime * MFLOPS);
        delete[] y;
        delete[] x;
        n = n + STEP;
    }
}
