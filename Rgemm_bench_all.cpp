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
#include <chrono>

#include <gmpxx.h>
#include <qd/qd_real.h>
#include <quadmath.h>

#include "mpblas.hpp"
#include "typetype.hpp"

#define MFLOPS 1e-6
#define NANOSECOND 1e-9

// cf. https://netlib.org/lapack/lawnspdf/lawn41.pdf p.120
double flops_gemm(int64_t k_i, int64_t m_i, int64_t n_i) {
  double adds, muls, flops;
  double k, m, n;
  m = (double)m_i;
  n = (double)n_i;
  k = (double)k_i;
  muls = m * (k + 2) * n;
  adds = m * k * n;
  flops = muls + adds;
  return flops;
}

template <typename REAL>
void bench(int argc, char *argv[]) {
  REAL alpha, beta, dummy;
  double elapsedtime;

  char transa, transb, normtype;
  int64_t N0, M0, K0, STEPN = 3, STEPM = 3, STEPK = 3, LOOP = 3, TOTALSTEPS = 400;
  int64_t lda, ldb, ldc;
  int64_t i, m, n, k, ka, kb, p;
  int64_t check_flag;

  std::cout << "Test for " << TypeName<REAL>() << "\n";
  
  // initialization
  N0 = M0 = K0 = 1;
  transa = transb = 'n';
  normtype = 'm';
  if (argc != 1) {
    for (i = 1; i < argc; i++) {
      if (strcmp("-N", argv[i]) == 0) {
	N0 = atoi(argv[++i]);
      } else if (strcmp("-M", argv[i]) == 0) {
	M0 = atoi(argv[++i]);
      } else if (strcmp("-K", argv[i]) == 0) {
	K0 = atoi(argv[++i]);
      } else if (strcmp("-STEPN", argv[i]) == 0) {
	STEPN = atoi(argv[++i]);
      } else if (strcmp("-STEPM", argv[i]) == 0) {
	STEPM = atoi(argv[++i]);
      } else if (strcmp("-STEPK", argv[i]) == 0) {
	STEPK = atoi(argv[++i]);
      } else if (strcmp("-NN", argv[i]) == 0) {
	transa = transb = 'n';
      } else if (strcmp("-TT", argv[i]) == 0) {
	transa = transb = 't';
      } else if (strcmp("-NT", argv[i]) == 0) {
	transa = 'n';
	transb = 't';
      } else if (strcmp("-TN", argv[i]) == 0) {
	transa = 't';
	transb = 'n';
      } else if (strcmp("-NOCHECK", argv[i]) == 0) {
	check_flag = 0;
      } else if (strcmp("-LOOP", argv[i]) == 0) {
	LOOP = atoi(argv[++i]);
      } else if (strcmp("-TOTALSTEPS", argv[i]) == 0) {
	TOTALSTEPS = atoi(argv[++i]);
      }
    }
  }

  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::uniform_real_distribution<> urdist(-1.0, 1.0);

  m = M0;
  n = N0;
  k = K0;
  for (p = 0; p < TOTALSTEPS; p++) {
    if (Mlsame(&transa, "n")) {
      ka = k;
      lda = m;
    } else {
      ka = m;
      lda = k;
    }
    if (Mlsame(&transb, "n")) {
      kb = n;
      ldb = k;
    } else {
      kb = k;
      ldb = n;
    }
    ldc = m;

    REAL *a = new REAL [lda * ka];
    REAL *b = new REAL [ldb * kb];
    REAL *c = new REAL [ldc * n];
    REAL mone = -1;
    alpha = urdist(engine);
    beta = urdist(engine);
    for (i = 0; i < lda * ka; i++) {
      a[i] = urdist(engine);
    }
    for (i = 0; i < ldb * kb; i++) {
      b[i] = urdist(engine);
    }
    for (i = 0; i < ldc * n; i++) {
      c[i] = urdist(engine);
    }
    elapsedtime = 0.0;
    for (int j = 0; j < LOOP; j++) {
      std::chrono::steady_clock::time_point time_before;
      std::chrono::steady_clock::time_point time_after;

      time_before = std::chrono::steady_clock::now();
      mpblas::Rgemm<REAL>(&transa, &transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
      time_after = std::chrono::steady_clock::now();
      double time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(time_after - time_before).count();

      elapsedtime += time_in_ns;
    }
    elapsedtime = elapsedtime * NANOSECOND / (double)LOOP;
    printf("    m     n     k     MFLOPS    transa   transb\n");
    printf("%5d %5d %5d %10.3f         %c        %c\n", (int)m, (int)n, (int)k, flops_gemm(k, m, n) / elapsedtime * MFLOPS, transa, transb);
    delete[] c;
    delete[] b;
    delete[] a;
    m = m + STEPM;
    n = n + STEPN;
    k = k + STEPK;
  }
}

int main(int argc, char *argv[]) {
  bench<_Float16>(argc, argv);
  bench<float>(argc, argv);
  bench<double>(argc, argv);
  bench<dd_real>(argc, argv);
  bench<qd_real>(argc, argv);
  bench<_Float128>(argc, argv);
  bench<mpf_class>(argc, argv);
}
