#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include "algo.h"

void benchmark_quasi_hess_vs_hessenberg(unsigned int N, unsigned int ITE);

void benchmark_MPFR_quasi_hess_vs_hessenberg(unsigned int N, unsigned int ITE);

void benchmark_qr_decomposition(unsigned int N, unsigned int ITE);

void benchmark_quasi_hess(unsigned int N, unsigned int ITE);

void benchmark_hessenberg(unsigned int N, unsigned int ITE);

#endif /* BENCHMARK_H_ */