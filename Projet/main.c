#include "algo.h"
#include "unit_test.h"
#include "benchmark.h"

int main(){
    int n = 6;
    srand(time(NULL));
    test_base();
    test_matrix(n);
    test_algo(n);

    // benchmark_quasi_hess_vs_hessenberg(10, 2);
    benchmark_MPFR_quasi_hess_vs_hessenberg(10, 2);

   
    return 0;
}