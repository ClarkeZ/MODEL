#include "algo.h"
#include "unit_test.h"

int main(){
    srand(time(NULL));
    test_base();
    test_matrix(3);
    test_algo(3);

    // mpfr n;
    // mpfr_init2(n, 200);
    // mpfr_set_d(n, 2.0, MPFR_RNDD);
    // Matrix *mat = init_matrix(n);
    // mpfr_out_str (stdout, 10, 0, mat->n, MPFR_RNDD);
    // mpfr_out_str (stdout, 10, 0, mat->m[0], MPFR_RNDD);
    // mpfr_out_str (stdout, 10, 0, mat->m[1], MPFR_RNDD);
    // mpfr_out_str (stdout, 10, 0, mat->m[2], MPFR_RNDD);
    // mpfr_out_str (stdout, 10, 0, mat->m[3], MPFR_RNDD);
    // free_matrix(mat);

    // mpfr_clear(n);
    return 0;
}