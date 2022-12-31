#include "unit_test.h"

/* ---------- TEST BASE.C ---------- */
void test_base(){

    // test_add();

    // test_sub();

    // test_mul();

    /* ***** MPFR ***** */

    // test_mpadd();

    // test_mpsub();

    // test_mpmul();

}

/* ---------- TEST MATRIX.C ---------- */
void test_matrix(int n){
    // Pour éviter les warnings de variables non utilisées
    n += 1;
    n -= 1;
    
    // test_init_matrix(n);

    // test_init_eye(n);

    // test_copy_matrix(n);

    // test_matrix_add(n);

    // test_matrix_mul(n);

    // test_matrix_sub(n);

    // test_matrix_transpose(n);

    /* ***** MPFR ***** */

    // test_init_MPFR_matrix(n);

    // test_init_MPFR_eye(n);

    // test_MPFR_copy_matrix(n);

    // test_MPFR_matrix_add(n);

    // test_MPFR_matrix_sub(n);

    // test_MPFR_matrix_mul(n);

    // test_MPFR_matrix_transpose(n);

}

 /* ---------- TEST ALGO.C ---------- */
void test_algo(int n){
    // Pour éviter les warnings de variables non utilisées
    n += 1;
    n -= 1;

    // test_qr_decomposition(n);

    // test_quasi_hess(n);

    // test_hessenberg(n);

    // test_eigenvalues(n);

    /* ***** MPFR ***** */

    // test_MPFR_qr_decomposition(n);

    // test_MPFR_quasi_hess(n);

    // test_MPFR_hessenberg(n);

    // test_MPFR_eigenvalues(n);
}


