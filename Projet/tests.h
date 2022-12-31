#ifndef TESTS_H_
#define TESTS_H_

#include "algo.h"

/* ---------- TEST BASE.C ---------- */

void test_add();

void test_sub();

void test_mul();

/* ***** MPFR ***** */

void test_mpadd();

void test_mpsub();

void test_mpmul();

/* ---------- TEST MATRIX.C ---------- */

void test_init_matrix(int n);

void test_init_eye(int n);

void test_copy_matrix(int n);

void test_matrix_add(int n);

void test_matrix_sub(int n);

void test_matrix_mul(int n);

void test_matrix_transpose(int n);

/* ***** MPFR ***** */

void test_init_MPFR_matrix(int n);

void test_init_MPFR_eye(int n);

void test_MPFR_copy_matrix(int n);

void test_MPFR_matrix_add(int n);

void test_MPFR_matrix_sub(int n);

void test_MPFR_matrix_mul(int n);

void test_MPFR_matrix_transpose(int n);


/* ---------- TEST ALGO.C ---------- */

void test_qr_decomposition(int n);

void test_quasi_hess(int n);

void test_hessenberg(int n);

void test_eigenvalues(int n);

/* ***** MPFR ***** */

void test_MPFR_qr_decomposition(int n);

void test_MPFR_quasi_hess(int n);

void test_MPFR_hessenberg(int n);

void test_MPFR_eigenvalues(int n);


#endif // TESTS_H