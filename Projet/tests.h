#ifndef TESTS_H_
#define TESTS_H_

#include "algo.h"

/* ---------- TEST BASE.C ---------- */

void test_add();

void test_sub();

void test_mul();

void test_mpadd();

void test_mpsub();

void test_mpmul();

/* ---------- TEST MATRIX.C ---------- */

void test_init_matrix(u64 n);

void test_init_eye(u64 n);

void test_copy_matrix(u64 n);

void test_matrix_add(u64 n);

void test_matrix_sub(u64 n);

void test_matrix_mul(u64 n);

void test_matrix_mul_coef(u64 n);





/* ---------- TEST ALGO.C ---------- */


#endif // TESTS_H