#include "tests.h"

/* ---------- TEST BASE.C ---------- */

void test_add(){
    printf(" --- TEST ADD --- \n");
    u64 a = rand();
    u64 b = rand();
    printf("%llu + %llu = %llu\n", a, b, add(a, b));
}

void test_sub(){
    printf(" --- TEST SUB --- \n");
    u64 a = rand();
    u64 b = rand();
    printf("%llu - %llu = %llu\n", a, b, sub(a, b));
}

void test_mul(){
    printf(" --- TEST MUL --- \n");
    u64 a = rand();
    u64 b = rand();
    printf("%llu * %llu = %llu\n", a, b, mul(a, b));
}

// void test_mpadd(){
//     printf(" --- TEST MPADD --- \n");
//     mpfr a, b, c;
//     mpinit(a);
//     mpinit(b);
//     mpinit(c);

//     mpset(a, 1.0);
//     mpset(b, 2.0);

//     mpadd(c, a, b);

//     mpfr_out_str (stdout, 10, 0, c, MPFR_RNDD);
//     mpfr_clear(a);
//     mpfr_clear(b);
//     mpfr_clear(c);
// }

/* ---------- TEST MATRIX.C ---------- */

void test_init_matrix(u64 n){
    printf("--- TEST INIT_MATRIX ---\n");
    Matrix *M = init_matrix(n);
    print_matrix(M);

    free_matrix(M);
}

void test_init_eye(u64 n){
    printf("--- TEST INIT_EYE ---\n");
    Matrix *M = init_eye(n);
    print_matrix(M);

    free_matrix(M);
}

void test_copy_matrix(u64 n){
    printf("--- TEST COPY_MATRIX ---\n");
    Matrix *M = init_matrix(n);
    Matrix *copy;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = rand();

    printf("--- M AVANT ---\n");
    print_matrix(M);

    tic = wtime();
    copy = copy_matrix(M);
    toc = wtime();

    printf("--- M APRES ---\n");
    print_matrix(M);
    printf("--- COPIE de M ---\n");
    print_matrix(copy);

    printf("time = %f\n", toc - tic);

    free_matrix(M);
    free_matrix(copy);
}

void test_matrix_add(u64 n){
    printf("--- TEST MATRIX_ADD ---\n");
    Matrix *A = init_matrix(n);
    Matrix *B = init_matrix(n);
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand();
        B->m[i] = rand();
    }
    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- B ---\n");
    // print_matrix(B);

    tic = wtime();
    res = matrix_add(A, B);
    toc = wtime();

    // printf("--- A + B ---\n");
    // print_matrix(res);

    printf("time = %f\n", toc - tic);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_matrix_sub(u64 n){
    printf("--- TEST MATRIX_SUB ---\n");
    Matrix *A = init_matrix(n);
    Matrix *B = init_matrix(n);
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand();
        B->m[i] = rand();
    }
    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- B ---\n");
    // print_matrix(B);

    tic = wtime();
    res = matrix_sub(A, B);
    toc = wtime();

    // printf("--- A - B ---\n");
    // print_matrix(res);

    printf("time = %f\n", toc - tic);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_matrix_mul(u64 n){
    printf("--- TEST MATRIX_MUL ---\n");
    Matrix *A = init_matrix(n);
    Matrix *B = init_matrix(n);
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = rand();
        B->m[i] = rand();
    }
    printf("--- A ---\n");
    print_matrix(A);
    printf("--- B ---\n");
    print_matrix(B);

    tic = wtime();
    res = matrix_mul(A, B);
    toc = wtime();

    printf("--- A * B ---\n");
    print_matrix(res);

    printf("time = %f\n", toc - tic);

    free_matrix(A);
    free_matrix(B);
    free_matrix(res);
}

void test_matrix_mul_coef(u64 n){
    printf("--- TEST MATRIX_MUL_COEF ---\n");
    Matrix *A = init_matrix(n);
    u64 c = rand();
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i)
        A->m[i] = rand();

    // printf("--- A ---\n");
    // print_matrix(A);
    // printf("--- c ---\n");
    // printf("%u\n", c);

    tic = wtime();
    res = matrix_mul_coef(A, c);
    toc = wtime();

    // printf("--- A * c ---\n");
    // print_matrix(res);

    printf("time = %f\n", toc - tic);

    free_matrix(A);
    free_matrix(res);
}