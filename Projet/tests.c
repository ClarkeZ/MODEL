#include "tests.h"

/* ---------- TEST BASE.C ---------- */

void test_add(){
    printf(" --- TEST ADD --- \n");
    double a = randfrom(0 - rand(), rand());
    double b = randfrom(0 - rand(), rand());
    printf("%f + %f = %f\n", a, b, add(a, b));
}

void test_sub(){
    printf(" --- TEST SUB --- \n");
    double a = randfrom(0 - rand(), rand());
    double b = randfrom(0 - rand(), rand());
    printf("%f - %f = %f\n", a, b, sub(a, b));
}

void test_mul(){
    printf(" --- TEST MUL --- \n");
    double a = randfrom(0 - rand(), rand());
    double b = randfrom(0 - rand(), rand());
    printf("%f * %f = %f\n", a, b, mul(a, b));
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

void test_init_matrix(int n){
    printf("--- TEST INIT_MATRIX ---\n");
    Matrix *M = init_matrix(n);
    print_matrix(M);

    free_matrix(M);
}

void test_init_eye(int n){
    printf("--- TEST INIT_EYE ---\n");
    Matrix *M = init_eye(n);
    print_matrix(M);

    free_matrix(M);
}

void test_copy_matrix(int n){
    printf("--- TEST COPY_MATRIX ---\n");
    Matrix *M = init_matrix(n);
    Matrix *copy;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        M->m[i] = randfrom(0 - rand(), rand());

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

void test_matrix_add(int n){
    printf("--- TEST MATRIX_ADD ---\n");
    Matrix *A = init_matrix(n);
    Matrix *B = init_matrix(n);
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = randfrom(0 - rand(), rand());
        B->m[i] = randfrom(0 - rand(), rand());
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

void test_matrix_sub(int n){
    printf("--- TEST MATRIX_SUB ---\n");
    Matrix *A = init_matrix(n);
    Matrix *B = init_matrix(n);
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = randfrom(0 - rand(), rand());
        B->m[i] = randfrom(0 - rand(), rand());
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

void test_matrix_mul(int n){
    printf("--- TEST MATRIX_MUL ---\n");
    Matrix *A = init_matrix(n);
    Matrix *B = init_matrix(n);
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        A->m[i] = randfrom(0 - rand(), rand());
        B->m[i] = randfrom(0 - rand(), rand());
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

void test_matrix_mul_coef(int n){
    printf("--- TEST MATRIX_MUL_COEF ---\n");
    Matrix *A = init_matrix(n);
    double c = randfrom(0 - rand(), rand());
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i)
        A->m[i] = randfrom(0 - rand(), rand());

    printf("--- A ---\n");
    print_matrix(A);
    printf("--- c ---\n");
    printf("%f\n", c);

    tic = wtime();
    res = matrix_mul_coef(A, c);
    toc = wtime();

    printf("--- A * c ---\n");
    print_matrix(res);

    printf("time = %f\n", toc - tic);

    free_matrix(A);
    free_matrix(res);
}


void test_qr_decomposition(int n){
    printf("--- TEST QR_DECOMPOSITION ---\n");
    Matrix *M = init_matrix(n);
    QR *qr = (QR *) malloc(sizeof(QR));
    int i;
    double tic, toc;

    for (i = 0 ; i < n * n ; ++i) 
        // M->m[i] = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN));
        M->m[i] = rand() % (MAX-MIN); // Pour des entiers

    printf("--- AVANT ---\n");
    print_matrix(M);
    
    tic = wtime();
    qr = qr_decomposition(M);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- Q ---\n");
    print_matrix(qr->Q);
    printf("--- R ---\n");
    print_matrix(qr->R);

    printf("--- VERIF = AVANT ---\n");
    Matrix *tmp_qr = matrix_mul(qr->Q, qr->R);
    print_matrix(tmp_qr);

    free_matrix(tmp_qr);

    printf("\ntime = %f\n", toc - tic);

    free_matrix(M);
    free_qr(qr);
}

void test_quasi_hess(int n){
    printf("--- TEST QUASI HESSENBERG ---\n");
    Matrix *M = init_matrix(n);
    Matrix *hess = init_matrix(n);
    int i;
    double tic, toc;

    for (i = 0 ; i < n * n ; ++i) 
        // M->m[i] = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN));
        M->m[i] = rand() % (MAX-MIN); // Pour des entiers

    printf("--- AVANT ---\n");
    print_matrix(M);

    tic = wtime();
    hess = quasi_hess(M);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- QUASI HESSENBERG ---\n");
    print_matrix(hess);

    printf("\ntime = %f\n", toc - tic);

    free_matrix(M);
    free_matrix(hess);
}