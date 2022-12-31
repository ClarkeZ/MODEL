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

/* ***** MPFR ***** */

void test_mpadd(){
    printf(" --- TEST MPADD --- \n");
    mpfr a, b, c;
    mpinit(a);
    mpinit(b);
    mpinit(c);

    double x = randfrom(0 - rand(), rand());
    double y = randfrom(0 - rand(), rand());

    mpset(a, x);
    mpset(b, y);

    mpadd(c, a, b);

    printf("%f + %f = ", x, y);
    mpfr_out_str(stdout, 10, 0, c, MPFR_RNDD);
    printf("\n");
    mpfr_printf("%f + %f = %.10Rf\n",x, y, c); // 10 digits

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(c);
}

void test_mpsub(){
    printf(" --- TEST MPSUB --- \n");
    mpfr a, b, c;
    mpinit(a);
    mpinit(b);
    mpinit(c);

    double x = randfrom(0 - rand(), rand());
    double y = randfrom(0 - rand(), rand());

    mpset(a, x);
    mpset(b, y);

    mpsub(c, a, b);

    printf("%f - %f = ", x, y);
    mpfr_out_str(stdout, 10, 0, c, MPFR_RNDD);
    printf("\n");
    mpfr_printf("%f - %f = %.10Rf\n",x, y, c); // 10 digits

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(c);
}

void test_mpmul(){
    printf(" --- TEST MPMUL --- \n");
    mpfr a, b, c;
    mpinit(a);
    mpinit(b);
    mpinit(c);

    double x = randfrom(0 - rand(), rand());
    double y = randfrom(0 - rand(), rand());

    mpset(a, x);
    mpset(b, y);

    mpmul(c, a, b);

    printf("%f * %f = ", x, y);
    mpfr_out_str(stdout, 10, 0, c, MPFR_RNDD);
    printf("\n");
    mpfr_printf("%f * %f = %.10Rf\n",x, y, c); // 10 digits

    mpfr_clear(a);
    mpfr_clear(b);
    mpfr_clear(c);
}

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
    printf("--- A ---\n");
    print_matrix(A);
    printf("--- B ---\n");
    print_matrix(B);

    tic = wtime();
    res = matrix_add(A, B);
    toc = wtime();

    printf("--- A + B ---\n");
    print_matrix(res);

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
    printf("--- A ---\n");
    print_matrix(A);
    printf("--- B ---\n");
    print_matrix(B);

    tic = wtime();
    res = matrix_sub(A, B);
    toc = wtime();

    printf("--- A - B ---\n");
    print_matrix(res);

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

void test_matrix_transpose(int n){
    printf("--- TEST MATRIX_TRANSPOSE ---\n");
    Matrix *A = init_matrix(n);
    Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i)
        A->m[i] = randfrom(0 - rand(), rand());

    printf("--- A ---\n");
    print_matrix(A);

    tic = wtime();
    res = matrix_transpose(A);
    toc = wtime();

    printf("--- A^T ---\n");
    print_matrix(res);

    printf("time = %f\n", toc - tic);

    free_matrix(A);
    free_matrix(res);
}

/* ***** MPFR ***** */

void test_init_MPFR_matrix(int n){
    printf("--- TEST INIT_MPFR_MATRIX ---\n");
    MPFR_Matrix *M = init_MPFR_matrix(n);
    print_MPFR_matrix(M);

    free_MPFR_matrix(M);
}

void test_init_MPFR_eye(int n){
    printf("--- TEST INIT_MPFR_EYE ---\n");
    MPFR_Matrix *M = init_MPFR_eye(n);
    print_MPFR_matrix(M);

    free_MPFR_matrix(M);
}

void test_MPFR_copy_matrix(int n){
    printf("--- TEST COPY_MPFR_MATRIX ---\n");
    MPFR_Matrix *M = init_MPFR_matrix(n);
    MPFR_Matrix *copy;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < M->n * M->n ; ++i) 
        mpset(M->m[i], randfrom(0 - rand(), rand()));

    printf("--- M AVANT ---\n");
    print_MPFR_matrix(M);

    tic = wtime();
    copy = copy_MPFR_matrix(M);
    toc = wtime();

    printf("--- M APRES ---\n");
    print_MPFR_matrix(M);
    printf("--- COPIE de M ---\n");
    print_MPFR_matrix(copy);

    printf("time = %f\n", toc - tic);

    free_MPFR_matrix(M);
    free_MPFR_matrix(copy);
}

void test_MPFR_matrix_add(int n){
    printf("--- TEST MPFR_MATRIX_ADD ---\n");
    MPFR_Matrix *A = init_MPFR_matrix(n);
    MPFR_Matrix *B = init_MPFR_matrix(n);
    MPFR_Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        mpset(A->m[i], randfrom(0 - rand(), rand()));
        mpset(B->m[i], randfrom(0 - rand(), rand()));
    }
    printf("--- A ---\n");
    print_MPFR_matrix(A);
    printf("--- B ---\n");
    print_MPFR_matrix(B);

    tic = wtime();
    res = MPFR_matrix_add(A, B);
    toc = wtime();

    printf("--- A + B ---\n");
    print_MPFR_matrix(res);

    printf("time = %f\n", toc - tic);

    free_MPFR_matrix(A);
    free_MPFR_matrix(B);
    free_MPFR_matrix(res);
}

void test_MPFR_matrix_sub(int n){
    printf("--- TEST MPFR_MATRIX_SUB ---\n");
    MPFR_Matrix *A = init_MPFR_matrix(n);
    MPFR_Matrix *B = init_MPFR_matrix(n);
    MPFR_Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        mpset(A->m[i], randfrom(0 - rand(), rand()));
        mpset(B->m[i], randfrom(0 - rand(), rand()));
    }
    printf("--- A ---\n");
    print_MPFR_matrix(A);
    printf("--- B ---\n");
    print_MPFR_matrix(B);

    tic = wtime();
    res = MPFR_matrix_sub(A, B);
    toc = wtime();

    printf("--- A - B ---\n");
    print_MPFR_matrix(res);

    printf("time = %f\n", toc - tic);

    free_MPFR_matrix(A);
    free_MPFR_matrix(B);
    free_MPFR_matrix(res);
}

void test_MPFR_matrix_mul(int n){
    printf("--- TEST MPFR_MATRIX_MUL ---\n");
    MPFR_Matrix *A = init_MPFR_matrix(n);
    MPFR_Matrix *B = init_MPFR_matrix(n);
    MPFR_Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        // mpset(A->m[i], randfrom(0 - rand(), rand()));
        // mpset(B->m[i], randfrom(0 - rand(), rand()));
        // mpset(A->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)));
        // mpset(B->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)));
        mpset(A->m[i], rand() % (MAX-MIN));
        mpset(B->m[i], rand() % (MAX-MIN));
    }
    printf("--- A ---\n");
    print_MPFR_matrix(A);
    printf("--- B ---\n");
    print_MPFR_matrix(B);

    tic = wtime();
    res = MPFR_matrix_mul(A, B);
    toc = wtime();

    printf("--- A * B ---\n");
    print_MPFR_matrix(res);

    printf("time = %f\n", toc - tic);

    free_MPFR_matrix(A);
    free_MPFR_matrix(B);
    free_MPFR_matrix(res);
}

void test_MPFR_matrix_transpose(int n){
    printf("--- TEST MPFR_MATRIX_TRANSPOSE ---\n");
    MPFR_Matrix *A = init_MPFR_matrix(n);
    MPFR_Matrix *res;
    double tic, toc;
    unsigned int i;

    for (i = 0 ; i < A->n * A->n ; ++i){
        // mpset(A->m[i], randfrom(0 - rand(), rand()));
        // mpset(A->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)));
        mpset(A->m[i], rand() % (MAX-MIN));
    }
    printf("--- A ---\n");
    print_MPFR_matrix(A);

    tic = wtime();
    res = MPFR_matrix_transpose(A);
    toc = wtime();

    printf("--- A^T ---\n");
    print_MPFR_matrix(res);

    printf("time = %f\n", toc - tic);

    free_MPFR_matrix(A);
    free_MPFR_matrix(res);
}

/* ---------- TEST ALGO.C ---------- */

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

void test_hessenberg(int n){
    printf("--- TEST HESSENBERG ---\n");
    Matrix *M = init_matrix(n);
    Matrix *hess = init_matrix(n);
    int i;
    double tic, toc;

    for (i = 0 ; i < n * n ; ++i) 
        M->m[i] = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)); // Pour des doubles
        // M->m[i] = rand() % (MAX-MIN); // Pour des entiers


    printf("--- AVANT ---\n");
    print_matrix(M);

    tic = wtime();
    hess = hessenberg(M);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- HESSENBERG ---\n");
    print_matrix(hess);

    printf("\ntime = %f\n", toc - tic);

    free_matrix(M);
    free_matrix(hess);
}

void test_eigenvalues(int n){
    printf("--- TEST EIGENVALUE ---\n");
    Matrix *M = init_matrix(n);
    int i;
    double tic, toc;

    for (i = 0 ; i < n * n ; ++i) 
        // M->m[i] = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN));
        M->m[i] = rand() % (MAX-MIN); // Pour des entiers

    printf("--- MATRICE A ---\n");
    print_matrix(M);
    Matrix *hess = hessenberg(M);
    Matrix *quasi_h = quasi_hess(M);

    tic = wtime();
    double *eigen = eigenvalues(M);
    toc = wtime();
    free_matrix(M);
    
    printf("\n--- EIGENVALUE DE LA MATRICE A ---\n");
    for(i = 0 ; i < n ; ++i)
        printf("%f\t", eigen[i]);
    printf("\n");
    printf("\ntime = %f\n", toc - tic);


    printf("\n--- MATRICE QUASI HESSENBERG ---\n");
    print_matrix(quasi_h);
    tic = wtime();
    double *eigen_quasi_h = eigenvalues(quasi_h);
    toc = wtime();
    free_matrix(quasi_h);

    printf("\n--- EIGENVALUE DE LA MATRICE QUASI HESSENBERG ---\n");
    for(i = 0 ; i < n ; ++i)
        printf("%f\t", eigen_quasi_h[i]);
    printf("\n");
    printf("\ntime = %f\n", toc - tic);


    printf("\n--- MATRICE HESSENBERG ---\n");
    print_matrix(hess);
    tic = wtime();
    double *eigen_hess = eigenvalues(hess);
    toc = wtime();
    free_matrix(hess);

    printf("\n--- EIGENVALUE DE LA MATRICE HESSENBERG ---\n");
    for(i = 0 ; i < n ; ++i)
        printf("%f\t", eigen_hess[i]);
    printf("\n");
    printf("\ntime = %f\n", toc - tic);


    free(eigen);
    free(eigen_quasi_h);
    free(eigen_hess);
}

/* ***** MPFR ***** */

void test_MPFR_qr_decomposition(int n){
    printf("--- TEST MPFR QR_DECOMPOSITION ---\n");
    MPFR_Matrix *M = init_MPFR_matrix(n);
    MPFR_QR *qr = (MPFR_QR *) malloc(sizeof(MPFR_QR));
    int i;
    double tic, toc;

    for (i = 0 ; i < n * n ; ++i) 
        // mpset(M->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)));
        mpset(M->m[i], rand() % (MAX-MIN)); // Pour des entiers

    printf("--- AVANT ---\n");
    print_MPFR_matrix(M);
    
    tic = wtime();
    qr = MPFR_qr_decomposition(M);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- Q ---\n");
    print_MPFR_matrix(qr->Q);
    printf("--- R ---\n");
    print_MPFR_matrix(qr->R);

    printf("--- VERIF = AVANT ---\n");
    MPFR_Matrix *tmp_qr = MPFR_matrix_mul(qr->Q, qr->R);
    print_MPFR_matrix(tmp_qr);

    free_MPFR_matrix(tmp_qr);

    printf("\ntime = %f\n", toc - tic);

    free_MPFR_matrix(M);
    free_MPFR_qr(qr);
}

void test_MPFR_quasi_hess(int n){
    printf("--- TEST MPFR QUASI HESSENBERG ---\n");
    MPFR_Matrix *M = init_MPFR_matrix(n);
    MPFR_Matrix *hess = init_MPFR_matrix(n);
    int i;
    double tic, toc;

    for (i = 0 ; i < n * n ; ++i) 
        // mpset(M->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)));
        mpset(M->m[i], rand() % (MAX-MIN)); // Pour des entiers

    printf("--- AVANT ---\n");
    print_MPFR_matrix(M);

    tic = wtime();
    hess = MPFR_quasi_hess(M);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- QUASI HESSENBERG ---\n");
    print_MPFR_matrix(hess);

    printf("\ntime = %f\n", toc - tic);

    free_MPFR_matrix(M);
    free_MPFR_matrix(hess);
}

void test_MPFR_hessenberg(int n){
    printf("--- TEST MPFR HESSENBERG ---\n");
    MPFR_Matrix *M = init_MPFR_matrix(n);
    MPFR_Matrix *hess = init_MPFR_matrix(n);
    int i;
    double tic, toc;

    for (i = 0 ; i < n * n ; ++i) 
        // mpset(M->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)));
        mpset(M->m[i], rand() % (MAX-MIN)); // Pour des entiers

    printf("--- AVANT ---\n");
    print_MPFR_matrix(M);

    tic = wtime();
    hess = MPFR_hessenberg(M);
    toc = wtime();

    printf("--- APRES ---\n");
    printf("--- HESSENBERG ---\n");
    print_MPFR_matrix(hess);

    printf("\ntime = %f\n", toc - tic);

    free_MPFR_matrix(M);
    free_MPFR_matrix(hess);
}


void test_MPFR_eigenvalues(int n){
    printf("--- TEST MPFR EIGENVALUES ---\n");
    MPFR_Matrix *M = init_MPFR_matrix(n);
    int i;
    double tic, toc;

    mpfr *eigen = (mpfr *) malloc(n * sizeof(mpfr));
    mpfr *eigen_quasi_h = (mpfr *) malloc(n * sizeof(mpfr));
    mpfr *eigen_h = (mpfr *) malloc(n * sizeof(mpfr));
    for(int i = 0 ; i < n ; i++){
        mpinit(eigen[i]);
        mpinit(eigen_quasi_h[i]);
        mpinit(eigen_h[i]);
    }

    for (i = 0 ; i < n * n ; ++i) 
        // mpset(M->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)));
        mpset(M->m[i], rand() % (MAX-MIN)); // Pour des entiers

    printf("--- MATRICE A ---\n");
    print_MPFR_matrix(M);
    MPFR_Matrix *hess = MPFR_hessenberg(M);
    MPFR_Matrix *quasi_h = MPFR_quasi_hess(M);

    tic = wtime();
    MPFR_eigenvalues(M, eigen);
    toc = wtime();
    free_MPFR_matrix(M);
    
    printf("\n--- EIGENVALUE DE LA MATRICE A ---\n");
    for(i = 0 ; i < n ; ++i)
        mpfr_printf("%.10Rf\t", eigen[i]);
    printf("\n");
    printf("\ntime = %f\n", toc - tic);

    printf("\n--- MATRICE QUASI HESSENBERG ---\n");
    print_MPFR_matrix(quasi_h);
    tic = wtime();
    MPFR_eigenvalues(quasi_h, eigen_quasi_h);
    toc = wtime();
    free_MPFR_matrix(quasi_h);

    printf("\n--- EIGENVALUE DE LA MATRICE QUASI HESSENBERG ---\n");
    for(i = 0 ; i < n ; ++i)
        mpfr_printf("%.10Rf\t", eigen_quasi_h[i]);
    printf("\n");
    printf("\ntime = %f\n", toc - tic);

    printf("\n--- MATRICE HESSENBERG ---\n");
    print_MPFR_matrix(hess);
    tic = wtime();
    MPFR_eigenvalues(hess, eigen_h);
    toc = wtime();
    free_MPFR_matrix(hess);

    printf("\n--- EIGENVALUE DE LA MATRICE HESSENBERG ---\n");
    for(i = 0 ; i < n ; ++i)
        mpfr_printf("%.10Rf\t", eigen_h[i]);
    printf("\n");
    printf("\ntime = %f\n", toc - tic);
    
    free(eigen);
    free(eigen_quasi_h);
    free(eigen_h);
}