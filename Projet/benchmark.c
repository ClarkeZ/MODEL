#include "benchmark.h"

void benchmark_quasi_hess_vs_hessenberg(unsigned int N, unsigned int ITE) {
    printf("\n--- Quasi Hessenberg vs Hessenberg (double precision) ---\n");
    fflush(stdout);
    unsigned int i, k;

    double time_quasiH = 0;
    double time_H = 0;
    double tic, toc, tic1, toc1;

    Matrix *A;
    Matrix *tmp;

    Matrix *res_quasiH;
    Matrix *res_H;
    for(k = 0 ; k < ITE ; k++){
        A = init_matrix(N);
        for (i = 0 ; i < A->n * A->n ; ++i) 
            A->m[i] = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)); // Pour des doubles
            // M->m[i] = rand() % (MAX-MIN); // Pour des entiers

        // printf("--- Matrice A ---\n");
        // print_matrix(A);
        tmp = copy_matrix(A);

        // printf("--- Quasi Hessenberg ---\n");
        tic1 = wtime();
        res_quasiH = quasi_hess(tmp);
        toc1 = wtime();
        time_quasiH += toc1 - tic1;
        // print_matrix(res_quasiH);
        free_matrix(tmp);
        free_matrix(res_quasiH);
        
        // printf("--- Hessenberg ---\n");
        tic = wtime();
        res_H = hessenberg(A);
        toc = wtime();
        time_H += toc - tic;
        // print_matrix(res_H);
        free_matrix(A);
        free_matrix(res_H);
    }

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH);
    printf("--- Hessenberg          = %f ---\n", time_H);
    printf("##### Temps moyen #####\n");
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH / ITE);
    printf("--- Hessenberg          = %f ---\n", time_H / ITE);

    printf("%s\n", (time_quasiH < time_H) ? "Quasi Hessenberg plus rapide" : "Hessenberg plus rapide");
}

void benchmark_MPFR_quasi_hess_vs_hessenberg(unsigned int N, unsigned int ITE) {
    printf("\n--- Quasi Hessenberg vs Hessenberg (MPFR) ---\n");
    fflush(stdout);
    unsigned int i, k;

    double time_quasiH = 0;
    double time_H = 0;
    double tic, toc, tic1, toc1;

    MPFR_Matrix *A;
    MPFR_Matrix *tmp;

    MPFR_Matrix *res_quasiH;
    MPFR_Matrix *res_H;
    for(k = 0 ; k < ITE ; k++){
        A = init_MPFR_matrix(N);
        for (i = 0 ; i < A->n * A->n ; ++i){
            mpset(A->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN))); // Pour des doubles
            // mpset(A->m[i], rand() % (MAX-MIN)); // Pour des entiers
        }

        // printf("--- Matrice A ---\n");
        // print_matrix(A);
        tmp = copy_MPFR_matrix(A);

        // printf("--- Quasi Hessenberg ---\n");
        tic1 = wtime();
        res_quasiH = MPFR_quasi_hess(tmp);
        toc1 = wtime();
        time_quasiH += toc1 - tic1;
        // print_matrix(res_quasiH);
        free_MPFR_matrix(tmp);
        free_MPFR_matrix(res_quasiH);
        
        // printf("--- Hessenberg ---\n");
        tic = wtime();
        res_H = MPFR_hessenberg(A);
        toc = wtime();
        time_H += toc - tic;
        // print_matrix(res_H);
        free_MPFR_matrix(A);
        free_MPFR_matrix(res_H);
    }

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH);
    printf("--- Hessenberg          = %f ---\n", time_H);
    printf("##### Temps moyen #####\n");
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH / ITE);
    printf("--- Hessenberg          = %f ---\n", time_H / ITE);

    printf("%s\n", (time_quasiH < time_H)? "Quasi Hessenberg plus rapide" : "Hessenberg plus rapide");
}

void benchmark_qr_decomposition(unsigned int N, unsigned int ITE){
    printf("\n--- QR decomposition (double précision vs MPFR) ---\n");
    fflush(stdout);
    unsigned int i, k;

    double time_qr = 0;
    double time_MPFR_qr = 0;
    double tic, toc, tic1, toc1;

    double r;

    Matrix *A;
    MPFR_Matrix *MPFR_A;
    // Matrix *tmp = init_matrix(N);
    // MPFR_Matrix *MPFR_tmp = init_MPFR_matrix(N);

    QR *res_qr;
    MPFR_QR *res_MPFR_qr;

    for(k = 0 ; k < ITE ; k++){
        A = init_matrix(N);
        MPFR_A = init_MPFR_matrix(N);
        for (i = 0 ; i < N * N ; ++i){
            r = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)); // Pour des doubles
            // r = rand() % (MAX-MIN); // Pour des entiers
            A->m[i] = r; 
            mpset(MPFR_A->m[i], r); 
        }

        // printf("--- Matrice A ---\n");
        // print_matrix(A);
        // printf("--- MPFR A ---\n");
        // print_MPFR_matrix(MPFR_A);
        // tmp = copy_matrix(A);
        // MPFR_tmp = copy_MPFR_matrix(MPFR_A);

        // printf("--- QR (double precision) ---\n");
        tic1 = wtime();
        res_qr = qr_decomposition(A);
        toc1 = wtime();
        time_qr += toc1 - tic1;
        // printf("Q = \n");
        // print_matrix(res_qr->Q);
        // printf("R = \n");
        // print_matrix(res_qr->R);
        free_matrix(A);
        free_qr(res_qr);
        
        // printf("--- QR (MPFR) ---\n");
        tic = wtime();
        res_MPFR_qr = MPFR_qr_decomposition(MPFR_A);
        toc = wtime();
        time_MPFR_qr += toc - tic;
        // printf("MPFR Q = \n");
        // print_MPFR_matrix(res_MPFR_qr->Q);
        // printf("MPFR R = \n");
        // print_MPFR_matrix(res_MPFR_qr->R);
        free_MPFR_matrix(MPFR_A);
        free_MPFR_qr(res_MPFR_qr);
    }

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- QR decomposition (double precision) = %f ---\n", time_qr);
    printf("--- QR decomposition (MPFR)             = %f ---\n", time_MPFR_qr);
    printf("##### Temps moyen #####\n");
    printf("--- QR decomposition (double precision) = %f ---\n", time_qr / ITE);
    printf("--- QR decomposition (MPFR)             = %f ---\n", time_MPFR_qr / ITE);

    printf("%s\n", (time_qr < time_MPFR_qr)? "double precision plus rapide" : "MPFR plus rapide");
}

void benchmark_quasi_hess(unsigned int N, unsigned int ITE){
    printf("\n--- Quasi Hessenberg (double precision vs MPFR) ---\n");
    fflush(stdout);
    unsigned int i, k;

    double time_H = 0;
    double time_MPFR_H = 0;
    double tic, toc, tic1, toc1;

    double r;

    Matrix *A;
    MPFR_Matrix *MPFR_A;
    // Matrix *tmp = init_matrix(N);
    // MPFR_Matrix *MPFR_tmp = init_MPFR_matrix(N);

    Matrix *res_H;
    MPFR_Matrix *res_MPFR_H;

    for(k = 0 ; k < ITE ; k++){
        A = init_matrix(N);
        MPFR_A = init_MPFR_matrix(N);
        for (i = 0 ; i < N * N ; ++i){
            r = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)); // Pour des doubles
            // r = rand() % (MAX-MIN); // Pour des entiers
            A->m[i] = r; 
            mpset(MPFR_A->m[i], r); 
        }

        // printf("--- Matrice A ---\n");
        // print_matrix(A);
        // printf("--- MPFR A ---\n");
        // print_MPFR_matrix(MPFR_A);
        // tmp = copy_matrix(A);
        // MPFR_tmp = copy_MPFR_matrix(MPFR_A);

        // printf("--- Quasi Hessenberg (double precision) ---\n");
        tic1 = wtime();
        res_H = quasi_hess(A);
        toc1 = wtime();
        time_H += toc1 - tic1;
        // printf("H = \n");
        // print_matrix(res_H);
        free_matrix(A);
        free_matrix(res_H);
        
        // printf("--- Quasi Hessenberg (MPFR) ---\n");
        tic = wtime();
        res_MPFR_H = MPFR_quasi_hess(MPFR_A);
        toc = wtime();
        time_MPFR_H += toc - tic;
        // printf("MPFR H = \n");
        // print_MPFR_matrix(res_MPFR_H);
        free_MPFR_matrix(MPFR_A);
        free_MPFR_matrix(res_MPFR_H);
    }

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- Quasi Hessenberg (double precision)  = %f ---\n", time_H);
    printf("--- Quasi Hessenberg (MPFR)              = %f ---\n", time_MPFR_H);
    printf("##### Temps moyen #####\n");
    printf("--- Quasi Hessenberg (double precision)  = %f ---\n", time_H / ITE);
    printf("--- Quasi Hessenberg (MPFR)              = %f ---\n", time_MPFR_H / ITE);

    printf("%s\n", (time_H < time_MPFR_H)? "double precision plus rapide" : "MPFR plus rapide");
}

void benchmark_hessenberg(unsigned int N, unsigned int ITE){
    printf("\n--- Hessenberg (double precision vs MPFR)  ---\n");
    fflush(stdout);
    unsigned int i, k;

    double time_H = 0;
    double time_MPFR_H = 0;
    double tic, toc, tic1, toc1;

    double r;

    Matrix *A;
    MPFR_Matrix *MPFR_A;
    // Matrix *tmp = init_matrix(N);
    // MPFR_Matrix *MPFR_tmp = init_MPFR_matrix(N);

    Matrix *res_H;
    MPFR_Matrix *res_MPFR_H;

    for(k = 0 ; k < ITE ; k++){
        A = init_matrix(N);
        MPFR_A = init_MPFR_matrix(N);
        for (i = 0 ; i < N * N ; ++i){
            r = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)); // Pour des doubles
            // r = rand() % (MAX-MIN); // Pour des entiers
            A->m[i] = r; 
            mpset(MPFR_A->m[i], r); 
        }

        // printf("--- Matrice A ---\n");
        // print_matrix(A);
        // printf("--- MPFR A ---\n");
        // print_MPFR_matrix(MPFR_A);
        // tmp = copy_matrix(A);
        // MPFR_tmp = copy_MPFR_matrix(MPFR_A);

        // printf("--- Hessenberg (double precision) ---\n");
        tic1 = wtime();
        res_H = hessenberg(A);
        toc1 = wtime();
        time_H += toc1 - tic1;
        // printf("H = \n");
        // print_matrix(res_H);
        free_matrix(A);
        free_matrix(res_H);
        
        // printf("--- MPFR Hessenberg (MPFR) ---\n");
        tic = wtime();
        res_MPFR_H = MPFR_hessenberg(MPFR_A);
        toc = wtime();
        time_MPFR_H += toc - tic;
        // printf("MPFR H = \n");
        // print_MPFR_matrix(res_MPFR_H);
        free_MPFR_matrix(MPFR_A);
        free_MPFR_matrix(res_MPFR_H);
    }

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- Hessenberg (double precision) = %f ---\n", time_H);
    printf("--- Hessenberg (MPFR)             = %f ---\n", time_MPFR_H);
    printf("##### Temps moyen #####\n");
    printf("--- Hessenberg (double precision) = %f ---\n", time_H / ITE);
    printf("--- Hessenberg (MPFR)             = %f ---\n", time_MPFR_H / ITE);

    printf("%s\n", (time_H < time_MPFR_H)? "double precision plus rapide" : "MPFR plus rapide");
}