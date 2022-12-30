#include "benchmark.h"

void benchmark_quasi_hess_vs_hessenberg(unsigned int N, unsigned int ITE) {
    printf("\n--- Quasi Hessenberg vs Hessenberg ---\n");
    unsigned int i, k;

    double time_quasiH = 0;
    double time_H = 0;
    double tic, toc, tic1, toc1;

    Matrix *A = init_matrix(N);
    Matrix *tmp = init_matrix(N);
    Matrix *tmp2 = init_matrix(N);

    Matrix *res_quasiH;
    Matrix *res_H;
    for(k = 0 ; k < ITE ; k++){
        for (i = 0 ; i < A->n * A->n ; ++i) 
            A->m[i] = randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN)); // Pour des doubles
            // M->m[i] = rand() % (MAX-MIN); // Pour des entiers

        // printf("--- Matrice A ---\n");
        // print_matrix(A);
        tmp = copy_matrix(A);
        tmp2 = copy_matrix(A);

        // printf("--- Quasi Hessenberg ---\n");
        tic1 = wtime();
        res_quasiH = quasi_hess(tmp);
        toc1 = wtime();
        free_matrix(res_quasiH);
        time_quasiH += toc1 - tic1;
        // print_matrix(res_quasiH);
        
        // printf("--- Hessenberg ---\n");
        tic = wtime();
        res_H = hessenberg(tmp2);
        toc = wtime();
        free_matrix(res_H);
        time_H += toc - tic;
        // print_matrix(res_H);

        free_matrix(A);
        free_matrix(tmp);
        free_matrix(tmp2);
    }

    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH);
    printf("--- Hessenberg          = %f ---\n", time_H);
    printf("##### Temps moyen #####\n");
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH / ITE);
    printf("--- Hessenberg          = %f ---\n", time_H / ITE);

    printf("%s\n", (time_quasiH-time_H >0)? "Quasi Hessenberg plus rapide" : "Hessenberg plus rapide");
}

void benchmark_MPFR_quasi_hess_vs_hessenberg(unsigned int N, unsigned int ITE) {
    printf("\n--- MPFR Quasi Hessenberg vs Hessenberg ---\n");
    unsigned int i, k;

    double time_quasiH = 0;
    double time_H = 0;
    double tic, toc, tic1, toc1;

    MPFR_Matrix *A;
    MPFR_Matrix *tmp = init_MPFR_matrix(N);
    MPFR_Matrix *tmp2 = init_MPFR_matrix(N);

    MPFR_Matrix *res_quasiH;
    MPFR_Matrix *res_H;
    for(k = 0 ; k < ITE ; k++){
        A = init_MPFR_matrix(N);
        for (i = 0 ; i < A->n * A->n ; ++i){
            // mpinit(A->m[i]);
            mpset(A->m[i], randfrom(0 - rand() % (MAX - MIN), rand() % (MAX - MIN))); // Pour des doubles
        }
            // M->m[i] = rand() % (MAX-MIN); // Pour des entiers

        // printf("--- Matrice A ---\n");
        // print_matrix(A);
        tmp = copy_MPFR_matrix(A);
        tmp2 = copy_MPFR_matrix(A);

        // printf("--- Quasi Hessenberg ---\n");
        tic1 = wtime();
        res_quasiH = MPFR_quasi_hess(tmp);
        toc1 = wtime();
        free_MPFR_matrix(res_quasiH);
        time_quasiH += toc1 - tic1;
        // print_matrix(res_quasiH);
        
        // printf("--- Hessenberg ---\n");
        tic = wtime();
        res_H = MPFR_hessenberg(tmp2);
        toc = wtime();
        free_MPFR_matrix(res_H);
        time_H += toc - tic;
        // print_matrix(res_H);

        free_MPFR_matrix(A);
        free_MPFR_matrix(tmp);
        free_MPFR_matrix(tmp2);
    }



    printf("##### Temps pour %d itérations #####\n", ITE);
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH);
    printf("--- Hessenberg          = %f ---\n", time_H);
    printf("##### Temps moyen #####\n");
    printf("--- Quasi Hessenberg    = %f ---\n", time_quasiH / ITE);
    printf("--- Hessenberg          = %f ---\n", time_H / ITE);

    printf("%s\n", (time_quasiH-time_H >0)? "Quasi Hessenberg plus rapide" : "Hessenberg plus rapide");
}
