#include "matrix.h"

/*
Initialise une matrice carree de taille n et remplie de 0
@param n : la taille de la matrice

@return une matrice de coefficient initialise a 0 de taille n*n 
*/
Matrix *init_matrix(int n){
    Matrix *mat = malloc(sizeof(Matrix));

    if(mat == NULL){
        printf("init_matrix : Erreur allocation mémoire struct Matrix !\n");
        return NULL;
    }

    mat->n = n;
    mat->m = calloc(n*n,sizeof(double));

    if(mat->m == NULL){
        printf("init_matrix : Erreur allocation mémoire Matrix !\n");
        return NULL;
    }
    return mat;
}

/*
Initialise une matrice identite de taille n*n
@param n : la taille de la matrice

@return une matrice identite de taille n*n 
*/
Matrix *init_eye(int n) {
    int i;
    
    Matrix *mat = malloc(sizeof(Matrix));

    if(mat == NULL){
        printf("init_matrix : Erreur allocation mémoire struct Matrix !\n");
        return NULL;
    }

    mat->n = n;
    mat->m = calloc(n*n,sizeof(double));

    if(mat->m == NULL){
        printf("init_matrix : Erreur allocation mémoire Matrix !\n");
        return NULL;
    }
    for (i = 0; i < n; ++i) {
        mat->m[i*n + i] = 1;
    }
    return mat;
}

/*
Libere la memoire allouee a une structure Matrix
@param mat : la matrice a liberer
*/
void free_matrix(Matrix *mat){
    free(mat->m);
    free(mat);
}

/*
Libere la memoire allouee a une structure 
@param qr : la structure qr a liberer
*/
void free_qr(QR *qr){
    free_matrix(qr->Q);
    free_matrix(qr->R);
    free(qr);
}

/*
Copie une matrice dans une nouvelle matrice
@param M : la matrice a copier

@return une copie de la matrice M
*/
Matrix *copy_matrix(Matrix *M){
    Matrix *copy = init_matrix(M->n);

    memcpy(copy->m, M->m, sizeof(double) * M->n * M->n);

    return copy;
}

/*
Affiche la matrice
@param M : la matrice a afficher
*/
void print_matrix(Matrix *mat){
    int i, j;
    for (i = 0; i < mat->n; ++i) {
        for (j = 0; j < mat->n; ++j) {
            printf("%f ", mat->m[i * mat->n + j]);
            fflush(stdout);
        }
        printf("\n");
    }
}

/*
Additionne deux matrices (A+B)
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A+B
*/
Matrix *matrix_add(Matrix *A, Matrix *B){
    int i;

    if(A->n == B->n){
        Matrix *mat_add = init_matrix(A->n);
        for (i = 0; i < A->n * A->n; ++i) {
            mat_add->m[i] = add(A->m[i], B->m[i]);
        }
        return mat_add;
    }
    else{
        printf("matrix.c/matrix_add : Les matrices n'ont pas la meme taille !\n");
        return NULL;
    }
}

/*
Soustraction de deux matrices (A-B)
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A-B
*/
Matrix *matrix_sub(Matrix *A, Matrix *B){
    int i;

    if (A->n == B->n) {
        Matrix *mat_sub = init_matrix(A->n);
        for (i = 0 ; i < A->n * A->n ; ++i) {
            mat_sub->m[i] = sub(A->m[i], B->m[i]);
        }
        return mat_sub;   
    }
    else {
        printf("Soustraction de matrice impossible, dimensions : %d != %d", A->n, B->n);
        return NULL;
    }
}


/*
Multiplication naive de deux matrices (A*B)
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A*B
*/
Matrix *matrix_mul(Matrix *A, Matrix *B) {
    int i, j, k;

    if (A->n == B->n) {
        int n = A->n;
        Matrix *mat_mult = init_matrix(A->n);
        for (i = 0 ; i < n ; i++)
            for (j = 0 ; j < n ; j++){

                for (k = 0 ; k < n ; k++) {
                    double x = mat_mult->m[i * n + j];
                    double y = A->m[k + n * i];
                    double z = B->m[k * n + j];
                    mat_mult->m[i * n + j] = add(x, mul(y, z));
                }
                // printf("res[%d %d] = %f\n", i*n, j, mat_mult->m[i*n + j]);
            }
        return mat_mult;   
    }
    else {
        printf("matrix.c/matrix_mul : Multiplication de matrice impossible, dimensions : %d != %d", A->n, B->n);
        return NULL;
    }
}

/*
Transpose une matrice A
@param A : la matrice a transposer

@return la matrice A transposee
*/
Matrix *matrix_transpose(Matrix *A){
    int i, j;

    Matrix *mat_transp = init_matrix(A->n);
    for(i = 0 ; i < A->n ; i++){
        for(j = 0 ; j < A->n ; j++){
            mat_transp->m[i * A->n + j] = A->m[j * A->n + i];
        }
    }
    return mat_transp;
}

/* ***** MPFR ***** */

/*
Initialise une matrice carree de taille n et remplie de 0
En utilisant la librairie MPFR
@param n : la taille de la matrice

@return une matrice de coefficient initialise a 0 de taille n*n 
*/
MPFR_Matrix *init_MPFR_matrix(int n){
    MPFR_Matrix *mat = malloc(sizeof(MPFR_Matrix));

    if(mat == NULL){
        printf("init_MPFR_matrix : Erreur allocation mémoire struct MPFR_Matrix !\n");
        return NULL;
    }

    mat->n = n;
    mat->m = malloc(n*n * sizeof(mpfr_t));

    if(mat->m == NULL){
        printf("init_MPFR_matrix : Erreur allocation mémoire MPFR_Matrix !\n");
        return NULL;
    }

    // Initialisation des elements de la matrice
    for (int i = 0; i < n*n; i++) {
        mpinit(mat->m[i]);
        mpset(mat->m[i], 0);
    }
    return mat;
}

/*
Initialise une matrice identite de taille n*n
En utilisant la librairie MPFR
@param n : la taille de la matrice

@return une matrice identite de taille n*n 
*/
MPFR_Matrix *init_MPFR_eye(int n){
    MPFR_Matrix *mat = init_MPFR_matrix(n);

    for (int i = 0; i < n; i++) 
        mpset(mat->m[i*n + i], 1);
    
    return mat;
}

/*
Libere la memoire allouee a une structure Matrix
En utilisant la librairie MPFR
@param mat : la matrice a liberer
*/
void free_MPFR_matrix(MPFR_Matrix *mat){
    for (int i = 0; i < mat->n * mat->n; i++) {
        mpfr_clear(mat->m[i]);
    }
    free(mat->m);
    free(mat);
}

/*
Libere la memoire allouee a une structure 
@param qr : la structure qr a liberer
*/
void free_MPFR_qr(MPFR_QR *qr){
    free_MPFR_matrix(qr->Q);
    free_MPFR_matrix(qr->R);
    free(qr);
}

/*
Affiche la matrice
En utilisant la librairie MPFR
@param M : la matrice a afficher
*/
void print_MPFR_matrix(MPFR_Matrix *mat){
    int i, j;
    for (i = 0; i < mat->n; ++i) {
        for (j = 0; j < mat->n; ++j) {
            // mpfr_out_str(stdout, 10, 0, mat->m[i * mat->n + j], MPFR_RNDD);
            // printf(" ");
            mpfr_printf("%.10Rf ", mat->m[i * mat->n + j]); // 10 digits
            fflush(stdout); // Permet de vider le buffer de sortie standard, sinon le printf n'affiche pas tout
        }
        printf("\n");
    }
    fflush(stdout);
}

/*
Copie une matrice dans une nouvelle matrice
En utilisant la librairie MPFR
@param M : la matrice a copier

@return une copie de la matrice M
*/
MPFR_Matrix *copy_MPFR_matrix(MPFR_Matrix *M){
    MPFR_Matrix *copy = init_MPFR_matrix(M->n);

    for (int i = 0; i < M->n * M->n; i++) {
        mpfr_set(copy->m[i], M->m[i], MPFR_RNDD);
    }

    return copy;
}

/*
Additionne deux matrices (A+B)
En utilisant la librairie MPFR
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A+B
*/
MPFR_Matrix *MPFR_matrix_add(MPFR_Matrix *A, MPFR_Matrix *B){
    int i;

    if(A->n == B->n){
        MPFR_Matrix *mat_add = init_MPFR_matrix(A->n);
        for (i = 0; i < A->n * A->n; ++i) {
            mpadd(mat_add->m[i], A->m[i], B->m[i]);
        }
        return mat_add;
    }
    else{
        printf("matrix.c/MPFR_matrix_add : Les matrices n'ont pas la meme taille !\n");
        return NULL;
    }
}

/*
Soustraction de deux matrices (A-B)
En utilisant la librairie MPFR
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A-B
*/
MPFR_Matrix *MPFR_matrix_sub(MPFR_Matrix *A, MPFR_Matrix *B){
    int i;

    if(A->n == B->n){
        MPFR_Matrix *mat_sub = init_MPFR_matrix(A->n);
        for (i = 0; i < A->n * A->n; ++i) {
            mpsub(mat_sub->m[i], A->m[i], B->m[i]);
        }
        return mat_sub;
    }
    else{
        printf("matrix.c/MPFR_matrix_sub : Les matrices n'ont pas la meme taille !\n");
        return NULL;
    }
}

/*
Multiplication de deux matrices (A*B)
En utilisant la librairie MPFR
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A*B
*/
MPFR_Matrix *MPFR_matrix_mul(MPFR_Matrix *A, MPFR_Matrix *B){
    int i, j, k;

    if (A->n == B->n) {
        int n = A->n;
        MPFR_Matrix *mat_mult = init_MPFR_matrix(A->n);
        for (i = 0 ; i < n ; i++)
            for (j = 0 ; j < n ; j++){

                for (k = 0 ; k < n ; k++) {
                    mpfr_t x, y, z;
                    mpinit(x);
                    mpinit(y);
                    mpinit(z);

                    mpfr_set(x, mat_mult->m[i * n + j], MPFR_RNDD);
                    mpfr_set(y, A->m[k + n * i], MPFR_RNDD);
                    mpfr_set(z, B->m[k * n + j], MPFR_RNDD);

                    mpmul(y, y, z);
                    mpadd(mat_mult->m[i * n + j], x, y);

                    mpfr_clear(x);
                    mpfr_clear(y);
                    mpfr_clear(z);
                }
            }
        return mat_mult;
    }
    else {
        printf("matrix.c/MPFR_matrix_mul : Multiplication de matrice impossible, dimensions : %d != %d", A->n, B->n);
        return NULL;
    }                
}

/*
Transpose une matrice A
@param A : la matrice a transposer

@return la matrice A transposee
*/
MPFR_Matrix *MPFR_matrix_transpose(MPFR_Matrix *A){
    int i, j;

    MPFR_Matrix *mat_transp = init_MPFR_matrix(A->n);
    for(i = 0 ; i < A->n ; i++){
        for(j = 0 ; j < A->n ; j++){
            mpfr_set(mat_transp->m[i * A->n + j], A->m[j * A->n + i], MPFR_RNDD);
        }
    }
    return mat_transp;
}