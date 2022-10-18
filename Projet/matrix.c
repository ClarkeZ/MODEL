#include "matrix.h"


/*
Initialise une matrice carree de taille n et remplie de 0
@param n : la taille de la matrice

@return une matrice de coefficient initialise a 0 de taille n*n 
*/
*/
Matrix *init_matrix(mpfr n){
    Matrix *mat = malloc(sizeof(Matrix));

    if(mat == NULL){
        printf("init_matrix : Erreur allocation mémoire struct Matrix !\n");
        return NULL;
    }

    mat->n = n;
    mat->m = calloc(n*n,sizeof(mpfr));

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
*/
Matrix *init_eye(mpfr n) {
    unsigned int i;
    
    Matrix *mat = malloc(sizeof(Matrix));

    if(mat == NULL){
        printf("init_matrix : Erreur allocation mémoire struct Matrix !\n");
        return NULL;
    }

    mat->n = n;
    mat->m = calloc(n*n,sizeof(mpfr));

    if(mat->m == NULL){
        printf("init_matrix : Erreur allocation mémoire Matrix !\n");
        return NULL;
    }
    for (i = 0; i < n; ++i) {
        mat->m[i*n+i] = 1;
    }
    return mat;
}

/*
Libere la memoire allouee a une structure Matrix
*/
void free_matrix(Matrix *mat){
    free(mat->m);
    free(mat);
}

/*
Libere la memoire allouee a une structure LU
*/
void free_lu(LU *lu){
    free_matrix(lu->L);
    free_matrix(lu->U);
    free(lu);
}

/*
Libere la memoire allouee a une structure PLUQ
*/
void free_pluq(PLUQ* pluq){
    free_matrix(pluq->P);
    free_matrix(pluq->L);
    free_matrix(pluq->U);
    free_matrix(pluq->Q);
    free(pluq);
}

/*
Copie une matrice dans une nouvelle matrice
@param M : la matrice a copier

@return une copie de la matrice M
*/
Matrix *copy_matrix(Matrix *M){
    Matrix *copy = init_matrix(M->n);

    memcpy(copy->m, M->m, sizeof(mpfr) * M->n * M->n);

    return copy;
}

/*
Affiche une matrice
@param M : la matrice a afficher
*/
void print_matrix(Matrix *M){
    unsigned int i, j;
    for (i = 0; i < M->n; ++i) {
        for (j = 0; j < M->n; ++j) {
            printf("%d ", M->m[i*M->n+j]);
            // mpfr_out_str (stdout, 10, 0, M->m[i*M->n+j], MPFR_RNDD);
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
    unsigned int i;

    if(A->n == B->n){
        Matrix *mat_add = init_matrix(A->n);
        for (i = 0; i < A->n * A->n; ++i) {
            mat_add->m[i] = mpadd(A->m[i], B->m[i]);
        }
        return mat_add;
    }
    else{
        printf("matrix.c/matrix_add : Les matrices n'ont pas la meme taille !\n");
        return NULL;
    }
}

/*
Soustrait deux matrices (A-B)
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A-B
*/
Matrix *matrix_sub(Matrix *A, Matrix *B){
    unsigned int i;

    if(A->n == B->n){
        Matrix *mat_sub = init_matrix(A->n);
        for (i = 0; i < A->n * A->n; ++i) {
            mat_sub->m[i] = mpsub(A->m[i], B->m[i]);
        }
        return mat_sub;
    }
    else{
        printf("matrix.c/matrix_sub : Les matrices n'ont pas la meme taille !\n");
        return NULL;
    }
}

/*
Multiplication naive de deux matrices (A*B)
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A*B
*/
Matrix *matrix_mul(Matrix *A, Matrix *B){
    unsigned int i, j, k;

    if(A->n == B->n){
        mpfr n = A->n;
        Matrix *mat_mul = init_matrix(A->n);
        for (i = 0; i < A->n; ++i) {
            for (j = 0; j < A->n; ++j) {
                for (k = 0; k < A->n; ++k) {
                    mpfr x = mat_mul->m[i * n + j];
                    mpfr y = A->m[i * n + k];
                    mpfr z = B->m[k * n + j];
                    mat_mul->m[i * n + j] = mpadd(x, mpmul(y, z));
                }
            }
        }
        return mat_mul;
    }
    else{
        printf("matrix.c/matrix_mul : Les matrices n'ont pas la meme taille !\n");
        return NULL;
    }
}

/*
Multiplication d'une matrice par un scalaire (A*c)
@param A : la matrice
@param c : le scalaire

@return la matrice A*c
*/
Matrix *matrix_mul_scalar(Matrix *A, mpfr scalar){ 
    unsigned int i, j;
    Matrix *mat_mul = init_matrix(A->n);
    
    for (i = 0; i < A->n ; ++i) {
        for(j = 0; j < A->n; ++j){
            mat_mul->m[i * A->n + j] = mpmul(A->m[i * A->n + j], scalar);
        }
    }
    return mat_mul;
}

/*
Multiplication d'une matrice par un vecteur (A*v)
@param A : la matrice
@param v : le vecteur

@return le vecteur A*v
*/
mpfr *matrix_mul_vector(Matrix *A, mpfr *v){
    unsigned int i, j;
    mpfr *res = calloc(A->n, sizeof(mpfr));
    
    for (i = 0; i < A->n ; ++i) {
        for(j = 0; j < A->n; ++j){
            res[i] = mpadd(res[i], mpmul(A->m[i * A->n + j], v[j]));
        }
    }
    return res;
}

/*
Transpose une matrice A
@param A : la matrice a transposer

@return la matrice A transposee
*/
Matrix *matrix_transpose(Matrix *A){
    unsigned int i, j;

    Matrix *mat_transp = init_matrix(A->n);
    for(i = 0 ; i < A->n ; i++){
        for(j = 0 ; j < A->n ; j++){
            mat_transp->m[i * A->n + j] = A->m[j * A->n + i];
        }
    }
    return mat_transp;
}