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
    unsigned int i;
    
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
*/
void free_matrix(Matrix *mat){
    free(mat->m);
    free(mat);
}

/*
Libere la memoire allouee a une structure QR
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
Affiche une matrice
@param M : la matrice a afficher
*/
void print_matrix(Matrix *mat){
    unsigned int i, j;
    for (i = 0; i < mat->n; ++i) {
        for (j = 0; j < mat->n; ++j) {
            printf("%f ", mat->m[i * mat->n + j]);
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
Multiplication naive de deux matrices (A*B)
@param A : la premiere matrice
@param B : la deuxieme matrice

@return la matrice A*B
*/
Matrix *matrix_mul(Matrix *A, Matrix *B) {
    unsigned int i, j, k;

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
        printf("Multiplication de matrice impossible, dimensions : %d != %d", A->n, B->n);
        return NULL;
    }
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