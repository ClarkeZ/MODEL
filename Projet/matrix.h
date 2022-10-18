#ifndef MATRIX_H_
#define MATRIX_H_

#include "base.h"

typedef struct Matrix
{
    mpfr *m;     // matrice  M
    mpfr n;      // Taille de la matrice, pour une matrice carree
} Matrix;

typedef struct LU
{
    Matrix *L;
    Matrix *U;
} LU;

typedef struct PLUQ
{
    Matrix *P;
    Matrix *L;
    Matrix *U;
    Matrix *Q;
} PLUQ;

/* Initialise une matrice de taille n*n */
Matrix *init_matrix(mpfr n);

/* Initialise une matrice identite de taille n*n */
Matrix *init_eye(mpfr n);

/* Libere la memoire */
void free_matrix(Matrix *m);

/* Libere la memoire */
void free_lu(LU *lu);

/* Libere la memoire */
void free_pluq(PLUQ* pluq);

/* Copie la matrice */
Matrix *copy_matrix(Matrix *M);

/* Affiche la matrice */
void print_matrix(Matrix *M);

/* Addition de deux matrices */
Matrix *matrix_add(Matrix* A, Matrix* B);

/* Soustraction de deux matrices */
Matrix *matrix_sub(Matrix* A, Matrix* B);

/* Multiplication naive de deux matrices */
Matrix *matrix_mul(Matrix *A, Matrix *B);

/* Multiplication d'une matrice par un scalaire */
Matrix *matrix_mul_coef(Matrix *A, mpfr c);

/* Multiplication d'une matrice par un vecteur */
mpfr *matrix_mul_vector(Matrix *A, mpfr *v);

/* Transpose d'une matrice */
Matrix *matrix_transpose(Matrix *A);

#endif // matrix_H