#ifndef MATRIX_H_
#define MATRIX_H_

#include "base.h"

typedef struct Matrix {
    double *m;
    int n;
} Matrix;

typedef struct QR {
    Matrix *Q;
    Matrix *R;
} QR;

Matrix *init_matrix(int n);

Matrix *init_eye(int n);

void free_matrix(Matrix *mat);

void free_qr(QR *qr);

void print_matrix(Matrix *mat);

Matrix *copy_matrix(Matrix *M);

Matrix *matrix_mul(Matrix *A, Matrix *B);

Matrix *matrix_add(Matrix *A, Matrix *B);

Matrix *matrix_transpose(Matrix *A);


/* ===== Fonctions MPFR ===== */

typedef struct MPFR_Matrix
{
    mpfr *m;     // matrice  M
    int n;      // Taille de la matrice, pour une matrice carree
} MPFR_Matrix;

typedef struct MPFR_QR
{
    MPFR_Matrix *Q;
    MPFR_Matrix *R;
} MPFR_QR;

/* Initialise une matrice de taille n*n */
MPFR_Matrix *init_MPFR_matrix(int n);

/* Initialise une matrice identite de taille n*n */
MPFR_Matrix *init_MPFR_eye(int n);

/* Libere la memoire */
void free_MPFR_matrix(MPFR_Matrix *m);

/* Libere la memoire */
void free_MPFR_qr(MPFR_QR *qr);

/* Affiche la matrice */
void print_MPFR_matrix(MPFR_Matrix *M);

/* Copie la matrice */
MPFR_Matrix *copy_MPFR_matrix(MPFR_Matrix *M);

/* Addition de deux matrices */
MPFR_Matrix *MPFR_matrix_add(MPFR_Matrix* A, MPFR_Matrix* B);

/* Multiplication naive de deux matrices */
MPFR_Matrix *MPFR_matrix_mul(MPFR_Matrix *A, MPFR_Matrix *B);

/* Transpose d'une matrice */
MPFR_Matrix *MPFR_matrix_transpose(MPFR_Matrix *A);

#endif // matrix_H