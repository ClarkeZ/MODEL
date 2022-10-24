#ifndef MATRIX_H_
#define MATRIX_H_

#include "base.h"

typedef struct Matrix {
    u64 *m;
    u64 n;
} Matrix;

typedef struct QR {
    Matrix *Q;
    Matrix *R;
} QR;

Matrix *init_matrix(u64 n);

Matrix *init_eye(u64 n);

void free_matrix(Matrix *mat);

void free_qr(QR *qr);

void print_matrix(Matrix *mat);

Matrix *copy_matrix(Matrix *M);

Matrix *matrix_mul(Matrix *A, Matrix *B);

Matrix *mult_matrix_scalar(Matrix *A, mpfr scalar);

Matrix *matrix_add(Matrix *A, Matrix *B);

Matrix *matrix_sub(Matrix *A, Matrix *B);

Matrix *matrix_transpose(Matrix *A);

Matrix *matrix_inverse(Matrix *A);

Matrix *matrix_mul_coef(Matrix *A, u64 c);


/* ===== Fonctions MPFR ===== */

typedef struct MPFR_Matrix
{
    mpfr *m;     // matrice  M
    mpfr n;      // Taille de la matrice, pour une matrice carree
} MPFR_Matrix;

typedef struct MPFR_QR
{
    MPFR_Matrix *Q;
    MPFR_Matrix *R;
} MPFR_QR;

/* Initialise une matrice de taille n*n */
MPFR_Matrix *init_MPFR_matrix(mpfr n);

/* Initialise une matrice identite de taille n*n */
MPFR_Matrix *init_MPFR_eye(mpfr n);

/* Libere la memoire */
void free_MPFR_matrix(Matrix *m);

/* Libere la memoire */
void free_MPFR_qr(QR *qr);

/* Copie la matrice */
MPFR_Matrix *copy_MPFR_matrix(Matrix *M);

/* Affiche la matrice */
void print_MPFR_matrix(Matrix *M);

/* Addition de deux matrices */
MPFR_Matrix *matrix_MPFR_add(Matrix* A, Matrix* B);

/* Soustraction de deux matrices */
MPFR_Matrix *matrix_MPFR_sub(Matrix* A, Matrix* B);

/* Multiplication naive de deux matrices */
MPFR_Matrix *matrix_MPFR_mul(Matrix *A, Matrix *B);

/* Multiplication d'une matrice par un scalaire */
MPFR_Matrix *matrix_mul_MPFR_coef(Matrix *A, mpfr c);

/* Multiplication d'une matrice par un vecteur */
mpfr *matrix_mul_MPFR_vector(Matrix *A, mpfr *v);

#endif // matrix_H