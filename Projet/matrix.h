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

Matrix *matrix_sub(Matrix *A, Matrix *B);

Matrix *matrix_transpose(Matrix *A);

Matrix *matrix_inverse(Matrix *A);

Matrix *matrix_mul_coef(Matrix *A, double c);


/* ===== Fonctions MPFR ===== */

// typedef struct MPFR_Matrix
// {
//     mpfr *m;     // matrice  M
//     int n;      // Taille de la matrice, pour une matrice carree
// } MPFR_Matrix;

// typedef struct MPFR_QR
// {
//     MPFR_Matrix *Q;
//     MPFR_Matrix *R;
// } MPFR_QR;

/* Initialise une matrice de taille n*n */
// MPFR_Matrix *init_MPFR_matrix(int n);

/* Initialise une matrice identite de taille n*n */
// MPFR_Matrix *init_MPFR_eye(int n);

/* Libere la memoire */
// void free_MPFR_matrix(Matrix *m);

/* Libere la memoire */
// void free_MPFR_qr(QR *qr);

/* Copie la matrice */
// MPFR_Matrix *copy_MPFR_matrix(Matrix *M);

/* Affiche la matrice */
// void print_MPFR_matrix(Matrix *M);

/* Addition de deux matrices */
// MPFR_Matrix *matrix_MPFR_add(Matrix* A, Matrix* B);

/* Soustraction de deux matrices */
// MPFR_Matrix *matrix_MPFR_sub(Matrix* A, Matrix* B);

/* Multiplication naive de deux matrices */
// MPFR_Matrix *matrix_MPFR_mul(Matrix *A, Matrix *B);

/* Multiplication d'une matrice par un scalaire */
// MPFR_Matrix *matrix_mul_MPFR_coef(Matrix *A, mpfr c);

/* Multiplication d'une matrice par un vecteur */
// mpfr *matrix_mul_MPFR_vector(Matrix *A, mpfr *v);

#endif // matrix_H