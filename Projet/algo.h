#ifndef ALGO_H_
#define ALGO_H_

#include <stdbool.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include "matrix.h"

#define NB_ITE 30

double wtime();

double randfrom(double min, double max);

double func_sqrt(double n);

/* ***** ALGO ***** */

double find_cos(int i, int j, Matrix *A);

double find_cos_matrix(int i, int j, Matrix *A);

double find_sin(int i, int j, Matrix *A);

double find_sin_matrix(int i, int j, Matrix *A);

Matrix *givens(int i, int j, Matrix *A);

void givens_matrix(int i, int j, Matrix *G, Matrix *A);

QR *qr_decomposition(Matrix *A);

Matrix *quasi_hess(Matrix *A);

Matrix *hessenberg(Matrix *A);

double *eigenvalues(Matrix *A, int k);

/* ***** MPFR ***** */

void MPFR_find_cos(int i, int j, MPFR_Matrix *A, mpfr cos);

void MPFR_find_cos_matrix(int i, int j, MPFR_Matrix *A, mpfr cos);

void MPFR_find_sin(int i, int j, MPFR_Matrix *A, mpfr sin);

void MPFR_find_sin_matrix(int i, int j, MPFR_Matrix *A, mpfr sin);

MPFR_Matrix *MPFR_givens(int i, int j, MPFR_Matrix *A);

void MPFR_givens_matrix(int i, int j, MPFR_Matrix *G, MPFR_Matrix *A);

MPFR_QR *MPFR_qr_decomposition(MPFR_Matrix *A);

MPFR_Matrix *MPFR_quasi_hess(MPFR_Matrix *A);

MPFR_Matrix *MPFR_hessenberg(MPFR_Matrix *A);

void MPFR_eigenvalues(MPFR_Matrix *A, mpfr *eigen, int k);

#endif /* ALGO_H_ */