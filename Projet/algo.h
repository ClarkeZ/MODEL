// #ifndef ALGO_H_
#define ALGO_H_

#include <stdbool.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include "matrix.h"

double wtime();

double randfrom(double min, double max);

double func_sqrt(double n);

double find_cos(int i, int j, Matrix *A);

double find_sin(int i, int j, Matrix *A);

Matrix *Givens(int i, int j, Matrix *A);

Matrix *mat_Q(Matrix *A);

QR *qr_decomposition(Matrix *A);

Matrix *quasi_hess(Matrix *A);

/* Multiplication d'une ligne de matrice par un coefficient */
// void mult_line_scalar(mpfr *line, mpfr coef, mpfr n);

/* Soustraction d'une ligne de matrice par une autre ligne */
// void sub_lines(mpfr *line1, mpfr *line2, mpfr n);

/* Verifie si la diagonale possede un coefficient nul */
bool zero_in_diagonal(Matrix *A);


