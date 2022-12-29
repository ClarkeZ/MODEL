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

double find_cos_matrix(int i, int j, Matrix *A);

double find_sin(int i, int j, Matrix *A);

double find_sin_matrix(int i, int j, Matrix *A);

Matrix *givens(int i, int j, Matrix *A);

void givens_matrix(int i, int j, Matrix *G, Matrix *A);

QR *qr_decomposition(Matrix *A);

Matrix *quasi_hess(Matrix *A);

Matrix *hessenberg(Matrix *A);





