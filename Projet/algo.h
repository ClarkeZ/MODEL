#ifndef ALGO_H_
#define ALGO_H_

#include <stdbool.h>
#include <sys/time.h>
#include <time.h>

#include "matrix.h"

double wtime();

/* Multiplication d'une ligne de matrice par un coefficient */
void mult_line_scalar(mpfr *line, mpfr coef, mpfr n);

/* Soustraction d'une ligne de matrice par une autre ligne */
void sub_lines(mpfr *line1, mpfr *line2, mpfr n);

/* Verifie si la diagonale possede un coefficient nul */
bool zero_in_diagonal(Matrix *A);

/* Decomposition LU */
LU *lu(Matrix *A);

/* Decomposition PLUQ */
PLUQ *pluq(Matrix *A);