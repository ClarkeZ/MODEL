#ifndef BASE_H_
#define BASE_H_

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

// GNU MPFR Library
#include <gmp.h>
#include <mpfr.h>

#define MIN 0
#define MAX 10

typedef mpfr_t mpfr;

double add(double a, double b);

double sub(double a, double b);

double mul(double a, double b);

void mpinit(mpfr a);

void mpset(mpfr a, double b);

/* Addition de deux entiers */
void mpadd(mpfr res, mpfr a, mpfr b);

/* Soustraction de deux entiers */
void mpsub(mpfr res, mpfr a, mpfr b);

/* Multiplication de deux entiers */
void mpmul(mpfr res, mpfr a, mpfr b);

#endif /* BASE_H_ */