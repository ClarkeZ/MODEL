/*
Tout les calculs sont faites avec la librairie MPFR (Multiple Precision Floating-point Reliable) 
qui permet de faire des calculs avec des nombres flottants de grande precision.
*/

#ifndef BASE_H_
#define BASE_H_

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

// GNU MPFR Library
#include <gmp.h>
#include <mpfr.h>

typedef mpfr_t mpfr;
// typedef uint64_t u64;

double add(double a, double b);

double sub(double a, double b);

double mul(double a, double b);

// u64 div(u64 a, u64 b);

// void mpinit(mpfr a);

// void mpset(mpfr a, double b);

/* Addition de deux entiers */
// void mpadd(mpfr res, mpfr a, mpfr b);

/* Soustraction de deux entiers */
// void mpsub(mpfr res, mpfr a, mpfr b);

/* Multiplication de deux entiers */
// void mpmul(mpfr res, mpfr a, mpfr b);

/* Division de deux entiers */
// void mpdiv(mpfr res, mpfr a, mpfr b);

#endif /* BASE_H_ */