#include "base.h"

#define PRECISION 200

double add(double a, double b) {
    return a + b;
}

double sub(double a, double b) {
    return a - b;
}

double mul(double a, double b) {
    return a * b;
}

/*
Initialisation d'un nombre flottant de grande precision
@param a : le nombre a initialiser
*/
void mpinit(mpfr a){
    mpfr_init2(a, PRECISION);
}

/*
Set d'un nombre flottant de grande precision
@param a : le nombre à initialiser
@param b : le nombre à mettre dans a
*/
void mpset(mpfr a, double b){
    mpfr_set_d(a, b, MPFR_RNDD);
}

/*
Addition de deux nombres flottants de grande precision 
@param res : resultat de l'addition
@param a : le premier entier
@param b : le deuxieme entier
*/
void mpadd(mpfr res, mpfr a, mpfr b){
    mpfr_add(res, a, b, MPFR_RNDN);
}

/*
Soustraction de deux nombres flottants de grande precision
@param res : resultat de la soustraction
@param a : le premier entier
@param b : le deuxieme entier
*/
void mpsub(mpfr res, mpfr a, mpfr b){
    mpfr_sub(res, a, b, MPFR_RNDN);
}

/*
Multiplication de deux nombres flottants de grande precision
@param res : resultat de la multiplication
@param a : le premier entier
@param b : le deuxieme entier
*/
void mpmul(mpfr res, mpfr a, mpfr b){
    mpfr_mul(res, a, b, MPFR_RNDN);
}