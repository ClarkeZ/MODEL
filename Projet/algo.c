#include "algo.h"

double wtime(){
    struct timeval ts;
    gettimeofday(&ts, NULL);
    return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

/* generate a random floating point number from min to max */
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

/*
Calcul la racine carrée d'un nombre n
@param n : un double

@return un double qui est la racine carrée de n
*/
double func_sqrt(double n){
    double sqrt = n / 2;
    double tmp = 0;

    while(sqrt != tmp){
        tmp = sqrt;
        sqrt = (n / tmp + tmp) / 2;
    }

    return sqrt;
}

/*
Trouve le cosinus pour la matrice de Givens Gij de la matrice A de taille n*n 
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le cosinus de A pour la matrice de Givens Gij
*/
double find_cos(int i, int j, Matrix *A){
    // printf("Aij = %f\n",A->m[i*A->n + j]);
    // printf("Ajj = %f\n",A->m[j*A->n + j]);
    if(A->m[i*A->n + j] == 0 && A->m[j*A->n + j] == 0)
        return 1;
    else 
        return A->m[j*A->n + j] / sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
        // si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[j*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

/*
Trouve le cosinus pour la matrice de Givens Gij pour la fonction Hessenberg
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le cosinus de A pour la matrice de Givens Gij
*/
double find_cos_matrix(int i, int j, Matrix *A){
    if(A->m[(i-1)*A->n + j] == 0 && A->m[i*A->n + j] == 0)
        return 1;
    else 
        return A->m[(i-1)*A->n + j] / sqrt(A->m[(i-1)*A->n + j]*A->m[(i-1)*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
        // si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[j*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

/*
Trouve le sinus pour la matrice de Givens Gij de la matrice A de taille n*n 
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le sinus de A pour la matrice de Givens Gij
*/
double find_sin(int i, int j, Matrix *A){
    // printf("Aij = %f\n",A->m[i*A->n + j]);
    // printf("Ajj = %f\n",A->m[j*A->n + j]);
    if(A->m[i*A->n + j] == 0 && A->m[j*A->n + j] == 0)
        return 1;
    else   
        return A->m[i*A->n + j] / sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
        //si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[i*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

/*
Trouve le cosinus pour la matrice de Givens Gij pour la fonction Hessenberg
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le cosinus de A pour la matrice de Givens Gij
*/
double find_sin_matrix(int i, int j, Matrix *A){
    if(A->m[i*A->n + j] == 0 && A->m[(i-1)*A->n + j] == 0)
        return 1;
    else   
        return A->m[i*A->n + j] / sqrt(A->m[(i-1)*A->n + j]*A->m[(i-1)*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
        //si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[i*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

/*
Calcul la matrice de Givens Gij de la matrice A de taille n*n 
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return une matrice de taille n*n qui est la matrice de Givens Gij
*/
Matrix *givens(int i, int j, Matrix *A){
    Matrix *G = init_eye(A->n);

    double c = find_cos(i, j, A);
    double s = find_sin(i, j, A);
    // printf("c = %f\n", c);
    // printf("s = %f\n", s);

    // Verification de la formule de Givens (c^2 + s^2 = 1) 
    // if(c*c + s*s != 1){ 
    //     printf("Error: c^2 + s^2 != 1\n");
    //     // exit(EXIT_FAILURE);
    // }    

    G->m[j*A->n + j] = c; // En haut a gauche
    G->m[j*A->n + i] = s; // En haut a droite
    G->m[i*A->n + j] = -s; // En bas a gauche
    G->m[i*A->n + i] = c; // En bas a droite

    // printf("Givens(%d, %d)\n", i+1, j+1);
    // print_matrix(G);
    return G;
}

/*
Calcul la matrice de Givens Gij de la matrice A de taille n*n pour la fonction Hessenberg
@param i : un entier
@param j : un entier
@param *G : la matrice G, qui est la matrice de Givens Gij
@param *A : la matrice A
*/
void givens_matrix(int i, int j, Matrix *G, Matrix *A){
    double c = find_cos_matrix(i, j, A);
    double s = find_sin_matrix(i, j, A);
    // printf("c = %f\n", c);
    // printf("s = %f\n", s);

    // Verification de la formule de Givens (c^2 + s^2 = 1) 
    // if(c*c + s*s != 1){ 
    //     printf("Error: c^2 + s^2 != 1\n");
    //     // exit(EXIT_FAILURE);
    // }    

    G->m[(i-1)*A->n + (i-1)] = c; // En haut a gauche
    G->m[(i-1)*A->n + i] = s; // En haut a droite
    G->m[i*A->n + (i-1)] = -s; // En bas a gauche
    G->m[i*A->n + i] = c; // En bas a droite

    // printf("Givens(%d, %d)\n", i+1, j+1);
    // print_matrix(G);
}

/*
QR decomposition de la matrice A de taille n*n 
@param *A : la matrice A

@return une structure de type QR qui contient la matrice Q et la matrice R
*/
QR *qr_decomposition(Matrix *A){
    QR *qr = (QR *) malloc(sizeof(QR));
    qr->Q = init_eye(A->n);
    qr->R = copy_matrix(A);

    for(int j = 0 ; j < A->n ; j++){ // From j = 1 to n
        for(int i = j + 1 ; i < A->n ; i++){ // From i = j+1 to m
            Matrix *G = givens(i, j, qr->R);

            qr->R = matrix_mul(G, qr->R);
            Matrix *Gt = matrix_transpose(G);
            qr->Q = matrix_mul(qr->Q, Gt);

            free_matrix(G);
            free_matrix(Gt);
        }
    }
    return qr;
}

/*
Calcul une matrice quasi Hessenberg de la matrice A de taille n*n
@param *A : la matrice A

@return une matrice de taille n*n qui est une matrice quasi Hessenberg
*/
Matrix *quasi_hess(Matrix *A){
    int i = 0;
    unsigned int cpt = 0;
    
    // Nombre d'iterations max pour la boucle while (pour eviter les boucles infinies)
    unsigned int it = 0;
    
    double threshold = 1e-4;
    // double threshold = 0.5;
    
    Matrix *A_prime = copy_matrix(A);
    QR *qr = qr_decomposition(A_prime);

    // Recupere les coefficients de la subdiagonale inferieure
    double lower_diag[A->n - 1];
    for(i = 0 ; i < A->n - 1 ; i++){
        lower_diag[i] = A->m[(i + 1)*A->n + i];
    }

    while(it != NB_ITE){
        // Juste pour afficher si le nombre max d'iteration est atteint
        // if(it == nb_it-1) 
        //     printf("Trop d'iterations\n");
        
        for(i = 0 ; i < A->n -1 ; i++){
            if(lower_diag[i] > threshold)
                cpt++;
        }

        if(cpt == 0)
            break;

        Matrix *Qt = matrix_transpose(qr->Q);
        Matrix *tmp = matrix_mul(Qt, A_prime);

        A_prime = matrix_mul(tmp, qr->Q);

        free_qr(qr);
        free_matrix(Qt);
        free_matrix(tmp);

        qr = qr_decomposition(A_prime);

        // Recupere les coefficients de la subdiagonale inferieure
        for(i = 0 ; i < A->n - 1 ; i++)
            lower_diag[i] = A_prime->m[(i+1)*A->n + i];
        
        it++;
    }

    free_qr(qr);
    return A_prime;
}

/*
Calcul la matrice de Hessenberg de la matrice A de taille n*n
@param *A : la matrice A

@return une matrice de taille n*n qui est une matrice de Hessenberg
*/
Matrix *hessenberg(Matrix *A){
    int i = 0;
    int j = 0;
    int n = A->n;
    Matrix *hess = copy_matrix(A);
    Matrix *G = init_eye(n);
    Matrix *Gt = init_eye(n);

    // Parcours les coeff sous la sub diagonal
    for(j = 0 ; j < n - 2 ; j++){
        for(i = n-1 ; i > j + 1 ; i--){
            if(hess->m[i*n + j] != 0){
                G = init_eye(n);
                givens_matrix(i, j, G, hess);
                Gt = matrix_transpose(G);

                hess = matrix_mul(G, hess);
                hess = matrix_mul(hess, Gt);

                free_matrix(G);
                free_matrix(Gt);
            }
        }
    }
    return hess;
}

/*
Calcul les valeurs propres de la matrice A de taille n*n avec la methode QR 
@param *A : la matrice A
@param k : nombre d'iterations

@return un tableau de taille n qui contient les valeurs propres de la matrice A
*/
double *eigenvalues(Matrix *A, int k){
    double *eigen = (double *) malloc(A->n * sizeof(double));

    QR *qr;
    Matrix *A_prime = copy_matrix(A);
    // Matrix *tmp;

    int i = 0;
    
    // k Nombre d'iterations, plus k est grand, plus la precision est bonne
    // Mais plus la matrice est grande, plus k doit etre grand (encore plus que pour une petite matrice)
    // Car on fait k fois la decomposition QR sur la matrice A
    // Donc plus A, a de chance de converger vers une matrice diagonale 
    while(i < k){
        qr = qr_decomposition(A_prime);
        // Matrix *Qt = matrix_transpose(qr->Q);
        // tmp = matrix_mul(A_prime, qr->Q);
        // A_prime = matrix_mul(Qt, tmp);
        // free_matrix(Qt);
        // free_matrix(tmp);

        // Faire RQ revient a faire Q*AR
        A_prime = matrix_mul(qr->R, qr->Q);
        free_qr(qr);
        i++;
    }

    for(int i = 0 ; i < A->n ; i++){
        eigen[i] = A_prime->m[i*A->n + i];
    }

    return eigen;
}


/* ***** MPFR ***** */

/*
Trouve le cosinus pour la matrice de Givens Gij de la matrice A de taille n*n 
En utilisant la librairie MPFR
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le cosinus de A pour la matrice de Givens Gij
*/
void MPFR_find_cos(int i, int j, MPFR_Matrix *A, mpfr cos){
    if(mpfr_cmp_d(A->m[i*A->n + j], 0) == 0 && mpfr_cmp_d(A->m[j*A->n + j], 0) == 0)
        mpset(cos, 1);
    else {
        mpfr tmp, tmp2, tmp3, tmp4;
        mpinit(tmp);
        mpinit(tmp2);
        mpinit(tmp3);
        mpinit(tmp4);

        mpfr_mul(tmp, A->m[i*A->n + j], A->m[i*A->n + j], MPFR_RNDN);
        mpfr_mul(tmp2, A->m[j*A->n + j], A->m[j*A->n + j], MPFR_RNDN);
        mpfr_add(tmp3, tmp, tmp2, MPFR_RNDN);

        mpfr_sqrt(tmp4, tmp3, MPFR_RNDN);

        mpfr_div(cos, A->m[j*A->n + j], tmp4, MPFR_RNDN);
        
        mpfr_clear(tmp);
        mpfr_clear(tmp2);
        mpfr_clear(tmp3);
        mpfr_clear(tmp4);
    }
}

/*
Trouve le cosinus pour la matrice de Givens Gij pour la fonction Hessenberg
En utilisant la librairie MPFR
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le cosinus de A pour la matrice de Givens Gij
*/
void MPFR_find_cos_matrix(int i, int j, MPFR_Matrix *A, mpfr cos){
    if(mpfr_cmp_d(A->m[(i-1)*A->n + j], 0) == 0 && mpfr_cmp_d(A->m[i*A->n + j], 0) == 0)
        mpset(cos, 1);
    else {
        mpfr a, b, sqrt;
        mpinit(a);
        mpinit(b);
        mpinit(sqrt);

        mpfr_set(a, A->m[(i-1)*A->n + j], MPFR_RNDN);
        mpfr_set(b, A->m[i*A->n + j], MPFR_RNDN);
        mpfr_mul(sqrt, a, a, MPFR_RNDN); // sqrt = a*a
        mpfr_fma(sqrt, b, b, sqrt, MPFR_RNDN); // sqrt = b*b + (a*a) 
        mpfr_sqrt(sqrt, sqrt, MPFR_RNDN); // sqrt = sqrt((b*b + a*a))
        mpfr_div(cos, a, sqrt, MPFR_RNDN); // cos = a / sqrt((b*b + a*a))
        
        mpfr_clear(a);
        mpfr_clear(b);
        mpfr_clear(sqrt);
    }
}

/*
Trouve le sinus pour la matrice de Givens Gij de la matrice A de taille n*n 
En utilisant la librairie MPFR
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le sinus de A pour la matrice de Givens Gij
*/
void MPFR_find_sin(int i, int j, MPFR_Matrix *A, mpfr sin){
    if(mpfr_cmp_d(A->m[i*A->n + j], 0) == 0 && mpfr_cmp_d(A->m[j*A->n + j], 0) == 0)
        mpset(sin, 1);
    else {
        mpfr tmp, tmp2, tmp3, tmp4;
        mpinit(tmp);
        mpinit(tmp2);
        mpinit(tmp3);
        mpinit(tmp4);

        mpfr_mul(tmp, A->m[i*A->n + j], A->m[i*A->n + j], MPFR_RNDN);
        mpfr_mul(tmp2, A->m[j*A->n + j], A->m[j*A->n + j], MPFR_RNDN);
        mpfr_add(tmp3, tmp, tmp2, MPFR_RNDN);

        mpfr_sqrt(tmp4, tmp3, MPFR_RNDN);

        mpfr_div(sin, A->m[i*A->n + j], tmp4, MPFR_RNDN);

        mpfr_clear(tmp);
        mpfr_clear(tmp2);
        mpfr_clear(tmp3);
        mpfr_clear(tmp4);
    }
}

/*
Trouve le cosinus pour la matrice de Givens Gij pour la fonction Hessenberg
En utilisant la librairie MPFR
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return un double qui est le cosinus de A pour la matrice de Givens Gij
*/
void MPFR_find_sin_matrix(int i, int j, MPFR_Matrix *A, mpfr sin){
    if(mpfr_cmp_d(A->m[(i-1)*A->n + j], 0) == 0 && mpfr_cmp_d(A->m[i*A->n + j], 0) == 0)
        mpset(sin, 1);
    else {
        mpfr a, b, sqrt;
        mpinit(a);
        mpinit(b);
        mpinit(sqrt);

        mpfr_set(a, A->m[(i-1)*A->n + j], MPFR_RNDN);
        mpfr_set(b, A->m[i*A->n + j], MPFR_RNDN);
        mpfr_mul(sqrt, a, a, MPFR_RNDN); // sqrt = a*a
        mpfr_fma(sqrt, b, b, sqrt, MPFR_RNDN); // sqrt = b*b + (a*a) 
        mpfr_sqrt(sqrt, sqrt, MPFR_RNDN); // sqrt = sqrt((b*b + a*a))
        mpfr_div(sin, b, sqrt, MPFR_RNDN); // sin = b / sqrt((b*b + a*a))
        
        mpfr_clear(a);
        mpfr_clear(b);
        mpfr_clear(sqrt);
    }
}

/*
Calcul la matrice de Givens Gij de la matrice A de taille n*n 
En utilisant la librairie MPFR
@param i : un entier
@param j : un entier
@param *A : la matrice A

@return une matrice de taille n*n qui est la matrice de Givens Gij
*/
MPFR_Matrix *MPFR_givens(int i, int j, MPFR_Matrix *A){
    MPFR_Matrix *G = init_MPFR_eye(A->n);
    mpfr cos, sin, negsin;

    mpinit(cos);
    mpinit(sin);
    mpinit(negsin);

    MPFR_find_cos(i, j, A, cos);
    MPFR_find_sin(i, j, A, sin);
    // mpfr_printf("cos = %.10Rf\n", cos);
    // mpfr_printf("sin = %.10Rf\n", sin);

    mpfr_neg(negsin, sin, MPFR_RNDN);

    mpfr_set(G->m[j*G->n + j], cos, MPFR_RNDN);
    mpfr_set(G->m[j*G->n + i], sin, MPFR_RNDN);
    mpfr_set(G->m[i*G->n + j], negsin, MPFR_RNDN);
    mpfr_set(G->m[i*G->n + i], cos, MPFR_RNDN);

    // printf("Givens G = \n");
    // print_MPFR_matrix(G);

    mpfr_clear(cos);
    mpfr_clear(sin);
    mpfr_clear(negsin);

    return G;
}

/*
Calcul la matrice de Givens Gij de la matrice A de taille n*n pour la fonction Hessenberg
En utilisant la librairie MPFR
@param i : un entier
@param j : un entier
@param *G : la matrice G, qui est la matrice de Givens Gij
@param *A : la matrice A
*/
void MPFR_givens_matrix(int i, int j, MPFR_Matrix *G, MPFR_Matrix *A){
    mpfr cos, sin, negsin; 

    mpinit(cos);
    mpinit(sin);
    mpinit(negsin);

    MPFR_find_cos_matrix(i, j, A, cos);
    MPFR_find_sin_matrix(i, j, A, sin);
    // mpfr_printf("cos = %.10Rf\n", cos);
    // mpfr_printf("sin = %.10Rf\n", sin);
    
    mpfr_neg(negsin, sin, MPFR_RNDN);

    mpfr_set(G->m[(i-1)*G->n + (i-1)], cos, MPFR_RNDN); // En haut a gauche
    mpfr_set(G->m[(i-1)*G->n + i], sin, MPFR_RNDN); // En haut a droite
    mpfr_set(G->m[i*G->n + (i-1)], negsin, MPFR_RNDN); // En bas a gauche
    mpfr_set(G->m[i*G->n + i], cos, MPFR_RNDN); // En bas a droite

    // printf("Givens G = \n");
    // print_MPFR_matrix(G);

    mpfr_clear(cos);
    mpfr_clear(sin);
    mpfr_clear(negsin);
}

/*
QR decomposition de la matrice A de taille n*n 
En utilisant la librairie MPFR
@param *A : la matrice A

@return une structure de type QR qui contient la matrice Q et la matrice R
*/
MPFR_QR *MPFR_qr_decomposition(MPFR_Matrix *A){
    MPFR_QR *qr = (MPFR_QR *) malloc(sizeof(MPFR_QR));
    MPFR_Matrix *G, *Gt;
    qr->Q = init_MPFR_eye(A->n);
    qr->R = copy_MPFR_matrix(A);

    for(int j = 0 ; j < A->n ; j++){ // From j = 1 to n
        for(int i = j + 1 ; i < A->n ; i++){ // From i = j+1 to m
            G = MPFR_givens(i, j, qr->R);
            // printf("Givens(%d, %d)\n", i+1, j+1);
            // print_MPFR_matrix(G);

            qr->R = MPFR_matrix_mul(G, qr->R);
            // printf("R = \n");
            // print_MPFR_matrix(qr->R);
            Gt = MPFR_matrix_transpose(G);
            qr->Q = MPFR_matrix_mul(qr->Q, Gt);
            // printf("Q = \n");
            // print_MPFR_matrix(qr->Q);

            free_MPFR_matrix(G);
            free_MPFR_matrix(Gt);
        }
    }
    
    return qr;
}

/*
Calcul une matrice quasi Hessenberg de la matrice A de taille n*n
En utilisant la librairie MPFR
@param *A : la matrice A

@return une matrice de taille n*n qui est une matrice quasi Hessenberg
*/
MPFR_Matrix *MPFR_quasi_hess(MPFR_Matrix *A){
    int i = 0;
    unsigned int cpt = 0;
    
    // Nombre d'iterations max pour la boucle while (pour eviter les boucles infinies)
    unsigned int it = 0;

    mpfr threshold;
    mpinit(threshold);
    mpset(threshold, 1e-4);
    
    MPFR_Matrix *A_prime = copy_MPFR_matrix(A);
    MPFR_QR *qr = MPFR_qr_decomposition(A_prime);

    // Recupere les coefficients de la subdiagonale inferieure
    mpfr lower_diag[A->n - 1];
    for(i = 0 ; i < A->n - 1 ; i++){
        mpinit(lower_diag[i]);
        mpfr_set(lower_diag[i], A->m[(i + 1)*A->n + i], MPFR_RNDN);
    }

    while(it != NB_ITE){
        // if(it == nb_it-1) // Juste pour afficher le message
        //     printf("Trop d'iterations\n");
        
        for(i = 0 ; i < A->n -1 ; i++){
            if(mpfr_cmp(lower_diag[i], threshold) > 0)
                cpt++;
        }

        if(cpt == 0)
            break;
        
        MPFR_Matrix *Qt = MPFR_matrix_transpose(qr->Q);
        MPFR_Matrix *tmp = MPFR_matrix_mul(Qt, A_prime);

        A_prime = MPFR_matrix_mul(tmp, qr->Q);

        free_MPFR_qr(qr);
        free_MPFR_matrix(Qt);
        free_MPFR_matrix(tmp);

        qr = MPFR_qr_decomposition(A_prime);

        // Recupere les coefficients de la subdiagonale inferieure
        for(i = 0 ; i < A->n - 1 ; i++)
            mpfr_set(lower_diag[i], A_prime->m[(i+1)*A->n + i], MPFR_RNDN);

        it++;

    }

    mpfr_clear(threshold);
    for(i = 0 ; i < A->n - 1 ; i++)
        mpfr_clear(lower_diag[i]);

    free_MPFR_qr(qr);
    return A_prime;
}

/*
Calcul la matrice de Hessenberg de la matrice A de taille n*n
En utilisant la librairie MPFR
@param *A : la matrice A

@return une matrice de taille n*n qui est une matrice de Hessenberg
*/
MPFR_Matrix *MPFR_hessenberg(MPFR_Matrix *A){
    int i = 0;
    int j = 0;
    int n = A->n;

    MPFR_Matrix *hess = copy_MPFR_matrix(A);
    MPFR_Matrix *G = init_MPFR_eye(n);
    MPFR_Matrix *Gt = init_MPFR_eye(n);

    // Parcours les coeff sous la sub diagonal
    for(j = 0 ; j < n - 2 ; j++){
        for(i = n-1 ; i > j + 1 ; i--){
            // printf("i = %d, j = %d\n", i, j);
            if(mpfr_cmp_ui(hess->m[i*n + j], 0) != 0){ 
                // Ne pas utiliser mpfr_cmp() : Segmentation fault, ui marche car unsigned int
                G = init_MPFR_eye(n);
                MPFR_givens_matrix(i, j, G, hess);
                Gt = MPFR_matrix_transpose(G);
                hess = MPFR_matrix_mul(G, hess);
                hess = MPFR_matrix_mul(hess, Gt);
                free_MPFR_matrix(G);
                free_MPFR_matrix(Gt);
            }
        }
    }

    return hess;
}

/*
Calcul les valeurs propres de la matrice A de taille n*n avec la methode QR 
En utilisant la librairie MPFR
@param *A : la matrice A
@param k : nombre d'iterations

@return un tableau de taille n qui contient les valeurs propres de la matrice A
*/
void MPFR_eigenvalues(MPFR_Matrix *A, mpfr *eigen, int k){
    MPFR_QR *qr;
    MPFR_Matrix *A_prime = copy_MPFR_matrix(A);

    int i = 0;
    while(i < k){
        qr = MPFR_qr_decomposition(A_prime);

        A_prime = MPFR_matrix_mul(qr->R, qr->Q);

        free_MPFR_qr(qr);
        i++;
    }
    for(int i = 0 ; i < A->n ; i++){
        mpfr_set(eigen[i], A_prime->m[i*A->n + i], MPFR_RNDN);
    }
}