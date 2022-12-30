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

double func_sqrt(double n){
    double sqrt = n / 2;
    double tmp = 0;

    while(sqrt != tmp){
        tmp = sqrt;
        sqrt = (n / tmp + tmp) / 2;
    }

    return sqrt;
}

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

double find_cos_matrix(int i, int j, Matrix *A){
    // printf("Aij = %f\n",A->m[i*A->n + j]);
    // printf("Ajj = %f\n",A->m[j*A->n + j]);
    if(A->m[(i-1)*A->n + j] == 0 && A->m[i*A->n + j] == 0)
        return 1;
    else 
        return A->m[(i-1)*A->n + j] / sqrt(A->m[(i-1)*A->n + j]*A->m[(i-1)*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
        // si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[j*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

double find_sin(int i, int j, Matrix *A){
    if(A->m[i*A->n + j] == 0 && A->m[j*A->n + j] == 0)
        return 1;
    else   
        return A->m[i*A->n + j] / sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
        //si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[i*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

double find_sin_matrix(int i, int j, Matrix *A){
    if(A->m[i*A->n + j] == 0 && A->m[(i-1)*A->n + j] == 0)
        return 1;
    else   
        return A->m[i*A->n + j] / sqrt(A->m[(i-1)*A->n + j]*A->m[(i-1)*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
        //si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[i*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

Matrix *givens(int i, int j, Matrix *A){
    Matrix *G = init_eye(A->n);
    // printf("n%d-j%d-1 = %d\n",A->n, j, k);
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

void givens_matrix(int i, int j, Matrix *G, Matrix *A){
    // printf("n%d-j%d-1 = %d\n",A->n, j, k);
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

QR *qr_decomposition(Matrix *A){
    QR *qr = (QR *) malloc(sizeof(QR));
    qr->Q = init_eye(A->n);
    qr->R = copy_matrix(A);

    for(int j = 0 ; j < A->n ; j++){ // From j = 1 to n
        for(int i = j + 1 ; i < A->n ; i++){ // From i = j+1 to m
            Matrix *G = givens(i, j, qr->R);
            // printf("Givens(%d, %d)\n", i+1, j+1);
            // print_matrix(G);

            qr->R = matrix_mul(G, qr->R);
            // printf("R = \n");
            // print_matrix(qr->R);
            Matrix *Gt = matrix_transpose(G);
            qr->Q = matrix_mul(qr->Q, Gt);
            // printf("Q = \n");
            // print_matrix(qr->Q);

            free_matrix(G);
            free_matrix(Gt);
        }
    }

    return qr;
}

Matrix *quasi_hess(Matrix *A){
    unsigned int i = 0;
    // unsigned int j = 0;
    // unsigned int k = 0;
    unsigned int cpt = 0;
    
    // Nombre d'iterations max pour la boucle while (pour eviter les boucles infinies)
    unsigned int it = 0;
    unsigned int nb_it = 100; 

    // unsigned int nb_zero = 0;
    // unsigned int cpt_zero = 0;
    
    double threshold = 1e-4;
    // double threshold = 0.5;
    // printf("threshold = %f\n", threshold);


    // QR *qr = (QR *) malloc(sizeof(QR));
    
    Matrix *A_prime = copy_matrix(A);
    QR *qr = qr_decomposition(A_prime);

    // for(i = 1 ; i < A->n - 1 ; i++){
    //     nb_zero += i;
    // }
    // printf("nb zero : %d\n", nb_zero);
    // Recupere les zeros
    // double zeros[nb_zero]; // +1 ici car sinon il y a un segmentation fault
    // for(i = 0 ; i < A->n - 2 ; i++){
    //     for(j = 0 ; j < i + 1 ; j++)
    //         zeros[k++] = A->m[(i + 2)*A->n + j];
    // }

    // Recupere les coefficients de la subdiagonale inferieure
    double lower_diag[A->n - 1];
    for(i = 0 ; i < A->n - 1 ; i++){
        lower_diag[i] = A->m[(i + 1)*A->n + i];
    }

    while(it != nb_it){
        // if(it == nb_it-1) // Juste pour afficher le message
        //     printf("Trop d'iterations\n");
        
        for(i = 0 ; i < A->n -1 ; i++){
            if(lower_diag[i] > threshold)
                cpt++;
            // printf("coef = %f\n", lower_diag[i]);
        }

        // for(i = 0 ; i < nb_zero ; i++){
        //     // Si le coefficient est plus petit que le threshold car on peut pas faire de comparaison entre un double avec un 0
        //     if(zeros[i] < threshold) 
        //         cpt_zero++;
        //     // printf("coef zero = %f\ncpt_zero = %d\n", zeros[i], cpt_zero);
        // }
        if(cpt == 0 /*||cpt_zero == nb_zero*/){
            // printf("cpt = %d\ncpt_zero = %d\nbreak\n", cpt, cpt_zero);
            break;
        }
        // if((A->n == 4) && (cpt_zero >= nb_zero - 1)){
        //     printf("cpt_zero >= nb_zero - 1 : %d >= %d\n", cpt_zero, nb_zero-1);
        //     break;
        // }
        // cpt_zero = 0;
        // cpt = 0;

        Matrix *Qt = matrix_transpose(qr->Q);
        Matrix *tmp = matrix_mul(Qt, A_prime);

        A_prime = matrix_mul(tmp, qr->Q);

        free_qr(qr);
        free_matrix(Qt);
        free_matrix(tmp);

        qr = qr_decomposition(A_prime);
        // print_matrix(A_prime);

        // Recupere les coefficients de la subdiagonale inferieure
        for(i = 0 ; i < A->n - 1 ; i++)
            lower_diag[i] = A_prime->m[(i+1)*A->n + i];
        // Recupere les zeros 
        // k = 0;
        // for(i = 0 ; i < A->n - 2 ; i++)
        //     for(j = 0 ; j < i+1 ; j++) 
        //         zeros[k++] = A_prime->m[(i + 2)*A->n + j];
        
        it++;
    }

    free_qr(qr);
    return A_prime;
}


Matrix *hessenberg(Matrix *A){
    int i = 0;
    int j = 0;
    int n = A->n;
    Matrix *hess = copy_matrix(A);
    Matrix *G = init_eye(n);
    Matrix *Gt = init_eye(n);
    // Matrix *tmp = init_matrix(n);
    // Si le dernier element de la premiere ligne est different de 0 alors on fait la rotation
    // if(hess->m[(n-1)*n] != 0) 
    // G = givens(n-1, 0, hess);
    // printf("G = \n");
    // print_matrix(G);

    // Gt = matrix_transpose(G);

    // hess = matrix_mul(G, hess);

    // hess = matrix_mul(hess, Gt);
    // free_matrix(G);
    // free_matrix(Gt);

    // printf("Hess : G_(1, n) A G*_(1, n)= \n");
    // print_matrix(hess);

    // Parcours les coeff sous la sub diagonal
    for(j = 0 ; j < n - 2 ; j++){
        for(i = n-1 ; i > j + 1 ; i--){
            // if(i == n-1 && j == 0)
            //     continue;
            printf("i = %d, j = %d\n", i, j);
            if(hess->m[i*n + j] != 0){
                G = init_eye(n);
                givens_matrix(i, j, G, hess);
                // printf("G = \n");
                // print_matrix(G);
                // printf("Hess = \n");
                // print_matrix(hess);
                Gt = matrix_transpose(G);
                hess = matrix_mul(G, hess);
                hess = matrix_mul(hess, Gt);
                free_matrix(G);
                free_matrix(Gt);

                // givens_matrix(i, j, G, hess);
            }
        }
    }
    // printf("G = \n");
    // print_matrix(G);
    // hess = matrix_mul(G, hess);

    return hess;
}


/* ***** MPFR ***** */

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

void MPFR_find_cos_matrix(int i, int j, MPFR_Matrix *A, mpfr cos){
    if(mpfr_cmp_d(A->m[(i-1)*A->n + j], 0) == 0 && mpfr_cmp_d(A->m[i*A->n + j], 0) == 0)
        mpset(cos, 1);
    else {
        mpfr tmp, tmp2, tmp3, tmp4;
        mpinit(tmp);
        mpinit(tmp2);
        mpinit(tmp3);
        mpinit(tmp4);

        mpfr_mul(tmp, A->m[(i-1)*A->n + j], A->m[(i-1)*A->n + j], MPFR_RNDN);
        mpfr_mul(tmp2, A->m[i*A->n + j], A->m[i*A->n + j], MPFR_RNDN);
        mpfr_add(tmp3, tmp, tmp2, MPFR_RNDN); // tmp3 = (a^2 + b^2)

        mpfr_sqrt(tmp4, tmp3, MPFR_RNDN); // tmp4 = sqrt(a^2 + b^2)

        mpfr_div(cos, A->m[(i-1)*A->n + j], tmp4, MPFR_RNDN); // cos = a / sqrt(a^2 + b^2)

        mpfr_clear(tmp);
        mpfr_clear(tmp2);
        mpfr_clear(tmp3);
        mpfr_clear(tmp4);

        // mpfr a, b, sqrt;
        // mpinit(a);
        // mpinit(b);
        // mpinit(sqrt);

        // mpfr_set(a, A->m[(i-1)*A->n + j], MPFR_RNDN);
        // mpfr_set(b, A->m[i*A->n + j], MPFR_RNDN);
        // mpfr_mul(sqrt, a, a, MPFR_RNDN); // sqrt = a*a
        // mpfr_fma(sqrt, b, b, sqrt, MPFR_RNDN); // sqrt = b*b + (a*a) 
        // mpfr_sqrt(sqrt, sqrt, MPFR_RNDN); // sqrt = sqrt((b*b + a*a))
        // mpfr_div(cos, a, sqrt, MPFR_RNDN); // cos = a / sqrt((b*b + a*a))
        
        // mpfr_clear(a);
        // mpfr_clear(b);
        // mpfr_clear(sqrt);
    }
}

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

void MPFR_find_sin_matrix(int i, int j, MPFR_Matrix *A, mpfr sin){
    if(mpfr_cmp_d(A->m[(i-1)*A->n + j], 0) == 0 && mpfr_cmp_d(A->m[i*A->n + j], 0) == 0)
        mpset(sin, 1);
    else {
        mpfr tmp, tmp2, tmp3, tmp4;
        mpinit(tmp);
        mpinit(tmp2);
        mpinit(tmp3);
        mpinit(tmp4);

        mpfr_mul(tmp, A->m[(i-1)*A->n + j], A->m[(i-1)*A->n + j], MPFR_RNDN);
        mpfr_mul(tmp2, A->m[i*A->n + j], A->m[i*A->n + j], MPFR_RNDN);
        mpfr_add(tmp3, tmp, tmp2, MPFR_RNDN);

        mpfr_sqrt(tmp4, tmp3, MPFR_RNDN);

        mpfr_div(sin, A->m[i*A->n + j], tmp4, MPFR_RNDN);

        mpfr_clear(tmp);
        mpfr_clear(tmp2);
        mpfr_clear(tmp3);
        mpfr_clear(tmp4);

        // mpfr a, b, sqrt;
        // mpinit(a);
        // mpinit(b);
        // mpinit(sqrt);

        // mpfr_set(a, A->m[(i-1)*A->n + j], MPFR_RNDN);
        // mpfr_set(b, A->m[i*A->n + j], MPFR_RNDN);
        // mpfr_mul(sqrt, a, a, MPFR_RNDN); // sqrt = a*a
        // mpfr_fma(sqrt, b, b, sqrt, MPFR_RNDN); // sqrt = b*b + (a*a) 
        // mpfr_sqrt(sqrt, sqrt, MPFR_RNDN); // sqrt = sqrt((b*b + a*a))
        // mpfr_div(sin, b, sqrt, MPFR_RNDN); // sin = b / sqrt((b*b + a*a))
        
        // mpfr_clear(a);
        // mpfr_clear(b);
        // mpfr_clear(sqrt);
    }
}

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

MPFR_Matrix *MPFR_quasi_hess(MPFR_Matrix *A){
    unsigned int i = 0;
    unsigned int cpt = 0;
    
    // Nombre d'iterations max pour la boucle while (pour eviter les boucles infinies)
    unsigned int it = 0;
    unsigned int nb_it = 100; 

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

    while(it != nb_it){
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
            if(mpfr_cmp_ui(hess->m[i*n + j], 0) != 0){ // Ne pas utiliser mpfr_cmp() : Segmentation fault, ui marche car unsigned int
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