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
    if(A->m[i*A->n + j] == 0 && A->m[j*A->n + j] == 0)
        return 1;
    else 
        // return A->m[j*A->n + j] / sqrt(A->m[j*A->n + j] + A->m[i*A->n + i]);
        // si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        return A->m[j*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

double find_sin(int i, int j, Matrix *A){
    if(A->m[i*A->n + j] == 0 && A->m[j*A->n + j] == 0)
        return 1;
    else   
        // return A->m[i*A->n + j] / sqrt(A->m[j*A->n + j] + A->m[i*A->n + i]);
        //si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        return A->m[i*A->n + j] / func_sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
}

Matrix *Givens(int i, int j, Matrix *A){
    Matrix *G = init_eye(A->n);
    double c = find_cos(i, j, A);
    double s = find_sin(i, j, A);
    // printf("c = %f\n", c);
    // printf("s = %f\n", s);
    G->m[i*A->n + i] = c; // En haut a gauche 
    G->m[i*A->n + j] = -s; // En haut a droite
    G->m[j*A->n + i] = s; // En bas a gauche
    G->m[j*A->n + j] = c; // En bas a droite

    return G;
}

QR *qr_decomposition(Matrix *A){
    QR *qr = (QR *) malloc(sizeof(QR));
    // Matrix *Q = mat_Q(A);
    // Matrix *R = matrix_mul(Q, A);
    qr->Q = init_eye(A->n);
    qr->R = copy_matrix(A);

    for(int j = 0 ; j < A->n ; j++){ // From j = 1 to n
        for(int i = j + 1 ; i < A->n ; i++){ // From i = j+1 to m
            Matrix *G = Givens(i, j, qr->R);
            // printf("Givens(%d, %d)\n", i+1, j+1);
            // print_matrix(G);

            qr->R = matrix_mul(G, qr->R);
            // printf("R = \n");
            // print_matrix(qr->R);
            qr->Q = matrix_mul(qr->Q, matrix_transpose(G));
            // printf("Q = \n");
            // print_matrix(qr->Q);

        }
    }

    return qr;
}

Matrix *quasi_hess(Matrix *A){
    QR *qr = (QR *) malloc(sizeof(QR));
    
    Matrix *A_prime = copy_matrix(A);
    qr = qr_decomposition(A_prime);

    double threshold = 0.00001;
    // double threshold = -1e-4;
    printf("threshold = %f\n", threshold);
    int size_n = (A->n / 2);
    double coef1, coef2;
    // printf("coef sub : \n");
    coef1 = A->m[size_n * A->n + 0];
    coef2 = A->m[(size_n + 1) * A->n + 1];
    // printf("coef 1 = %f\n", coef1);
    // printf("coef 2 = %f\n", coef2);
    // while((coef1 > threshold) || (coef2 > threshold)){
    //     coef1 -= 0.1;
    //     coef2 -= 0.1;
    //     printf("coef 1 = %f\n", coef1);
    //     printf("coef 2 = %f\n", coef2);
    // }
    int i = 0;
    // Attention CAS OU les coef = 0, boucle inf
    while(abs(coef1) > threshold && abs(coef2) > threshold){
        coef1 = A_prime->m[size_n * A->n + 0];
        coef2 = A_prime->m[(size_n + 1) * A->n + 1];
        printf("coef 1 = %f\n", coef1);
        printf("coef 2 = %f\n", coef2);
        if(abs(coef1) < threshold)
            printf("***** coef 1 < threshold : %d *****\n", i);
        if(abs(coef2) < threshold)
            printf("***** coef 2 < threshold : %d *****\n", i);
        if(abs(coef1) < threshold && abs(coef2) < threshold){
            printf("***** les 2 coef < threshold : %d *****\n", i);
            print_matrix(A_prime);
            return A_prime;
        }
            

        // if(coef1 <= threshold && coef2 <= threshold)
        //     return A_prime;
        // else{
            A_prime = matrix_mul(matrix_mul(qr->Q, A_prime), matrix_transpose(qr->Q));
            free_qr(qr);
            qr = qr_decomposition(A_prime);
        // }
        // printf("Matrix Quasi hessenberg \n");
        // print_matrix(A_prime);
        i++;
    }

    free_qr(qr);
    return A_prime;
}
