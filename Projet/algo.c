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
        return A->m[j*A->n + j] / sqrt(A->m[j*A->n + j] + A->m[i*A->n + i]);
        // si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[j*A->n + j] / func_sqrt(A->m[j*A->n + j] + A->m[i*A->n + i]);
}

double find_sin(int i, int j, Matrix *A){
    if(A->m[i*A->n + j] == 0 && A->m[j*A->n + j] == 0)
        return 1;
    else   
        return A->m[i*A->n + j] / sqrt(A->m[j*A->n + j] + A->m[i*A->n + i]);
        //si sqrt de math.h ne marche pas utiliser la fonction ci-dessous
        // return A->m[i*A->n + j] / func_sqrt(A->m[j*A->n + j] + A->m[i*A->n + i]);
}

Matrix *Givens(int i, int j, Matrix *A){
    Matrix *G = init_eye(A->n);
    
    double c = find_cos(i, j, A);
    double s = find_sin(i, j, A);
    
    G->m[j*A->n + j] = c;
    G->m[j*A->n + i] = s;
    G->m[i*A->n + j] = 0 - s;
    G->m[i*A->n + i] = c;

    return G;
}

Matrix *mat_Q(Matrix *A){
    Matrix *Q = init_eye(A->n);
    Matrix *G = init_matrix(A->n);

    for(int i = 0 ; i < A->n ; i++){
        for(int j = 0 ; j < A->n ; j++){
            G = Givens(i, j, A);
            Q = matrix_mul(G, Q);
            A = matrix_mul(G, A);
        }
    }

    return Q;
}

QR *qr_decomposition(Matrix *A){
    QR *qr = (QR *) malloc(sizeof(QR));
    Matrix *Q = mat_Q(A);
    Matrix *R = matrix_mul(Q, A);

    for(int i = 0 ; i < A->n ; i++){
        for(int j = 0; j < i ; j++){
            R->m[i*A->n + j] = 0;
        }
    }

    Q = matrix_transpose(Q);

    qr->Q = Q;
    qr->R = R;

    return qr;
}
