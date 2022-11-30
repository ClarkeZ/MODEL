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
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;
    unsigned int cpt = 0;

    unsigned int nb_zero = 0;
    unsigned int cpt_zero = 0;
    
    double threshold = 0.0001;
    // double threshold = 1e-4;
    // double threshold = 2;
    printf("threshold = %f\n", threshold);


    QR *qr = (QR *) malloc(sizeof(QR));
    
    Matrix *A_prime = copy_matrix(A);
    qr = qr_decomposition(A_prime);

    for(i = 1 ; i < A->n - 1 ; i++){
        nb_zero += i;
    }
    printf("nb zero : %d\n", nb_zero);
    double zeros[nb_zero]; // +1 ici car sinon il y a un segmentation fault
    // Recupere les zeros
    for(i = 0 ; i < A->n - 2 ; i++){
        for(j = 0 ; j < i + 1 ; j++)
            zeros[k++] = A->m[(i + 2)*A->n + j];
    }

    double lower_diag[A->n - 1];
    // Recupere les coefficients de la subdiagonale inferieure
    for(i = 0 ; i < A->n - 1 ; i++){
        lower_diag[i] = A->m[(i + 1)*A->n + i];
    }

    while(1){
        for(i = 0 ; i < A->n -1 ; i++){
            if(lower_diag[i] > threshold)
                cpt++;
            printf("coef = %f\n", lower_diag[i]);
        }
        // Probleme avec les zeros
        for(i = 0 ; i < nb_zero ; i++){
            if(zeros[i] == 0)
                cpt_zero++;
            printf("coef zero = %f\n", zeros[i]);
        }
        if(cpt == 0 || cpt_zero == nb_zero){
            printf("cpt = %d\ncpt_zero = %d\n", cpt, cpt_zero);
            printf("break\n");
            break;
        }
        // if((A->n == 4) && (cpt_zero >= nb_zero - 1)){
        //     printf("cpt_zero >= nb_zero - 1 : %d >= %d\n", cpt_zero, nb_zero-1);
        //     break;
        // }
        cpt_zero = 0;
        cpt = 0;
        A_prime = matrix_mul(matrix_mul(qr->Q, A_prime), matrix_transpose(qr->Q));
        free_qr(qr);
        qr = qr_decomposition(A_prime);
        print_matrix(A_prime);
        // Recupere les coefficients de la subdiagonale inferieure
        for(i = 0 ; i < A->n - 1 ; i++)
            lower_diag[i] = A_prime->m[(i+1)*A->n + i];
        // Recupere les zeros 
        k = 0;
        for(i = 0 ; i < A->n - 2 ; i++)
            for(j = 0 ; j < i+1 ; j++) 
                zeros[k++] = A_prime->m[(i + 2)*A->n + j];
        
    }

    free_qr(qr);
    return A_prime;
}
