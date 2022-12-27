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
        return A->m[j*A->n + j] / sqrt(A->m[j*A->n + j]*A->m[j*A->n + j] + A->m[i*A->n + j]*A->m[i*A->n + j]);
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

void givens_matrix(int i, int j, Matrix *G, Matrix *A){
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
    unsigned int nb_it = 100000; 

    // unsigned int nb_zero = 0;
    // unsigned int cpt_zero = 0;
    
    double threshold = 1e-6;
    // double threshold = 0.5;
    printf("threshold = %f\n", threshold);


    QR *qr = (QR *) malloc(sizeof(QR));
    
    Matrix *A_prime = copy_matrix(A);
    qr = qr_decomposition(A_prime);

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
        if(it == nb_it-1) // Juste pour afficher le message
            printf("Trop d'iterations\n");
        
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


// Matrix *hessenberg(Matrix *A){
//     int n;

//     if((n = A->n) < 3){
//         perror("Taille de la matrice doit etre superieur a 3 pour la matrice Hessenberg");
//         exit(EXIT_FAILURE);
//     }

//     Matrix *hess = init_matrix(n);
//     for(int i = 0 ; i < n ; i++)
//         for(int j = 0 ; j < n ; j++){
//             printf("i = %d, j = %d\n", i, j);
//             // Matrix *G = givens(n - 2, n-1, A); // -1 car on commence a 0
//             Matrix *G = givens(i, j, A);
//             Matrix *tmp = matrix_mul(G,A);
//             printf("G * A = \n");
//             print_matrix(tmp);
//             Matrix *G2 = givens(i, j, tmp);
//             Matrix *Gt = matrix_transpose(G2);
//             // Matrix *Gt = matrix_transpose(G);  
//             Matrix *tmp2 = matrix_mul(tmp, Gt);
//             printf("G A G* = \n");
//             print_matrix(tmp2);
//             printf("====================\n");
//         }
//     // printf("G = \n");
//     // print_matrix(G);
//     // printf("A = \n");
//     // print_matrix(A);

//     // Matrix *Gt = matrix_transpose(G);
//     // hess = matrix_mul(tmp, Gt);
//     // free_matrix(tmp);
//     // free_matrix(Gt);
//     // printf("Hess : G_(n-1, n) A G*_(n-1, n)= \n");
//     // print_matrix(hess);
    
//     // G = givens(n - 2, n - 1, A); // A ou A' ?
//     // hess = matrix_mul(matrix_mul(G, hess), matrix_transpose(G));
//     // print_matrix(hess);

//     // // G = givens(n - 1, n, hess); // (?) G de A'' ?
//     // hess = matrix_mul(matrix_mul(G, hess), matrix_transpose(G));

//     // // G = givens(n - 2, n - 1, hess); // Logiquement ? comme en haut
//     // print_matrix(hess);

//     return hess;
// }

Matrix *hessenberg(Matrix *A){
    int i = 0;
    int j = 0;
    int n = A->n;
    Matrix *hess = copy_matrix(A);
    Matrix *G = init_matrix(n);
    Matrix *Gt = init_matrix(n);
    // Matrix *tmp = init_matrix(n);
    G = givens(0, 0, hess);
    Gt = matrix_transpose(G);
    hess = matrix_mul(G, hess);
    hess = matrix_mul(hess, Gt);
    free_matrix(G);
    free_matrix(Gt);

    printf("Hess : G_(1, n) A G*_(1, n)= \n");
    print_matrix(hess);

    G = init_eye(n);
    for(i = 0 ; i < A->n - 2 ; i++){
        for(j = 0 ; j < i + 1 ; j++){
            // if(i == j)
            //     continue;
            if(hess->m[i*n + j] != 0){
                // G = givens(i, j, hess);
                // Gt = matrix_transpose(G);
                // tmp = matrix_mul(G, hess);
                // hess = matrix_mul(tmp, Gt);
                // hess = matrix_mul(G, hess);
                // A = matrix_mul(G, A);
                // free_matrix(G);
                // free_matrix(Gt);
                // free_matrix(tmp);
                givens_matrix(i, j, G, hess);
                printf("i = %d, j = %d\n", i, j);
            }
        }
    }
    hess = matrix_mul(G, hess);
    // parcours the under the subdiagonal elements


    return hess;


}