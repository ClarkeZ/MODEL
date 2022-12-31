// #include "algo.h"
#include "unit_test.h"
#include "benchmark.h"

void choix_bench(){
    printf("\n--- Choix du benchmark ---\n");
    printf(" 1 - Comparaison entre quasi Hessenberg et Hessenberg (double précision)\n");
    printf(" 2 - Comparaison entre quasi Hessenberg et Hessenberg (MPFR)\n");
    printf(" 3 - Comparaison des différents QR decomposition\n");
    printf("        => double precision vs MPFR\n");
    printf(" 4 - Comparaison des différents algorithmes rendant une matrice quasi Hessenberg supérieure\n");
    printf("        => double precision vs MPFR\n");
    printf(" 5 - Comparaison des différents algorithmes rendant une matrice Hessenberg supérieure\n");
    printf("        => double precision vs MPFR\n");
    printf(" 6 - Quitter\n");
    printf("\n Votre choix :\n");
}

int main(int argc, char *argv[]){
    srand(time(NULL));
    
    int n = 12;
    unsigned int ite = 1;

    /* ----- AFFECTATION DES VALEURS -----*/
    
    if(argc > 1){
        if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0){
            printf("Commande : \n./main [-size n] [-iteration i]\n");
            printf("Commande : \n./main [-n n] [-i it]\n");
            return EXIT_SUCCESS;
        }
        else{
            while(argc > 1){
                if(strcmp(argv[argc-2], "-size") == 0 || strcmp(argv[argc-2], "-n") == 0){
                    n = atoi(argv[argc-1]);
                    printf("n = %d\n", n);
                }
                if(strcmp(argv[argc-2], "-iteration") == 0 || strcmp(argv[argc-2], "-i") == 0){
                    ite = atoi(argv[argc-1]);
                    printf("ite = %u\n", ite);
                }
                argc -= 2;
            }
        }
    }

    printf("\n --- Initialisation ---\n");
    printf("Taille de la matrice : %d x %d\n", n, n);
    printf("Nombre d'itération(s) : %u\n\n", ite);

    /* ----- TESTS UNITAIRES AVEC LE TEMPS D'EXECUTION -----*/
    
    test_base();
    
    test_matrix(n);
    
    test_algo(n);

    /* ----- BENCHMARKS -----*/

    int choix = 0;

    choix_bench();
    fflush(stdout);
    scanf("%d", &choix);

    while(choix >= 1 || choix < 6){
        switch(choix){
            case 1:
                benchmark_quasi_hess_vs_hessenberg(n, ite);
                break;
            case 2:
                benchmark_MPFR_quasi_hess_vs_hessenberg(n, ite);
                break;
            case 3:
                benchmark_qr_decomposition(n, ite);
                break;
            case 4:
                benchmark_quasi_hess(n, ite);
                break;
            case 5:
                benchmark_hessenberg(n, ite);
                break;
            case 6:
                return EXIT_SUCCESS;
            default:
                printf("Veuillez choisir un nombre entre 1 et 6\n");
                break;
        }
        if(choix != 6){
            choix_bench();
            fflush(stdout);
            scanf("%d", &choix);
        }
    }
   
    return 0;
}