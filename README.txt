Ce code marche pour les versions MPFR et GMP suivantes :
MPFR library: 4.1.0-p13
MPFR header: 4.1.0-p13 (based on 4.1.0)
GMP : 6.2.1


Pour compiler le programme, mettez vous dans le dossier /Projet et faites :
$ make

Si tout ce passe bien vous aurez :
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o algo.o -c algo.c
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o base.o -c base.c
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o benchmark.o -c benchmark.c
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o main.o -c main.c
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o matrix.o -c matrix.c
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o tests.o -c tests.c
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o unit_test.o -c unit_test.c
gcc -g -Wall -O2 -pg -lm -lgmp -lmpfr -Wextra -Werror -o main algo.o base.o benchmark.o main.o matrix.o tests.o unit_test.o -lm -lmpfr -lgmp


Pour lancer l'exécutable vous pouvez utiliser les différents options :
$ ./main [-size n] [-iteration i]
ou
$ ./main [-n n] [-i i]

Les arguments peuvent être rentré dans l'ordre souhaité.

Si aucun argument n'est mis, ou si des arguments sont manquants alors les valeurs par défaut seront :
	- size : 4
	- iteration : 1

Exemple :
$ ./main
$ ./main -n 10
$ ./main -n 4 -i 5

Si vous avez un doute, vous pouvez utiliser la commande :
$ ./main -h
ou
$ ./main -help


Une fois exécuté, le programme vous demandera de choisir quel benchmark lancer. 




Pour exécuter les tests unitaires il faut décommenter les fonctions dans le fichier unit_test.c et exécuter le programme.

Par exemple vous pouvez tester l'algorithme qui rend une matrice A Hessenberg supérieure,
il suffit alors de décommenter la ligne 60, test_hessenberg(n).