# Travailleurs MPI sur les loci du chromosome 1

Le travail ROMEO 839351 compare le pool OpenMP au pool OpenMP de Kokkos 4.5.01 dans le même exécutable MPI. Les deux pools appellent le noyau `nvpd-indexed`. Chaque rang construit ses tables de dernières occurrences et les partage entre huit travailleurs. Le calcul porte sur 768 loci de CHM13, soit 294 528 distances. L’empreinte SHA-256 de chaque matrice doit correspondre à celle du [calcul edlib](../839322/edlib-oracle.json).

Cinq blocs permutent les huit configurations. Un placement occupe 24 rangs sur un nœud ou 48 rangs sur deux nœuds. Chaque processus exécute un échauffement puis une mesure. Le chronomètre démarre après la diffusion des entrées. Il couvre la préparation des tâches et des index, les scores et la réduction de la matrice. Les deux politiques répartissent les paires selon l’ordre de Hilbert pondéré par les produits de longueurs ou par vol de paquets. Le solveur calcule toutes les distances sans activer le raccourci d’inclusion.

Le fichier `trials.csv` distingue les cœurs échantillonnés au placement des travailleurs ayant traité au moins une paire. Le fichier `summary.csv` donne les moyennes, les écarts-types et les intervalles de Student à 95 % sur les cinq blocs. Pour deux nœuds, le gain apparié vaut \(T_1/T_2\) et l’efficacité vaut \(T_1/(2T_2)\). Le rapport entre les deux pools vaut \(T_{\rm OpenMP}/T_{\rm Kokkos}\).

Le [script Slurm](../../../../bench/mpi/chromosome/workers.sbatch) prend le FASTA et le JSON de référence en arguments. La variable `KOKKOS_SOURCE` désigne un dépôt Kokkos contenant la révision `175257a51ff29a0059ec48bcd233ee096b2c0438`.
