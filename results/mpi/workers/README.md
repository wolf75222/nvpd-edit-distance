# Travailleurs Kokkos dans les rangs MPI

L’allocation 839197 compare deux façons de lancer les huit travailleurs d’un rang. OpenMP ouvre une région parallèle. Kokkos distribue huit itérations statiques sur son espace hôte OpenMP. Chaque travailleur appelle le même noyau NvPD par paire. Le sélecteur `--workers` s’applique aux matrices avec une politique `--schedule`.

La matrice humaine triée par longueur contient 128 fenêtres et 8 128 paires. Chaque moteur exécute les cinq politiques sur un puis deux nœuds EPYC, avec 24 rangs par nœud. Les 100 lancements ont retrouvé les scores d’edlib et de WFA2. Les cinq blocs mélangent l’ordre des lancements. Le chronomètre couvre la construction de l’ordre, les alignements et la réduction après un premier appel de chauffe.

Pour Hilbert pondéré sur deux nœuds, le rapport apparié des temps OpenMP/Kokkos vaut 0,998 avec un intervalle de Student à 95 % de 0,969 à 1,027. Les durées moyennes sont respectivement 0,3960 et 0,3968 s. Les mesures d’ordonnancement sur une autre paire de nœuds forment une campagne distincte. [Observations](839197/trials.csv) · [Moyennes et comparaisons appariées](839197/summary.csv) · [Panel et références](../schedules/README.md)

Pour répéter la comparaison sur ROMEO, définir `KOKKOS_SOURCE` comme le chemin du dépôt Kokkos contenant le commit de la version 4.5.01. Le premier argument du script est le FASTA trié par longueur, le second le TSV des distances de référence. Le binaire exige les options CMake `CHPS_MPI=ON` et `CHPS_KOKKOS=ON`, ainsi qu’un Kokkos compilé avec `Kokkos_ENABLE_OPENMP=ON`.

```sh
sbatch -A "$ROMEO_PROJECT" bench/mpi/workers/run.sbatch \
  "$PWD/data/mpi-schedules/length.fa" "$PWD/results/mpi/schedules/839194/scores.tsv"
python3 bench/mpi/workers/summarize.py results/mpi/workers/839197
```
