# Reconstruction de l’alignement

Les deux méthodes utilisent un fil sur Apple M1 Pro avec Apple Clang 21.0.0 en mode Release. Les entrées sont les préfixes de 2 000, 8 000 et 25 000 bases de la première paire synthétique de 100 000 bases à 2 %. Chaque configuration démarre cinq processus. Le chronométrage de l'appel inclut les allocations et l'écriture du CIGAR en mémoire. Le RSS couvre le processus. Le script rejoue chaque chemin et compare son coût au score global d'edlib. [Essais individuels](trials.csv)

```sh
python3 bench/cpu/traceback/traceback.py build/editdist data/synthetic/pairs/n100000_r0.02_p1_*.fa build/traceback.csv
```

## Limite de la matrice dense

Le [job ROMEO 839287](boundary-accounting.psv) compare les préfixes de la même paire sur quatre cœurs Grace et un GH200, avec GCC 11.5, CUDA 12.6.85 et Kokkos 4.5.01. À 20 000 × 20 000 bases, la référence séquentielle, NvPD/OpenMP, les deux variantes CUDA, Kokkos/CUDA et Hirschberg produisent un CIGAR de coût 454. Le contrôle rejoue chaque opération et compare le score à celui d’edlib. Le calcul par anti-diagonales obtient le même score.

En ajoutant une base au premier préfixe, la distance passe à 453. En l’ajoutant au second, elle passe à 455. Les cinq variantes denses signalent la limite dans les deux cas. Le calcul du coût et la reconstruction par Hirschberg aboutissent. Pour le carré de 20 000 bases, les tampons GPU occupent 480 024 octets en coût seul et 1 600 480 020 octets avec la matrice complète. [Résultats par appel](boundary-gh200.csv)

```sh
python3 tests/trace_boundary.py build/editdist seq nvpd-omp nvpd-cuda nvpd-cuda-local nvpd-kokkos
```
