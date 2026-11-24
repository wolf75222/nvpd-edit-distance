# Reconstruction de l’alignement

Les deux méthodes utilisent un fil sur Apple M1 Pro avec Apple Clang 21.0.0 en mode Release. Les entrées sont les préfixes de 2 000, 8 000 et 25 000 bases de la première paire synthétique de 100 000 bases à 2 %. Chaque configuration démarre cinq processus. Le chronométrage de l'appel inclut les allocations et l'écriture du CIGAR en mémoire. Le RSS couvre le processus. Le script rejoue chaque chemin et compare son coût au score global d'edlib. [Essais individuels](trials.csv)

```sh
python3 bench/traceback.py build/editdist data/synthetic/pairs/n100000_r0.02_p1_*.fa build/traceback.csv
```
