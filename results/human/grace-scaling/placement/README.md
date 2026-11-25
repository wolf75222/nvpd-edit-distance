# Placement des séquences sur les processeurs Grace

Quatre cœurs alignent les chromosomes 1 à 4 de CHM13 entre les assemblages v1.0 et v1.1. Chaque cœur traite le même chromosome pendant seize appels à `nvpd_sparse`. Les 64 scores portent sur 28,3445 Gnt en additionnant les longueurs des deux assemblages. Les quatre distances valent respectivement 176, 16, 659 et 488. [Entrées et scores de référence](../inputs.tsv)

Dans le mode `master`, le premier cœur lit les quatre FASTA à la suite. Dans le mode `owner`, chaque cœur lit son chromosome en parallèle puis calcule son score. Chaque essai démarre dans un processus neuf. Le chronomètre de lecture englobe le décodage des quatre FASTA. Celui du calcul démarre après un alignement d’échauffement par chromosome. Les six blocs comprennent chacun les deux modes, qui passent chacun trois fois en première position. Les écarts-types et les intervalles de Student à 95 % décrivent les six blocs sur ce nœud. [Mesures](trials.csv) · [Rapports appariés](paired.csv)

Slurm a attribué les cœurs logiques 0 à 3 d’un même processeur Grace au job 840884 sur `romeo-a002`. Les douze processus ont utilisé ces quatre cœurs dans le même ordre. Après l’échauffement, Linux situait sur le nœud NUMA 0 la totalité des 1 689,75 Mio résidents des régions anonymes contenant les huit séquences. Le relevé porte sur les pages des régions virtuelles qui recouvrent les chaînes, y compris les octets de ces régions réservés par l’allocateur. [Placement des cœurs](environment.txt) · [Pages par nœud](residency-summary.csv)

Avec la lecture parallèle, le temps de lecture et de décodage passe de 1,826 ± 0,007 s à 0,619 ± 0,026 s. Le rapport moyen par bloc atteint 2,955 (IC à 95 % de 2,809 à 3,102). Le calcul des 64 scores passe de 339,31 ± 1,45 ms à 337,90 ± 1,04 ms, soit un écart apparié de 1,42 ms (IC à 95 % de 0,58 à 2,26 ms). La lecture et le calcul réunis prennent respectivement 2,165 s et 0,957 s. L’essentiel du gain provient de la lecture simultanée des fichiers sur le même nœud NUMA. [Durées par phase](summary.csv) · [Comptabilité Slurm](accounting.psv)

La commande suivante recalcule les statistiques à partir des observations.

```sh
python3 bench/repeated/divergence/launch/placement_analysis.py results/human/grace-scaling/placement
```
