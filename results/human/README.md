# Deux références humaines sur le chromosome 21

Les fenêtres comparent la séquence primaire du chromosome 21 de GRCh38.p14 et celle de T2T-CHM13v2.0. Les deux assemblages représentent des matériels biologiques différents. La distance mesure les différences entre les chaînes déposées. [GRCh38, NC_000021.9](https://www.ncbi.nlm.nih.gov/nuccore/NC_000021.9/) · [T2T-CHM13, NC_060945.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_060945.1/)

| Cas | Intervalle GRCh38 | Intervalle CHM13 | Distance attendue ou vérifiée |
|---|---|---|---:|
| chr21-2m | [34 000 000, 36 000 000) | [32 381 795, 34 383 058) | 5 599 |
| chr21-4m | [34 000 000, 38 000 000) | [32 381 795, 36 383 559) | 12 482 |
| chr21-full-insert64 | — | 45 090 682 bases contre la même séquence augmentée de 64 bases | 64 |

Le script vérifie l’empreinte de chaque chromosome et de la chaîne de correspondance UCSC. Chaque borne des fenêtres se trouve dans un bloc aligné de la chaîne 23, sur le brin positif. Les positions commencent à zéro et la borne droite est exclue. Les fenêtres de deux et quatre millions de bases se recouvrent. Le script conserve les insertions et les délétions à l’intérieur des fenêtres. Les 41 bases à chaque extrémité n’apparaissent qu’une fois sur leur chromosome d’origine. Les deux fenêtres contiennent uniquement A, C, G et T. [Chaîne de correspondance hg38 vers hs1](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz)

Pour le contrôle sur chromosome entier, le générateur insère un G à chacune des 64 frontières `floor(k × 45 090 682 / 65)`, pour `k` de 1 à 64. La différence de longueur impose au moins 64 opérations. Les 64 insertions atteignent cette borne. Les accessions désignent les séquences sources, y compris pour la copie modifiée. [Manifeste des FASTA](manifest.tsv)

## Exécution

```sh
python3 scripts/human_reference.py
sbatch --account="$ROMEO_PROJECT" bench/human.sbatch
```

Le job réserve un GH200, quatre cœurs Arm, 16 Gio et 45 minutes. Chaque ligne de `trials.csv` correspond à un nouveau processus. Le champ `seconds` mesure l’appel de l’algorithme. Le champ `wall_seconds` couvre en plus la lecture du FASTA et le lancement du processus. Les temps CUDA incluent l’initialisation du contexte et les copies. `device_bytes` compte les allocations CUDA explicites. Une valeur nulle de `peak_rss_kib` renvoyée par le système pour un processus ayant effectivement chargé des données indique une mesure indisponible.

Le moteur CUDA compare les deux fenêtres homologues à coût unitaire, avec contrôle par edlib et WFA2. La dépendance maximale à gauche atteint 290 cellules sur la fenêtre de deux millions et 622 sur celle de quatre millions. La limite de 64 cellules de la tuile déclenche le noyau à une ligne dans les deux cas. Pour le chromosome entier, edlib, WFA2 et la méthode en bande avec suppression des extrémités communes traitent le contrôle à 64 insertions. Un calcul dense sur cette paire demanderait plus de deux millions de milliards de cellules.

La variable facultative `KOKKOS_SOURCE` désigne un checkout de Kokkos qui contient la révision `175257a51ff29a0059ec48bcd233ee096b2c0438`. Le script de reproduction conserve CUDA 13.0.2 pour les essais natifs et charge CUDA 12.6.3 dans un sous-processus pour compiler cette version de Kokkos. Le script soumis au job 838220 a utilisé CUDA 13 pour les deux étapes. [Script soumis](838220/submitted.sbatch)

## Mesures du job 838220

Les dix-huit appels natifs ont tous donné la distance de référence. Les temps ci-dessous sont les médianes de deux processus, en secondes. Les colonnes de temps couvrent l’algorithme. Pour la paire de chromosomes entiers, le temps médian de bout en bout atteint 0,476 s avec WFA2, dont 0,0113 s dans l’alignement. [Essais individuels](838220/trials.csv)

| Paire | edlib | WFA2 | NvPD CUDA | Bande avec extrémités communes retirées |
|---|---:|---:|---:|---:|
| Fenêtre de 2 M | 0,790 | 0,0538 | 18,35 | — |
| Fenêtre de 4 M | 2,991 | 0,348 | 90,89 | — |
| CHM13 entier + 64 insertions | 0,643 | 0,0113 | — | 6,804 |

NvPD atteint respectivement 218,1 et 176,1 GCUPS sur les deux fenêtres. Il alloue 48 030 336 et 96 042 360 octets sur le GPU. Les pics RSS des processus CUDA sont de 216 et 252 Mio côté hôte. Sur chromosome entier, les pics RSS sont de 180 Mio pour edlib, 144 Mio pour WFA2 et 108 Mio pour la méthode en bande.

Les étapes Slurm 838220.0 et 838220.1 ont terminé avec le code 0. L’étape facultative 838220.2 a échoué pendant la compilation de Kokkos 4.5.01 avec CUDA 13.0.2. Cette combinaison utilise des signatures CUDA devenues incompatibles pour les graphes et la gestion de mémoire. Le job global porte donc le code 2. [Comptabilité Slurm](838220/accounting.psv) · [Diagnostic compilateur](838220/kokkos-build-error.txt)
