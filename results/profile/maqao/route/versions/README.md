# Profil MAQAO du répartiteur Edlib et WFA2

Sur quatre cœurs de l’EPYC 9654, MAQAO 2026.1.1 profile l’exécutable `editdist-route` pendant un calcul d’échauffement suivi de six calculs de la matrice des 64 génomes. La politique `hybrid-counts` confie 610 paires à Edlib et 1 406 à WFA2. La matrice obtenue coïncide avec la référence. [Empreintes de l’exécutable et de la matrice](environment.txt) · [Allocation Slurm](accounting.tsv)

Le calcul de Myers dans Edlib concentre 84,96 % des échantillons de cycles, boucles internes comprises. L’extension AVX-512 de WFA2 en reçoit 11,09 %, contre 0,145 % pour la construction des esquisses. LProf attribue 0,238 % des échantillons à OpenMP. [Fonctions dont la couverture atteint 0,1 %](functions.csv) · [Répartition des échantillons](profiles.csv)

## Boucle de Myers et extension des fronts

La boucle interne de Myers représente à elle seule 81,17 % des échantillons. Le modèle CQA estime une itération assembleur à 11 cycles. Il signale des dépendances entre itérations et suggère d’examiner les accès aux structures et le déroulage de la boucle. Dans la boucle externe d’extension de WFA2, une instruction au moins échappe au modèle Zen 4. CQA ne fournit donc pas d’estimation en cycles pour cette boucle. [Couverture et coût des boucles](loops.csv)

Le calcul des bornes du front (`wavefront_compute_trim_ends`) représente 0,0173 % des échantillons. CQA y relève des conversions et des possibilités de vectorisation. Dans l’allocation des fronts (`wavefront_slab_allocate`, 0,0189 %), il signale aussi les appels de fonctions. [Profils des fonctions](function_cqa.csv)

## Paramètres de collecte

ONEVIEW détecte les symboles de débogage. Il recommande de compiler avec `-fno-omit-frame-pointer` pour reconstruire les piles d’appels. Le journal signale aussi un risque de perte d’échantillons lié au tampon Linux de `perf`. La commande `lstopo-no-graphics` manque dans l’environnement de cette collecte, et l’accès aux compteurs de l’espace noyau est restreint. [Indicateurs ONEVIEW](oneview.csv) · [Effet de la fréquence d’échantillonnage sur NvPD](../../versions/README.md)

## Reproduction

Depuis la racine du dépôt sur ROMEO, préparer le FASTA de la cohorte puis lancer la commande suivante avec un répertoire de sortie neuf.

```sh
sbatch -A "$ROMEO_PROJECT" bench/cpu/maqao/route.sbatch \
  build/hybrid-matrix/viral/surveillance.fa results/phylogeny/sars-edit/matrix.tsv \
  build/revue-maqao/route-versions 2026.1.1
```

Dans l’archive produite, la page `RESULTS/editdist-route_one_html/index.html` présente le profil. Les diagnostics de chaque boucle se trouvent dans `static_data/cqa/`.
