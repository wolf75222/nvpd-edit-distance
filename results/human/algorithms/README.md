# BiWFA sur deux assemblages du chromosome 21

`editdist pair --algorithm biwfa` calcule la distance d’édition globale à coût unitaire. L’option `--trace` renvoie un CIGAR étendu orienté de la première séquence vers la seconde. `wfa2 --trace` fournit un témoin unidirectionnel avec le même format de sortie. Les deux moteurs comparent littéralement les symboles et désactivent les heuristiques. [Interface C](../../../src/libraries/wfa2/wfa_adapter.c) · [Interface C++](../../../src/libraries/baselines.cpp)

## Algorithme et chemin exécuté

Pour un score entier donné et chaque diagonale, le front retient le point le plus avancé accessible puis prolonge les égalités. Une substitution, une insertion ou une suppression coûte une unité. BiWFA fait progresser deux fronts depuis les extrémités opposées. Dans la métrique `edit`, il détecte leur chevauchement sur la même diagonale dans la composante M. La somme des deux coûts donne le score optimal à la jonction. Le mode score rend ce coût ou le coût d’arrivée à l’autre extrémité. Le mode trace coupe les séquences à la jonction et reconstruit récursivement les deux morceaux. Pour une longueur maximale $n$ et une distance $s\geq1$, l’article donne un temps $O(ns)$ et un espace de travail $O(s)$. La lecture et la restitution d’un alignement identique prennent un temps linéaire même pour $s=0$. Dans le processus mesuré, les copies de séquences et le CIGAR de la bibliothèque occupent aussi de la mémoire. [Marco-Sola et al., 2023](https://doi.org/10.1093/bioinformatics/btad074)

La révision WFA2 `bcf473a6` construit le moteur bidirectionnel lorsque `memory_mode=wavefront_memory_ultralow`. Elle choisit `wavefront_compute_edit` pour les deux fronts. Elle sait restituer le chemin global avec `compute_alignment`. La récursion passe en WFA auxiliaire unidirectionnel pour une sous-distance restante au plus égale à 250 et dès l’appel initial si les deux séquences font au plus 100 bases. Si un front atteint l’extrémité pendant la recherche de la jonction, la récupération du chemin par WFA auxiliaire admet un score allant jusqu’à 500. Le projet vérifie le statut de fin de l’alignement. Le score seul en mode `high` utilise des fronts circulaires, tandis que le chemin `high` garde les fronts antérieurs nécessaires à la remontée. [Allocation](https://github.com/Marco-Sola-Lab/WFA2-lib/blob/bcf473a6561fa297934a80eaa7c04b4ee412360c/wavefront/wavefront_aligner.c#L424-L456) · [Fonctions de jonction et de récursion](https://github.com/Marco-Sola-Lab/WFA2-lib/blob/bcf473a6561fa297934a80eaa7c04b4ee412360c/wavefront/wavefront_bialign.c)

## Entrées et mesures

Les fenêtres GRCh38.p14 / T2T-CHM13v2.0 de 2 et 4 Mb partent de la même coordonnée GRCh38. Les deux loci de 250 kb sont disjoints et entièrement contenus dans la fenêtre de 4 Mb. Leurs bornes homologues proviennent de la chaîne UCSC 23, brin positif. Le pilote vérifie sur la paire de 4 Mb que les 41 bases entourant chaque borne sont identiques entre assemblages et que cette chaîne y apparaît une seule fois de chaque côté. Toutes les longueurs ci-dessous sont celles des chaînes effectivement alignées. Les coordonnées des deux assemblages sont en base zéro, borne droite exclue. [Origine des assemblages](../README.md) · [Coordonnées et SHA-256](838910/cases.csv)

| Cas | Longueur GRCh38 | Longueur CHM13 | Distance |
|---|---:|---:|---:|
| Locus 34,5 Mb | 250 000 | 250 017 | 455 |
| Locus 37,0 Mb | 250 000 | 249 997 | 538 |
| Fenêtre 2 Mb | 2 000 000 | 2 001 263 | 5 599 |
| Fenêtre 4 Mb | 4 000 000 | 4 001 764 | 12 482 |

ROMEO, job 838910. Chaque bloc permute les 24 configurations (quatre paires, trois algorithmes, score ou trace). Chacun des cinq blocs démarre un processus neuf par configuration et exécute les processus en série sur un cœur fixé. Chaque processus charge la paire, chauffe le backend par un appel complet, libère le résultat, puis chronomètre un second appel complet avec `steady_clock`. L’allocation de l’aligneur, son calcul et la conversion du CIGAR font partie de la mesure. Les lectures du FASTA, l’écriture de la sortie et la validation du chemin ont lieu hors chronomètre. Le pic du processus provient du champ Linux `VmHWM` après le calcul. Il inclut le runtime, les entrées, la chauffe et l’appel mesuré. Chaque processus a une limite d’espace d’adressage de 4 Gio et un délai de 150 s. La réservation Slurm plafonne le job à 8 Gio et 30 min.

Le lecteur FASTA de cette allocation correspond à la révision `a728923` (`src/io/fasta.cpp`, SHA-256 `28242b98ce45a329bdc94a54e7fec7cdde53a030ab407a171408466710e528f0`). Pour reproduire le périmètre mémoire historique dans une copie de travail dédiée, y restaurer ce fichier avec `git restore --source=a728923 -- src/io/fasta.cpp`.

Le pilote compare chaque score à un calcul edlib séparé, rejoue chaque CIGAR, vérifie les égalités et substitutions, la consommation des deux chaînes et le coût total. Pour estimer la dispersion de l’exécution sur ce nœud, le résumé donne la moyenne, l’écart-type échantillonnal et la demi-largeur de l’IC à 95 % de Student sur cinq processus ($t_4=2{,}776$). Les quatre fenêtres constituent un jeu de cas déterministe et se recouvrent en partie. [Mesures brutes](838910/trials.csv) · [Résumé](838910/summary.csv)

## Temps et mémoire observés

Pour le chemin complet, chaque cellule indique le temps moyen en millisecondes avec la demi-largeur de l’IC à 95 %, suivi du pic RSS moyen en Mio. L’absence de désaccord avec edlib porte sur les 120 scores. Les 60 chemins chronométrés ont passé le rejeu complet.

| Cas | BiWFA, ms / Mio | WFA2 unidirectionnel, ms / Mio | edlib, ms / Mio |
|---|---:|---:|---:|
| Locus 34,5 Mb | 2,028 ± 0,164 / 9,81 | 1,254 ± 0,299 / 8,24 | 99,169 ± 0,081 / 4,50 |
| Locus 37,0 Mb | 2,281 ± 0,071 / 10,70 | 1,452 ± 0,202 / 9,17 | 114,286 ± 0,096 / 4,50 |
| Fenêtre 2 Mb | 46,726 ± 0,774 / 31,30 | 71,498 ± 0,865 / 139,73 | 2 501,080 ± 1,005 / 22,89 |
| Fenêtre 4 Mb | 218,716 ± 1,238 / 54,02 | 375,899 ± 2,912 / 627,02 | 8 670,735 ± 1,645 / 41,63 |

Sur la fenêtre de 4 Mb avec trace, la moyenne géométrique des rapports appariés WFA2/BiWFA est de 1,719 pour le temps (IC à 95 % de 1,707 à 1,730) et de 11,61 pour le RSS (de 11,36 à 11,86). Le calcul prend le logarithme de chaque rapport par bloc, applique l’IC de Student à la moyenne des cinq logarithmes, puis reprend l’exponentielle. Sur les deux loci plus courts, WFA2 unidirectionnel termine avant BiWFA. Le calcul du score seul présente un autre compromis à 4 Mb. BiWFA y demande 113,300 ± 1,283 ms et 41,46 Mio, WFA2 unidirectionnel 211,943 ± 0,806 ms et 23,87 Mio. Le score seul de WFA2 garde les fronts nécessaires à l’étape suivante. Le contrôle de dispatch à 2 Mb compte 62 changements de sous-intervalles par chemin BiWFA et aucun pour son calcul du score seul. [Série complète](838910/summary.csv) · [Contrôle de dispatch](838910/dispatch.csv)

## Extension au bras long du chromosome 21

Le job 838918 compare quatre fenêtres emboîtées du bras q des mêmes assemblages. GRCh38 commence à 14 000 100 et CHM13 à 12 344 175. Chaque borne GRCh38 est couverte par un seul bloc de chaîne à cet emplacement. La chaîne 23 la relie à chr21 CHM13 sur le brin positif. Les 41 bases centrées sur chacune des cinq bornes sont identiques entre les deux séquences et uniques dans chacun des deux chromosomes entiers. Le script contrôle les SHA-256 des chromosomes, de la chaîne de correspondance et des FASTA extraits. Les `N` sont conservés et traités comme des symboles littéraux. [Manifeste du bras q](arm/838918/cases.csv) · [Extraction et contrôles](../../../bench/human/algorithms/arm.py)

| Taille GRCh38 | Fin GRCh38 exclue | Fin CHM13 exclue | Longueur CHM13 | N GRCh38 / CHM13 |
|---:|---:|---:|---:|---:|
| 8 Mb | 22 000 100 | 20 367 714 | 8 023 539 | 0 / 0 |
| 16 Mb | 30 000 100 | 28 366 012 | 16 021 837 | 0 / 0 |
| 24 Mb | 38 000 100 | 36 383 659 | 24 039 484 | 30 / 0 |
| 32 Mb | 46 000 100 | 44 375 522 | 32 031 347 | 50 160 / 0 |

Sur les 32 Mb, sept blocs contiennent des `N` dans GRCh38. Le plus long occupe [43 212 462, 43 262 462). Chacun des 50 160 `N` requiert une substitution ou une suppression puisque CHM13 n'en contient aucun. La différence de longueur impose en plus au moins 31 347 insertions. Les deux contraintes donnent une borne inférieure de 81 507 opérations sur les chaînes publiées.

Les quatre tailles démarrent simultanément sur quatre cœurs fixes d’un même EPYC. Chaque cœur traite une seule taille, et ses variantes se succèdent dans un ordre mélangé par bloc. Son travail se termine après le cinquième bloc. Les derniers essais à 32 Mb se déroulent donc avec moins de calculs concurrents que les premiers. Cinq processus neufs par variante effectuent une chauffe et un appel mesuré, avec le même chronomètre et le même contrôle du RSS que l’étude précédente. Le score BiWFA, le chemin BiWFA et le score WFA2 unidirectionnel couvrent tous les cas. edlib vérifie séparément le score et restitue aussi un chemin à 8 Mb. Les deux variantes WFA partagent la bibliothèque WFA2 et utilisent des directions d’exploration différentes. Les CIGAR sont entièrement rejoués et chaque coût est confronté à l’oracle WFA2 unidirectionnel. Le délai de 720 s et le plafond de 4 Gio d’espace d’adressage s’appliquent à chaque processus, chauffe et appel mesuré compris. Les limites Slurm sont 24 Gio et 90 min pour l’allocation. [Mesures du bras q](arm/838918/trials.csv) · [Moyennes, écarts-types et IC](arm/838918/summary.csv)

Sur les 70 processus mesurés, les 25 chemins ont passé le rejeu intégral. À 32 millions de bases, BiWFA calcule le score en 45,623 ± 0,224 s et restitue le chemin en 89,775 ± 0,108 s (moyenne ± écart-type). Le RSS moyen atteint respectivement 272,83 et 311,17 Mio. WFA2 unidirectionnel calcule le score en 102,038 ± 0,906 s avec un RSS moyen de 143,81 Mio. Les temps portent sur l’appel complet, la lecture du FASTA ayant lieu avant l’échauffement. [Mesures et limites Slurm](arm/838918/accounting.psv) · [Tous les cas](arm/838918/summary.csv)

## Reproduction

```sh
python3 scripts/datasets/human_reference.py --out data/human --manifest results/human/manifest.tsv
sbatch --account=r260105 --partition=instant bench/human/algorithms/human.sbatch
```

Pour un nœud de calcul hors ligne, placer les versions épinglées de WFA2 et edlib dans `deps/wfa2.tar.gz` et `deps/edlib.tar.gz`. Le script compile le projet en mode Release avec GCC et `-march=native`, lance CTest, puis effectue un contrôle de dispatch par interposition GNU ld sur un binaire séparé. La série temporelle utilise le binaire sans interposition. [Allocation et compilation](../../../bench/human/algorithms/human.sbatch) · [Pilote des cas et calcul statistique](../../../bench/human/algorithms/study.py) · [Empreintes sources](838910/sources.sha256) · [Compteurs d’appels](838910/dispatch.csv)

Pour reproduire l’extension du bras q, conserver la chaîne UCSC dans le cache local.

```sh
mkdir -p data/human/reference-cache
python3 scripts/datasets/human_reference.py --whole-assemblies --cache data/human/reference-cache
CHPS_CHAIN=data/human/reference-cache/hg38ToHs1.over.chain.gz sbatch --account=r260105 --partition=short bench/human/algorithms/arm.sbatch
```

Après récupération des dossiers de résultats, `python3 bench/human/algorithms/combine.py` produit le [tableau commun](summary.csv), avec le numéro de job et le protocole d’exécution de chaque série.
