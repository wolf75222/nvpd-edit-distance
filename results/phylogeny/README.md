# Distances d’édition et arbres

Chaque matrice d’édition contient les 2 016 distances de Levenshtein globales entre 64 entrées. La substitution, l’insertion et la suppression coûtent une unité. Le calcul compare les caractères IUPAC littéralement, dans l’ordre déposé, extrémités comprises. Une distance dépend donc aussi des bases ambiguës et des extrémités couvertes.

## SARS-CoV-2

Le manifeste provient de NCBI Datasets. Les identifiants versionnés et les longueurs des 64 enregistrements publics lus dans NCBI Nucleotide correspondaient à ce manifeste le 24 novembre 2026. Les entrées sont classées complètes par NCBI et totalisent 1 905 372 caractères, dont 65 588 hors A/C/G/T (3,44 %). Les entrées ont alimenté en mémoire le comparateur générique du projet. Les matrices calculées par Edlib, WFA2 et NvPD concordent octet par octet. La vérification s’est déroulée sur un Apple M1 Pro à huit cœurs, avec GCC 16.1.0 en configuration Release (`-O3`). Le fichier de vérification contient les empreintes SHA-256 de la matrice et du binaire ainsi que les durées murales des passages de vérification. [Manifeste NCBI](../ncbi/sars_metadata.tsv) · [Matrice](sars-edit/matrix.tsv) · [Vérification](sars-edit/verification.csv).

Le passage dense de NvPD couvre 1 786 858 049 040 cellules de programmation dynamique. Une première exécution a réparti les colonnes d’une paire à la fois sur huit threads. Ce parcours prévoit 59 899 758 synchronisations de fin de ligne pour la cohorte entière. Le processus a atteint la limite de 1 200 s avant d’émettre sa matrice. Le second passage a distribué les paires entre six travailleurs. Chacun a lancé le même binaire en mode paire avec un thread NvPD et vérifié son score. Les 2 016 invocations et la collecte des scores ont pris 268,828 s au total. Edlib et WFA2 ont utilisé un processus et quatre threads chacun. Dans le fichier de vérification, `threads` indique le réglage par processus et `processes` le nombre maximal de processus de calcul simultanés. [Tentative plafonnée](sars-edit/censored.csv) · [Trois vérifications achevées](sars-edit/verification.csv).

La référence Wuhan-Hu-1 et les neuf entrées de chacun des sept lignages participent à la construction des deux arbres. Sur les 63 entrées avec un groupe de neuf, la médiane vaut 112 au sein d’un lignage (252 paires, de 1 à 8 793) et 402 entre lignages (1 701 paires, de 77 à 12 229). Les 63 comparaisons à la référence figurent dans la matrice. L’échantillonnage de taille fixe et le traitement littéral des ambiguïtés limitent la portée biologique des écarts. [Distances par groupe](sars-edit/group_distances.csv) · [Statistiques regroupées](sars-edit/pooled_distances.csv).

Une arête de longueur strictement positive isole les neuf B.1.617.2 et une autre les neuf P.1 dans chaque arbre. NJ isole aussi les neuf B.1.1.529. Le critère est un clade descendant pour l’arbre enraciné UPGMA et l’un ou l’autre côté d’une bipartition pour NJ. L’ordre des feuilles dans la matrice départage les égalités de l’algorithme. Les groupes désignent ici les étiquettes Pango exactes et excluent leurs sous-lignages. Des étiquettes ancestrales telles que B.1 peuvent donc correspondre à des groupes paraphylétiques dans cet échantillon. [Arbre UPGMA](sars-edit/upgma.nwk) · [Arbre NJ](sars-edit/nj.nwk) · [Groupes et tailles des clades englobants](sars-edit/group_clades.csv).

### Comparaison avec une phylogénie publique

Le Newick public UCSC du 23 novembre 2026 contient 8 472 515 feuilles. La restriction conserve 52 accessions identifiées par un champ exact, dont la référence sous son synonyme GenBank MN908947.3. Elle contracte les nœuds à un seul descendant en additionnant leurs longueurs. Six entrées du groupe A, chacune des neuf entrées des groupes B.1, B.1.1.7, B.1.351, B.1.617.2 et P.1, ainsi que la référence ont un correspondant. Le groupe B.1.1.529 en a zéro. L’empreinte SHA-256 du Newick source compressé est `7c8e49e7a8ddcdc887dacbf82428a8bd75d2dcd093f07eea7ade047109c14d95`. [Index UCSC et conventions](https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) · [Version datée](https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/2026/11/23/public-2026-11-23.version.txt) · [Couverture](sars-surveillance/coverage.csv).

Dans cet arbre restreint, cinq des six groupes représentés par au moins deux feuilles forment un clade d’arête positive. Le sixième, B.1, compte neuf entrées dispersées dans un clade de 45 feuilles. Pour comparer les bipartitions, on reconstruit UPGMA et NJ sur la sous-matrice d’édition des mêmes 52 feuilles, puis on considère chaque arbre comme non enraciné. Sur les 27 bipartitions internes positives de l’arbre public, UPGMA en retrouve 6 et NJ 8. L’arbre de surveillance contient des polytomies. Une bipartition supplémentaire du nouvel arbre peut résoudre une polytomie sans la contredire. [Newick restreint](sars-surveillance/usher.nwk) · [Groupes](sars-surveillance/group_clades.csv) · [Bipartitions](sars-surveillance/shared_splits.csv).

La distance patristique publique additionne les longueurs des branches sur le chemin entre deux feuilles. Les branches UShER s’expriment par défaut en nombres d’événements assignés par parcimonie. Sur les 1 326 paires communes, le coefficient de Pearson entre ce chemin et la distance d’édition vaut 0,124. Ces paires partagent leurs feuilles et ce coefficient est descriptif. [Méthode UShER](https://doi.org/10.1038/s41588-021-00862-7) · [Distances patristiques](sars-surveillance/patristic.tsv) · [Comparaison](sars-surveillance/overlap.csv).

## Forme des distances

Pour chaque triplet, on calcule l’écart entre les deux plus grandes distances. Pour chaque quartet, on applique la même opération aux trois sommes de deux distances entre paires disjointes. Une matrice ultramétrique donne un écart nul pour chaque triplet. Une distance d’arbre additive respecte la condition des quatre points pour chaque quartet. Les tests portent sur toutes les combinaisons. L’inégalité triangulaire ordinaire est satisfaite pour les 41 664 triplets de chaque matrice. [Comptages SARS-CoV-2](sars-edit/gaps.csv) · [Comptages φX174](phix174/gaps.csv).

| Cohorte | Triplets d’écart positif / 41 664 | Quartets d’écart positif / 635 376 | RMSE UPGMA | RMSE NJ |
|---|---:|---:|---:|---:|
| SARS-CoV-2 | 40 747 | 547 882 | 429,30 | 281,62 |
| Phage φX174 | 40 003 | 389 529 | 14,89 | 3,56 |

La RMSE compare la distance d’édition au chemin inféré entre chaque paire, en unités du score. Les fichiers `pair_fit.csv` donnent les deux chemins pour chaque paire, les fichiers `fit.csv` les résidus globaux. Dans `branches.csv`, la RMSE d’une branche porte sur les paires dont le chemin traverse cette branche. Une même paire contribue ainsi à plusieurs branches. L’ajustement en moyenne laisse subsister des écarts importants, même pour NJ. [Résidus SARS-CoV-2](sars-edit/fit.csv) · [Résidus φX174](phix174/fit.csv).

Chez φX174, une seule paire a une distance nulle entre accessions distinctes, HM753660.1 et HM753661.1. Leurs empreintes de séquence sont identiques. L’ensemble présente trois arêtes nulles dans chaque arbre. En retenant une feuille par classe de distance nulle, 38 275 des 39 711 triplets et 368 769 des 595 665 quartets violent encore les conditions correspondantes. Les écarts maximaux restent à 60 pour les triplets et 61 pour les quartets. [Manifeste φX174](../ncbi/phix174.tsv) · [Matrice mesurée sur ROMEO](../romeo/cohort-837271/matrix.tsv) · [Sensibilité aux doublons](phix174/duplicate_sensitivity.csv).

Sur ce phage, la RMSE des paires traversant une branche varie de 9,62 à 26,63 pour UPGMA et de 0,98 à 11,98 pour NJ. Les maxima correspondent aux arêtes terminales de AF299312.1 pour UPGMA et de AF299307.1 pour NJ, avec 63 paires dans chaque cas. [Diagnostics par branche](phix174/branches.csv).

## Recalcul à partir des scores publiés

```sh
python3 scripts/phylo_fit.py results/romeo/cohort-837271/matrix.tsv --out results/phylogeny/phix174
python3 scripts/phylo_fit.py results/phylogeny/sars-edit/matrix.tsv --out results/phylogeny/sars-edit
python3 scripts/phylo_groups.py results/phylogeny/sars-edit/matrix.tsv results/ncbi/sars_metadata.tsv \
  --reference-group reference --out results/phylogeny/sars-edit \
  --comparison-tree results/phylogeny/sars-surveillance/usher.nwk \
  --comparison-out results/phylogeny/sars-surveillance
```

Les distances de chemin de chaque arbre publié ont aussi été recalculées avec Bio.Phylo. Les 1 326 chemins du Newick UCSC concordent exactement. Pour les quatre arbres reconstruits, l’écart d’arrondi maximal sur les chemins exportés est inférieur à 0,00000005 unité.

Pour reproduire la restriction depuis le Newick public `public-2026-11-23.all.nwk.gz` indiqué dans l’index UCSC, utiliser la commande suivante à la racine du projet.

```sh
python3 scripts/restrict_newick.py public-2026-11-23.all.nwk.gz results/ncbi/sars_metadata.tsv \
  --alias MN908947.3=NC_045512.2 --out results/phylogeny/sars-surveillance/usher.nwk
```
