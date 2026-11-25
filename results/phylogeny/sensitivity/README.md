# Largeur de bande et stabilité des arbres

Les tables comparent UPGMA et Neighbor-Joining à la topologie produite par la même méthode sur les distances exactes. Les 64 accessions versionnées du panel public ont été relues en mémoire. Les longueurs correspondent au manifeste et les 2 016 distances indépendantes d’edlib à la matrice publiée. Le calcul conserve les symboles ambigus et compare les caractères à coût unitaire. Seuls des dénombrements et des distances entre topologies sont conservés. [Manifeste](../../ncbi/sars_metadata.tsv) · [Matrice de référence](../sars-edit/matrix.tsv) · [Calcul des arbres](../../../scripts/phylogeny/tree.py)

## Lecture des résultats

Dans `scores.csv`, `width` est la demi-largeur de base. La règle `fixed` retient les cellules pour lesquelles l’écart absolu entre indices vaut au plus cette largeur. Une paire dont l’écart de longueur dépasse la largeur est inaccessible. Les colonnes `overestimated`, `mean_excess` et `max_excess` portent seulement sur les paires accessibles. La moyenne inclut les erreurs nulles. Avec `endpoint_floor`, la largeur effective d’une paire est le maximum de la largeur de base et de son écart de longueur. La colonne `widened` compte les élargissements. Le témoin `equal_length` utilise la longueur la plus fréquente du panel, 29 580 bases (neuf feuilles). [Scores du panel](panel/scores.csv)

Dans `trees.csv`, un arbre n’est construit que si toutes les paires de son ensemble de feuilles ont un score. Le RF compte les partitions internes présentes dans un seul des deux arbres. UPGMA conserve l’enracinement. Neighbor-Joining utilise les bipartitions sans racine. Le dénominateur de `normalized_rf` est la somme des nombres de partitions. L’ordre des feuilles est identique pour les deux calculs et la comparaison conserve les arêtes de longueur nulle. Les RF concernent les topologies effectivement rendues, y compris le départage des ex æquo. Sur les 64 génomes, une demi-largeur fixe de 323 laisse deux surestimations et un RF nul pour les deux méthodes. Avec les demi-largeurs ajustées par paire, le RF brut de Neighbor-Joining passe de 102 à 104 entre les demi-largeurs de base 0 et 4. Une demi-largeur fixe de 512 restitue tous les scores. [Arbres du panel](panel/trees.csv) · [Définition de Robinson et Foulds](https://doi.org/10.1016/0025-5564%2881%2990043-2)

Les témoins comprennent huit jeux de 32 chaînes de 384 chiffres de l’alphabet 0 à 3, issus des graines 2026 à 2033. Une construction binaire de profondeur cinq introduit, à chaque branche, quatre substitutions et un couple suppression/insertion de même longueur, tirée entre 1 et 24. L’ordre des chaînes est ensuite mélangé. Pour chaque jeu, l’arbre de référence est reconstruit à partir de ses propres scores exacts. L’unité statistique est le jeu entier. La cohorte publique constitue un panel fixe. [Scores synthétiques](synthetic/scores.csv) · [RF synthétiques](synthetic/trees.csv)

## Reproduction

Depuis la racine d’un dépôt configuré avec CMake dans `build`, compiler la petite bibliothèque de calcul avec la version d’edlib déjà téléchargée. La commande suivante fonctionne sur Linux et macOS avec un compilateur C++20.

```sh
mkdir -p build/sensitivity
c++ -std=c++20 -O3 -fPIC -shared scripts/phylogeny/band_cost.cpp \
  build/_deps/edlib-src/edlib/src/edlib.cpp -Ibuild/_deps/edlib-src/edlib/include \
  -o build/sensitivity/band_cost.so
python3 scripts/phylogeny/band_sensitivity.py --native build/sensitivity/band_cost.so \
  --synthetic 2026 --replicates 8 --study neutral --widths 0 2 4 8 16 32 64 \
  --workers 3 --out build/sensitivity/synthetic
```

Le même programme accepte `--fasta fichier.fa` ou `--fasta -` pour l’entrée standard, et `--reference matrice.tsv` pour vérifier chaque distance exacte. Le panel a utilisé les demi-largeurs 0, 4, 16, 32, 64, 128, 256, 323, 512 et 1 024. Un chemin optimal d’edlib qui reste entièrement dans une bande suffit à prouver l’égalité des scores. Les autres cases du tableau sont calculées sans saturation par programmation dynamique. L’élargissement ne peut qu’abaisser le coût d’un chemin contraint, jusqu’au coût exact. [Calcul en bande](../../../scripts/phylogeny/band_cost.cpp) · [Comparaison des topologies](../../../scripts/phylogeny/band_sensitivity.py)
