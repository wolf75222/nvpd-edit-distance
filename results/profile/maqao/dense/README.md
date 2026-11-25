# Profil des lignes de NvPD sur AMD EPYC

Les calculs utilisent deux paires de fenêtres du chromosome 21. Chaque paire comprend une séquence de GRCh38 et une de CHM13. Les premières coordonnées des fenêtres GRCh38 sont 32 000 000 et 22 500 000, en base zéro. À 100 000 bases, leurs distances d’édition valent respectivement 3 345 et 51 331. Le manifeste identifie les séquences, leur longueur et le score calculé par edlib. [Manifeste](inputs.csv)

## Protocole

La sélection des noyaux utilise le nœud `romeo-c022` et la qualification indépendante le nœud `romeo-c023`. Tous deux sont équipés de processeurs AMD EPYC 9654. Le pilote lit le FASTA puis chronomètre l’appel à `nvpd_omp` avec une horloge monotone. Le temps inclut la construction de l’index et le calcul des lignes. Le pilote compare chaque score à celui d’edlib.

Chaque processus effectue un échauffement avant la mesure. Un bloc mesure toutes les variantes sur les mêmes entrées, avec un ou quatre cœurs. L’ordre des variantes tourne entre les blocs et l’ordre des configurations change à chaque bloc.

Le dossier `screen/` contient dix blocs par paire et par effectif pour les longueurs de 16 384 et 32 768 bases. La sélection maximise la moyenne des logarithmes des gains appariés à 32 768 bases, avec le même poids pour les deux paires. Huit blocs supplémentaires comparent le candidat retenu à la référence sur les préfixes de 100 000 bases des mêmes fenêtres.

Le dossier `production/` contient les mesures du point d’entrée public avec la sélection du noyau à l’exécution. Pour les préfixes de 32 à 4 096 bases, chaque temps conservé est la moyenne de douze appels consécutifs après l’échauffement.

ONEVIEW échantillonne les cycles de la fenêtre proche à 32 768 bases. CQA prévoit le nombre de cycles par itération de la boucle assembleur pour des données en L1. La boucle de référence met à jour huit cellules à chaque itération et la boucle AVX-512 en met à jour seize. Les temps de `screen/` et `production/` proviennent d’exécutions séparées du profileur. [Guide ONEVIEW](https://www.maqao.org/documentation.html)

Les horaires de début et de fin figurent dans le relevé Slurm. Les jobs `840152` et `840171` se sont arrêtés avant les mesures. [Relevé Slurm](accounting.tsv)

Le premier profil (job `840153`) contient seize appels par effectif. Ses boucles portent les noms `threads-1` et `threads-4` dans l’export. Le job `840172` répète 136 appels sur un cœur pour la référence (`baseline-1`) et pour le candidat (`blocked-1`). Le job `840174` reprend les 136 appels sur la version de production (`production-1`). Le job `840180` mesure les deux noyaux à quatre cœurs sur 437 appels (`original-4` et `production-4`). [Boucles et diagnostics](hotloops.csv)

## Reproduction

Depuis la racine du dépôt, construire `editdist` avec edlib comme indiqué dans le README du projet. Les deux premières commandes préparent les fenêtres et leur score. La dernière commande reconstruit les variantes à partir du noyau antérieur à l’optimisation.

```sh
python3 bench/hybrid/matrix/prepare.py
python3 bench/cpu/maqao/inputs.py build/hybrid-matrix/inputs build/editdist build/maqao/inputs
git worktree add --detach build/maqao/base ba97a3d
python3 bench/cpu/maqao/prepare_variants.py build/maqao/base build/maqao/candidates
```

Les profils de cette campagne utilisent MAQAO 2025.1.4 sur EPYC. Pour les reproduire sur ROMEO, charger cette version avant la soumission. Les scripts utilisent le programme `maqao` du `PATH`. La variable `MAQAO` permet de fournir un chemin absolu. Le premier job produit un classement des cinq variantes, puis un rapport HTML pour la référence et le candidat retenu. Le second compile le noyau de production en modes natif et générique. Il teste les scores et la reconstruction de l’alignement, puis profile et chronomètre le noyau.

```sh
romeo_load_x64cpu_env
spack load maqao@2025.1.4
export MAQAO="$(command -v maqao)"
ab=$(sbatch --parsable -A "$ROMEO_PROJECT" -p instant bench/cpu/maqao/compare.sbatch \
  build/maqao/base build/maqao/candidates build/maqao/inputs build/maqao/screen)
sbatch -A "$ROMEO_PROJECT" -p instant --dependency=afterok:"$ab" bench/cpu/maqao/confirm.sbatch \
  build/maqao/base . build/maqao/inputs build/maqao/screen build/maqao/candidates build/maqao/final
```

Chaque archive `VARIANTE-FILS.tar.gz` contient le rapport HTML sous `RESULTS/`, les exports LProf dans `shared/run_0/` et les données CQA dans `static_data/`. La commande suivante extrait le diagnostic principal d’un rapport décompressé.

```sh
mkdir -p build/maqao/production-4
tar -xf build/maqao/final/production-4.tar.gz -C build/maqao/production-4
python3 bench/cpu/maqao/hotloops.py build/maqao/production-4
```

## Indications du second profil

CQA relève la conversion des indices sur un octet et des transferts non alignés dans la boucle par blocs. Il suggère aussi de dérouler la boucle. La seconde campagne évalue chaque indication sur le noyau du commit `50d9bd9`.

La variante `index32` stocke seize indices de 32 bits avant le calcul. Chaque descripteur occupe 68 octets contre 20 dans le noyau de départ. La variante `aligned` aligne l’écriture de la ligne et la lecture du terme supérieur sur 64 octets. Le terme diagonal commence à l’entier précédent. Les variantes `unroll2` et `unroll4` déroulent respectivement deux et quatre blocs par itération.

```sh
git worktree add --detach build/maqao/first-blocks 50d9bd9
python3 bench/cpu/maqao/prepare_followups.py build/maqao/first-blocks/src/cpu build/maqao/followups
sbatch -A "$ROMEO_PROJECT" -p instant bench/cpu/maqao/follow.sbatch \
  build/maqao/base build/maqao/first-blocks build/maqao/followups \
  build/maqao/inputs build/maqao/candidates build/maqao/follow
```

Dans les mesures de cette campagne, `baseline` désigne le premier noyau par blocs et `original` la référence antérieure à MAQAO.

Une comparaison supplémentaire porte sur des longueurs impaires de 24 001 à 80 009 bases. Les deux versions du noyau travaillent sur les mêmes préfixes à un et quatre cœurs. Changer la longueur déplace le début de la seconde ligne et modifie le nombre de cellules dans le dernier bloc. Cette série mesure la sensibilité à la longueur sans isoler l’effet des ensembles de cache.

## Qualification de la disposition des lignes

Le job `840193` utilise un second nœud EPYC 9654, `romeo-c023`. NvPD place la séquence la plus longue en colonnes quand seul le score est demandé. Les rectangles et leurs transposés reprennent les mêmes enregistrements dans l’ordre inverse. Le manifeste donne les deux longueurs et les coordonnées de chaque séquence. [Matrices de qualification](qualification/inputs.csv)

La variante `layout` aligne les lignes dès que les deux séquences comptent au moins 40 000 bases. Le noyau SIMD reste commun aux deux dispositions. La variante `index-pitch` ajoute du remplissage entre les tables de prédécesseurs associées aux bases. Leur pas contient alors un nombre impair de lignes de cache de 64 octets. La comparaison des largeurs voisines de 65 536 bases teste l’effet de ce remplissage à pas de matrice fixé. ONEVIEW profile à 95 003 bases le noyau `50d9bd9` et le candidat `layout`, sur 40 appels à un cœur et 120 appels à quatre cœurs.

```sh
python3 bench/cpu/maqao/prepare_followups.py --layout build/maqao/first-blocks/src/cpu build/maqao/layout-candidates
python3 bench/cpu/maqao/qualify_inputs.py build/hybrid-matrix/cache build/hybrid-matrix/inputs/manifest.tsv \
  build/editdist build/maqao/qualification
sbatch -A "$ROMEO_PROJECT" -p instant bench/cpu/maqao/qualify.sbatch \
  build/maqao/base build/maqao/first-blocks build/maqao/layout-candidates/layout \
  build/maqao/qualification build/maqao/candidates build/maqao/qualified build/maqao/layout-candidates/index-pitch
```

Sur le second nœud, le remplissage de l’index donne un gain de 0,833 à 65 536 bases et de 1,350 à 65 537 bases, sur quatre cœurs. Ces valeurs sont les moyennes géométriques des rapports de temps appariés, avec le même poids pour chaque paire. La variation d’une base suffit ici à inverser le classement. Le noyau public conserve les indices sur un octet et la disposition des lignes du commit `50d9bd9`. [Mesures des largeurs](qualification/widths.csv)

À 131 072 bases et quatre cœurs, la paire proche passe de 1,940 à 0,497 s, pour une distance de 3 391 éditions. Le gain moyen des neuf rapports appariés vaut 3,901, avec un intervalle de Student à 95 % de [3,850 ; 3,952]. Sur la paire distante de 66 919 éditions, le temps passe de 1,908 à 0,498 s. [Mesures à 131 072 bases](qualification/large-four.csv)
