# NvPD · Distance d’édition parallèle

`editdist` calcule la distance de Levenshtein globale entre deux séquences d’ADN. Chaque insertion, suppression ou substitution coûte une unité. NvPD repère la dernière occurrence de chaque symbole dans la séquence cible. Les colonnes d’une même ligne se calculent alors en parallèle. Le dépôt réunit les implémentations CPU et GPU, leurs mesures sur ROMEO et les outils qui construisent un arbre à partir des distances. [Article de Sadiq et coll.](https://doi.org/10.1007/s10586-019-02962-w)

[Rapport complet](report/rapport.pdf) · [Synthèse courte](report/tldr.pdf) · [Carnet exécuté](notebooks/profil_nvpd.ipynb) · [Bibliographie](report/tex/references.bib) · [Planche des figures](figures/montage.png) · [Résultats par exigence du sujet](results/README.md)

Sur les douze paires synthétiques chronométrées à chaud, NvPD sur GPU devance les deux bibliothèques sur une seule paire. Elle compte un million de bases de référence et 10 % d’éditions injectées. Sur cinq répétitions, NvPD prend 3,591 ± 0,011 s sur le GPU du GH200. Sur un cœur Grace du même nœud, Edlib prend 5,903 ± 0,060 s et WFA2-lib 9,835 ± 0,023 s (moyenne ± écart-type). [Mesures sur GH200](results/repeated/gh200/summary.csv)

## Compiler et tester sur CPU

La compilation requiert CMake 3.22 ou ultérieur, un compilateur C et un compilateur C++20 compatible avec OpenMP. Les exemples de données et d’arbres utilisent Python 3.

Sur Debian ou Ubuntu, installer les outils de compilation et Python.

```sh
sudo apt update
sudo apt install git build-essential cmake python3
```

Sur macOS avec Homebrew, installer les outils puis sélectionner GCC avant de configurer CMake. Les noms ci-dessous correspondent à GCC 16. Adapter le suffixe si une autre version est installée.

```sh
brew install git cmake gcc python
export CC=gcc-16 CXX=g++-16
```

Cloner le dépôt, compiler en mode Release, puis lancer les tests.

```sh
git clone https://github.com/wolf75222/nvpd-edit-distance.git
cd nvpd-edit-distance
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

Pour contrôler les sources du projet, ajouter `-DCHPS_WARNINGS_AS_ERRORS=ON` à la commande de configuration. La CI active cette option pour C, C++ et les backends OpenMP, MPI et Kokkos.

CMake télécharge Edlib et WFA2-lib à la première configuration et vérifie leurs empreintes SHA-256. Pour compiler hors ligne, placer les archives dans `deps/edlib.tar.gz` et `deps/wfa2.tar.gz`. [Révisions et empreintes](CMakeLists.txt)

## Aligner deux séquences

Les deux séquences de cet exemple diffèrent par une substitution à la deuxième position.

```sh
printf '>alpha\nACGT\n>beta\nAGGT\n' > build/paire.fa
build/editdist pair --algorithm nvpd-omp --threads 4 --input build/paire.fa
build/editdist pair --algorithm hirschberg --trace --input build/paire.fa
```

Les deux commandes écrivent un objet JSON avec `"distance":1`. La seconde ajoute `"cigar":"1=1X2="`, soit une correspondance, une substitution et deux correspondances. Le champ `seconds` mesure l’alignement après lecture du FASTA. `peak_rss_kib` donne le pic de mémoire du processus en Kio et `device_bytes` la taille des tampons explicites sur le GPU.

Le mode `pair` attend deux entrées FASTA. Le premier mot de chaque en-tête sert d’identifiant et doit être unique dans le fichier. Le programme ignore la casse et compare les symboles IUPAC littéralement. Ainsi, `N` contre `N` coûte zéro et `N` contre `A` coûte un.

Pour un FASTA à lignes de largeur régulière, installer `samtools` avec APT ou Homebrew et créer l’index avec `faidx`. L’option `--fai` utilise les longueurs indexées pour réserver la mémoire des séquences avant la lecture.

```sh
samtools faidx build/paire.fa
build/editdist pair --algorithm biwfa --trace --input build/paire.fa --fai build/paire.fa.fai
```

## Construire une matrice et deux arbres

Ajouter une troisième séquence au lot, puis construire les arbres UPGMA et Neighbor-Joining à partir de la même matrice.

```sh
cp build/paire.fa build/lot.fa
printf '>gamma\nTTTT\n' >> build/lot.fa
build/editdist matrix --algorithm edlib --input build/lot.fa > build/distances.tsv
python3 scripts/phylogeny/tree.py build/distances.tsv --method upgma --out build/upgma.nwk
python3 scripts/phylogeny/tree.py build/distances.tsv --method nj --out build/nj.nwk
```

Le TSV conserve les identifiants en ligne et en colonne. Ici, `gamma` est à trois éditions de chacune des deux autres séquences. UPGMA regroupe d’abord `alpha` et `beta`.

La cohorte SARS-CoV-2 réunit Wuhan-Hu-1 et neuf génomes complets de chacun des sept lignages retenus, soit 64 accessions versionnées. Avec la commande NCBI `datasets` installée, `python3 scripts/datasets/ncbi_sars.py` reconstitue leur FASTA dans `data/sars/cohort.fa`. Le dépôt contient les 2 016 distances et les deux arbres. Une bande fixe de demi-largeur 323 surestime deux distances de 5 et 19 éditions. La distance de Robinson-Foulds entre l’arbre contraint et l’arbre de référence vaut zéro pour UPGMA comme pour Neighbor-Joining. [Données et manifestes](data/README.md) · [Matrices des cohortes SARS-CoV-2 et φX174](results/phylogeny/README.md) · [Effet de la bande](results/phylogeny/sensitivity/README.md)

Sur 32 cœurs EPYC, subdiviser les rangées de cette matrice ramène Edlib de 1,645 à 0,343 s et WFA2 de 3,616 à 0,775 s. Pour l’hybride NvPD/BiWFA, un filtre sur les effectifs de lettres réduit le temps du noyau de 584 à 448 ms. WFA-GPU, l’aligneur affine sur CUDA d’Aguado-Puig et coll., confie 15 paires de cette cohorte au GPU et 2 001 à WFA2 sur CPU. Avec quatre travailleurs CPU, ce partage CPU/GPU prend 6,844 s, contre 6,874 s pour WFA2 sur CPU. L’intervalle à 95 % du rapport des temps couvre 1. [Temps et dispersion sur ROMEO](results/phylogeny/sars-romeo/README.md) · [Ablation du filtre](results/phylogeny/sars-ablation/README.md) · [WFA-GPU](results/wfa-gpu/README.md)

## Générer des paires et mesurer le temps

Générer une référence de 5 000 bases, puis appliquer 5 % d’éditions. Le manifeste TSV conserve le chemin du FASTA, la graine, les longueurs et les effectifs de chaque opération.

```sh
python3 scripts/datasets/synthetic.py --lengths 5000 --rates 0.05 --pairs 1 --out build/synthetique
python3 bench/campaign/run.py --binary build/editdist --manifest build/synthetique/manifest.tsv \
  --target cpu --algorithms seq nvpd-omp --threads 1 4 --trials 3 --output build/mesures.csv
python3 bench/campaign/summarize.py build/mesures.csv > build/resume.csv
```

Avec la graine par défaut (2026), les neuf mesures donnent une distance de 245 pour 250 éditions injectées. Pour répéter le banc, choisir un nouveau nom pour `--output`.

Sans option, le générateur écrit 60 paires dans `data/synthetic/`. Les références de 5 000, 20 000, 100 000 et 1 000 000 nt suivent une loi uniforme sur A, C, G et T. Pour chaque taux (2, 5 et 10 %), il produit cinq paires par longueur. Il tire les positions sans remise et choisit indépendamment chaque type d’édition avec probabilité 1/3. À 1 Mnt, le calcul séquentiel utilise la première paire de chaque taux. La commande `python3 bench/campaign/check_course_corpus.py` régénère les FASTA et vérifie leurs empreintes et les scores archivés. [Manifeste synthétique](data/synthetic/manifest.tsv) · [Couverture des campagnes ROMEO](results/romeo/README.md)

Le résumé de cette campagne donne le temps médian et l’accélération par rapport au séquentiel sur la même paire. Pour une méthode dense, le débit en GCUPS vaut $mn/(10^9 T)$, avec $T$ en secondes. Les mesures répétées du rapport utilisent des processus préchauffés et publient les moyennes, les écarts-types et les intervalles de confiance. [Soumission Slurm sur ROMEO](bench/README.md) · [Mesures répétées](results/repeated/README.md)

## Aligner des chromosomes humains

Les assemblages HG002 v1.1 et v1.2 fournissent quatre paires de chromosomes entiers, de 44 à 252 millions de bases. Le script extrait les chromosomes du FASTA et applique les variants du VCF à la première version pour vérifier qu’il retrouve la seconde. Les commandes ci-dessous calculent la distance entre les deux versions de chaque haplotype.

```sh
python3 bench/cpu/chromosome-dispatch/prepare.py data/human/hg002
for chrom in chr1_MATERNAL chr1_PATERNAL chr21_MATERNAL chr21_PATERNAL; do
  cat "data/human/hg002/${chrom}.v1.1.fa" "data/human/hg002/${chrom}.v1.2.fa" > "build/${chrom}.fa"
  build/editdist pair --algorithm nvpd-sparse --threads 8 --input "build/${chrom}.fa"
done
```

Sur le chromosome 1, les distances maternelle et paternelle valent respectivement 54 et 111 éditions. Sur le chromosome 21, elles valent 6 et 15. Sur un cœur EPYC 9654, `nvpd-sparse` traite ces révisions 3,03 à 4,05 fois plus vite que WFA2. Sur le chromosome 1 maternel, le score prend 28,19 ± 0,06 ms et le pic de mémoire résidente du processus vaut en moyenne 477 Mio. Le chromosome 21 paternel contient 2 300 690 symboles N dans chaque version. [Données HG002](https://github.com/marbl/HG002) · [Cinq blocs par configuration](results/human/chromosome-dispatch/)

La comparaison des révisions v1.0 et v1.1 de CHM13 porte sur 18 chromosomes entiers. Sur le chromosome 1, le chemin d’alignement compte 176 éditions. Le calcul par fragments prend 11,81 ms à huit cœurs, contre 183,11 ms pour WFA2 sur un cœur. Le banc de pages mémoire distingue le coût d’un alignement isolé de celui d’une série qui réutilise les séquences. [CIGAR chromosomiques](results/trace/chromosomes/) · [Pages de 4 Kio et 2 Mio](results/profile/huge-pages/839316/)

```sh
python3 scripts/datasets/chm13_panel.py --cache build/t2t-cache --out data/human/chm13-panel
build/editdist pair --algorithm nvpd-split --threads 8 --trace --input data/human/chm13-panel/chr1.fa
```

## Comparer une cohorte de 1000 Genomes

Le script reconstruit les deux haplotypes de 512 individus de la phase 3 de 1000 Genomes dans une fenêtre de 32 768 bases du chromosome 21 (GRCh37). Les 1 024 haplotypes forment 523 776 paires. Sur cette fenêtre, 386 séquences sont distinctes. L’option `--unique` aligne les 74 305 paires distinctes et restitue la matrice dans l’ordre initial, avec les 1 024 identifiants. Pour $u$ séquences à comparer, la commande `editdist matrix` stocke $u(u-1)/2$ scores et écrit les deux moitiés du tableau TSV. [Manifeste et protocole](results/human/packed-cohort/)

```sh
python3.12 -m venv build/phase3-env
build/phase3-env/bin/python -m pip install pysam==0.24.1
build/phase3-env/bin/python bench/cpu/packed-cohort/prepare.py data/human/phase3-window
build/editdist matrix --algorithm packed-wfa --unique --threads 8 \
  --input data/human/phase3-window/cohort.fa > build/phase3.tsv
```

Sur 64 cœurs EPYC, les fronts bidirectionnels à deux bits calculent la matrice complète en 56,31 ± 0,32 ms. Les mêmes fronts sur octets prennent 54,51 ± 0,03 ms, WFA2 prend 74,24 ± 0,26 ms et Edlib 5,712 ± 0,002 s. Après regroupement, WFA2 prend 17,61 ± 0,12 ms et le front à deux bits 17,94 ± 0,33 ms, restitution des doublons comprise. Ces fronts parcourent les diagonales par coût croissant. NvPD calcule les lignes par indices de dernière occurrence. [Mesures à charge fixe](results/human/packed-cohort/trials-64.csv) · [Compteurs matériels](results/human/packed-pmu/counters.csv)

## Distribuer une matrice humaine

Préparer le chromosome 1 de CHM13 comme indiqué plus haut, puis en extraire 768 régions disjointes de 512 à 3 072 bases. Le script Slurm réserve deux nœuds EPYC de 192 cœurs chacun et mesure la même matrice sur un puis deux nœuds. Placer les archives de dépendances dans `deps/` et définir `ROMEO_PROJECT` avec le nom du compte Slurm. [Campagne et mesures](results/mpi/chromosome/)

```sh
python3 bench/mpi/chromosome/cohort.py data/human/chm13-panel/chr1.fa data/human/mpi-chr1
sbatch -A "$ROMEO_PROJECT" bench/mpi/chromosome/run.sbatch data/human/mpi-chr1
```

Dans une allocation de deux nœuds, les commandes suivantes placent 24 rangs Open MPI par nœud et huit travailleurs OpenMP par rang, soit 384 cœurs physiques. Elles exécutent successivement la répartition statique, le vol de travail et le parcours de Hilbert pondéré par le produit des longueurs. Avec `--length-bound 1`, le programme teste d’abord l’inclusion si la séquence la plus longue mesure au moins trois fois la plus courte. Si le test réussit, la distance vaut la différence des longueurs.

```sh
cmake -S . -B build-mpi -DCMAKE_BUILD_TYPE=Release -DCHPS_MPI=ON
cmake --build build-mpi --parallel 8 --target editdist-mpi
for policy in static steal hilbert-weighted; do
  mpirun -n 48 --map-by ppr:24:node:PE=8 --bind-to core build-mpi/editdist-mpi matrix \
    --input data/human/mpi-chr1/chr1.fa --algorithm nvpd-indexed --threads 8 \
    --workers openmp --schedule "$policy" --grain 4 --length-bound 1 > "build/${policy}.tsv"
done
```

Le test de sous-séquence résout 39 679 des 294 528 paires. En ordre de Hilbert pondéré sur 384 cœurs, il ramène WFA2 de 1,2931 ± 0,0021 à 1,0742 ± 0,0053 s et Edlib de 0,1680 ± 0,0006 à 0,1645 ± 0,0007 s. Les options `--schedule hilbert` et `--schedule cyclic` donnent le parcours de Hilbert sans pondération et la distribution cyclique. Dans une compilation avec `-DCHPS_KOKKOS=ON` et Kokkos/OpenMP, `--workers kokkos` lance les travailleurs via Kokkos. [Mesures sur les 768 régions](results/mpi/chromosome/839322/trials.csv)

Avec le noyau NvPD indexé et le vol de tâches, les travailleurs Kokkos traitent les 294 528 paires en 2,015 s sur un nœud et 1,183 s sur deux, soit 85,2 % d’efficacité forte. Sur deux nœuds, le rapport de temps OpenMP/Kokkos vaut 1,0009 (IC à 95 % de 0,9962 à 1,0056). [Comparaison des équipes OpenMP et Kokkos](results/mpi/chromosome/workers/)

## Choisir une implémentation

| Méthode ou option CMake | Utilisation |
| --- | --- |
| `seq`, `diagonal`, `nvpd-omp` | Programmation dynamique séquentielle, parcours par anti-diagonales OpenMP ou calcul NvPD par ligne. `--threads` fixe la taille de l’équipe OpenMP. |
| `edlib`, `wfa2`, `biwfa` | Edlib de Šošić et Šikić utilise le parallélisme sur bits. WFA2 de Marco-Sola et coll. suit des fronts de coût croissant. BiWFA les propage depuis les deux extrémités. Les deux parcours WFA acceptent `--trace`. |
| `nvpd-adaptive` | Score seul. Une borne sur le nombre d’éditions limite les colonnes NvPD. La bande s’élargit jusqu’à couvrir un chemin optimal. |
| `nvpd-split`, `nvpd-sparse` | Des coupes vérifiées répartissent les fragments entre cœurs. Pour le score, `nvpd-sparse` essaie un front borné à 128 éditions par grand fragment, puis passe à BiWFA. Les deux commandes reconstruisent aussi le CIGAR par le découpage initial. |
| `packed-wfa` | Pour le score sur A, C, G et T, le front bidirectionnel compare 32 bases par mot de 64 bits. BiWFA prend le relais si la paire dépasse le budget de 96 éditions ou contient un autre symbole. |
| `banded`, `banded-trim` | Le calcul élargit la bande jusqu’à couvrir la distance. `banded-trim` retire d’abord le préfixe et le suffixe identiques. |
| `hirschberg` | Reconstruit le chemin en mémoire linéaire. `--trace` fonctionne aussi avec Edlib et, jusqu’à 20 000 bases par séquence, avec `seq` et NvPD. |
| `-DCHPS_CUDA=ON` | Active `nvpd-cuda` et `nvpd-cuda-local` avec un compilateur CUDA C++20. La variante `nvpd-cuda-local` calcule le halo nécessaire à chaque bloc de colonnes. Sur GH200, ajouter `-DCMAKE_CUDA_ARCHITECTURES=90`. |
| `-DCHPS_KOKKOS=ON` | Active `nvpd-kokkos` avec Kokkos ≥ 4. Le mode score garde une équipe sur CPU et fusionne deux lignes sur GPU. |
| `-DCHPS_MPI=ON` | Construit `editdist-mpi`. `pair` répartit un alignement entre les rangs, `pairs` traite les paires consécutives du FASTA et `matrix` calcule toutes les distances par paires. |
| `-DCHPS_QUICKED=ON` | Construit `editdist-quicked`, avec contrôle du score par Edlib. |

`build/editdist --help` affiche les modes et les options. A*PA2 de Groot Koerkamp associe un parcours en blocs SIMD à une heuristique de coût restant. QuickEd de Doblas et coll. estime une borne pour reconstruire l’alignement en bande par parallélisme sur bits. Si tous les chemins optimaux sortent de la bande, le meilleur alignement accessible peut coûter davantage. Une erreur distincte de lecture au bord d’un mot de bits affecte deux témoins sur EPYC. L’adaptateur protège ce bord, rejoue le chemin et élargit la recherche si son coût diffère de celui d’Edlib. [Comparaison sur le chromosome 21](results/alternatives/human/README.md) · [Comparaison avec QuickEd d’origine](results/alternatives/README.md)

## Profiler les tuiles CUDA

Sur GH200, `CHPS_CUDA_TILE_ROWS` fixe le nombre de lignes traitées par tuile (de 1 à 32, avec 8 par défaut). L’option `CHPS_CUDA_PACKED_MATCHES` active la mise en mémoire partagée de quatre décalages sur huit bits par colonne. Les décalages trop grands utilisent l’index entier.

```sh
cmake -S . -B build-gpu -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90 -DCHPS_CUDA_TILE_ROWS=16 -DCHPS_CUDA_PACKED_MATCHES=ON
cmake --build build-gpu --parallel 8 --target editdist
build-gpu/editdist pair --algorithm nvpd-cuda-local --input build/paire.fa
```

À une fenêtre GRCh38 de 524 288 bases correspond une région de CHM13 qui contient un segment supplémentaire de 9 540 bases entre deux blocs homologues. À seize lignes par tuile, le compactage ramène le temps moyen de 1,878 à 1,461 s et l’énergie du GPU de 534,0 à 427,7 J par alignement. L’énergie provient du compteur cumulé du GPU. Le modèle Roofline combine le débit mesuré par une copie étalon et le trafic estimé à partir de la géométrie des tuiles. [Relevés énergétiques et profils](results/energy/gh200/) · [Modèle de trafic](results/roofline/gpu/)

## Se repérer dans le dépôt

Les noyaux et la ligne de commande se trouvent dans `src/`, les interfaces dans `include/chps/` et les tests dans `tests/`. Les outils de `scripts/` préparent les données et les arbres. Les bancs de `bench/` produisent les tables de `results/`. Le carnet dans `notebooks/` utilise les fonctions de tracé de `python/`. Les deux documents LaTeX de `report/` utilisent les graphiques de `figures/`. [Sources](src/) · [Tests](tests/) · [Résultats](results/)

Depuis la racine, les commandes suivantes régénèrent le carnet, les figures et les deux PDF avec Python 3.11 et Tectonic. La conversion du schéma TikZ en SVG utilise `pdftocairo` (Poppler).

```sh
python3.11 -m venv .venv-figures
. .venv-figures/bin/activate
python -m pip install -r notebooks/requirements.txt
python report/figures/phylograms.py
python report/figures/render_decision.py
python -m jupyter nbconvert --to notebook --execute --inplace notebooks/profil_nvpd.ipynb
(cd report && tectonic rapport.tex && tectonic tldr.tex)
```
