# Géométrie locale du chromosome 20

Les en-têtes du FASTA donnent les intervalles semi-ouverts à partir de zéro dans les versions CP068258.1 et CP068258.2 de CHM13. Les 33 bases de chaque fenêtre proviennent du FASTA de la comparaison entre assemblages. Pour les extraire avec samtools, les bornes ci-dessous sont inclusives et comptées à partir de un. [Assemblages et téléchargement](../../human/whole/README.md) · [Entrées de la figure](chr20.fa)

```sh
samtools faidx data/human/chr20-assemblies.fa chr20-assemblies_A:66210139-66210171 chr20-assemblies_B:66210147-66210179
```

Depuis la racine, le pilote du banc adaptatif fournit le score et les compteurs de travail. Sous Linux avec GCC, après la construction de `build/editdist`, on peut le compiler et rejouer les deux mesures ainsi. [Pilote](../../../bench/cpu/adaptive/driver.cpp)

```sh
g++ -O3 -std=c++20 -fopenmp -Iinclude -Isrc/cpu bench/cpu/adaptive/driver.cpp src/io/fasta.cpp \
    build/libdistances.a build/libedlib_static.a build/libwfa_static.a -lm -o build/adaptive-driver
build/adaptive-driver align results/adaptive/geometry/chr20.fa nvpd-adaptive
build/editdist pair --algorithm edlib --trace --input results/adaptive/geometry/chr20.fa
```

Le pilote écrit le temps, le score, le nombre de mises à jour, le nombre de passages, le dernier budget et les longueurs retirées aux deux extrémités. Le relevé conserve les compteurs et le CIGAR. La figure utilise le travail de la dernière bande complète et le compteur total pour retrouver la ligne d’arrêt du premier passage. [Relevé](check.csv) · [Construction de la figure](../../../python/studies/band_geometry.py)
