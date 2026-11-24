# Comparaison sur Apple M1 Pro

Cinq processus par couple méthode–entrée, dans un ordre mélangé avec la graine 2026. GCC 16.1.0 compile C et C++ en mode Release avec `-mcpu=native`. Les cinq méthodes utilisent un fil, sans affinité imposée. Chaque essai lit d'abord le FASTA, puis chronomètre l'appel, allocations comprises. Le pic RSS couvre tout le processus. L'exécutable QuickEd 1.2 utilise le chemin exact par défaut, avec reconstruction du CIGAR. Les autres méthodes demandent le coût seul.

Les cas 20k utilisent la première paire de chaque taux du générateur de référence. Les deux autres cas utilisent les premières bases des fenêtres humaines homologues du chromosome 21. Le contrôle de 100 000 bases conserve les deux extrémités et permute A→C, C→G, G→T, T→A dans les 128 bases centrales. Les cinq méthodes retrouvent les mêmes distances. [Manifeste des fenêtres](../human/manifest.tsv) · [Mesures individuelles](trials.csv) · [Méthode QuickEd](https://doi.org/10.1093/bioinformatics/btaf112)

```sh
cmake -S . -B build-alt -DCMAKE_BUILD_TYPE=Release -DCHPS_QUICKED=ON \
  -DCMAKE_C_FLAGS=-mcpu=native -DCMAKE_CXX_FLAGS=-mcpu=native
cmake --build build-alt --parallel
python3 scripts/synthetic.py
python3 scripts/human_reference.py
python3 bench/alternatives.py --build build-alt --out build-alt/trials.csv --machine apple-m1-pro
```
