# Comparaison sur Apple M1 Pro

Cinq processus par couple méthode–entrée, dans un ordre mélangé avec la graine 2026. GCC 16.1.0 compile C et C++ en mode Release avec `-mcpu=native`. Les cinq méthodes utilisent un fil, sans affinité imposée. Chaque essai lit d'abord le FASTA, puis chronomètre l'appel, allocations comprises. Le pic RSS couvre tout le processus. L'exécutable QuickEd 1.2 utilise le chemin exact par défaut, avec reconstruction du CIGAR. Les autres méthodes demandent le coût seul.

Les cas 20k utilisent la première paire de chaque taux du générateur de référence. Les deux autres cas utilisent les premières bases des fenêtres humaines homologues du chromosome 21. Le contrôle de 100 000 bases conserve les deux extrémités et permute A→C, C→G, G→T, T→A dans les 128 bases centrales. Les cinq méthodes retrouvent les mêmes distances sur ces entrées. Deux tests supplémentaires de 512 bases à motifs répétés exposent une erreur de QuickEd 1.2. L’exécutable compare désormais chaque score à edlib après la mesure et retourne une erreur si les valeurs diffèrent. [Manifeste des fenêtres](../human/manifest.tsv) · [Mesures individuelles](trials.csv) · [Méthode QuickEd](https://doi.org/10.1093/bioinformatics/btaf112)

```sh
cmake -S . -B build-alt -DCMAKE_BUILD_TYPE=Release -DCHPS_QUICKED=ON \
  -DCMAKE_C_FLAGS=-mcpu=native -DCMAKE_CXX_FLAGS=-mcpu=native
cmake --build build-alt --parallel
python3 scripts/synthetic.py
python3 scripts/human_reference.py
python3 bench/alternatives.py --build build-alt --out build-alt/trials.csv --machine apple-m1-pro
```

## A*PA2 sur les mêmes paires

L'expérience [astar.csv](astar.csv) rassemble 175 processus neufs sur Apple M1 Pro, soit cinq entrées, sept variantes et cinq répétitions. Les deux modes publiés, A*PA2-simple et A*PA2-full, calculent la distance globale unitaire **avec chemin**. Les fonctions Rust de l'interface C construisent aussi le CIGAR sous forme de chaîne. Le temps de l'adaptateur couvre l'appel et la libération de cette chaîne. Les deux modes donnent les mêmes distances qu'edlib sur les cinq entrées. [Article WABI 2024](https://doi.org/10.4230/LIPIcs.WABI.2024.17) · [Code amont utilisé](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner/tree/bf2e14e0cbc3a9a03600dcda0641d7f89e401e63) · [Interface C](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner/blob/bf2e14e0cbc3a9a03600dcda0641d7f89e401e63/astarpa-c/src/lib.rs) · [Adaptateur](../../bench/astar_trial.cpp)

| Entrée | Distance | Médiane simple (ms) | Médiane full (ms) |
| --- | ---: | ---: | ---: |
| 20k, 2% | 396 | 0,841 | 1,587 |
| 20k, 5% | 980 | 1,167 | 1,748 |
| 20k, 10% | 1 908 | 1,593 | 2,323 |
| Fenêtre chr21, 50k | 698 | 2,160 | 3,416 |
| Contrôle central chr21, 100k | 72 | 2,002 | 6,196 |

Reproduction depuis la racine du projet après la préparation des entrées et de `build-alt` ci-dessus. Les fichiers amont fixent le compilateur Rust à `nightly-2026-05-21`, les dépendances via `Cargo.lock` et les instructions natives via `-C target-cpu=native`. La révision date du 22 mai 2026. Sur la machine mesurée, `g++-16` désigne GCC 16.1.0. L'adaptateur inclut l'en-tête C publié par l'amont. [Collecteur](../../bench/alternatives.py)

```sh
up=build/astar-pairwise-aligner
git clone https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner "$up"
git -C "$up" checkout --detach bf2e14e0cbc3a9a03600dcda0641d7f89e401e63
(cd "$up" && cargo build --release --locked -p astarpa-c)
g++-16 -std=c++20 -O3 -mcpu=native -Iinclude -I"$up/astarpa-c" \
  bench/astar_trial.cpp src/fasta.cpp -L"$up/target/release" -lastarpa_c \
  -Wl,-rpath,"$PWD/$up/target/release" -o build-alt/astar-trial
OMP_NUM_THREADS=1 python3 bench/alternatives.py --build build-alt \
  --astar-binary build-alt/astar-trial --out build-alt/astar.csv --machine apple-m1-pro
```
