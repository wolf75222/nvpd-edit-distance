# Comparaison des bibliothèques exactes

## ROMEO, mesures à chaud sur AMD EPYC 9654

Le job 838762 a exécuté six méthodes sur dix paires. Chaque bloc tire au hasard l’ordre des 60 combinaisons entrée–méthode (graine 216814). Chaque combinaison démarre un processus, lit ses deux séquences, réalise un alignement d’échauffement, puis chronomètre dix alignements avec `steady_clock`. Les huit blocs donnent 80 appels mesurés par cellule. Le tableau fournit la moyenne des huit moyennes de processus, leur écart-type échantillonnal et un intervalle bilatéral de Student à 95 % (7 degrés de liberté). L’incertitude décrit les répétitions sur ce nœud. La compilation utilisait huit cœurs et la mesure un seul cœur physique, sans chevauchement d’alignements. Les intervalles temporels du stock comprennent aussi les appels ayant rendu un mauvais score. Les champs `exact`, `observed` et `correct_calls` permettent de les repérer. [Mesures individuelles](romeo/838762/trials.csv) · [Moyennes, SD, IC](romeo/838762/summary.csv) · [Entrées et empreintes](romeo/838762/cases.csv) · [Placement et compilateur](romeo/838762/environment.txt) · [Empreintes des sources chronométrées](romeo/838762/sources.sha256)

QuickEd stock inclut la construction et la libération du CIGAR amont. QuickEd certifié inclut le correctif de bord, la relecture du CIGAR, l’appel exact edlib et toute reprise. WFA2 et edlib calculent le coût. A*PA2-simple et full construisent puis libèrent le CIGAR. Les trois cas 20k et les deux contrôles humains ont les mêmes SHA-256 que la campagne M1. Les deux cas recherchés sur EPYC sont des contre-exemples sélectionnés. Ils renseignent le mécanisme d’erreur, sans estimer sa fréquence sur des données naturelles. [Pilote chronométré](../../bench/quicked/warm.cpp) · [Sélection des témoins](../../bench/quicked/find_witness.cpp) · [FASTA conservés](romeo/838762/witnesses)

### Borne adaptative et traceback de QuickEd

L’article distingue l’estimation rapide de la borne et la reconstruction de l’alignement. Dans la révision utilisée, `run_quicked` commence par une fenêtre de deux mots de 64 bits et un mot de recouvrement. Selon le nombre de fenêtres à forte erreur, il affine la borne avec des fenêtres plus larges dans les deux sens, puis éventuellement avec un calcul bandé de score. Le paramètre `bandwidth` de la voie par défaut intervient dans cette étape facultative. La valeur de coût finalement estimée est transmise à Hirschberg. Modifier ce seul paramètre ne force donc pas la largeur de la reconstruction finale. [Article](https://doi.org/10.1093/bioinformatics/btaf112) · [Code de la révision](https://github.com/maxdoblas/QuickEd/blob/2b1d44b8204a17fa0752fbe0b69b081c86a56e10/quicked/src/quicked.c)

Sur M1, les deux témoins à 512 bases donnent une borne 128, égale à leur distance. Au passage entre deux blocs de 64 colonnes, la routine de traceback consulte un prédécesseur dont l’indice relatif à la colonne suivante vaut −1. Le test d’origine décale alors un bit d’un nombre négatif de positions (comportement indéfini en C). Le premier passage se situe en (511,447) ou (255,191) en coordonnées (colonne,ligne) indexées dès zéro. Les scores amont sont 954 et 450. Les deux gardes interdisent de tester un bit extérieur aux mots stockés et rétablissent le coût 128 sur les deux entrées. Sur EPYC natif, l’estimateur donne 130 pour ces mêmes entrées et la reconstruction rend 128. Deux autres séquences à motifs répétés, de longueurs 256 et 768, y révèlent les erreurs amont 506 au lieu de 128 et 760 au lieu de 256. Chacune des 80 répétitions reproduit ces scores. Le résultat corrigé est exact aux 800 appels de la campagne. [Trace M1](diagnosis/trace-m1.csv) · [Justesse M1](diagnosis/correctness-m1.csv) · [Trace EPYC](romeo/838762/trace.csv) · [Correctif C appliqué par CMake](../../cmake/patches/quicked-trace-band.cmake)

Le programme certifié rejoue tous les mouvements, les égalités et les substitutions du CIGAR, vérifie les deux longueurs consommées et compare le coût à la distance globale edlib. En cas d’échec, il appelle explicitement le mode Hirschberg avec une borne de coût exprimée en pourcentage. Le premier pourcentage est arrondi vers le haut à partir de la distance certifiée, puis double jusqu’à 100 % si nécessaire. Si toutes les reprises échouent, le JSON identifie la distance de secours par `source=edlib-fallback`. Dans les dix cas chronométrés, le premier chemin corrigé a suffi. Un test distinct force la première bande à 1 % sur deux séquences de 1 152 bases. Il produit d’abord 896, puis le bon score 384 en une reprise à 34 % sur M1. [Adaptateur de production](../../src/libraries/quicked/quicked_adapter.c) · [Test de reprise](../../tests/quicked_widen.c)

Une expérience indépendante force la DP classique à rester dans `|i−j|≤r`. Sur les témoins de 512 bases, les coûts restreints sont 133 et 129 à r=63, puis 128 sur les deux entrées à r=64. Elle illustre l’effet mathématique d’une bande trop étroite. L’instrumentation de la voie amont sur M1 localise l’autre cause, au bord du tableau de bits. [Valeurs par rayon](diagnosis/fixed-band.csv) · [Calcul indépendant](../../bench/quicked/fixed_band.py)

### Reproduire la campagne à partir d’un clone

Depuis la racine du projet, préparer les données et l’archive A*PA2. Python 3, Git, CMake ≥3.22, un compilateur C/C++20 et OpenMP sont requis. Le téléchargement des références humaines et de la chaîne UCSC est vérifié par le générateur.

```sh
python3 scripts/datasets/synthetic.py
python3 scripts/datasets/human_reference.py
python3 bench/quicked/cases.py data/controls/quicked
mkdir -p deps build
git clone https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner build/quicked-astar
git -C build/quicked-astar archive --format=tar.gz \
    -o "$PWD/deps/astar-src.tar.gz" bf2e14e0cbc3a9a03600dcda0641d7f89e401e63
```

Sur ROMEO, soumettre depuis le clone contenant `data/controls/quicked/` et `deps/astar-src.tar.gz`. CMake vérifie les archives SHA-256 de QuickEd, edlib et WFA2. Il les télécharge quand elles sont absentes de `deps/`. Le job vérifie la révision intégrée dans l’archive Git A*PA, l’empreinte du Cargo.lock (`ab826558…7e12c`) et celle du manifeste Rust nightly du 21 mai 2026 (`6c7ed921…ec5ab`). Rustup contrôle les composants du manifeste et Cargo compile avec `--locked`. Le cache de téléchargement Cargo peut être fourni facultativement sous `deps/cargo-cache.tar.gz` avec `registry/cache`, `registry/index` et `git/db`. Le cache n’altère pas les versions. Les outils et la compilation Rust utilisent un dossier projet temporaire supprimé en fin de job. L’accès HTTPS au registre et à static.rust-lang.org est nécessaire pour les composants manquants. Le répertoire de sortie est `results/alternatives/romeo/$SLURM_JOB_ID`. [Script Slurm](../../bench/quicked/run.sbatch) · [Compilateur Rust mesuré](romeo/838762/rust-version.txt)

```sh
sbatch -A r260105 -p instant bench/quicked/run.sbatch
```

Pour refaire les diagnostics sans Slurm, utiliser un répertoire CMake compilé avec `CHPS_QUICKED=ON`, les cibles `measure-quicked`, `measure-quicked-stock` et `measure-quicked-native`, puis lancer `bench/quicked/diagnose.py BUILD data/controls/quicked trace.csv`, `bench/quicked/check_witnesses.py BUILD data/controls/quicked MACHINE correctness.csv` et `bench/quicked/fixed_band.py data/controls/quicked fixed-band.csv` avec Python. Pour l’instrumentation, choisir les mêmes compilateurs C et C++ que CMake avec `--cc` et `--cxx`.

## Série historique sur Apple M1 Pro

Cinq processus par couple méthode–entrée, dans un ordre mélangé avec la graine 2026. GCC 16.1.0 compile C et C++ en mode Release avec `-mcpu=native`. Les cinq méthodes utilisent un fil, sans affinité imposée. Chaque essai lit d'abord le FASTA, puis chronomètre l'appel, allocations comprises. Le pic RSS couvre tout le processus. L’exécutable historique QuickEd 1.2 utilise le chemin amont par défaut, avec reconstruction du CIGAR. Les autres méthodes demandent le coût seul.

Les cas 20k utilisent la première paire de chaque taux du générateur de référence. Les deux autres cas utilisent les premières bases des fenêtres humaines homologues du chromosome 21. Le contrôle de 100 000 bases conserve les deux extrémités et permute A→C, C→G, G→T, T→A dans les 128 bases centrales. Les cinq méthodes retrouvent les mêmes distances sur ces entrées. Deux tests supplémentaires de 512 bases à motifs répétés exposent une erreur de QuickEd 1.2. La révision historique de l’exécutable comparait chaque score à edlib après la mesure et signalait tout désaccord. La série ROMEO ci-dessus mesure le chemin corrigé avec certification incluse. [Manifeste des fenêtres](../human/manifest.tsv) · [Mesures individuelles](trials.csv) · [Méthode QuickEd](https://doi.org/10.1093/bioinformatics/btaf112)

```sh
# À la révision historique ba18624.
cmake -S . -B build-alt -DCMAKE_BUILD_TYPE=Release -DCHPS_QUICKED=ON \
  -DCMAKE_C_FLAGS=-mcpu=native -DCMAKE_CXX_FLAGS=-mcpu=native
cmake --build build-alt --parallel
python3 scripts/datasets/synthetic.py
python3 scripts/datasets/human_reference.py
python3 bench/cpu/alternatives/alternatives.py --build build-alt --out build-alt/trials.csv --machine apple-m1-pro
```

## A*PA2 sur les mêmes paires

L'expérience [astar.csv](astar.csv) rassemble 175 processus neufs sur Apple M1 Pro, soit cinq entrées, sept variantes et cinq répétitions. Les deux modes publiés, A*PA2-simple et A*PA2-full, calculent la distance globale unitaire **avec chemin**. Les fonctions Rust de l'interface C construisent aussi le CIGAR sous forme de chaîne. Le temps de l'adaptateur couvre l'appel et la libération de cette chaîne. Les deux modes donnent les mêmes distances qu'edlib sur les cinq entrées. [Article WABI 2024](https://doi.org/10.4230/LIPIcs.WABI.2024.17) · [Code amont utilisé](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner/tree/bf2e14e0cbc3a9a03600dcda0641d7f89e401e63) · [Interface C](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner/blob/bf2e14e0cbc3a9a03600dcda0641d7f89e401e63/astarpa-c/src/lib.rs) · [Adaptateur](../../bench/cpu/alternatives/astar_trial.cpp)

| Entrée | Distance | Médiane simple (ms) | Médiane full (ms) |
| --- | ---: | ---: | ---: |
| 20k, 2% | 396 | 0,841 | 1,587 |
| 20k, 5% | 980 | 1,167 | 1,748 |
| 20k, 10% | 1 908 | 1,593 | 2,323 |
| Fenêtre chr21, 50k | 698 | 2,160 | 3,416 |
| Contrôle central chr21, 100k | 72 | 2,002 | 6,196 |

Reproduction depuis la racine du projet après la préparation des entrées et de `build-alt` ci-dessus. Les fichiers amont fixent le compilateur Rust à `nightly-2026-05-21`, les dépendances via `Cargo.lock` et les instructions natives via `-C target-cpu=native`. La révision date du 22 mai 2026. Sur la machine mesurée, `g++-16` désigne GCC 16.1.0. L'adaptateur inclut l'en-tête C publié par l'amont. [Collecteur](../../bench/cpu/alternatives/alternatives.py)

```sh
up=build/astar-pairwise-aligner
git clone https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner "$up"
git -C "$up" checkout --detach bf2e14e0cbc3a9a03600dcda0641d7f89e401e63
(cd "$up" && cargo build --release --locked -p astarpa-c)
g++-16 -std=c++20 -O3 -mcpu=native -Iinclude -I"$up/astarpa-c" \
  bench/cpu/alternatives/astar_trial.cpp src/io/fasta.cpp -L"$up/target/release" -lastarpa_c \
  -Wl,-rpath,"$PWD/$up/target/release" -o build-alt/astar-trial
OMP_NUM_THREADS=1 python3 bench/cpu/alternatives/alternatives.py --build build-alt \
  --astar-binary build-alt/astar-trial --out build-alt/astar.csv --machine apple-m1-pro
```
