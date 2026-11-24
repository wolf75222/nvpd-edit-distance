# Chemins d’édition sur trois fenêtres du chromosome 21 humain

Trois méthodes calculent un alignement global à coût unitaire entre GRCh38 et T2T-CHM13. Les fenêtres de 2 et 4 millions de bases commencent à la position 34 000 000 de GRCh38. La fenêtre de 8 millions porte sur un autre locus du bras q, à partir de 14 000 100. Leurs distances exactes valent respectivement 5 599, 12 482 et 47 493. Les manifestes donnent les coordonnées des deux assemblages, les longueurs et les empreintes SHA-256 des FASTA. [Fenêtres de 2 et 4 Mb](../../human/algorithms/838910/cases.csv) · [Fenêtre de 8 Mb](../../human/algorithms/arm/838918/cases.csv)

## Résultats sur AMD EPYC 9654

Les temps incluent la construction du chemin. Les valeurs sont des moyennes ± écarts-types entre cinq processus. Le maximum de RSS inclut le lecteur FASTA, le pilote et les bibliothèques chargées. Les processeurs logiques 133, 134 et 135 du nœud c022 correspondent à trois cœurs physiques, avec un seul fil visible par cœur. [Essais individuels](839083/trials.csv) · [Moyennes, écarts-types et IC à 95 %](839083/summary.csv) · [Environnement](839083/environment.txt) · [Compte rendu Slurm](839083/accounting.psv)

| Fenêtre | Méthode | Durée (s) | Pic RSS (Mio) |
| --- | --- | ---: | ---: |
| 2 Mb | BiWFA | 0,04782 ± 0,00060 | 29,70 ± 0,99 |
| 2 Mb | A*PA2-full | 0,14134 ± 0,00398 | 125,25 ± 0,82 |
| 2 Mb | QuickEd avec contrôle du chemin | 1,63906 ± 0,00020 | 59,24 ± 0,55 |
| 4 Mb | BiWFA | 0,22115 ± 0,00121 | 48,82 ± 1,47 |
| 4 Mb | A*PA2-full | 0,48407 ± 0,00513 | 246,63 ± 1,58 |
| 4 Mb | QuickEd avec contrôle du chemin | limite du processus | — |
| 8 Mb | BiWFA | 2,8451 ± 0,0079 | 91,48 ± 3,98 |
| 8 Mb | A*PA2-full | 10,040 ± 0,628 | 482,92 ± 0,54 |
| 8 Mb | QuickEd avec contrôle du chemin | limite du processus | — |

Les 35 alignements chronométrés ont retrouvé les distances de référence. Le rejeu indépendant couvre 15 CIGAR de BiWFA et 15 d’A*PA2-full. Chacune des cinq exécutions terminées de QuickEd a validé le premier chemin proposé. À 4 et à 8 Mb, le premier processus QuickEd a atteint sa limite de 120 secondes. Le champ de temps reste vide pour ces deux entrées. La dispersion la plus forte concerne A*PA2-full à 8 Mb, avec un IC à 95 % de [9,260 ; 10,820] secondes. Cette fenêtre contient aussi davantage d’éditions que les deux autres. [Contrôle edlib et WFA2 à 2 et 4 Mb](../../human/algorithms/838910/trials.csv) · [Contrôle à 8 Mb](../../human/algorithms/arm/838918/trials.csv)

## Mesure et contrôle

Chaque fenêtre occupe un cœur physique fixe. Les trois fenêtres peuvent s’exécuter simultanément. Pour une même fenêtre, les méthodes se succèdent dans un ordre qui change d’un bloc au suivant. Un processus neuf lit le FASTA, effectue un alignement complet d’échauffement, puis chronomètre le second avec une horloge monotone. Le pic RSS provient du champ Linux `VmHWM` et couvre toute la vie du processus jusqu’à la fin du second appel. Le fichier est lu et son empreinte contrôlée hors chronomètre. Chaque processus dispose de 4 Gio d’espace d’adressage et de 120 secondes pour la lecture et les deux appels réunis.

A*PA2-full renvoie un CIGAR étendu par son interface C. BiWFA renvoie son CIGAR par le même pilote C++. Après le chronométrage, un lecteur indépendant parcourt chaque chemin contre les bases des deux entrées, contrôle les opérations, les longueurs consommées et le coût final. Le CSV conserve la taille et le SHA-256 de la chaîne dans son encodage d’origine. QuickEd rejoue son chemin en mémoire et vérifie son coût avec edlib, dans la période mesurée. Son interface libère ensuite le chemin et restitue uniquement le score, sa provenance et le nombre de reprises. Les empreintes de CIGAR ne sont donc renseignées que pour A*PA2-full et BiWFA. [Pilote C++](../../../bench/quicked/human/trial.cpp) · [Vérification du chemin QuickEd](../../../src/libraries/quicked/quicked_adapter.c) · [Collecteur](../../../bench/quicked/human/measure.py)

Une moyenne et son écart-type échantillonnal décrivent cinq processus à entrée fixée. L’intervalle bilatéral de Student à 95 % a quatre degrés de liberté. Une limite atteinte est enregistrée au premier essai de la méthode sur cette entrée. Les blocs suivants sont marqués comme non exécutés. La limite s’applique à l’ensemble du processus et ne donne pas de borne inférieure de 120 secondes au seul alignement à chaud.

## Reproduction

Depuis la racine du dépôt, générer les trois entrées avec les mêmes en-têtes, coordonnées et retours à la ligne que sur ROMEO. La chaîne UCSC est mise en cache. Les fenêtres auxiliaires de 16, 24 et 32 Mb produites par le générateur du bras q disparaissent à la fermeture du répertoire temporaire. [Générateur des assemblages publics](../../../scripts/datasets/human_reference.py) · [Découpage du bras q](../../../bench/human/algorithms/arm.py)

```sh
python3 scripts/datasets/human_reference.py --whole-assemblies \
    --cache data/human/reference-cache --manifest data/human/manifest.tsv
PYTHONPATH=bench/human/algorithms python3 - <<'PYTHON'
from pathlib import Path
from tempfile import TemporaryDirectory
from arm import windows
root = Path("data/human")
with TemporaryDirectory(dir=root) as temporary:
    work = Path(temporary)
    windows(work, root / "reference-cache/hg38ToHs1.over.chain.gz")
    (work / "chr21q-8m.fa").replace(root / "chr21q-8m.fa")
PYTHON
```

CMake vérifie les archives de QuickEd, edlib et WFA2. Le banc demande l’archive du commit A*PA2 et un cache Cargo. Depuis un clone neuf, les commandes suivantes préparent les deux fichiers dans `deps/`. Un cache initialement vide convient avec un accès HTTPS à GitHub, au registre Cargo et aux distributions Rust. Un cache existant est conservé. [Versions des dépendances](../../../CMakeLists.txt) · [Campagne sur les entrées courtes](../README.md)

```sh
mkdir -p build deps
git clone https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner build/astar-source
git -C build/astar-source archive --format=tar.gz \
    -o "$PWD/deps/astar-src.tar.gz" bf2e14e0cbc3a9a03600dcda0641d7f89e401e63
test -f deps/cargo-cache.tar.gz || tar -czf deps/cargo-cache.tar.gz -T /dev/null
sbatch --account=r260105 bench/quicked/human/run.sbatch
```

Le script vérifie le commit, le verrou Cargo et le manifeste Rust. Rust et les compilations C/C++ utilisent un répertoire créé sous `/project/$SLURM_JOB_ACCOUNT/$USER/`, supprimé à la fin du job. Les sorties restent dans `results/alternatives/human/$SLURM_JOB_ID/`. Les sept fichiers de sources utilisés ont été contrôlés par SHA-256. Le lecteur FASTA de cette mesure correspond à la révision `945a36be46fae7599fff32f53338c05ad0ee77b3` de `src/io/fasta.cpp`. Trois tentatives de préparation ont précédé 839083. Elles se sont arrêtées avant les essais, lors de la création du répertoire temporaire Rust (838930), du chargement hors ligne d’une dépendance Linux (839080) et de l’extraction de QuickEd sur scratch (839082). [Script Slurm](../../../bench/quicked/human/run.sbatch) · [Empreintes des sources](839083/sources.sha256) · [Empreintes des binaires](839083/binaries.sha256)
