# Fenêtre de population sur le chromosome 21

La matrice comprend les deux haplotypes de 512 personnes tirées avec la graine 1001 parmi les 2 504 individus du VCF de phase 3. La fenêtre couvre l’intervalle GRCh37 `chr21:[20000000,20032768)` en coordonnées de base zéro. Le générateur conserve les appels PASS bialléliques, phasés et complets, dont les allèles contiennent uniquement A, C, G et T. Il écarte les variants qui chevauchent un intervalle déjà retenu ou une borne de la fenêtre. Les URL publiques et les empreintes figurent dans [dataset.json](dataset.json). [Collection IGSR](https://www.internationalgenome.org/data-portal/data-collection/phase3)

Depuis la racine du dépôt, Python 3.12 et pysam 0.24.1 reconstruisent le FASTA. Sur ROMEO, les deux tâches Slurm utilisent le fichier préparé et les archives locales des bibliothèques.

```sh
python -m venv build/cohort-env
build/cohort-env/bin/pip install pysam==0.24.1
build/cohort-env/bin/python bench/cpu/packed-cohort/prepare.py data/human/phase3-window
sbatch -A "$ROMEO_PROJECT" bench/cpu/packed-cohort/run.sbatch data/human/phase3-window/cohort.fa
sbatch -A "$ROMEO_PROJECT" bench/cpu/packed-cohort/profile.sbatch data/human/phase3-window/cohort.fa
```

Les fichiers `trials-*.csv` contiennent les temps à charge fixe. Le bloc zéro échauffe chaque configuration. Les cinq blocs suivants changent l’ordre des solveurs selon une permutation reproductible. Le chronomètre couvre la préparation des représentations, le calcul triangulaire et, si les séquences identiques sont regroupées, la restitution de chaque score. La lecture du FASTA précède le chronomètre. La matrice Edlib sert de référence pour chaque répétition.

Les fichiers `weak-*.csv` suivent des préfixes emboîtés du FASTA. Le nombre de paires par cœur reste proche de celui des 1 024 haplotypes sur 64 cœurs. La ligne à 64 cœurs est commune aux deux plans. Les distances de la matrice complète sont regroupées dans [distances.csv](distances.csv).

Le profil matériel mesure séparément la matrice des 128 premiers haplotypes sur un cœur. Chaque processus échauffe le solveur sur ces mêmes 8 128 paires, puis active les compteurs avec la FIFO de `perf stat`. Le compteur Zen 4 `rF045` mesure les échecs du TLB de données de premier niveau qui manquent aussi le second niveau. Les données brutes figurent dans [counters.csv](../packed-pmu/counters.csv).

Le champ `unique_rows` compte les entrées passées au solveur. Lorsque `collapse=0`, il comprend les haplotypes répétés. La colonne `encoding_bytes` donne le volume des mots compactés, en plus du FASTA en mémoire.

Les entrées publiques `packed-wfa --unique` et `edlib --unique` produisent le même fichier TSV, dont l’empreinte figure dans [matrix.sha256](matrix.sha256). Le profil des préfixes suit la charge des fronts en prenant le carré du score comme indicateur du nombre d’états. NumPy et pandas suffisent pour le recalcul.

```sh
build/editdist matrix --algorithm packed-wfa --unique --threads 8 \
  --input data/human/phase3-window/cohort.fa > build/cohort.tsv
python bench/cpu/packed-cohort/distance_profile.py build/cohort.tsv results/human/packed-cohort
```
