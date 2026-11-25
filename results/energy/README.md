# Énergie et tuiles CUDA

Les fenêtres du chromosome 21 proviennent des assemblages GRCh38 et CHM13. Les coordonnées et les empreintes SHA-256 figurent dans [le manifeste](../kokkos/human/references.csv). Les tâches Slurm utilisent un GH200 par allocation. Les temps de lecture du FASTA sont exclus.

| Allocation | Expérience | Données |
|---|---|---|
| 839310 | Profondeur de tuile, étalonnage de bande passante et durée des noyaux | [Lots](gh200/839310/trials.csv), [statistiques](gh200/839310/summary.csv), [étalonnage](gh200/839310/calibration.csv), [activité CUDA](gh200/839310/activity.csv) |
| 839317 | Lecture directe ou compacte des dernières occurrences à profondeur 16 | [Lots appariés](gh200/839317/packing.csv), [cas de bord](gh200/839317/validation.csv) |
| 839319 | Durée d’intégration du compteur, géométrie des halos | [Énergie](gh200/839319/window.csv), [chemin de calcul](gh200/839319/geometry.csv) |
| 839323 | Candidats homologues et inter-locus avec les trois solveurs publics | [Entrées](gh200/839323/inputs.csv), [mesures](gh200/839323/loci.csv) |

Chaque configuration énergétique utilise cinq processus. Le premier alignement initialise CUDA. Le deuxième détermine le nombre d'appels à grouper. Le chronomètre et le compteur couvrent le même lot. La colonne `joules` donne la différence entre deux lectures du compteur d’énergie du GPU. Le quotient `joules / iterations` donne l’énergie moyenne d’un alignement dans le lot. Dans le lot `insertion-524288`, la fenêtre de GRCh38 mesure 524 288 bases et son homologue de CHM13 en mesure 537 308. Sur cinq blocs, l’énergie moyenne passe de 534,0 J par alignement avec lecture directe à 427,7 J avec les indices compactés. [Mesures appariées](gh200/839317/packing.csv)

Les essais inter-locus chronomètrent un appel par processus, initialisation CUDA comprise. Le bloc zéro sert à l’échauffement.

Les scripts attendent le dossier FASTA produit par [inputs.py](../../bench/kokkos/human/inputs.py). Depuis la racine du dépôt sur ROMEO, transmettre ce dossier à `sbatch -A "$ROMEO_PROJECT" bench/gpu/energy/run.sbatch "$inputs"`. Les scripts `packing.sbatch`, `sensor.sbatch` et `loci.sbatch` acceptent le même argument. La commande `python bench/gpu/energy/summarize.py results/energy/gh200/839310` recalcule les moyennes, les écarts-types et les intervalles de Student à 95 %.

Le réglage par défaut de CUDA est `CHPS_CUDA_TILE_ROWS=8` avec `CHPS_CUDA_PACKED_MATCHES=OFF`. Pour reproduire le balayage en lecture directe, régler séparément la profondeur à 4, 8, 16 et 32 lignes. Le compactage de l’index utilise quatre symboles et trois tableaux en mémoire partagée. Si ces tableaux dépassent la limite du bloc, la tuile utilise les indices entiers. Si le halo dépasse aussi la capacité des deux tableaux de scores, le calcul avance ligne par ligne.

L’exemple suivant compile la variante compacte chronométrée sur GH200. Il s’exécute depuis la racine du dépôt.

```sh
cmake -S . -B build/gh200-compact -DCMAKE_BUILD_TYPE=Release -DCHPS_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90 -DCHPS_CUDA_TILE_ROWS=16 -DCHPS_CUDA_PACKED_MATCHES=ON
cmake --build build/gh200-compact --target editdist -j4
build/gh200-compact/editdist pair --algorithm nvpd-cuda-local --input pair.fa
```

La [source soumise](gh200/submitted.tar.gz) conserve la variante expérimentale des trois premières allocations. Sur les processus courts, le système de comptabilité peut renvoyer un RSS nul. Les comparaisons inter-locus utilisent la durée de l’appel.
