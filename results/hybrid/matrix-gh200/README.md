# Temps de calcul des matrices sur GH200

Sur ROMEO, chaque essai réserve quatre cœurs Grace et un GPU GH200. Chaque processus lit le FASTA, calcule une matrice d’échauffement, puis chronomètre une nouvelle matrice. Le chronomètre couvre l’allocation de la matrice, la construction des esquisses, l’ordonnancement et le remplissage des deux triangles. Le programme écrit le TSV après avoir rapatrié le dernier score. Les blocs mesurent la variation du temps sur des entrées fixes.

## Fenêtres du chromosome 21

Les trois lots utilisent des fenêtres de 100 000 bases de GRCh38 et CHM13v2.0. Dans le lot à un locus, les fenêtres sont décalées de 0 à 384 bases au sein de la région homologue. Le lot à deux loci couvre deux régions distantes de 2 Mb sur GRCh38 et le troisième réunit huit loci disjoints. Le [manifeste](manifest.tsv) donne l’accession, les coordonnées sur le brin fourni (origine zéro, fin exclue) et l’empreinte SHA-256 de chaque fenêtre. Les paires d’un même lot partagent des séquences.

Le job 839949 calcule chaque matrice avec sept politiques. L’ordre de départ est tiré une fois par lot, puis décalé d’une position à chaque bloc. Le fichier [trials.csv](trials.csv) contient les matrices chronométrées et [warmups.csv](warmups.csv) les échauffements des mêmes processus. La colonne `matrix_sha256` identifie dans les deux fichiers la matrice écrite en sortie du processus. Chaque matrice chronométrée concorde avec la matrice de référence calculée par Edlib sur un cœur.

La politique `static` alterne entre NvPD sur GPU et Edlib sur CPU dans l’ordre des paires du triangle supérieur. La politique `cpu-sketch` choisit WFA2 ou Edlib sur CPU selon les mots partagés. Avec `hybrid`, WFA2 traite les paires proches et le GPU traite les autres. Avec `hybrid-share`, les fils CPU libérés peuvent traiter des paires restantes avec Edlib. La colonne `sketch_seconds` inclut l’estimation et la mise en file pour les politiques qui utilisent l’esquisse. Le diagnostic par paire donne la décision prise avant l’alignement. Les colonnes `cpu_pairs` et `gpu_pairs` comptent les paires effectivement alignées par chaque dispositif.

## Matrice des 64 génomes de surveillance

Le [premier essai](surveillance/trials.csv) compare Edlib, WFA2 et le routage strict sur six blocs (job 839971). Le [second essai](composition/trials.csv) répète la comparaison sur huit blocs et ajoute le choix du solveur selon la composition (job 839973). Chaque gain se calcule par rapport au témoin du même bloc et du même nœud. Les trois réservations et leurs environnements sont consignés dans les fichiers `accounting.tsv` et `environment.txt` de leurs répertoires respectifs. La [matrice de référence](../../phylogeny/sars-edit/matrix.tsv) contient les 2 016 distances entre génomes.

## Reproduction sur ROMEO

Depuis la racine du dépôt, préparer les entrées puis soumettre le protocole voulu. Le compte Slurm doit être défini dans `ROMEO_PROJECT`.

```sh
python3 bench/hybrid/matrix/prepare.py
sbatch --account="$ROMEO_PROJECT" bench/hybrid/matrix/run.sbatch build/hybrid-matrix/inputs

python3 scripts/datasets/ncbi_sars.py --fasta build/hybrid-matrix/viral/surveillance.fa
sbatch --account="$ROMEO_PROJECT" bench/hybrid/matrix/run.sbatch build/hybrid-matrix/viral surveillance
sbatch --account="$ROMEO_PROJECT" bench/hybrid/matrix/run.sbatch build/hybrid-matrix/viral composition
```

Le [banc de mesure](../../../bench/hybrid/matrix/measure.py) vérifie les identifiants, les scores et la symétrie avant de conserver les temps. Les [graphiques](../../../python/studies/hybrid_matrix.py) utilisent les blocs appariés pour calculer les moyennes et les intervalles de Student à 95 %.
