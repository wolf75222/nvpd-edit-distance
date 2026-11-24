# Fronts d’onde affines sur GH200

WFA-GPU d’Aguado-Puig et coll. calcule un alignement global par fronts d’onde sur GPU. Avec `-g 1,0,1`, une lacune de longueur ℓ coûte ℓ et une substitution coûte un. Le programme fournit ici le score seul en mode sans bande. La comparaison de chemins avec la version eWFA-GPU de 2022 figure dans une expérience séparée. [Article de 2023](https://doi.org/10.1093/bioinformatics/btad701) · [Révision du programme](https://github.com/Marco-Sola-Lab/WFA-GPU/tree/85acbcf4224811b3e52814b1bbe9364b5e2ba5ed) · [Expérience eWFA-GPU](../ewfa/README.md)

## Calcul des temps

Pour les lots synthétiques, φX174 et le chromosome humain, les quatre programmes lisent le même fichier et écrivent un score par paire. Le chronomètre du processus parent inclut le démarrage, l’initialisation CUDA pour les méthodes GPU et les entrées-sorties. Edlib et WFA2 utilisent un fil sur Grace. NvPD réutilise un espace de travail CUDA et traite les paires en série. WFA-GPU envoie jusqu’à 512 paires indépendantes au GPU. Le programme utilise un fil hôte pour les reprises. Le mode affine de sa dépendance WFA2 reconstruit en interne le chemin des paires reprises.

Un bloc de chauffe précède cinq blocs mesurés. L’ordre des méthodes change dans chaque bloc selon une permutation pseudo-aléatoire. La dispersion correspond à l’écart-type entre processus. L’intervalle de Student à 95 % a pour demi-largeur $2{,}776\,s/\sqrt{5}$. Les rapports utilisent deux temps du même bloc. Les deux premières allocations partagent l’hôte `romeo-a058` pendant une partie des mesures, avec des GPU et des cœurs distincts. [Placement du premier job](839291/placement.txt) · [Placement du second](839294/placement.txt) · [Comptabilité Slurm](accounting.psv)

## Lots synthétiques et panel φX174

Les jeux synthétiques croisent 5 et 20 knt avec 2 et 10 % d’éditions équiprobables. Chaque lot de cinq paires reprend les cinq premières du lot de 128 de la même configuration. Le plafond du GPU dépasse de huit le nombre d’éditions injectées. Un contrôle séparé à `-e 32` force les cinq reprises sur CPU et retrouve leurs scores edlib. Les 2 016 paires du panel φX174 utilisent `-e 256`. Les séquences proviennent du manifeste NCBI et les scores de référence de la matrice du job 837271. [Essais bruts](839291/trials.csv) · [Temps et rapports](839291/summary.csv) · [Contrôle des reprises](839291/fallback-control.txt)

À 20 knt et 10 %, le lot de 128 paires prend 0,395 ± 0,006 s avec WFA-GPU, 0,518 ± 0,001 s avec WFA2 et 0,250 ± 0,001 s avec edlib. Le rapport apparié WFA-GPU/WFA2 vaut 0,762 [0,748 à 0,776]. Tous les scores synthétiques concordent avec edlib. Parmi les paires de φX174, 186 comportent au moins un N. Le GPU les renvoie toutes au CPU. Le temps du lot atteint 0,487 ± 0,047 s avec WFA-GPU et 0,0246 ± 0,0012 s avec WFA2. [Manifeste φX174](../ncbi/phix174.tsv) · [Matrice de référence](../romeo/cohort-837271/matrix.tsv)

## Fenêtres homologues du chromosome humain

Les 768 fenêtres apparient les chromosomes 21 de GRCh38 et de T2T-CHM13 sur la chaîne colinéaire 23. Elles couvrent 12 582 912 bases de référence sur des intervalles disjoints. Le programme contrôle chaque empreinte de séquence et chaque score contre les résultats du calcul MPI. Les distances vont de 0 à 11 107. Le format signé sur 16 bits de la révision amont réserve au CPU les chaînes de 32 768 bases ou plus. Les fenêtres de ce lot restent sous cette limite. [Coordonnées et empreintes](../mpi/human/windows.tsv) · [Scores de référence](../mpi/human/scores.tsv) · [Mesures](839294/trials.csv) · [Statistiques](839294/summary.csv)

Les plafonds de 64, 256 et 1 024 entraînent respectivement 152, 41 et 10 reprises. Avec la limite 64, le noyau vérifie les scores strictement inférieurs à 63. Les 152 fenêtres de score 63 ou plus passent sur CPU. La durée moyenne varie de 1,069 ± 0,036 s à 1,023 ± 0,050 s entre les plafonds extrêmes. Le rapport apparié vaut 0,958 [0,899 à 1,017]. WFA2 termine ce lot sur un cœur Grace en 0,0725 ± 0,0009 s. [Noyau de distance amont](https://github.com/Marco-Sola-Lab/WFA-GPU/blob/85acbcf4224811b3e52814b1bbe9364b5e2ba5ed/lib/kernels/sequence_distance_kernel.cu#L322)

## Répartition de la matrice SARS-CoV-2

Le pilote lit le FASTA de 64 accessions et compare les symboles IUPAC à l’identique. Le GPU reçoit les paires dont les deux séquences contiennent exclusivement A, C, G et T et comptent moins de 32 768 bases. Le plafond vaut 256. Le pilote écrit les autres paires dans quatre fichiers pour les travailleurs WFA2. Les numéros de paire déterminent l’affectation cyclique à ces quatre travailleurs. Le travail sur le GPU et le calcul des quatre fichiers s’exécutent simultanément.

Les témoins edlib et WFA2 utilisent quatre travailleurs sur le même hôte. Le chronomètre couvre le programme du FASTA au TSV, y compris la répartition des entrées. Après chaque chronométrage, le pilote compare les 2 016 identifiants de paire et les scores du TSV à la matrice de référence. Les résultats publiés contiennent uniquement les mesures et les empreintes. [Manifeste](../ncbi/sars_content.tsv) · [Matrice attendue](../phylogeny/sars-edit/matrix.tsv) · [Pilote](../../bench/gpu/wfa_gpu/sars_cohort.py) · [Temps individuels](839298/trials.csv) · [Résultats par méthode](839298/summary.csv)

Le GPU calcule les quinze paires entre les six séquences ACGT, toutes de distance inférieure ou égale à 216. Le CPU calcule les 2 001 autres. Le temps de bout en bout vaut 6,844 ± 0,015 s pour le calcul hybride, 6,874 ± 0,030 s pour WFA2 et 2,628 ± 0,014 s pour edlib. Le rapport hybride/WFA2 vaut 0,9956 [0,9908 à 1,0004]. Le lot GPU dure en moyenne 0,498 s. Le plus lent des quatre lots CPU prend 6,684 s.

Le front de score $s$ admet au plus $2s+1$ diagonales à coût unitaire. Pour une paire de distance $d$, la somme $H(d)=(d+1)^2$ compte les positions possibles avant troncature aux bords. Sur la matrice de référence, $\sum H(d)$ vaut 223 071 pour les quinze paires GPU et 26 084 927 678 pour le panel. La fraction vaut 0,000855 %. Ce modèle ne compte pas les extensions le long de chaque diagonale. [Scores par paire](../phylogeny/sars-edit/matrix.tsv)

## Reproduction

Depuis la racine du dépôt sur ROMEO, préparer le FASTA de φX174 et réserver un GPU pour chaque étude. La génération du panel humain télécharge les références et la chaîne de correspondance dans `build/human-reference`. Le répertoire `deps` peut contenir les archives `wfa-gpu.tar.gz` et `wfa-gpu-wfa2.tar.gz` aux révisions fixées dans le lanceur. Le lanceur vérifie leurs empreintes, les extrait dans le répertoire temporaire du nœud et compile pour `sm_90`.

```sh
python3 scripts/datasets/ncbi_phix174.py
sbatch -A "$ROMEO_PROJECT" -p instant bench/gpu/wfa_gpu/run.sbatch
python3 bench/mpi/human/windows.py build/human-reference build/human-mpi/inputs build/human-mpi/windows.tsv
sbatch -A "$ROMEO_PROJECT" -p instant bench/gpu/wfa_gpu/run.sbatch build/human-mpi/inputs/matched-768.fa
sbatch -A "$ROMEO_PROJECT" -p instant --cpus-per-task=5 bench/gpu/wfa_gpu/run.sbatch --sars "$SARS_FASTA"
```

La variable `SARS_FASTA` désigne le FASTA des 64 accessions du manifeste. Le compteur de reprises du code amont s’écrit sur la sortie de diagnostic dès le mode optimisé. Les étapes de compilation et les révisions sont dans le lanceur. Le premier job utilisait le même programme avec le seul point d’entrée synthétique et φX174. [Lanceur](../../bench/gpu/wfa_gpu/run.sbatch) · [Script du premier job](839291/submitted-launcher.patch) · [Script du second job](839294/submitted-launcher.patch)
