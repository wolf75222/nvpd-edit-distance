# Échelle MPI sur les fenêtres homologues du chromosome humain 21

Les enregistrements adjacents du FASTA représentent le même locus dans GRCh38.p14 et T2T-CHM13v2.0. Le mode `editdist-mpi pairs` calcule une distance globale entre les deux enregistrements successifs, puis écrit les scores dans l'ordre du FASTA. Les rangs MPI se répartissent les paires de façon cyclique. Dans chaque rang, les fils OpenMP se partagent les paires attribuées. Une seule réduction rassemble leurs scores. [Implémentation](../../../src/mpi/mpi.cpp) · [Test de distribution et de mémoire](../../../bench/mpi/human/checks.py)

## Origine des fenêtres

Le générateur contrôle les empreintes des deux références chromosomiques et de la chaîne de correspondance UCSC `hg38ToHs1`. Il retient la chaîne 23 du chromosome 21, orientée sur le brin positif. Les deux bornes de chaque fenêtre doivent se situer dans un bloc aligné et leurs voisinages de 41 bases doivent être identiques dans les deux assemblages. Les fenêtres contenant un N sont écartées. À l'intérieur, les insertions et délétions restent présentes dans les séquences comparées. Les coordonnées du manifeste commencent à zéro et excluent la borne de fin. La sélection privilégie des régions alignables et porte sur deux assemblages de référence. [GRCh38, NC_000021.9](https://www.ncbi.nlm.nih.gov/nuccore/NC_000021.9) · [CHM13, NC_060945.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_060945.1) · [Chaînes de correspondance UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/) · [Manifeste](windows.tsv)

Les 768 fenêtres disjointes mesurent chacune 16 384 bases sur GRCh38 et de 5 297 à 19 840 bases sur CHM13. Leurs scores vont de 0 à 11 107 éditions, avec une médiane de 28 et 15 scores nuls. La première moitié du FASTA reprend une fenêtre sur deux de la grille géographique. Les 384 fenêtres couvrent donc la même région du chromosome que le lot complet. Leur travail dense totalise 103 269 613 568 cellules, celui des 768 fenêtres 206 353 809 408. Deux fenêtres plus longues examinent séparément le partage des colonnes d'une paire entre rangs. Sur 65 536 × 65 544 et 131 072 × 131 048 bases, les scores valent 120 et 2 710. Ces deux comparaisons ne s'ajoutent pas au travail du lot. Les 770 distances de référence proviennent d'edlib en mode global. Un second calcul WFA2-lib les a reproduites individuellement. [Scores par locus](scores.tsv) · [Vérification indépendante](verification.txt)

## Mesure

Le job ROMEO 838908 a réservé deux nœuds EPYC 9654. Sur chaque nœud actif, 24 rangs MPI occupent chacun huit cœurs physiques. Chaque bloc exécute les 15 configurations dans un ordre pseudo-aléatoire. Chaque lancement MPI effectue un alignement complet d'échauffement puis une mesure. L'horloge démarre après la lecture du FASTA, sa diffusion et une barrière. Elle inclut la construction de l'index NvPD, le calcul et le retour du score aux rangs concernés. La valeur publiée est le maximum des temps de rang. Le débit est le nombre exact de cellules divisé par ce temps. La mémoire est le maximum parmi les rangs des pics RSS de processus mesurés par `getrusage`. Elle comprend les entrées, l'échauffement et les bibliothèques. [Pilote de mesure](../../../bench/mpi/human/measure.py) · [Mesures brutes](trials.csv) · [Environnement](environment.txt)

Les moyennes, écarts-types d'échantillon et intervalles bilatéraux de Student à 95 % portent sur cinq blocs. Chaque accélération se calcule à l'intérieur du bloc avant l'agrégation. Pour un travail fixe, l'efficacité à deux nœuds vaut (T₁/T₂)/2. Pour le doublement de 384 à 768 loci, l'efficacité en débit par cœur vaut (W₂/W₁) × (T₁/T₂)/2 avec W₂/W₁ = 1,9982045. Le tableau fournit aussi le rapport brut T₁/T₂. Dans le calcul par bandes sur une paire, H=0 désigne l'échange collectif à chaque ligne, H=32 et H=128 les hauteurs des paquets de bord entre voisins. [Durées, débits et RSS](summary.csv) · [Comparaisons appariées](scaling.csv)

## Résultats sur ROMEO

Sur le même lot de 768 paires, un puis deux nœuds prennent 652,385 ± 19,337 et 404,446 ± 24,161 ms. Le débit moyen passe de 316,5 ± 9,4 à 511,6 ± 29,4 GCUPS. L'accélération appariée vaut 1,617 (IC à 95 % de 1,503 à 1,731), soit une efficacité forte de 80,84 %. En doublant le nombre de loci de 384 à 768, l'efficacité faible corrigée du travail atteint 80,89 % (72,80 à 88,97 %). Le pic RSS maximal par rang pour le lot de 768 vaut 44,96 ± 1,00 Mio sur un nœud et 45,12 ± 0,28 Mio sur deux. La diffusion des entrées donne à chaque rang une copie du lot complet. [Durées, débits et mémoire](summary.csv) · [Rapports appariés](scaling.csv)

Pour la paire de 131 072 bases de référence, le collectif à chaque ligne demande 6,672 ± 0,182 s sur deux nœuds. Les paquets de 128 lignes ramènent ce temps à 116,892 ± 1,440 ms, contre 159,522 ± 3,911 ms avec ce même pipeline sur un nœud. Le gain du pipeline à ressources identiques vaut 57,09 (54,79 à 59,38). Son accélération de un à deux nœuds vaut 1,365 (1,317 à 1,412). Sur la paire de 65 536 bases, les deux comparaisons donnent respectivement 74,12 (73,15 à 75,10) et 1,193 (1,169 à 1,218). H=32 présente davantage de dispersion sur certaines configurations, notamment sur la plus longue paire à deux nœuds. Ces résultats concernent les deux loci mesurés, la même paire de nœuds et le temps de calcul distribué. [Résultats par configuration](summary.csv) · [Statut Slurm](accounting.psv)

## Reproduction

Depuis la racine du projet, le cache peut contenir les deux références et la chaîne déjà utilisées par les calculs sur chromosomes humains. Les séquences restent dans le répertoire de travail.

```sh
python3 bench/mpi/human/windows.py build/human-cache data/mpi-human results/mpi/human/windows.tsv
job=$(sbatch --parsable --wait -A "$ROMEO_PROJECT" bench/mpi/human/run.sbatch)
python3 bench/mpi/human/summarize.py results/mpi/human/"$job"
python3 bench/mpi/human/verify.py data/mpi-human results/mpi/human/"$job"/scores.tsv build/editdist
```
