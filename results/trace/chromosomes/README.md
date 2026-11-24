# Reconstruction par fragments sur EPYC

Le jeu compare deux révisions de CHM13 sur les chromosomes 1, 18 et 20. La quatrième paire couvre deux millions de bases du chromosome 21 entre GRCh38 et CHM13. Les quatre distances de référence sont 176, 5, 30 et 5 599 éditions. Les empreintes des entrées figurent dans chaque dossier de mesure.

La campagne 839187 compare le score global et le CIGAR complet sur un nœud EPYC 9654 de ROMEO. Chaque bloc mélange l’ordre des configurations. Un processus lit le FASTA, effectue un appel d’échauffement puis chronomètre un deuxième appel avec les entrées en mémoire. Le rejeu du CIGAR se déroule après le chronométrage. Chaque configuration dispose d’un à huit cœurs physiques et d’un délai de 120 secondes pour le processus entier, échauffement compris. Le RSS maximal inclut le FASTA et les allocations de l’échauffement. Le processus relève le RSS après le rejeu du chemin.

Le CIGAR d’edlib sur le chromosome 1 atteint le plafond de 120 s du protocole principal. L’essai 839200 reprend cette paire avec un plafond de 720 s par processus. Chacune des cinq tâches relance edlib, WFA2, BiWFA et le solveur à huit cœurs dans son propre bloc. Jusqu’à trois tâches ont partagé le nœud. Les rapports de temps du complément emploient ses propres témoins. [Comparaison longue sur le chromosome 1](839200/summary.csv)

Le pilote compare le rejet de toutes les coupes dès qu’une échoue (`split-all`) à la conservation des coupes vérifiées (`split`). La campagne 839193 isole le rôle du découpage. Son témoin `trim-only` retire les extrémités communes et utilise les mêmes noyaux sur un seul fragment. Le solveur NvPD traite les petits rectangles. BiWFA résout les autres en mémoire proportionnelle au score pour ses fronts. Le champ `wavefront_leaves` compte les appels à BiWFA dans les fragments, y compris les recalculs. Un repli direct vers BiWFA en mode `split-all` est signalé par `retried`.

Les temps sont exprimés en secondes et le RSS en Kio. Les moyennes, écarts-types et intervalles de Student à 95 % portent sur cinq blocs. Une accélération est le rapport des temps mesurés dans un même bloc. Les données brutes conservent le nombre de coupes proposées et retenues. [Comparaison des solveurs](839187/trials.csv) · [Ablation du découpage](839193/summary.csv)

Depuis la racine, préparer les FASTA avec les scripts indiqués dans le README du projet. Le pilote prend le dossier des chromosomes complets et le chemin absolu du FASTA de la fenêtre. Sur ROMEO, le nœud de calcul doit aussi pouvoir lire les archives de dépendances sous `deps/`.

```sh
sbatch -A "$ROMEO_PROJECT" bench/cpu/chromosome-trace/run.sbatch \
  "$PWD/data/human/chm13-panel" "$PWD/data/human/chr21-2m.fa"
CHPS_TRACE_ABLATION=1 sbatch -A "$ROMEO_PROJECT" bench/cpu/chromosome-trace/run.sbatch \
  "$PWD/data/human/chm13-panel" "$PWD/data/human/chr21-2m.fa"
sbatch --array=1-5%3 --time=00:15:00 --mem=8G -A "$ROMEO_PROJECT" \
  bench/cpu/chromosome-trace/run.sbatch "$PWD/data/human/chm13-panel" "$PWD/data/human/chr21-2m.fa"
python3 bench/cpu/chromosome-trace/summarize.py results/trace/chromosomes/839187
```
