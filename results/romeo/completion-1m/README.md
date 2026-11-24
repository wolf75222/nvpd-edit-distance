# Mesures OpenMP à un million de bases

La campagne initiale avait mesuré la première paire de chacun des trois taux d’erreur avec l’algorithme par anti-diagonales et NvPD OpenMP. Le complément reprend les quinze paires pour les deux algorithmes. Il fournit ainsi vingt-quatre temps auparavant absents et six témoins répétés. Chaque ligne de la nouvelle série correspond à une mesure sur une paire. Les cinq paires d’un taux sont des entrées distinctes.

## Protocole

Le job 838440 utilise la révision `f08b4f6fcdf81de3a6ffcd60a8934cd1aca31c83`, GCC 11.5 et les options Release `-O3 -DNDEBUG -std=c++20 -fopenmp`. Le SHA-256 de l’exécutable est identique à celui des jobs historiques 837264, 837266, 837267 et 837268. La version NvPD mesurée est donc celle de la campagne initiale, avant la vectorisation introduite ultérieurement.

ROMEO a alloué 24 cœurs sur chacun de deux nœuds AMD EPYC 9654. Chaque vague correspond à une réplique. Sur un nœud, trois processus traitent simultanément les taux 2 %, 5 % et 10 % avec le même algorithme. L’autre nœud calcule les trois mêmes distances avec l’autre algorithme. Les nœuds échangent les algorithmes à chaque nouvelle réplique. Pour une même paire, les deux temps sont donc mesurés simultanément sur deux nœuds distincts. Le taux reste associé au même rang local dans chaque vague. Les trois processus d’un nœud utilisent des ensembles disjoints de huit cœurs. Les nœuds sont partagés. Les trois jobs historiques témoins avaient eux aussi tourné simultanément sur un nœud, c022.

Le temps mesuré entoure l’appel complet au calcul de distance, après la lecture du FASTA. Il inclut la construction de l’index de NvPD. La référence a exactement 1 000 000 de bases. Les mutants vont de 999 542 à 1 000 699 bases. Le débit dense utilise les dimensions réelles de chaque paire, `m × n / (temps × 10⁹)`. La mémoire est le pic RSS du processus. Une étape distincte recalcule chaque distance avec edlib et WFA2. Le collecteur vérifie le SHA-256 du FASTA, les deux longueurs, la graine et la distance avec les deux résultats historiques et les deux nouveaux résultats de bibliothèque.

## Résultats

| Calcul, 8 threads | Taux | Temps médian (s) | Min.–max. (s) | Débit médian (GCUPS) | Pic RSS maximal (Mio) |
|---|---:|---:|---:|---:|---:|
| Anti-diagonales | 2 % | 190,29 | 190,11–190,46 | 5,256 | 16,36 |
| Anti-diagonales | 5 % | 190,07 | 189,96–190,98 | 5,262 | 16,36 |
| Anti-diagonales | 10 % | 190,04 | 189,81–191,03 | 5,260 | 16,36 |
| NvPD initial | 2 % | 310,91 | 308,83–311,53 | 3,217 | 28,18 |
| NvPD initial | 5 % | 309,96 | 309,20–312,15 | 3,227 | 27,13 |
| NvPD initial | 10 % | 310,25 | 308,67–311,01 | 3,224 | 27,97 |

Sur les quinze paires, le rapport médian des temps appariés NvPD initial / anti-diagonales atteint 1,628. Les six témoins se rapprochent des mesures historiques à 0,12–0,17 % près pour les anti-diagonales et à 1,41–2,90 % près pour NvPD. Les témoins ont tous pris moins de temps lors de la nouvelle exécution. Les cinq vagues de calcul, la compilation et les vérifications ont occupé l’allocation pendant 27 min 34 s.

## Fichiers

- `trials.csv` contient les trente mesures à huit threads de cette campagne.
- `references.csv` contient les trente contrôles de bibliothèque sur les quinze paires.
- `coverage.tsv` joint les mesures historiques et nouvelles par algorithme et paire. Les six lignes témoins donnent le rapport du nouveau temps à l’ancien. Pour la réplique 1, le temps séquentiel historique permet aussi de calculer une accélération sur la même entrée.
- `summary.csv` résume uniquement les cinq paires de la nouvelle série pour chaque algorithme et taux. L’écart type porte sur la variation entre entrées.
- `placement.tsv` et `environment.txt` enregistrent l’affinité de chaque mesure OpenMP et les paramètres de compilation.
- `accounting.tsv` conserve les états et durées Slurm de chaque étape. `source-check.txt` consigne le contrôle des sources et `benchmark.sha256` identifie les versions du lanceur et du collecteur utilisées par le job.

Depuis la racine du dépôt, générer le corpus synthétique avec sa graine par défaut puis soumettre le job. Il extrait la révision initiale dans un répertoire temporaire et vérifie son manifeste de sources.

```bash
python3 scripts/synthetic.py --out data/synthetic
sbatch --account="$ROMEO_PROJECT" bench/complete_1m.sbatch
```

Les fichiers d’une nouvelle exécution se trouvent dans `results/romeo/completion-1m/jobs/<job>/`.
