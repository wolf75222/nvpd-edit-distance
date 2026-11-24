# Campagne ROMEO du 24 novembre 2026

Les fichiers `cpu-*.csv` et `gpu-*.csv` contiennent les essais individuels. Une ligne correspond à un alignement dans un nouveau processus. Les mesures biologiques figurent dans les dossiers `cohort-*`. Les empreintes des entrées permettent de comparer les jobs portant sur la même paire.

| Entrée | Usage |
|---|---|
| [submissions.tsv](submissions.tsv) | Rejouer les dix commandes depuis la racine du projet, avec son compte Slurm. La colonne `expected_trials` fixe l’effectif attendu de chaque job. |
| [accounting.tsv](accounting.tsv) | Vérifier les états Slurm, les codes de sortie et les ressources allouées. |
| [summary.csv](summary.csv) | Lire le temps médian par configuration, le débit, l’accélération et les pics de mémoire. Les accélérations utilisent la référence séquentielle mesurée sur la même paire. |
| [sources.sha256](sources.sha256) | Vérifier les sources C/C++/CUDA, CMake et le test compilés sur ROMEO. |
| [cohort-837271](cohort-837271/) | Consulter la matrice φX174 commune à NvPD-CUDA, edlib et WFA2, les deux arbres et les temps sur le nœud Arm/GPU. |
| [cohort-837261](cohort-837261/) | Consulter les temps edlib et WFA2 pour la même matrice sur le nœud x86. |

Ces commandes recalculent la synthèse depuis la racine du projet et vérifient les effectifs.

```sh
python3 bench/campaign/summarize.py --plan results/romeo/submissions.tsv \
  results/romeo/cpu-*.csv results/romeo/gpu-*.csv > summary-check.csv
cmp summary-check.csv results/romeo/summary.csv
git archive f08b4f6 -o /tmp/nvpd-initial.tar
# Les empreintes sources.sha256 correspondent aux sources de cette archive.
```

Les trois taux de mutation sont 2, 5 et 10 %. À 5 000 et 20 000 bases, la première paire de chaque taux sert à mesurer la montée en charge de 1 à 8 threads. À un million de bases, les trois méthodes denses CPU traitent la première paire de chaque taux. Le GPU, edlib et WFA2 traitent les cinq paires de chaque taux. La campagne utilise des nœuds partagés. Les temps décrivent la [révision initiale](https://github.com/wolf75222/nvpd-edit-distance/commit/f08b4f6fcdf81de3a6ffcd60a8934cd1aca31c83).

La [campagne OpenMP sur quinze paires d’un million de bases](completion-1m/) complète les mesures de la révision initiale.
