# Répartition des paires entre cœurs

Sur ROMEO, le job 839293 a calculé la matrice des 64 accessions sur 32 cœurs physiques d’un EPYC 9654. Le calcul par lignes confie à une tâche toutes les paires \((i,j)\) avec \(j>i\). Le second parcours partage chaque ligne entre 32 groupes. Le groupe \(g\) parcourt les colonnes \(j=i+1+g+32k\). OpenMP attribue les tâches disponibles à la demande. Les deux exécutables utilisent les mêmes noyaux et GCC 11.5.0 avec `-O3 -march=native`. [Protocole](../../../bench/cpu/matrix-schedule/run.sbatch) · [Hôte et empreintes](839293/environment.txt)

Chaque couple méthode/parcours reçoit une chauffe. Les cinq blocs de mesure alternent l’ordre des deux parcours pour chaque méthode. Le chronométrage couvre le lancement du processus, la lecture du FASTA et l’émission de la matrice. Les 36 matrices reproduisent les 2 016 distances publiées. Le pic de mémoire est le maximum résident du processus, relevé par GNU time. [Mesures](839293/trials.csv) · [Matrice](../sars-edit/matrix.tsv)

| Bibliothèque | Lignes, moyenne ± écart-type [s] | Groupes, moyenne ± écart-type [s] | Accélération appariée [IC à 95 %] |
|---|---:|---:|---:|
| edlib | 1,6450 ± 0,0006 | 0,3432 ± 0,0008 | 4,79 [4,78 à 4,81] |
| WFA2-lib | 3,6155 ± 0,0023 | 0,7747 ± 0,0019 | 4,67 [4,65 à 4,68] |
| BiWFA | 2,3137 ± 0,0046 | 0,4637 ± 0,0040 | 4,99 [4,93 à 5,05] |

Pour chaque bloc \(b\), on calcule \(z_b=\log(t_{\mathrm{row},b}/t_{\mathrm{group},b})\). L’intervalle est \(\exp(\bar z\pm t_{0.975,4}s_z/\sqrt{5})\). Sa portée est celle des répétitions sur ce nœud et cette cohorte. [Ordre des passages et contrôles](../../../bench/cpu/matrix-schedule/compare.py)
