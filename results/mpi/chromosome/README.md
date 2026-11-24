# Matrices de fenêtres du chromosome 1

La cohorte comporte 768 intervalles disjoints du chromosome 1 de CHM13 v1.1. Les coordonnées du [manifeste](manifest.tsv) commencent à zéro et excluent la borne droite. Chacune des longueurs de 512, 1 024, 2 048 et 3 072 bases est représentée par 192 fenêtres. Le sous-ensemble de 544 loci conserve 136 fenêtres de chaque longueur. Il s’agit de régions d’un assemblage individuel. L’unité biologique du panel d’haplotypes est décrite [dans sa propre étude](../../human/packed-cohort/).

Dans chacun des cinq blocs, le banc permute les configurations, lance de nouveaux processus MPI, chauffe chaque configuration une fois et mesure la matrice suivante. Le temps inclut la construction des index, la répartition des paires et le retour des scores au rang racine. La lecture du FASTA et l’écriture du TSV encadrent le chronomètre. Le RSS est le maximum des pics par rang. Le compteur de cellules équivalentes est la somme des produits de longueurs. Il décrit la taille du problème également pour les solveurs qui visitent moins de cellules.

Les [essais individuels](839322/trials.csv) conservent l’empreinte de chaque matrice, le nombre de cœurs physiques et les paires traitées par le test de sous-séquence. Le [lanceur Slurm](../../../bench/mpi/chromosome/run.sbatch) réserve deux nœuds EPYC 9654.

Le [code exécuté](https://github.com/wolf75222/nvpd-edit-distance/tree/d83fed7) est associé aux empreintes du binaire et des sources. Edlib recalcule les 294 528 paires dans un processus indépendant et retrouve les 39 679 distances égales à la différence de longueur parmi les paires admissibles. [Contrôle de la matrice](839322/edlib-oracle.json)
