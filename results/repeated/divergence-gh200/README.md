# Comparaison sur cinq paires par configuration

Le job Slurm 839826 utilise un nœud GH200 de ROMEO. Les 60 paires du corpus synthétique couvrent quatre longueurs et trois taux d’éditions introduites. Pour chaque paire, trois blocs comparent NvPD-CUDA, Edlib et WFA2. Chaque méthode est exécutée dans un processus distinct, avec un échauffement suivi d’un appel chronométré. NvPD utilise le GPU et les bibliothèques un cœur Grace. La lecture du FASTA précède le chronométrage.

Le fichier `divergence.csv` contient les 1 080 appels, échauffements compris. Les scores ont été confrontés aux références de la campagne synthétique et les empreintes SHA-256 aux fichiers d’entrée. Dans `summary.csv`, chaque rapport par paire est la moyenne des trois rapports appariés. La moyenne, l’écart-type et l’intervalle de Student à 95 % portent sur les cinq paires de chaque configuration.

La collecte et la figure sont définies dans `bench/repeated/divergence/`. Le fichier `accounting.tsv` décrit la réservation et `environment.txt` les versions utilisées.
