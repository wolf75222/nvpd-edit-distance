# Comparaison sur cinq paires par configuration

Sur un nœud GH200 de ROMEO, le job Slurm 839826 compare le calcul du score sur 60 paires synthétiques. Le corpus couvre quatre longueurs et trois taux d’éditions introduites, avec cinq paires par configuration. Le noyau dense initial de NvPD utilise le GPU, Edlib et WFA2 utilisent chacun un cœur Grace.

Les trois méthodes traitent chaque paire dans trois blocs successifs. Leur ordre tourne d’un bloc au suivant. Dans chaque bloc, un processus par méthode lit le FASTA, effectue un alignement d’échauffement, puis chronomètre le deuxième alignement. Le fichier [divergence.csv](divergence.csv) contient 540 échauffements et 540 mesures. Les empreintes SHA-256 identifient les entrées et les scores concordent avec ceux de la campagne synthétique.

Pour chaque bibliothèque, le rapport par paire est la moyenne des trois quotients du temps de NvPD par celui de la bibliothèque dans le même bloc. Le fichier [summary.csv](summary.csv) donne, pour cinq paires, la moyenne de ces rapports, leur écart-type et la demi-largeur de l’intervalle de Student à 95 % (colonne `half95`).

Les [scripts du banc](../../../bench/repeated/divergence/) reproduisent la collecte et la figure. La [réservation Slurm](accounting.tsv) et les [versions des logiciels](environment.txt) décrivent les conditions de mesure.

La comparaison appariée du noyau dense et de la version en bande porte sur 25 paires. [Calcul en bande sur GH200](../nvpd-focus/)
