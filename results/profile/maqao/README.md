# Profils MAQAO sur ROMEO

Sur un nœud AMD EPYC 9654, ONEVIEW associe l’échantillonnage des cycles par LProf à l’analyse statique des boucles par CQA. Deux programmes sollicitent le processeur différemment.

Le premier calcule avec NvPD la distance entre deux fenêtres du chromosome 21. La comparaison mesure l’effet de la disposition des indices et de la largeur des vecteurs, avec la construction de l’index dans le temps de chaque appel. [Chronométrages, profils et reproduction](dense/README.md) · [Lecture indirecte](../../../report/figures/maqao/boucle-initiale.png) · [Permutation en registres](../../../report/figures/maqao/boucle-permutation.png)

Le second distribue les 2 016 paires d’une matrice entre Edlib et WFA2. Le profil couvre le processus entier sur quatre cœurs, depuis l’échauffement jusqu’à la dernière répétition. [Profil du répartiteur](route/README.md) · [Classement des fonctions](../../../report/figures/maqao/fonctions-routage.png)

Les profils initiaux ont été établis avec MAQAO 2025.1.4. La relance avec 2026.1.1 compare les diagnostics sur des exécutables identiques et vérifie les avertissements de collecte. [NvPD](versions/README.md) · [Répartiteur Edlib et WFA2](route/versions/README.md)
