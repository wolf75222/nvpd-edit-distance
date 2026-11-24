# NvPD dans une bande certifiée

`editdist pair --algorithm nvpd-adaptive --input paire.fa` renvoie le coût global à pénalités unitaires. Les caractères IUPAC sont littéraux. Le calcul de chaque ligne exploite l'indépendance des colonnes de NvPD et conserve deux lignes compactes. La commande `matrix` peut répartir les paires sur une équipe OpenMP. [Implémentation](../../src/cpu/nvpd_adaptive.cpp) · [Boucle SIMD](../../src/cpu/nvpd_band_row.hpp)

## Exactitude

On retire les extrémités communes et oriente les chaînes restantes pour avoir $m\le n$. Pour un budget $k$, une cellule de diagonale $e=j-i$ exige au moins $|e|$ insertions ou suppressions avant elle et $|n-m-e|$ après elle. Le corridor utilisé est donc $-h\le e\le n-m+h$, où $h=\lfloor(k-(n-m))/2\rfloor$. Deux sentinelles au coût $k+1$ marquent les bornes de la ligne précédente. Pour chaque colonne $j$, le prédécesseur issu du dernier match $t$ est la case $t-1$ de cette ligne. S’il manque dans le corridor ou si le symbole est absent, sa contribution est saturée. Les autres prédécesseurs sont ceux de la récurrence NvPD, sur la même colonne et la diagonale.

Un chemin global de coût au plus $k$ reste entier dans le corridor. Pour toute cellule de ce chemin, un préfixe optimal suivi de son suffixe respecte encore le budget. On peut donc choisir le prédécesseur optimal donné par la récurrence de NvPD et développer son saut en un match suivi d’insertions. Chaque cellule traversée reste admissible, y compris la source $t-1$. Par induction sur les lignes, les préfixes nécessaires au chemin optimal conservent leur coût. Réciproquement, une valeur terminale au plus égale à $k$ provient de prédécesseurs réels. Une sentinelle saturée ne peut la produire. Une ligne dont toutes les valeurs dépassent le budget permet l’arrêt anticipé de l’essai. [NvPD de Sadiq et al.](https://doi.org/10.1007/s10586-019-02962-w) · [Calcul par seuil d’Ukkonen](https://doi.org/10.1016/S0019-9958(85)80046-2)

Le déficit total de symboles en faveur de la chaîne la plus longue donne une borne inférieure. Un chemin qui substitue les désaccords positionnels puis insère le reste donne une borne supérieure. Si elles coïncident, le score est connu. Sinon, on construit une seule fois l’index de dernière occurrence, puis on double le budget jusqu’au succès en le plafonnant à la borne supérieure. Pour un alphabet de $s$ symboles présents du côté des lignes et une distance $d$, le travail est $O(sn+m(d+1))$ sur les chaînes résiduelles, auquel s’ajoute le balayage initial linéaire dans les longueurs d’origine. L’espace des tableaux est $O(sn+d+1)$, y compris l’index. La version dense et les variantes par fronts d’onde restent accessibles.

## Mesures sur ROMEO

Le profil antérieur de NvPD dense sur un carré de 20 000 bases et de distance 991 compte 853 millions d’instructions après vectorisation, pour 400 millions de cellules. Le nouveau banc mesure l’effet d’une réduction du nombre de cellules. [Profil dense antérieur](../profile/cpu/counters.csv)

Sept cas couvrent deux paires synthétiques du générateur du projet, une paire de chaînes indépendantes et quatre fenêtres humaines emboîtées. Les bornes humaines ont été relevées sur la chaîne UCSC 23 du chromosome 21, brin positif, dans les assemblages GRCh38 et CHM13. Le pilote confirme leurs contextes de 41 bases dans la paire source de 4 Mb. Chaque résultat est comparé à un appel edlib indépendant. La version dense est mesurée jusqu'à 250 000 bases. La grille à un million contiendrait plus de mille milliards de cellules. [Génération des entrées](../../bench/cpu/adaptive/cases.py) · [Coordonnées, empreintes et scores](839188/cases.csv)

Le job 839188 réserve un nœud EPYC, quatre cœurs pour la compilation et 16 Gio sur la partition `instant` pendant 55 min au maximum. Chaque alignement utilise un seul cœur physique fixé. Dans chaque bloc, l’ordre des couples méthode–entrée est tiré avec une graine fixe. Un processus neuf lit le FASTA, effectue un appel complet d’échauffement, puis chronomètre un second appel complet. Le chronomètre inclut donc le prétraitement, la construction d’index, les allocations et les élargissements. Le pilote impose 180 s et 6 Gio d’espace adressable par processus. Un gain apparié est le temps du témoin divisé par le temps de la bande SIMD sur le même cas et dans le même bloc. Les résumés donnent la moyenne, l’écart-type et l’IC de Student bilatéral à 95 % sur les cinq processus ou cinq rapports appariés. [Temps bruts](839188/trials.csv) · [Temps et gains](839188/summary.csv) · [Cellules et tableaux](839188/work.csv)

Les tableaux de travail principaux comprennent une entrée entière de dernière occurrence par symbole et par colonne, puis deux lignes de scores. Sur les cas mesurés, chacun des quatre symboles A, C, G et T est présent. Le tableau de travail indique les octets de ces tableaux d’entiers, calculés après le retrait des extrémités communes pour la bande. L’allocation de l’exécutable, les chaînes FASTA et les autres métadonnées relèvent du RSS du processus. À un million de bases, l’index reste le premier poste mémoire même avec des lignes de largeur voisine de 2 700 cellules.

Le contrôle scalaire utilise le même algorithme et visite les mêmes cellules. Le script remplace les minima imbriqués par l’équivalent `std::min({a,b,c})` avant compilation dans une copie temporaire. Sur le GCC 11.5.0 mesuré, cette écriture scalaire empêche la vectorisation de la lecture indirecte. L’assembleur de production contient `vpgatherdd` sur huit entiers de 32 bits. Les deux binaires sont comparés dans un même bloc et sur le même cœur. Le script de reproduction donne un nom distinct à la fonction temporaire du contrôle. Sur cet EPYC, le GCC utilisé produit un objet entier identique à celui obtenu avec le nom initial, car la fonction est incorporée à son appelant. [Empreintes des deux objets](839188/control-object.sha256) · [Script soumis](839188/submitted.sbatch) [Diagnostic GCC et assembleur](839188/vectorization.txt)

## Modèle mémoire et profilage

Les deux groupes de compteurs `perf stat` sont activés et désactivés par FIFOs acquittés autour d’un appel échauffé. Le groupe cycles, instructions et branches s’exécute séparément du groupe chargements L1 et remplissages. Les octets demandés par une mise à jour comprennent une lecture de l’index, deux prédécesseurs contigus, un prédécesseur indirect et une écriture, soit 20 octets d’entiers dans l’expression scalaire. Plusieurs de ces adresses se recouvrent entre colonnes voisines. Les octets utiles de la triade représentent deux entiers lus et un écrit par élément. La lecture linéaire consomme un entier par élément. [Événements bruts et fraction de temps active](839188/counters.csv) · [Étalons de débit](839188/calibration.csv)

Pour le modèle de lecture externe, le volume estimé vaut $Q=64(C_{0x4844}+C_{0xff72})$ octets. Le premier compteur observe les remplissages L1 classés DRAM ou E/S mappée en mémoire, locale ou distante. Le second compte les prélectures qui manquent le L2 et le L3. L’estimateur porte sur les lectures et dépend du chemin de prélecture. Nous le comparons aux octets connus d’une lecture de 192 Mio qui dépasse le L3 de 32 Mio du groupe de cœurs. Le contrôleur mémoire et les écritures n’entrent pas dans ces compteurs. Le repère de débit est $\min(P_{L1},B_\mathrm{lecture}\,U/Q)$, avec $U$ le nombre de cellules réellement évaluées, $P_{L1}$ le débit de la récurrence dans des tableaux de 12 Kio, et $B_\mathrm{lecture}$ le débit de la lecture linéaire de 192 Mio. Les performances incluent aussi les barrières éventuelles, les recherches d’index et les changements de budget. [Table des événements Zen 4 de Linux v6.17](https://github.com/torvalds/linux/blob/v6.17/tools/perf/pmu-events/arch/x86/amdzen4/cache.json) · [Guide AMD uProf](https://www.amd.com/content/dam/amd/en/documents/developer/version-4-1-documents/uprof/uprof-ug-rev-4.1.pdf) · [Modèle Roofline](https://doi.org/10.1145/1498765.1498785) · [Indicateurs PMU par cas](839188/pmu-summary.csv)

## Résultats

Le job 839188 se termine en 627 secondes. Les 240 observations de temps reproduisent le score de référence. Les temps du tableau sont des moyennes en millisecondes. Le gain est la moyenne des rapports appariés NvPD dense / NvPD adaptatif. Les écarts-types et les intervalles à 95 % de toutes les cellules se trouvent dans le résumé CSV.

| Paire | Score | NvPD dense | NvPD adaptatif | Bande classique | edlib | WFA2 | BiWFA | Gain dense/adaptatif |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Synthétique 2 % | 396 | 166,253 | 5,905 | 57,547 | 1,033 | 0,226 | 0,268 | 28,15× |
| Synthétique 10 % | 1908 | 164,340 | 21,713 | 291,944 | 2,132 | 3,810 | 2,289 | 7,57× |
| Indépendantes | 10373 | 165,350 | 116,383 | 979,250 | 18,042 | 86,410 | 57,249 | 1,42× |
| Chr21 25 k | 34 | 255,436 | 1,911 | 7,640 | 0,487 | 0,038 | 0,149 | 133,69× |
| Chr21 100 k | 310 | 4108,371 | 22,391 | 109,325 | 6,152 | 0,219 | 0,378 | 183,49× |
| Chr21 250 k | 455 | 25115,794 | 97,651 | 457,518 | 10,137 | 0,519 | 0,793 | 257,20× |
| Chr21 1 M | 2336 | — | 1608,383 | 10484,604 | 242,858 | 7,070 | 5,414 | — |

Sur la fenêtre de 250 k, l’IC à 95 % du gain de 257,20 va de 256,58 à 257,82. Le contrôle scalaire évalue les mêmes cellules et prend 3,15 fois plus de temps. Sur les chaînes indépendantes de 20 k, les 271,23 millions de cellules visitées représentent 67,81 % de la grille et le gain sur le dense descend à 1,42. La bande analytique et le retrait des extrémités communes modifient tous deux le travail. Sur la fenêtre de 100 k, le retrait porte notamment sur 16 822 bases identiques du suffixe.

En lecture contiguë, le cœur atteint 30,64 ± 0,21 Go/s sur 192 Mio. La triade délivre 35,46 ± 0,46 Go/s utiles. L’estimateur matériel reproduit 95,85 ± 0,66 % du volume de lecture connu. Pour les trois alignements instrumentés, l’intensité estimée va de 133 à 461 cellules par octet lu hors cache. Le débit varie de 1,94 à 2,33 milliards de cellules par seconde. L’étalon de la boucle en L1 donne 2,17 ± 0,34 milliards de cellules par seconde. Son niveau dépend de la longueur de la ligne et du motif des derniers matchs. Ces valeurs utilisent l’écart-type entre les cinq processus.

## Reproduction

Préparer les FASTA avec les générateurs synthétique et humain du dépôt. Les archives épinglées de WFA2 et edlib peuvent être déposées dans `deps/` avant l’exécution hors ligne. Sur une autre installation, `PERF_BIN` peut désigner le chemin de l’exécutable Linux perf. La compilation et l’exécution Slurm sont décrites dans le [script du job](../../bench/cpu/adaptive/run.sbatch).

```sh
sbatch --account=r260105 --partition=instant bench/cpu/adaptive/run.sbatch
python3 bench/cpu/adaptive/analyze.py results/adaptive/839188
```
