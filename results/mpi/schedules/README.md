# Ordonnancement des comparaisons par paires sur deux nœuds EPYC

Les entrées sont 128 fragments de deux assemblages publics du chromosome 21 humain, GRCh38.p14 et T2T-CHM13v2.0. Seize loci espacés de 1,5 million de bases fournissent quatre fenêtres emboîtées par assemblage. Le générateur cherche la première position d’une grille déterministe à laquelle le début et les quatre fins sont ancrés sur la chaîne UCSC 23, brin positif. À chaque borne, les 41 bases de contexte sont identiques dans les deux assemblages. Les longueurs nominales de GRCh38 sont 512, 2 048, 4 096 et 12 288 bases. Les insertions et délétions de CHM13 sont conservées. Le manifeste donne les coordonnées, en base zéro avec fin exclue, et l’empreinte de chaque séquence. [Manifeste](windows.tsv) · [Extraction](../../../bench/mpi/schedules/cohort.py)

L’ensemble comprend 8 128 distances globales unitaires et représente 181 418 301 166 cellules de programmation dynamique. Il contient 448 comparaisons à l’intérieur d’un locus, dont 64 entre assemblages à même longueur nominale. Les 7 680 autres portent sur des loci distincts. Les 64 scores homologues vont de 0 à 467 (médiane 4), contre 244 à 12 234 entre loci (médiane 3 584). Ces dernières paires imposent un travail de calcul hétérogène. Leur distance n’est pas une mesure de divergence entre individus. Les cinq politiques reçoivent exactement ces mêmes séquences, dans l’ordre des loci puis par longueur croissante, avec les identifiants pour départager les longueurs égales.

## Mécanisme

Avec l’option MPI `matrix --schedule`, le triangle numérote d’abord les paires dans l’ordre lexicographique. La réduction restitue chaque distance à cette position. La sortie MPI conserve les lignes et colonnes dans l’ordre du FASTA pour toutes les politiques. Pour comparer aussi les deux ordres d’entrée, le banc vérifie chaque cellule par identifiants et conserve une empreinte commune du triangle trié par identifiant. Un compteur de tâches et un contrôle de tous les scores empêchent une paire dupliquée de masquer une paire omise. La réduction utilise le maximum et chaque rang initialise ses entrées non calculées à −1.

| Politique | Attribution des paires aux rangs |
| --- | --- |
| `static` | Blocs lexicographiques contigus, de mêmes effectifs à une paire près. |
| `cyclic` | Le rang *r* prend les numéros *r*, *r + P*, *r + 2P* pour *P* rangs. |
| `hilbert` | Tri des coordonnées du triangle selon l’indice de Hilbert dans un carré dont le côté est la puissance de deux immédiatement supérieure ou égale au nombre de séquences. Partage en blocs contigus de mêmes effectifs. |
| `hilbert-weighted` | Même courbe. Pour chaque rang, choix de la coupure dont le coût cumulé est le plus proche du quantile du travail total arrondi à l’entier inférieur. Le nombre de cellules donne le poids exact de la paire. |
| `steal` | Chaque rang commence par son propre bloc lexicographique. Tant qu’un voisin a des tâches non réservées, le rang peut lui en demander. Chaque victime gère sa file et attribue le travail une seule fois. |

Pour toutes les politiques, les huit fils d’un rang traitent les tâches par quatre et le noyau NvPD exploite SIMD à l’intérieur d’une paire. Pour le vol, les demandes MPI et les réponses de deux indices sont asynchrones. Le propriétaire prend ses propres tâches en tête de file et accorde au plus 32 tâches en queue. Le receveur conserve au plus une demande en vol et divise la réponse entre ses fils. Une réponse vide conduit au rang suivant. Un rang qui a épuisé tous les donneurs continue de répondre aux autres pendant une barrière non bloquante. Les appels MPI des fils sont sérialisés. Le calcul des cellules reste parallèle. [Solveur](../../../src/mpi/scheduled_matrix.hpp) · [Options](../../../src/mpi/mpi.cpp)

## Protocole

Un même lot de deux nœuds AMD EPYC 9654 héberge tous les essais, avec 24 rangs de huit cœurs physiques par nœud. Les placements sont contrôlés à chaque lancement. Un bloc lance chacune des 20 combinaisons (deux ordres × cinq politiques × un ou deux nœuds) dans un ordre tiré au sort. Chaque lancement crée de nouveaux processus MPI, réalise un calcul complet de chauffe puis chronomètre le deuxième. Le chronomètre commence après le chargement et la diffusion du FASTA. Il inclut la construction de la liste, son tri éventuel, le cumul des poids et la recherche des coupures pour Hilbert pondéré, la répartition, le calcul et la réduction du triangle. La borne finale est le temps du rang le plus lent. La collecte des diagnostics et l’écriture de la matrice ont lieu ensuite.

Le champ `gcups` divise le produit cumulé des longueurs par la durée en secondes et par un milliard. Le déséquilibre `max_mean_rank_cells` divise le plus grand travail en cellules d’un rang par la moyenne des rangs. La fraction `busy_fraction` rapporte la somme des durées de traitement des paquets par les fils à la capacité disponible (nombre de cœurs × temps de calcul distribué). Elle inclut le coût du noyau et la comptabilité à l’intérieur d’un paquet. Elle ne mesure pas la charge CPU du système. La mémoire est le maximum des pics RSS des rangs du lancement. Les agrégats donnent la moyenne, l’écart-type échantillonnal et l’intervalle bilatéral de Student à 95 % entre cinq blocs. L’accélération et l’efficacité sont calculées séparément dans chaque bloc avant agrégation.

## Mesures

Le job 839194 a exécuté les 100 lancers en 4 min 17 s. Les 812 800 scores chronométrés retrouvent les valeurs d’edlib et de WFA2. L’ordre inverse de la courbe et une deuxième implémentation des coupures pondérées ont aussi permis de vérifier les affectations des quatre politiques statiques, rang par rang, dans les deux ordres. La table contient le débit moyen ± l’écart-type à deux nœuds. [Essais et contrôles de placement](839194/trials.csv) · [Travail des 3 600 rangs échantillonnés](839194/rank_work.csv) · [Scores de référence](839194/scores.tsv) · [Statut du calcul](839194/accounting.psv)

| Politique | Par locus (GCUPS) | Par longueur (GCUPS) |
| --- | ---: | ---: |
| Blocs | 268,2 ± 2,7 | 106,5 ± 0,5 |
| Cyclique | 451,7 ± 5,2 | 419,6 ± 6,5 |
| Hilbert | 400,2 ± 3,6 | 107,2 ± 0,4 |
| Hilbert pondéré | 420,8 ± 5,5 | 587,4 ± 11,4 |
| Vol distribué | 435,2 ± 13,9 | 185,2 ± 13,0 |

En ordre de longueur, les coupures pondérées réduisent le maximum de cellules sur un rang de 6,830 à 1,038 fois la moyenne. Le temps à deux nœuds passe de 1,692 ± 0,006 à 0,309 ± 0,006 s. Le rapport apparié des temps vaut 5,480 (IC à 95 % de 5,331 à 5,629). La fraction du temps de cœur occupé par les lots monte de 13,56 à 75,98 %. Le pic RSS maximal des rangs de Hilbert pondéré est de 16,17 ± 0,10 Mio. À un nœud, la même méthode prend 0,607 ± 0,005 s. Son efficacité forte vaut 98,28 % (95,70 à 100,86 %). L’intervalle quantifie la variabilité entre blocs sur ces nœuds. [Statistiques par configuration](839194/summary.csv) · [Échelle forte](839194/scaling.csv) · [Gains appariés](839194/paired.csv)

Le vol distribué déplace 1 192 ± 245 paires entre les deux nœuds si les longueurs sont triées, contre 185 ± 37 en ordre de locus. Les rangs échangent les indices de tâches et conservent une copie locale des séquences. Le déplacement d’une tâche ne transfère donc pas son fragment d’ADN. En ordre de locus, le cyclique donne le plus grand débit. L’équilibrage des cellules entre rangs, la répartition entre fils et le coût de l’ordonnanceur contribuent tous au temps. Aucune mesure isolée des latences de message ou des défauts de cache n’a été réalisée. [Compte des vols par rang](839194/rank_work.csv) · [Comparaisons appariées à placement fixe](839194/paired.csv)

## Reproduction des entrées et des essais

Le script de préparation vérifie les deux chromosomes et la chaîne UCSC par SHA-256. Si les assemblages sont déjà présents, commencer au second appel Python. Les FASTA des deux ordres restent dans le répertoire de travail. La préparation, la double vérification des distances et le lancement complet emploient les commandes suivantes depuis la racine du dépôt.

```sh
python3 scripts/datasets/human_reference.py --whole-assemblies \
  --cache data/human/reference-cache --manifest data/human/manifest.tsv
python3 bench/mpi/schedules/cohort.py data/human/chr21-assemblies.fa \
  data/human/reference-cache data/mpi-schedules results/mpi/schedules/windows.tsv
job=$(sbatch --parsable --wait -A "$ROMEO_PROJECT" bench/mpi/schedules/run.sbatch)
python3 bench/mpi/schedules/summarize.py "results/mpi/schedules/$job"
```

Le contrôle des rangs sans paire, des séquences vides et des identifiants distincts à contenu identique tourne sans réseau. Le petit programme C++ contrôle la bijection de Hilbert, la continuité du parcours sur la grille complète et deux découpages pondérés de référence.

```sh
python3 bench/mpi/schedules/checks.py build/editdist-mpi build/editdist
mpicxx -std=c++20 -fopenmp -Iinclude -Isrc/mpi bench/mpi/schedules/cuts.cpp -o build/check-cuts
build/check-cuts
```

Une première implémentation par opérations atomiques MPI sur des fenêtres RMA a terminé les essais sur un seul nœud. Le premier calcul avec vol entre deux nœuds n’a pas abouti dans le job 839179, interrompu avant l'achèvement du bloc. Les campagnes suivantes emploient des échanges explicites entre rangs. Le pilote 839185 a exécuté les quatre premières politiques sur c036/c037. L’étude principale 839194 reprend les quatre témoins avec Hilbert pondéré sur c031/c035. Les intervalles de confiance de chaque expérience utilisent les cinq blocs de sa propre allocation. Les empreintes des dépendances sont communes aux deux campagnes. Le solveur de 839194 est celui du commit `5c16e21de152a7a2bf0c37b6b6fa61ae64b9d920`. Le CMake de compilation correspond au commit `b5d08d1`. [Sources communes](dependencies.sha256) · [Sortie pilote](839185/trials.csv) · [Comptabilité pilote](839185/accounting.psv) · [Différence de sources pour reproduire 839185 depuis 839194](839185/sources-from-839194.patch)
