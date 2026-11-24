# NvPD distribué par bandes de colonnes

La commande `editdist-mpi pair --input PAIRE.fa --threads 8` utilise une opération collective par ligne. L’option `--pipeline-rows H` choisit le pipeline, pour tout entier positif H représentable par `int`. Le score et le nom de la méthode sont fournis en JSON. Les deux implémentations acceptent les mêmes chaînes et comparent littéralement leurs symboles IUPAC.

## Récurrence locale

Considérons une bande de colonnes globales de s+1 à s+w. À la ligne i, posons P[x]=D(i−1,s+x) et C[x]=D(i,s+x). Le voisin de gauche fournit C[0]. Dérouler la récurrence de Wagner–Fischer pour les x dernières colonnes donne trois familles de chemins.

1. Les x derniers mouvements sont des insertions depuis C[0], pour un coût C[0]+x.
2. La dernière opération qui consomme A[i] est une délétion en colonne locale t. Son coût suivi d’insertions est P[t]+1+x−t. Le minimum sur 1≤t≤x vaut P[x]+1, car P[x]≤P[t]+x−t.
3. La dernière opération qui consomme A[i] est diagonale en t. Son coût est P[t−1]+[A[i]≠B[s+t]]+x−t. Pour une égalité, le dernier t admissible, noté k, minimise P[t−1]+x−t par la même inégalité sur P. Les substitutions coûtent au moins P[x−1]+1. Si les caractères en x diffèrent, ce coût est atteint en x. S’ils sont égaux, le terme de correspondance en x vaut déjà P[x−1].

On obtient ainsi C[x]=min(P[x]+1, P[x−1]+1, C[0]+x, P[k−1]+x−k), en omettant le dernier terme si aucun caractère de cette bande ne convient. Toutes les colonnes locales peuvent alors avancer en parallèle. La première bande a C[0]=i et chaque bande démarre avec P[x]=s+x à la ligne zéro. Ces deux initialisations, la récurrence et la transmission des bords exacts établissent l’exactitude pour toutes les bandes et toutes les lignes.

## Communications et travail

Un message transporte h valeurs de bord successives, sauf le dernier paquet éventuellement plus court. La bande j reçoit du rang j−1 et envoie au rang j+1. Les valeurs précédentes de chaque bande restent locales. Les rangs sans colonne participent à la diffusion finale du score. Pour m lignes, r bandes non vides et B=⌈m/h⌉, le pipeline envoie (r−1)B messages contenant au total (r−1)m entiers, puis diffuse un score. Le témoin appelle m fois `MPI_Allgather`, ce qui dépose 2r entiers par rang et par ligne lorsque toutes les bandes sont non vides. Les deux méthodes calculent mn cases. Chaque rang non vide du pipeline conserve deux lignes locales de w+1 entiers, l’index des dernières correspondances locales et deux paquets de min(h,m) entiers, en plus des chaînes répliquées et du tableau commun des bornes de répartition (un entier par rang, plus la borne finale).

Pour des tuiles de même coût et des communications instantanées, le chemin critique traverse B+r−1 tuiles. La phase de remplissage et la phase de vidange limitent donc l’occupation lorsque B est petit. Le nombre de travailleurs attachés à des cœurs physiques décrit le placement. Plusieurs travailleurs peuvent attendre un paquet de leur voisin.

## Protocole ROMEO

Le job 838796 utilise deux nœuds EPYC exclusifs, 24 rangs par nœud et 8 threads OpenMP par rang. Les essais à un nœud se déroulent sur le premier de ces deux nœuds. Les deux paires sont les entrées synthétiques des campagnes MPI précédentes. Edlib vérifie la distance et le pilote contrôle l’empreinte SHA-256 et le placement. Le chronomètre couvre l’indexation locale, les communications, les mises à jour et la diffusion finale. Il exclut le chargement du FASTA et la distribution initiale des chaînes. Le temps publié est le maximum des durées mesurées sur les rangs.

Chaque bloc tire au hasard l’ordre des douze configurations (deux paires, un ou deux nœuds et trois modes). Un nouveau lancement MPI effectue deux alignements d’échauffement et un alignement mesuré. Les sept blocs fournissent la moyenne, l’écart-type échantillonnal et l’IC de Student à 95 % sur sept observations. Les accélérations sont calculées dans chaque bloc avant agrégation. Le contrôle PMPI reproduit les douze configurations dans des processus distincts des processus chronométrés.

Depuis un clone sur ROMEO, lancer `sbatch -A r260105 -p instant bench/mpi/repeated.sbatch pipeline`. Le script utilise CMake, GCC et Open MPI, compile les deux exécutables et teste le pipeline contre la DP séquentielle et edlib. Les fichiers se trouvent sous `results/mpi/pipeline/$SLURM_JOB_ID/`. [Noyau et CLI](../../../src/mpi/mpi.cpp) · [Pilote](../../../bench/mpi/pipeline/measure.py) · [Vérifications](../../../bench/mpi/mpi_checks.py)

Les observations individuelles incluent le bloc, l’ordre de lancement, le placement échantillonné, les deux longueurs et la distance. La sonde PMPI porte exclusivement sur les collectifs `Allgather`. Le nombre de messages point-à-point ci-dessus découle du protocole. [Essais](838796/trials.csv) · [Résumé et IC](838796/summary.csv) · [Nombre d’Allgather par rang](838796/allgather.csv) · [Validation](838796/validation.txt) · [Environnement](838796/environment.txt)

## Temps et rapports appariés

Le tableau donne les durées en millisecondes (moyenne ± écart-type des sept processus MPI). Le gain du pipeline par rapport au témoin est calculé à nombre de nœuds fixé. Le score vaut 431 pour le carré de 20 000 × 20 037 et 998 979 pour la paire large de 1 024 × 1 000 003.

| Paire | Nœuds | Allgather par ligne | Pipeline h=8 | Pipeline h=32 |
| --- | ---: | ---: | ---: | ---: |
| Carré | 1 | 94,988 ± 3,515 | 13,833 ± 0,091 | 11,381 ± 0,104 |
| Carré | 2 | 996,885 ± 14,981 | 19,533 ± 0,577 | 13,074 ± 0,425 |
| Large | 1 | 14,843 ± 0,287 | 7,979 ± 0,092 | 10,906 ± 0,097 |
| Large | 2 | 51,925 ± 0,570 | 5,815 ± 0,103 | 8,953 ± 0,440 |

Sur le carré à deux nœuds, h=32 accélère le témoin d’un facteur 76,30 (IC95 de 74,52 à 78,07). À h=32 fixé, le gain d’échelle de un à deux nœuds vaut 0,871. Sur la paire large à deux nœuds, h=8 accélère le témoin d’un facteur 8,93 (IC95 de 8,77 à 9,09), et le gain d’échelle à h=8 atteint 1,372. Cette paire ne compte que 1 024 lignes, soit 128 paquets de huit ou 32 paquets de 32 pour alimenter 48 bandes. La sonde indépendante enregistre 20 000 et 1 024 appels à Allgather par rang dans la version initiale, puis zéro dans les deux modes pipeline.
