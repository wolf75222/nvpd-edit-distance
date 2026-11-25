# Profil du répartiteur sur EPYC

Sur quatre cœurs EPYC 9654, le calcul de Myers dans Edlib reçoit 85,07 % des échantillons de cycles de l’exécutable `editdist-route`. L’extension des correspondances par blocs AVX-512 dans WFA2 en reçoit 10,97 %. La construction des esquisses en reçoit 0,12 %. Ces proportions décrivent un profil du processus entier, pour une chauffe et six répétitions dans le même processus. [Fonctions LProf](lprof_functions.csv) · [Boucles LProf](lprof_loops.csv)

Le job 840155 a exécuté la politique `hybrid-counts` sous MAQAO. Chaque matrice répartit 610 paires sur Edlib et 1 406 sur WFA2. Le GPU reste inutilisé sur ces 2 016 paires. Le calcul sur EPYC donne la même matrice entière que la référence publiée. [Décisions individuelles](../../../hybrid/matrix-gh200/composition/surveillance-decisions.tsv) · [Empreinte SHA-256](matrix.sha256) · [Allocation Slurm](accounting.tsv) · [Compilation et placement](environment.txt)

Le profil du noyau dense répond à une autre question. Sur le même nœud, le job 840153 consacre 87,60 % des échantillons à la boucle intérieure de NvPD avec quatre fils et une paire de 32 768 bases. Le profil du répartiteur porte sur une matrice dont les paires utilisent les bibliothèques CPU. Les proportions décrivent la répartition du travail propre à chaque entrée.

## Reproduction

Depuis la racine du dépôt sur ROMEO, utiliser le FASTA de la cohorte et sa matrice de référence. Le script compile le binaire en mode Release avec les symboles de débogage, exécute OneView et compare la sortie du processus à la référence avant d’archiver le rapport HTML.

```sh
sbatch --account="$ROMEO_PROJECT" bench/cpu/maqao/route.sbatch \
    build/hybrid-matrix/viral/surveillance.fa results/phylogeny/sars-edit/matrix.tsv
```

Le rapport se trouve dans `build/maqao/route-<job>/profile.tar.gz`. La page `RESULTS/editdist-route_one_html/fcts_and_loops.html` classe les fonctions et les boucles. Dans les exports LProf, `inc_time_p` inclut les boucles internes d’une fonction. Le fichier de catégorisation identifie l’événement échantillonné. [Échantillonnage](lprof_categorization.csv)
