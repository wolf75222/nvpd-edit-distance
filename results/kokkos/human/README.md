# NvPD sur deux régions du chromosome 21 humain

Sur le GH200, fusionner deux lignes NvPD ramène le temps moyen de 2,520 à 1,345 s pour la fenêtre colinéaire de 524 288 bases de GRCh38. À huit cœurs EPYC, garder l’équipe Kokkos active pendant le calcul ramène le temps de 1,080 à 0,652 s pour 65 536 bases. OpenMP natif calcule cette dernière distance en 0,284 s. Chaque temps porte sur cinq processus indépendants. [Mesures individuelles](trials.csv) · [Moyennes, écarts types et IC à 95 %](summary.csv)

## Séquences et distances

Les fenêtres comparent GRCh38.p14 et T2T-CHM13v2.0, brin direct. Leurs deux bornes proviennent de blocs alignés sur la chaîne UCSC 23. Les coordonnées du manifeste commencent à zéro, borne de fin exclue. Chaque fenêtre contient uniquement A, C, G et T. Les dix FASTA sont identifiés par SHA-256. [GRCh38, NC_000021.9](https://www.ncbi.nlm.nih.gov/nuccore/NC_000021.9) · [CHM13, NC_060945.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_060945.1) · [Chaîne UCSC hg38–hs1](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHs1.over.chain.gz) · [Coordonnées et empreintes](inputs.csv)

| Centre GRCh38 | Caractéristique | Distance à 16 384 bases GRCh38 | Distance à 524 288 bases GRCh38 |
|---|---|---:|---:|
| 37 400 000 | Région colinéaire | 14 / 16 394 (0,085 %) | 979 / 524 431 (0,187 %) |
| 33 474 775 | Segment supplémentaire de 9 540 bases dans CHM13 entre deux blocs consécutifs | 9 568 / 25 940 (36,89 %) | 13 706 / 537 308 (2,55 %) |

Le dénominateur est la longueur du plus long intervalle. Les cinq fenêtres de chaque locus sont emboîtées. Sur chaque architecture, edlib et WFA2 ont retrouvé indépendamment les dix scores. [Références par entrée](references.csv)

## Politiques Kokkos et calcul mesuré

- `native` appelle NvPD OpenMP sur EPYC ou NvPD CUDA sur GH200. La portée maximale des dépendances vaut entre 86 et 1 813 colonnes sur ces entrées. Le seuil des tuiles CUDA natives est 64, ce qui sélectionne le noyau natif ligne par ligne pour ces dix fenêtres. [Portée par entrée](reach.csv)
- `kokkos-rows` initialise une première ligne puis soumet un `RangePolicy` pour chaque ligne suivante.
- `kokkos-opt` est la politique utilisée par `editdist pair --algorithm nvpd-kokkos`. Sur hôte, une équipe Kokkos traite toutes les lignes. Chaque travailleur conserve un bloc contigu de colonnes et l’équipe se synchronise en fin de ligne. Si une seule équipe ne peut pas réunir les travailleurs demandés, la politique ligne par ligne prend le relais. Sur accélérateur, pour le score seul, un noyau calcule deux lignes logiques. Chaque cellule de la seconde reconstruit à la demande jusqu’à trois valeurs intermédiaires depuis la dernière ligne terminée. Le calcul accepte les dépendances lointaines et traite la dernière ligne impaire. Le chemin complet sur accélérateur conserve le tableau et les soumissions ligne par ligne. [Implémentation](../../../src/kokkos/kokkos.cpp)

L’EPYC 9654 couvre les longueurs GRCh38 16 384 et 65 536 à 1, 4 et 8 cœurs par locus. Trois tâches Slurm tournent simultanément sur des ensembles de cœurs disjoints. Le GH200 couvre aussi 131 072, 262 144 et 524 288 bases avec un accélérateur et un cœur pour piloter le calcul. Le binaire est commun aux trois méthodes d’une même architecture. GCC 11.5, mode Release, Kokkos 4.5.01 et CUDA 12.6 sur GH200. [Compilation EPYC](environment-epyc.txt) · [Compilation GH200](environment-gh200.txt) · [Allocations Slurm](accounting.tsv)

Chaque bloc alterne l’ordre des trois méthodes. Chaque méthode démarre un nouveau processus, lit son FASTA, chauffe avec un calcul complet, puis chronomètre un seul appel au score par horloge monotone. Le chronomètre inclut la construction de l’index, les allocations, les copies et la récupération du résultat. La lecture du FASTA et l’initialisation du runtime précèdent l’échauffement. Les deux scores du processus sont contrôlés. Les 330 processus forment 66 configurations à cinq répétitions chacune. L’écart type d’échantillon et l’IC de Student à quatre degrés de liberté décrivent la dispersion entre processus sur l’allocation utilisée. Les gains entre méthodes sont les moyennes des cinq rapports appariés par bloc. [Rapports appariés et IC](paired.csv) · [Tests EPYC](validation-epyc.txt) · [Tests GH200](validation-gh200.txt)

## Temps, débit et mémoire

| Cas | Natif (s) | Kokkos ligne (s) | Kokkos optimisé (s) | IC 95 % du temps optimisé (s) | Débit optimisé (GCUPS) |
|---|---:|---:|---:|---:|---:|
| EPYC 8 cœurs, colinéaire 65 536 | 0,2841 ± 0,0012 | 1,0801 ± 0,0050 | 0,6520 ± 0,0004 | [0,6515 ; 0,6525] | 6,59 |
| GH200, colinéaire 524 288 | 1,9900 ± 0,0346 | 2,5196 ± 0,0339 | 1,3450 ± 0,0184 | [1,3222 ; 1,3679] | 204,45 |
| GH200, insertion 524 288 | 1,9912 ± 0,0264 | 2,5051 ± 0,0405 | 1,3760 ± 0,0144 | [1,3581 ; 1,3938] | 204,75 |

Les temps sont donnés en moyenne ± écart type. Le débit moyen compte les cellules logiques `m × n` par seconde. Le noyau fusionné peut recalculer des valeurs intermédiaires. Sur la fenêtre colinéaire de 65 536 bases, le rapport des temps moyens à un et huit cœurs vaut 7,30 pour l’équipe persistante, soit 91,2 % d’efficacité. OpenMP natif atteint 6,38, soit 79,7 %. Ces rapports descriptifs comparent des ensembles de cœurs différents. [Échelle forte sur EPYC](scaling.csv)

Le RSS est le maximum du processus recueilli par `getrusage(RUSAGE_SELF)` après l’échauffement et le calcul mesuré. Le noyau retourne 184 320 Kio pour tous les processus GH200 de cette série. La colonne `device_bytes` donne séparément les octets demandés pour l’index de correspondance et les deux lignes. À 524 288 bases, les trois méthodes utilisent 12 586 368 octets sur le GPU pour le cas colinéaire et 12 895 416 pour le cas avec insertion. Le runtime, le contexte et les espaces internes du matériel sont hors de ce décompte. `host_workspace_array_bytes` est la taille calculée des tableaux principaux sur l’hôte. Pour quatre symboles et une largeur `W = max(m,n) + 1`, elle vaut `24W` octets en OpenMP natif, `40W` en Kokkos hôte (index et copie Kokkos), et `16W` pour l’index préparé sur Grace. [Statistiques mémoire par configuration](summary.csv)

## Vérification des lancements GPU

Nsight Systems 2024.5 a suivi séparément un appel chauffé de chaque méthode pour trois entrées. À 16 384 bases GRCh38, chaque locus passe de 16 385 lancements avec Kokkos ligne à 8 193 avec la fusion, initialisation comprise. À 65 536 bases colinéaires, le décompte passe de 65 537 à 32 769 et le temps cumulé dans les noyaux de 95,83 à 64,06 ms. Les deux variantes transfèrent 1 049 300 octets en quatre copies sur cette entrée. Les temps de la table précédente proviennent d’appels sans profileur. [Neuf profils mesurés](profile.csv) · [Décompte par noyau CUDA](kernel_groups.csv) · [Environnement](environment-profile.txt) · [Validation](validation-profile.txt)

## Reproduction

Depuis la racine du projet, préparer les entrées et un clone Git de Kokkos accessible depuis les nœuds. Le lanceur extrait la révision `175257a51ff29a0059ec48bcd233ee096b2c0438` (4.5.01). Les partitions et empreintes Spack du lanceur correspondent à ROMEO.

```sh
export KOKKOS_INPUTS="$PWD/build/human-kokkos/inputs"
python3 bench/kokkos/human/inputs.py build/human-kokkos/cache "$KOKKOS_INPUTS"
git clone --branch 4.5.01 --depth 1 https://github.com/kokkos/kokkos.git build/human-kokkos/kokkos
export KOKKOS_SOURCE="$PWD/build/human-kokkos/kokkos"
cpu_job=$(sbatch --parsable -A "$ROMEO_PROJECT" -p instant -C x64cpu -n3 -c8 bench/kokkos/human/run.sbatch epyc)
gpu_job=$(sbatch --parsable -A "$ROMEO_PROJECT" -p instant -C armgpu -n1 -c4 --gpus-per-node=1 bench/kokkos/human/run.sbatch gh200)
sbatch -A "$ROMEO_PROJECT" -p instant -C armgpu -n1 -c4 --gpus-per-node=1 --export=ALL,KOKKOS_PROFILE=1 bench/kokkos/human/run.sbatch gh200
# Après les deux jobs de mesure
python3 bench/kokkos/human/collect.py "results/kokkos/human/jobs/$cpu_job" "results/kokkos/human/jobs/$gpu_job" build/human-kokkos/summary
```

Les jobs 838913 et 838914 ont mesuré les temps. Le job 838917 a enregistré les compteurs de lancement. La base source est `a728923aef087ce35020d6f1568af496074f0c43`. L’archive contient les seuls fichiers à superposer à cette révision pour retrouver les octets soumis. Extraire `timing/` pour les temps, puis superposer `profile/` pour Nsight. Le générateur publié prépare aussi le cache vide et écrit la portée dans `reach.csv`. Les sources natives, les tests et le CMake principal coïncident avec ceux de la révision de base. [Sources soumises](submitted.tar.gz) · [Lanceur](../../../bench/kokkos/human/run.sbatch)
