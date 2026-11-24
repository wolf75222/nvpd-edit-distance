# Plafonds entiers du GH200 et tuiles NvPD à halo local

## Calibration indépendante

La copie coalescée lit 128 Mio et écrit 128 Mio par passage. Les noyaux de calcul utilisent la même expression entière générique, `min(1 + min(up, diagonal), match + (j-k))`. Elle contient deux minima, une soustraction et deux additions. Un compilateur peut les fusionner. Le binaire SM90 comporte notamment des instructions `VIMNMX` et `VIADDMNMX`. Les nombres d’opcodes archivés proviennent du désassemblage statique. [Durées CUDA Event](839177/calibration.csv) · [Désassemblage, nombre d’opcodes statiques](839177/register-sass-counts.txt) · [Environnement](839177/environment.txt)

Le plafond de calcul applique 512 fois l’expression à huit chaînes indépendantes de registres par thread. Le balayage lit trois mots de 32 bits et en écrit un, puis exécute 1, 4, 16, 64 ou 256 mises à jour par élément de sortie. Le volume cumulé des quatre tableaux est 512 Mio. Trois positions du résultat sont vérifiées pour chaque intensité. Un processus chauffe chaque noyau et chronomètre sept appels entre deux événements CUDA sur le même flux. Cinq processus indépendants forment les cinq unités statistiques d’une configuration.

| Grandeur | Moyenne ± écart type |
|---|---:|
| Copie, octets demandés | 3 487,67 ± 7,05 Go/s |
| Calcul entier en registres | 5 480,58 ± 13,08 milliards de mises à jour/s |
| Balayage, une mise à jour pour 16 octets | 213,97 ± 0,18 milliards de mises à jour/s |
| Balayage, 256 mises à jour pour 16 octets | 5 067,94 ± 0,86 milliards de mises à jour/s |

Le changement de pente de l’enveloppe empirique se situe à 1,57 mise à jour par octet demandé. Les deux plafonds caractérisent le GH200, le compilateur et ces motifs d’accès. Nsight Compute a demandé les octets DRAM lus et écrits, les instructions exécutées et la durée du noyau. Le pilote 590.48.01 réserve les compteurs aux administrateurs et renvoie `ERR_NVGPUCTRPERM`. Le trafic physique HBM et l’occupation réelle du solveur restent inconnus. Les activités CUPTI accessibles par Nsight Systems donnent les lancements, leurs durées et les copies explicites. [Réponse de NCU sur la calibration](839177/ncu-calibration.txt) · [Documentation NVIDIA](https://developer.nvidia.com/nvidia-development-tools-solutions-err_nvgpuctrperm-permission-issue-performance-counters)

## Variantes et périmètre

Les fenêtres de GRCh38 et CHM13 proviennent de deux régions homologues du chromosome 21. Les premières sont colinéaires. Les secondes incluent une différence de 9 540 bases entre blocs adjacents de la chaîne UCSC. Les longueurs GRCh38 sont 16 384, 65 536, 131 072 et 524 288 bases. Chaque FASTA est vérifié par SHA-256. La référence de score est commune à edlib et WFA2. [Coordonnées, empreintes et distances](../../kokkos/human/references.csv)

- `native` appelle le NvPD CUDA du projet. Le témoin passe à une ligne par lancement lorsque la portée globale dépasse 64, comme sur les huit entrées humaines.
- `fused` produit deux lignes par lancement. La seconde ligne recalcule jusqu’à trois valeurs intermédiaires en registres. `shuffle` partage les valeurs intermédiaires des 32 colonnes d’un warp, avec recalcul si le prédécesseur se trouve hors du warp. Ces deux variantes sont disponibles dans le banc. [Expérience de fusion](839181/trials.csv) · [Profil d’activité](839181/activity.csv)
- `local`, accessible sous `--algorithm nvpd-cuda-local`, conserve les tuiles du témoin sur les entrées à faible portée. Pour une grande portée globale, le plan calcule huit bords gauches par bloc de sortie. Au bord gauche `L > 0`, les colonnes nécessaires au rang précédent sont au moins `min(L-1, min(k_s(L)-1 pour les k_s(L)>0))`. Cette borne vaut zéro à l’origine. L’indice de dernière occurrence ne peut que croître après le premier symbole. Sept applications de la borne suffisent pour huit lignes, y compris les blocs tronqués au bord droit. Le dernier paquet utilise le suffixe du plan correspondant au nombre de lignes restantes. Si le double tampon ne tient pas en mémoire partagée pour tous les blocs, l’appel prend le chemin natif par ligne. Le mode trace conserve le tableau et le CIGAR habituels. [Code](../../../src/cuda/local_tiles.cuh) · [Comparaison sur la même allocation](839191/trials.csv) · [Profil d’activité](839191/activity.csv)

Chaque bloc répartit les méthodes dans un ordre pseudo-aléatoire. Chaque processus lit le FASTA, effectue un appel complet d’échauffement, puis chronomètre un unique appel au score. Le temps inclut l’index hôte, le plan éventuel, les allocations et les copies et synchronisations nécessaires au résultat. Les deux scores sont validés. Les temps des bibliothèques correspondent à un thread Grace dans la tâche liée aux quatre cœurs réservés. Le job 839181 mesure les deux fusions face au témoin, puis le job 839191 remesure le témoin, la fusion simple et la tuile locale sur son allocation. Le sondage Nsight utilise un processus séparé par entrée et méthode. [Moyennes, écarts types et IC95](summary.csv) · [Ratios appariés par bloc](paired.csv) · [Tests du mode public](839191/cli.txt)

Les GCUPS logiques divisent `m × n` par le temps complet de calcul. Pour le témoin natif, chaque cellule intérieure charge l’index et écrit un score. Suivant la branche de la récurrence, elle référence un, deux ou trois scores antérieurs, ce qui donne 12, 16 ou 20 octets logiquement demandés. Ces comptes par adresse proviennent des séquences. Les caches, la coalescence et les relectures changent les transactions produites. Le plan publié compte également les cellules calculées dans le halo et les colonnes de bord. `device_bytes` compte les tableaux explicitement réservés sur le GPU, dont les 32 octets de plan par bloc dans la variante locale. Le pic RSS est celui du processus Grace après les deux appels. Certains petits processus edlib et WFA2 reçoivent un RSS noyau nul sur ce système. Ces zéros figurent dans le brut et les résumés RSS requièrent cinq observations strictement positives. [Modèle par entrée](traffic.csv) · [Plafonds et intervalles](calibration-summary.csv)

## Halos locaux et regroupement par largeur

La variante expérimentale `classes` divise les blocs de sortie au 90ᵉ percentile de largeur. Tous les blocs à égalité au seuil restent dans le premier groupe. Chaque groupe non vide lance son propre noyau pour le même paquet de huit lignes. La mémoire partagée par bloc suit la largeur maximale de son groupe. Les deux variantes calculent les mêmes cellules et gardent le même plan de huit entiers par bloc. [Code de regroupement](../../../bench/gpu/roofline/width_classes.cuh) · [Durées individuelles](839198/trials.csv)

Chaque colonne du tableau donne une moyenne de cinq rapports de temps appariés par bloc. Les deux colonnes utilisent des allocations distinctes. Une valeur supérieure à un favorise la méthode au dénominateur.

| Région de chr21 | Bases GRCh38 | Par ligne / halo local (839191) | Halo local / deux classes (839198) |
|---|---:|---:|---:|
| Colinéaire | 16 384 | 2,693 | 0,510 |
| Colinéaire | 65 536 | 2,645 | 0,532 |
| Colinéaire | 131 072 | 2,349 | 0,572 |
| Colinéaire | 524 288 | 1,637 | 0,586 |
| Insertion | 16 384 | 2,747 | 0,506 |
| Insertion | 65 536 | 2,431 | 0,540 |
| Insertion | 131 072 | 2,160 | 0,617 |
| Insertion | 524 288 | 1,161 | 0,703 |

Sur la plus grande insertion, le regroupement ramène le double tampon de 22 624 à 8 960 octets pour 474 des 525 blocs. Les 51 blocs restants conservent 22 624 octets. Les deux variantes calculent 292 824 350 720 cellules, bords et halos compris. Les lancements passent de 65 537 à 131 073 et leur durée cumulée de 1,686 à 2,374 s dans les traces d’activité séparées. Le temps complet moyen passe de 1,751 ± 0,020 à 2,491 ± 0,022 s. Le rapport apparié vaut 0,703 avec un IC de Student à 95 % de 0,697 à 0,709. Le mode public `nvpd-cuda-local` conserve un lancement par paquet. [Activités Nsight Systems](839198/activity.csv) · [Contrôles des scores](839198/validation.txt) · [Tailles des deux groupes](traffic.csv)

## Rejouer sur ROMEO

Depuis une copie Git du projet contenant `b5d08d1e5bb6d6b8e744bc6680b404caf569e991`, préparer les séquences et soumettre dans l’ordre ci-dessous. Les lanceurs construisent localement au nœud avec CUDA 13.0.2 et SM90. L’archive de sources regroupe les fichiers soumis par numéro de job. Pour une reproduction historique, appliquer dans l’ordre les contenus des dossiers 839181, 839191 puis 839198 à la révision de base, en s’arrêtant au job voulu. Le calibrage 839177 est indépendant. La première soumission du mode local (839190) a échoué avant les mesures, car le répertoire des tests manquait dans la copie distante. Le job 839191 inclut ces fichiers de la révision de base. [Sources soumises](submitted.tar.gz) · [Allocations Slurm](accounting.psv)

```sh
python3 bench/kokkos/human/inputs.py build/roof/cache build/roof/inputs
base=bench/gpu/roofline
out=results/roofline/gpu
cal=$(sbatch --parsable -A "$ROMEO_PROJECT" "$base/calibrate.sbatch")
fusion=$(sbatch --parsable -A "$ROMEO_PROJECT" --dependency=afterok:$cal "$base/solve.sbatch" build/roof/inputs)
local_job=$(sbatch --parsable -A "$ROMEO_PROJECT" --dependency=afterok:$fusion "$base/solve.sbatch" build/roof/inputs local)
classes_job=$(sbatch --parsable -A "$ROMEO_PROJECT" --dependency=afterok:$local_job "$base/solve.sbatch" build/roof/inputs classes)
# Après les quatre jobs
python3 "$base/collect.py" "$out" "$cal" "$fusion" "$local_job" "$classes_job"
python3 "$base/traffic.py" build/roof/inputs results/kokkos/human/references.csv "$out/traffic.csv" --trials "$out/$local_job/trials.csv"
```
