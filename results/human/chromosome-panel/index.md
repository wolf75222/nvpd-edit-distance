# Préallocation des séquences à partir d’un index FAI

Le lecteur SIMD accepte un index `.fai` optionnel. Il réserve une fois la capacité de chaque séquence avant de décoder ses bases. La commande suivante s’applique à un FASTA à deux enregistrements. Le mode `matrix` prend le même argument.

```sh
samtools faidx paire.fa
editdist pair --algorithm edlib --input paire.fa --fai paire.fa.fai
```

L’index décrit un FASTA non compressé, avec un repliage régulier par séquence. Les cinq champs sont le nom, la longueur, l’offset de la première base, les bases par ligne et les octets par ligne. Le lecteur attend les entrées dans l’ordre des enregistrements et vérifie ces champs pendant le décodage. Une longueur hors fichier, un index absent explicitement demandé ou une discordance de mise en lignes déclenche une erreur. Sans `--fai`, le programme conserve sa lecture habituelle des fichiers à repliage libre. [Format HTSlib](https://www.htslib.org/doc/faidx.html) · [Lecteur](../../../src/io/fasta.cpp)

## Protocole

Les six FASTA contiennent chacun deux versions complètes d’un chromosome humain (CHM13v1.0 et CHM13v1.1). L’indexation a lieu avant les mesures. Samtools 1.21, via pysam 0.23.3, a reproduit indépendamment les octets des six index. Les SHA-256 du FASTA brut, du FASTA normalisé et de l’index figurent dans les observations. Le même binaire lit chaque paire avec et sans index, et les deux sorties décodées correspondent au SHA-256 calculé par un troisième décodeur. Un FAI encode la géométrie du fichier. Une substitution conservant cette géométrie laisse l’index valide. [Manifestes des chromosomes](manifest.tsv) · [Index et offsets](index-entries.tsv)

Le job 838927 a utilisé deux cœurs physiques disjoints d’un nœud EPYC 9654, un processus à la fois par cœur. Chaque chromosome reste sur le même cœur. Chacun des neuf blocs lance les deux modes dans un ordre pseudo-aléatoire. Un processus frais lit entièrement le fichier deux fois. La seconde lecture fournit le temps et le maximum RSS après les deux appels. Le chronomètre inclut l’ouverture et la validation de l’index, la normalisation des bases et la construction des chaînes en mémoire. Les séquences et l’index sont servis par le cache de pages après l’échauffement. [Placement, compilateur et sources](index-environment.txt) · [Durée et état Slurm](index-accounting.psv)

Un deuxième processus frais répète la chauffe puis le calcul instrumenté pour chaque observation. Les compteurs PMU couvrent les instructions utilisateur, les cycles, les chargements L1D et leurs échecs, avec le pourcentage de temps actif dans les lignes brutes. La surcharge de `operator new/new[]` compte les appels ordinaires et additionne les octets demandés pendant la seule lecture instrumentée. Le cumul inclut les capacités intermédiaires libérées. Le maximum RSS couvre le processus entier et s’exprime en Kio. [108 essais et compteurs](index-trials.csv)

Les synthèses donnent la moyenne arithmétique, l’écart type d’échantillon et l’IC bilatéral de Student à huit degrés de liberté. Pour chaque chromosome et chaque métrique, le rapport est calculé dans chacun des neuf blocs avant d’être moyenné. Les intervalles décrivent la dispersion entre démarrages au sein de cette allocation. [Statistiques par mode](index-summary.csv) · [Rapports blocs/indexé](index-paired.csv)

## Résultats

| Chromosome | Blocs, s (moy. ± ET) | Indexé, s (moy. ± ET) | Gain apparié [IC 95 %] | Appels à new, blocs → indexé | Pic RSS moyen, Mio, blocs → indexé |
|---|---:|---:|---:|---:|---:|
| 1 | 0,9445 ± 0,0056 | 0,6568 ± 0,0051 | 1,438 [1,427 à 1,449] | 53 → 15 | 480,38 → 477,00 |
| 4 | 0,8012 ± 0,0078 | 0,5152 ± 0,0062 | 1,555 [1,539 à 1,572] | 53 → 15 | 427,98 → 372,16 |
| 7 | 0,7084 ± 0,0063 | 0,4133 ± 0,0029 | 1,714 [1,700 à 1,729] | 53 → 15 | 396,50 → 308,96 |
| 12 | 0,6459 ± 0,0095 | 0,3552 ± 0,0018 | 1,819 [1,798 à 1,839] | 54 → 16 | 370,52 → 256,50 |
| 16 | 0,3887 ± 0,0021 | 0,2466 ± 0,0014 | 1,577 [1,568 à 1,586] | 52 → 16 | 214,77 → 186,37 |
| 20 | 0,3236 ± 0,0062 | 0,1763 ± 0,0014 | 1,836 [1,807 à 1,865] | 52 → 16 | 185,89 → 129,00 |

La réservation supprime les agrandissements successifs des deux chaînes. Pour le chromosome 1, le cumul des tailles demandées tombe de 1 006,64 à 496,79 Mo, avec une baisse de 3,38 Mio du pic résident. Sur le chromosome 12, le pic baisse de 114,02 Mio. Le nombre de cycles utilisateur diminue d’un facteur de 1,33 à 1,69 et les échecs des chargements L1D d’un facteur de 1,67 à 2,40, selon le chromosome. Les contrôles de géométrie ajoutent des instructions à chaque ligne. Le nombre d’instructions utilisateur augmente de 5,6 à 7,2 % dans ces mesures. [Temps, allocations et compteurs PMU](index-summary.csv)

## Reproduire

À la racine du projet, `PANEL` désigne le répertoire des FASTA du [générateur de chromosomes](../../../scripts/datasets/chm13_panel.py). Le lanceur construit les index dans le répertoire de sortie.

```sh
job=$(sbatch --parsable -A "$ROMEO_PROJECT" bench/human/chromosome-panel/blocks/index.sbatch "$PANEL")
# Après la fin du job
raw=results/human/chromosome-panel/index-$job
python3 bench/human/chromosome-panel/blocks/index_summary.py "$raw" build/fai-summary
```

Pour vérifier une seconde génération, passer `--htslib-indexes CHEMIN` au collecteur. Ce répertoire doit contenir les six fichiers `chr*.fa.fai` créés par `samtools faidx` sur les mêmes FASTA. Chaque ligne de `index-entries.tsv` fournit les cinq champs de l’index précédés du chromosome.
