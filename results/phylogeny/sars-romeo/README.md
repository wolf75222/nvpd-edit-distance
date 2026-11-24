# Calcul de la matrice de surveillance sur ROMEO

Les 64 génomes complets du manifeste NCBI donnent 2 016 paires. Chaque processus lit le FASTA, calcule les distances globales avec les symboles IUPAC littéraux et écrit les 64 lignes du TSV. Le chronomètre mural englobe ces opérations et le lancement du processus. Le numéro de bloc zéro désigne la chauffe. Chacun des cinq blocs suivants permute l’ordre des configurations avec une graine fixe. Toutes les sorties sont comparées octet par octet à la matrice publiée. [Accès aux données](../../../scripts/datasets/ncbi_sars.py) · [Matrice de référence](../sars-edit/matrix.tsv) · [Chronométrage](../../../bench/campaign/surveillance/measure.py)

L’allocation 839288 compare quatre méthodes sur un et 32 cœurs EPYC et chronomètre le noyau NvPD dense sur 32 cœurs. L’allocation 839289 compare NvPD à halo local sur un GH200 avec Edlib et WFA2 sur un cœur Grace. Slurm déclare un thread par cœur sur ces deux nœuds. Le champ `threads` compte les threads hôtes. Pour NvPD sur GPU, un seul thread hôte soumet les noyaux CUDA. Les temps portent sur une cohorte fixe et les intervalles de Student sur cinq processus, avec quatre degrés de liberté. [EPYC](839288/trials.csv) · [Grace et Hopper](839289/trials.csv) · [Topologie Slurm](topology.txt)

Ces deux mesures utilisent le CLI du commit `d9f394c`. Le [diff du lanceur](initial-runner.patch) restitue les deux scripts soumis à partir de la version publiée ici. Sur EPYC, ce CLI attribuait des rangées entières aux bibliothèques. L'essai apparié 839293 compare ce parcours avec la répartition des paires en groupes, intégrée depuis au CLI. Les résultats de ces campagnes gardent leur propre identifiant d’allocation. [Effet de la répartition](../sars-schedule/README.md) · [Coût des coupes de l’hybride](../sars-ablation/)

Depuis la racine du dépôt, les commandes suivantes reconstruisent la cohorte et lancent la version présente dans le répertoire de travail. Le script de soumission compile le binaire dans l’allocation. Les archives verrouillées d’Edlib et de WFA2 peuvent être placées dans `deps/` avant le transfert sur ROMEO.

```sh
python3 scripts/datasets/ncbi_sars.py
sbatch -A "$ROMEO_PROJECT" --constraint=x64cpu -c 32 bench/campaign/surveillance/run.sbatch \
    cpu data/sars/cohort.fa results/phylogeny/sars-edit/matrix.tsv
sbatch -A "$ROMEO_PROJECT" --constraint=armgpu --gpus-per-node=1 -c 4 bench/campaign/surveillance/run.sbatch \
    gpu data/sars/cohort.fa results/phylogeny/sars-edit/matrix.tsv
```

L’allocation 839297, antérieure au changement des écritures de traceback de WFA2, reprend Edlib, WFA2, BiWFA et l’hybride après subdivision des rangées et ajout du filtre de composition. Les huit configurations couvrent un et 32 cœurs. La figure de durée conserve les deux branches de l’essai apparié 839293 pour illustrer le changement de répartition des bibliothèques à 32 cœurs. [Premier essai des huit configurations](839297/trials.csv)

```sh
sbatch -A "$ROMEO_PROJECT" --constraint=x64cpu -c 32 -t 00:10:00 \
    bench/campaign/surveillance/run.sbatch release data/sars/cohort.fa results/phylogeny/sars-edit/matrix.tsv
```

L’allocation 839304 répète les huit configurations avec les écritures bornées du CIGAR de WFA2. Dans le CLI livré, Edlib calcule la matrice en 0,3420 ± 0,0001 s sur 32 cœurs EPYC. Les temps de BiWFA, de l’hybride NvPD/BiWFA et de WFA2 valent respectivement 0,4624 ± 0,0045, 0,4668 ± 0,0052, 0,7740 ± 0,0016 s. L’accélération appariée de l’hybride entre un et 32 cœurs vaut 29,71 (IC à 95 % de 29,31 à 30,12), soit une efficacité de 92,9 %. Sur GH200, NvPD dense à halo local prend 104,07 ± 1,34 s pour la même matrice. [Essais du CLI livré](839304/trials.csv) · [Temps et passage à l’échelle](scaling.csv) · [Calcul des rapports appariés](../../../bench/campaign/surveillance/scaling.py)
