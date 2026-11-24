# Chromosomes humains entiers

Les assemblages T2T-CHM13 v1.0 et v1.1 proviennent du même échantillon. La v1.0 comporte cinq lacunes dans les réseaux d’ADN ribosomique. La comparaison porte sur les 17 autres autosomes et sur X, de la première à la dernière base. La somme des longueurs de référence vaut 2 643 920 725 nt. Chaque paire conserve son orientation et ses extrémités. Le générateur vérifie les MD5 des deux archives du consortium puis calcule le SHA-256 de chaque FASTA. [Versions du consortium](https://github.com/marbl/CHM13/blob/master/Earlier_assembly_releases_and_associated_data.md) · [Manifeste des 18 paires](manifest.tsv)

## Score et durée

Le job ROMEO 838916 a réparti les chromosomes entre huit cœurs physiques de l’EPYC 9654. Chaque alignement utilise un cœur. Dans chacun des cinq blocs, un nouveau processus calcule deux fois le score sur les mêmes séquences en mémoire. Seul le second appel est chronométré. Le pic RSS du processus comprend aussi la lecture et l’échauffement. L’ordre des bibliothèques varie entre blocs. Le lot comporte 180 temps d’alignement. Les intervalles à 95 % utilisent la loi de Student sur les cinq observations d’une configuration. Un rapport apparié se calcule avant la moyenne entre blocs. [Chronomètres](trials.csv) · [Moyennes, écarts-types et intervalles](summary.csv) · [Rapports appariés](paired.csv)

Edlib et WFA2 donnent le même score sur les 18 paires. Le chromosome 18 demande 5 éditions et le chromosome 3 en demande 659. Le cumul atteint 2 587 éditions, soit 0,978 édition par million de bases de la v1.0. Il décrit les différences entre ces deux assemblages sur les chromosomes retenus. [Résultats par chromosome](trials.csv)

## Position des modifications de référence

Le VCF officiel donne les remplacements entre versions, en coordonnées de la v1.0. Sur les 18 chromosomes, ses 310 entrées ont passé le filtre du consortium. Le pilote vérifie chaque allèle de référence et applique le premier allèle du génotype déposé. Sur 292 entrées, tout l’intervalle de référence tient dans les mille premières ou dernières bases du chromosome. Le remplacement reproduit la v1.1 sur 14 chromosomes. Sur les chromosomes 1, 3, 4 et X, seule la première base du résultat diffère de la v1.1. BCFtools 1.21, utilisé avec le premier allèle de chaque génotype, reproduit les 18 empreintes du rejeu.

Le fichier [release-intervals.csv](release-intervals.csv) conserve, pour chaque entrée du VCF, le chromosome, les bornes de référence semi-ouvertes à partir de zéro et le coût du remplacement. Depuis la racine du dépôt, la commande suivante le reconstruit à partir du VCF compressé du consortium.

```sh
python python/studies/chromosome_landscape.py /chemin/vers/patch.vcf.gz results/human/chromosome-panel/release-intervals.csv
```

Dans ce VCF, l’un des deux allèles de chaque remplacement a une seule base. Pour un allèle long de longueur $k$, le coût local vaut $k-1$ si la base isolée y apparaît et $k$ sinon. Les coûts locaux totalisent 2 663. Quatre substitutions supplémentaires donnent une transformation réalisable de coût 2 667 entre les séquences publiées. L’optimisation globale économise 80 opérations par rapport à ce chemin. [Fichier du consortium](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.0_patch.vcf.gz) · [Rejeu par chromosome et empreintes](release-patch.csv)

## Lecture et mémoire

Le job 838916 a comparé trois lecteurs sur les chromosomes 1, 4, 7, 12, 16 et 20. Le premier valide chaque caractère, le deuxième copie par ligne, le troisième ajoute une directive SIMD au deuxième. Les 90 temps excluent les compteurs matériels. Un processus séparé relève les instructions, les cycles et les lectures L1 de données, dont celles qui manquent le cache. Les empreintes portent sur tous les identifiants et toutes les bases décodées. [Lectures](reader.csv) · [Compteurs](counters.csv) · [Empreintes](checks.csv) · [Environnement](environment.txt)

Le job 838919 compare le lecteur par ligne de la première expérience à un test explicite par blocs de 16 octets. Six processus occupent chacun un cœur et un chromosome. Le pilote alterne les variantes dans cinq blocs et distingue la lecture chronométrée du relevé de compteurs. Il compte les appels C++ à `new` et `new[]` et les octets cumulés demandés pendant le décodage. Le RSS est le pic résident du processus en Kio. Les trois quantités mesurent des aspects différents de la mémoire. [Mesures appariées](blocks.csv) · [Configuration de compilation](blocks-environment.txt)

Le job 838927 isole la préallocation avec l’index standard FAI. Les deux modes utilisent le même binaire, neuf blocs par chromosome, un cœur fixé et des processus distincts pour le chronométrage et les compteurs. Le lecteur vérifie les positions et les longueurs déduites du fichier. [Protocole et commandes de reproduction](index.md) · [Mesures sur les six chromosomes](index-trials.csv)

## Reproduction

Depuis la racine du dépôt, reconstruire les FASTA puis vérifier le manifeste.

```sh
python3 scripts/datasets/chm13_panel.py --cache build/t2t-cache --out data/human/chm13-panel
cmp data/human/chm13-panel/manifest.tsv results/human/chromosome-panel/manifest.tsv
curl -fL https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.0_patch.vcf.gz -o build/t2t-cache/patch.vcf.gz
python3 bench/human/chromosome-panel/release_patch.py data/human/chm13-panel build/t2t-cache/patch.vcf.gz build/t2t-cache/release-patch.csv
cmp build/t2t-cache/release-patch.csv results/human/chromosome-panel/release-patch.csv
```

Pour effectuer aussi le contrôle avec BCFtools, indexer une copie BGZF du VCF, produire la séquence consensus, puis passer ce FASTA au pilote. Les exécutables `bgzip` et `bcftools` sont disponibles avec HTSlib et BCFtools.

```sh
gzip -dc build/t2t-cache/patch.vcf.gz | bgzip -c > build/t2t-cache/edits.vcf.gz
bcftools index build/t2t-cache/edits.vcf.gz
bcftools consensus -H 1 -f build/t2t-cache/chm13.v1.0.fa.gz -o build/t2t-cache/consensus.fa build/t2t-cache/edits.vcf.gz
python3 bench/human/chromosome-panel/release_patch.py data/human/chm13-panel build/t2t-cache/patch.vcf.gz build/t2t-cache/release-patch.csv --consensus build/t2t-cache/consensus.fa
```

Sur ROMEO, mettre le chemin absolu des entrées dans `panel`. Les deux bancs historiques compilent le lecteur contenu dans leur copie de travail. Les révisions ci-dessous restituent celui effectivement mesuré. CMake utilise les archives vérifiées de WFA2 et d’edlib dans `deps/`, ou les télécharge à la configuration.

```sh
panel=$(realpath data/human/chm13-panel)
git worktree add --detach build/panel-lines 421927e9
git worktree add --detach build/panel-blocks 345c06e
(cd build/panel-lines && sbatch -A "$ROMEO_PROJECT" bench/human/chromosome-panel/run.sbatch "$panel")
(cd build/panel-blocks && sbatch -A "$ROMEO_PROJECT" bench/human/chromosome-panel/blocks/run.sbatch "$panel")
```

Chaque soumission crée un dossier numéroté sous `results/human/chromosome-panel/` de sa copie de travail. Les commandes suivantes agrègent les fichiers produits par chaque tâche. [Pilote de vérification et d’agrégation](../../../bench/human/chromosome-panel/collect.py)

```sh
python3 bench/human/chromosome-panel/collect.py build/panel-lines/results/human/chromosome-panel/"$JOB_LINES" build/reproduction-lines
python3 bench/human/chromosome-panel/collect.py --reader-blocks build/panel-blocks/results/human/chromosome-panel/"$JOB_BLOCKS" build/reproduction-blocks
```
