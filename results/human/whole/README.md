# Assemblages chromosomiques entiers

Le chromosome 20 provient de deux éditions du même assemblage CHM13. Le chromosome 21 compare GRCh38.p14 et T2T-CHM13v2.0. Le [manifeste](manifest.tsv) donne les accessions versionnées et les empreintes des paires. [Versions de CHM13](https://github.com/marbl/CHM13/blob/master/Earlier_assembly_releases_and_associated_data.md)

Le job 838443 a réservé un cœur EPYC. Chaque processus dispose de 120 s et de 4 Gio d'espace d'adressage. Les [mesures](838443/trials.csv) distinguent le temps de l'appel et celui du processus. Une interruption laisse le score vide. Les [bornes du chromosome 21](bounds.csv) utilisent les effectifs par lettre et le chemin global qui suit les blocs de la chaîne UCSC hg38ToHs1 numéro 23.

```sh
python3 scripts/human_reference.py --whole-assemblies --cache build/refs --manifest results/human/whole/manifest.tsv
python3 bench/chromosome_bounds.py data/human/chr21-assemblies.fa build/refs/hg38ToHs1.over.chain.gz build/chr21-bounds.csv
sbatch --account="$ROMEO_PROJECT" bench/whole_chromosomes.sbatch
```
