# Assemblages chromosomiques entiers

Le chromosome 20 provient des assemblages CHM13v1.0 et v2.0 du même échantillon. La version intermédiaire v1.1 porte déjà l’accession CP068258.2. [CHM13v1.0](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.2_T2T-CHM13v1.0/GCA_009914755.2_T2T-CHM13v1.0_assembly_report.txt) · [CHM13v1.1](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_T2T-CHM13v1.1/GCA_009914755.3_T2T-CHM13v1.1_assembly_report.txt) · [CHM13v2.0](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_assembly_report.txt)

Le chromosome 21 compare GRCh38.p14 et T2T-CHM13v2.0. Le [manifeste](manifest.tsv) donne les accessions versionnées et les empreintes des paires. [Versions de CHM13](https://github.com/marbl/CHM13/blob/master/Earlier_assembly_releases_and_associated_data.md)

Sur le chromosome 20, un chemin optimal de CP068258.1 vers CP068258.2 contient 18 substitutions, 10 insertions et deux délétions de bases. Dix-neuf éditions précèdent un bloc identique de 66 209 720 bases et onze le suivent. Les 27 segments d'édition sont repérés sur les deux références par des intervalles à base zéro, borne droite exclue. Le job ROMEO 838781 a rejoué ce chemin sur la paire du manifeste. Les accès désignent deux versions du même assemblage. [Segments d’édition](chr20-edits.tsv) · [Polissage des références T2T](https://doi.org/10.1038/s41592-022-01440-3).

Le job 838443 a réservé un cœur EPYC. Chaque processus dispose de 120 s et de 4 Gio d'espace d'adressage. Les [mesures](838443/trials.csv) distinguent le temps de l'appel et celui du processus. Une interruption laisse le score vide. Les [bornes du chromosome 21](bounds.csv) utilisent les effectifs par lettre et le chemin global qui suit les blocs de la chaîne UCSC hg38ToHs1 numéro 23.

```sh
python3 scripts/human_reference.py --whole-assemblies --cache build/refs --manifest results/human/whole/manifest.tsv
python3 bench/chromosome_bounds.py data/human/chr21-assemblies.fa build/refs/hg38ToHs1.over.chain.gz build/chr21-bounds.csv
python3 bench/human/assembly_edits.py build/editdist data/human/chr20-assemblies.fa build/chr20-edits.tsv
sbatch --account="$ROMEO_PROJECT" bench/whole_chromosomes.sbatch
```
