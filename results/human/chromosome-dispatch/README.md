# Révisions chromosomiques de HG002

Les FASTA v1.1 et v1.2 du consortium T2T fournissent les chromosomes 1 et 21 des haplotypes maternel et paternel de HG002. Le manifeste donne la longueur et le SHA-256 de chaque suite de bases. Le script applique le VCF de révision à chaque chromosome v1.1 et compare le condensat obtenu à celui du FASTA v1.2. Les deux assemblages conservent les lacunes représentées par N. [Assemblages et VCF](https://github.com/marbl/HG002) · [Hansen et coll., 2026](https://doi.org/10.1016/j.cell.2026.06.016).

Depuis la racine du dépôt, préparer les entrées puis soumettre la campagne sur ROMEO avec les archives de dépendances décrites dans le README principal.

```sh
python3 bench/cpu/chromosome-dispatch/prepare.py data/human/hg002
sbatch -A "$ROMEO_PROJECT" bench/cpu/chromosome-dispatch/run.sbatch data/human/hg002
```

Chaque ligne de `trials.csv` correspond à un processus. Le temps commence après la lecture des deux FASTA. Le processus appelle le solveur une première fois pour l'échauffement, puis mesure le second appel. Les cinq blocs permutent les configurations. Le RSS inclut les chaînes, leur lecture et les allocations du solveur. Une seconde série à un cœur mesure les instructions et les cycles du noyau, ainsi que les lectures et les défauts de lecture L1. Les entrées de `counters.csv` donnent le compte brut et le pourcentage de temps de programmation du compteur.
