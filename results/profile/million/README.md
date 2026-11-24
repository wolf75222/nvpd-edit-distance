# NvPD OpenMP sur trois paires d'un million de bases

Le job ROMEO 838538 compare les deux boucles NvPD sur la première paire du manifeste pour chaque taux d'éditions injectées. La référence compte 1 000 000 de bases. Les insertions et délétions du générateur modifient la longueur de la seconde séquence. La version initiale teste la dernière occurrence dans chaque colonne. La version SIMD utilise la formule de minimum et l'annotation OpenMP de vectorisation.

| Taux injecté | Longueur du mutant | Distance | Initiale (s) | SIMD (s) | Rapport initiale / SIMD |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 % | 1 000 033 | 19 801 | 324,413 | 53,043 | 6,116× |
| 5 % | 1 000 084 | 48 887 | 326,446 | 51,9074 | 6,289× |
| 10 % | 999 734 | 95 720 | 322,361 | 51,6137 | 6,246× |

Les trois rapports correspondent aux secondes mesurées sous les deux versions pour la même entrée et le même placement. [Six mesures, CPU et SHA-256](838538/trials.csv) · [Machine et compilation](838538/environment.txt) · [Compte rendu Slurm](838538/accounting.tsv)

Un nœud EPYC 9654 exécute simultanément trois processus de calcul. Chaque paire dispose d'un groupe de huit cœurs physiques. Ses deux variantes utilisent le même groupe, successivement, dans deux processus neufs. L'ordre est initiale puis SIMD à 2 et 10 %, et SIMD puis initiale à 5 %. La mesure couvre un appel au noyau avec construction de l'index, allocations, démarrage de l'équipe OpenMP et calcul du score. La lecture du FASTA et l'échantillonnage des CPU sont hors chronométrage. Il n'y a ni appel d'échauffement ni répétition par variante. Les rapports décrivent ces trois observations sous charge concurrente. Les bandes de variabilité obtenues sur les tailles plus courtes appartiennent à l'autre campagne. [Essais répétés à tailles plus courtes](../cpu/README.md)

Les deux sources de `src/cpu/dp.cpp` viennent des commits `f08b4f6` et `4724f68` et utilisent le même pilote, GCC et les options `-std=c++20 -O3 -DNDEBUG -march=native -fopenmp`. La validation compare le SHA-256, les deux longueurs et le score de chaque mesure aux sorties indépendantes d'edlib et WFA2. Elle regroupe les identifiants de CPU par socket et cœur physique via sysfs avant de vérifier le placement des équipes. [Pilote C++](../../../bench/cpu/openmp/million_pair.cpp) · [Script Slurm et validation](../../../bench/cpu/openmp/profile_million.sbatch)

Depuis la racine d'un clone qui contient les deux commits, soumettre `sbatch --account="$ROMEO_PROJECT" --partition=instant bench/cpu/openmp/profile_million.sbatch`. Si les trois FASTA manquent, les lignes suivantes régénèrent uniquement ces entrées en conservant le manifeste de la campagne principale.

```sh
tmp=$(mktemp -d)
python3 scripts/datasets/synthetic.py --out "$tmp" --lengths 1000000 --rates .02 .05 .1 --pairs 1 --seed 2026
mkdir -p data/synthetic/pairs
cp "$tmp"/pairs/*.fa data/synthetic/pairs/
```
