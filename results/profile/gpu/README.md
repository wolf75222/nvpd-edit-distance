# Profil CUDA sur ROMEO

Les quatre allocations ont utilisé un GH200, quatre cœurs Arm et CUDA 13.0.2. Chaque variante a utilisé les mêmes sources hôtes, les mêmes entrées et les options de compilation `nvcc -std=c++20 -O3 -arch=sm_90 -Xcompiler=-fopenmp`. Les fichiers `variants.sha256` identifient la traduction CUDA. Les sources hôtes de programmation dynamique correspondent au commit `4724f68`.

| Job | Comparaison |
|---|---|
| 838065 | Premier profil de la version à une ligne par lancement, commit `f08b4f6`. |
| 838161 | Version initiale contre tuiles de huit lignes et 256 colonnes, commit `9afad07`, quatre longueurs à 2 %. |
| 838170 | Version initiale contre tuiles de largeur adaptative, commit `21d7b77`, six longueurs à 2 %. |
| 838173 | Même comparaison sur quatre paires à 10 %. Contrôle des 2 016 distances de la cohorte φX174 dans le programme de production. |

Pour chaque paire et chaque variante, le banc lance trois processus. Chaque processus effectue un premier calcul à froid, puis trois calculs à chaud. L’ordre des variantes alterne entre les processus. Le chronomètre couvre la construction de la table des dernières correspondances, les allocations, les transferts, les noyaux et la libération des tampons. Il exclut la lecture du FASTA. Les essais à froid comprennent aussi le démarrage du contexte CUDA. Chaque appel vérifie sa distance.

## Temps à chaud

Pour chaque processus, le banc retient la médiane de ses trois appels à chaud. Le tableau donne la médiane des trois résultats par variante, en secondes. Il divise ces deux médianes pour obtenir le rapport. Les essais à 2 % et à 10 % proviennent d’allocations distinctes. Les résultats individuels et les temps à froid restent consultables dans `trials.csv` pour chaque job.

| Longueur | Err. | Initial | Adaptatif | Rapport initial/adaptatif |
|---:|---:|---:|---:|---:|
| 5 000 | 2 % | 0,01774 | 0,01808 | 0,98 |
| 20 000 | 2 % | 0,05156 | 0,02306 | 2,24 |
| 100 000 | 2 % | 0,24528 | 0,10630 | 2,31 |
| 250 000 | 2 % | 0,72782 | 0,73540 | 0,99 |
| 500 000 | 2 % | 2,04855 | 1,14900 | 1,78 |
| 1 000 000 | 2 % | 5,14159 | 3,61835 | 1,42 |
| 5 000 | 10 % | 0,01103 | 0,00457 | 2,42 |
| 20 000 | 10 % | 0,05193 | 0,01613 | 3,22 |
| 100 000 | 10 % | 0,28864 | 0,12128 | 2,38 |
| 1 000 000 | 10 % | 5,15160 | 3,64234 | 1,41 |

La tuile fixe avait ralenti la paire à un million de bases de 2,4 % dans le job 838161. La largeur adaptative passe de 256 à 1 024 colonnes lorsque le nombre de blocs larges atteint deux par multiprocesseur. Cette règle réduit le nombre de colonnes recalculées dans le chevauchement des tuiles larges. Sur le job 838170, les appels adaptatifs à chaud de 5 000 bases varient de 4,13 à 182,37 ms. Les mesures courtes sur nœud partagé dépendent donc fortement de l’essai.

## Mesure des noyaux

Nsight Systems 2025.3.2 capture un appel distinct après un premier appel complet. Pour la paire de 5 000 × 4 992 bases du job 838170, la version initiale lance 4 992 noyaux de récurrence et cumule 6,456 ms de calcul GPU. La version adaptative lance 624 noyaux et cumule 3,338 ms. Les deux variantes lancent aussi un noyau d’initialisation. Les transferts aller et retour cumulent 9,12 µs pour la version initiale et 8,03 µs pour la version adaptative. Les temps des appels API hôtes peuvent recouvrir le calcul des noyaux. Les essais `trials.csv` s’exécutent sans le profileur.

La tuile couvre huit lignes au maximum. Chaque bloc recalcule à gauche les cellules dont il dépend. Une ligne nécessite au plus `reach` colonnes supplémentaires, déterminées par la dernière occurrence de chaque symbole. Le programme utilise la version à une ligne si `reach > 64`, pour le traceback ou pour moins de huit lignes. Les octets `device_bytes` comptent les allocations explicites. Une tuile utilise en plus deux lignes de mémoire partagée, soit au plus 11 776 octets par bloc large. Les compteurs matériels étaient inaccessibles aux deux variantes sur cette allocation (`ERR_NVGPUCTRPERM`). Les débits de mémoire DRAM et le taux d’occupation mesuré restent inconnus.

## Reproduire

Générer les paires synthétiques avec la graine maître 2026. Les commandes suivantes utilisent la version de référence épinglée et les fichiers exacts de la première paire de chaque configuration. En dehors d’un checkout Git, exporter la source `f08b4f6:src/cuda.cu` dans un fichier, puis définir `GPU_BASELINE` sur son chemin. Le banc vérifie son empreinte.

```sh
python3 scripts/synthetic.py --out data/synthetic --seed 2026
sbatch --account="$ROMEO_PROJECT" bench/profile_gpu.sbatch profile
GPU_CANDIDATE_REV=9afad07 GPU_INTERMEDIATE=0 sbatch --account="$ROMEO_PROJECT" bench/profile_gpu.sbatch compare
GPU_CANDIDATE_REV=21d7b77 GPU_COUNTERS=1 sbatch --account="$ROMEO_PROJECT" bench/profile_gpu.sbatch compare
GPU_CANDIDATE_REV=21d7b77 sbatch --account="$ROMEO_PROJECT" bench/profile_gpu.sbatch holdout
```

Le mode `holdout` utilise aussi `data/phix174/cohort.fa`, préparé par `scripts/ncbi_phix174.py`. Les trois moteurs du job 838173 donnent une matrice identique, enregistrée par son empreinte dans `cohort.csv`. L’exécutable CUDA prend 43,51 s pour cette matrice, edlib 0,0543 s et WFA2 0,00986 s. [Mesures de la cohorte](838173/cohort.csv)
