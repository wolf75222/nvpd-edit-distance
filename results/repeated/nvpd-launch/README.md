# Lancements CUDA et choix du solveur à 100 knt

Sur le nœud `romeo-a002`, le job Slurm 840841 compare le NvPD CUDA du commit `9a9e23e`, un noyau CUDA persistant et un choix entre WFA2 sur CPU et NvPD sur GPU. WFA2 et Edlib utilisent chacun le cœur Grace 0, sur le nœud NUMA du GPU GH200. Les cinq méthodes calculent le score de Levenshtein des mêmes quinze paires, cinq par taux d’éditions introduites (2, 5 et 10 %). Edlib fournit les scores de contrôle. [Allocation](accounting.tsv) · [Affinité](timing/binding.tsv) · [Empreintes et entrées](timing/inputs.tsv)

Chaque paire passe dans cinq blocs. À la fin du cinquième, chaque méthode a occupé une fois chaque position dans l’ordre de passage. Un processus lit le FASTA et effectue un premier calcul d’échauffement. Pendant qu’il attend la commande suivante, le pilote relève son RSS et son PSS dans `/proc/PID/smaps_rollup`. Il demande ensuite un second calcul sur les mêmes séquences. Le temps mesuré inclut le choix du solveur, l’indexation et les transferts qu’effectue cet appel. Pour chaque paire, le rapport des temps est la moyenne des cinq quotients appariés par bloc. Les intervalles de Student à 95 % portent sur cinq paires indépendantes. [375 mesures et leurs chauffes](timing/trials.csv) · [Intervalles par paire](analysis/pairs.csv)

À 10 %, le noyau persistant prend 68,90 ms en moyenne, contre 73,82 ms pour le NvPD CUDA initial et 100,18 ms pour WFA2. Son rapport moyen à WFA2 vaut 0,688 (IC à 95 % de 0,679 à 0,697) et les cinq paires favorisent NvPD. À 2 %, WFA2 prend 5,31 ms, contre 66,99 ms pour le noyau persistant. Sur la première paire à ce taux, Nsight Systems compte 6 249 lancements du noyau initial et 25 du noyau persistant. La somme des durées de ces noyaux passe de 65,24 à 66,20 ms dans les traces instrumentées. [Comparaison des temps](analysis/summary.csv) · [Figure des rapports à WFA2](../../../figures/comparison/nvpd_launch.png) · [Trace initiale](activity-baseline.csv) · [Trace persistante](activity-persistent.csv)

Le choix CPU/GPU envoie les cinq paires à 2 % vers WFA2. Son temps moyen atteint 5,51 ms, soit un rapport de 1,039 au temps de WFA2 appelé directement (IC à 95 % de 0,983 à 1,094). Les dix paires à 5 et 10 % empruntent le chemin NvPD CUDA initial. À 2 %, le processus hybride occupe après la chauffe 5,09 Mio de RSS et 3,13 Mio de PSS, contre 155,86 et 151,63 Mio pour le processus NvPD CUDA. Ces valeurs décrivent la mémoire hôte au rendez-vous. La colonne `device_bytes` indique séparément le pic des tampons alloués par le solveur pendant l’appel. [Décision pour chaque paire](analysis/routes.csv) · [Temps et mémoires par méthode](analysis/resources.csv)

## Confirmation du routage

La seconde expérience utilise les graines 12 à 16 pour chaque taux. Six blocs comparent les deux politiques et quatre solveurs sur les quinze nouvelles paires. À 5 %, les cinq paires passent de NvPD CUDA à WFA2. Le temps moyen passe de 72,96 à 28,36 ms. À 10 %, le routage révisé appelle le noyau persistant et prend 68,82 ms, contre 99,30 ms pour WFA2. La figure donne les temps individuels des paires et l’intervalle de confiance de leur moyenne. [540 chronométrages](refinement/timing/trials.csv) · [Résultats par paire](refinement/analysis/pairs.csv) · [Figure](../../../figures/comparison/nvpd_dispatch_confirm.png)

L’allocation a aussi servi à mesurer la répartition de chromosomes entiers sur quatre cœurs Grace. Le pilote a terminé les trois expériences et quitté avec le code 3 avant la comparaison supplémentaire des placements NUMA. Slurm a enregistré l’état `FAILED` après 12 min 33 s. Le job 840884 a ensuite terminé les six blocs de lecture et de placement sur quatre cœurs d’une même puce. [Résultats de placement](../../human/grace-scaling/placement/) · [Relevé Slurm](accounting.tsv) · [Mesures sur chromosomes](../../human/grace-scaling/README.md)

## Reproduction

Depuis la racine d’un clone du dépôt sur ROMEO, préparer les deux arbres de sources. Le premier contient le NvPD initial. Le second reprend le commit `148ecdd` et rétablit cinq fichiers dans leur état au moment du job 840841. Les deux exécutables ciblent l’architecture CUDA `sm_90`. [Arbre candidat](https://github.com/wolf75222/nvpd-edit-distance/tree/148ecdd) · [Différence avec le commit de publication](job-source.patch)

```bash
project_root=$PWD
study_dir="$project_root/build/launch"
mkdir -p "$study_dir"
git worktree add --detach "$study_dir/base-src" 9a9e23e
git worktree add --detach "$study_dir/candidate-src" 148ecdd
git -C "$study_dir/candidate-src" apply "$project_root/results/repeated/nvpd-launch/job-source.patch"
if [ -d "$project_root/deps" ]; then
    ln -s "$project_root/deps" "$study_dir/base-src/deps"
    ln -s "$project_root/deps" "$study_dir/candidate-src/deps"
fi
cd "$study_dir/candidate-src"
sbatch -A "$ROMEO_PROJECT" bench/repeated/divergence/launch/run.sbatch \
    "$study_dir/base-src" "$study_dir/run"
```

Après les mesures, le job vérifie le noyau coopératif avec Compute Sanitizer sur une paire de 34 knt à 2 %, puis exécute les deux traces Nsight sur la première paire de 100 knt à 2 %. L’analyse calcule les rapports au solveur effectivement choisi dans le même bloc. [Résultat de Compute Sanitizer](memcheck.txt) · [Versions et matériel](environment.txt) · [Empreintes des trois exécutables](timing/binaries.tsv)

```bash
cd "$project_root"
python3 bench/repeated/divergence/launch/analyze.py "$study_dir/run/timing" \
    "$study_dir/run/analysis" --figure "$study_dir/run/comparison.png"
```

Pour reproduire la confirmation, le répertoire courant contient le routage révisé. Le script utilise l’arbre candidat de la première expérience comme référence et lit les mêmes archives de bibliothèques pour les deux exécutables. [Empreintes de la confirmation](refinement/environment.txt)

```bash
sbatch -A "$ROMEO_PROJECT" bench/repeated/divergence/launch/refine.sbatch \
    "$study_dir/candidate-src" "$study_dir/confirmation"
# Après la fin du job
python3 bench/repeated/divergence/launch/refine_figure.py \
    "$study_dir/confirmation/analysis" "$study_dir/confirmation/comparison.png"
```
