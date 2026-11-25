# NvPD et WFA2 sur Grace Hopper

L’étude compare le temps du score de Levenshtein sur un GPU GH200 pour NvPD et sur un cœur Grace pour Edlib et WFA2. Le point de départ est le rapport moyen de 1,044 entre NvPD et WFA2 sur cinq paires de 100 knt à 10 % d’éditions introduites. Le solveur initial correspond au commit `0c916d5`. [Mesures d’origine](../divergence-gh200/README.md)

## Réglage du noyau et confirmation

Le premier essai modifie les tuiles du noyau dense (job Slurm 840625). Les jobs 840636 et 840637 examinent ensuite le calcul dans un corridor, à forte puis à faible divergence. Ces exécutions laissent chaque processus monothread s’exécuter sur les cœurs de son allocation. [Tuiles denses](pilot/trials.csv) · [Corridor à 10 et 5 %](gh200/timing/trials.csv) · [Corridor à 2 %](low-gh200/timing/trials.csv)

Le réglage des tuiles du corridor utilise la paire 11 à 1 Mnt et 2 %, distincte des cinq paires de cette taille retenues pour la confirmation. La paire 1 à 100 knt et 10 % sert de contrôle. Parmi les géométries dont le temps moyen sur le contrôle ne dépasse pas celui de la tuile de 32 × 128, le script choisit la plus rapide sur la paire 11. Un essai séparé mesure ensuite l’effet des ancres sur les deux paires. [Critère de sélection](../../../bench/repeated/divergence/corridor_sweep.py) · [Temps par géométrie](selection/trials.csv) · [Effet des ancres](selection/anchor-trials.csv) · [Paramètres retenus](selection/selected-flags.txt)

La confirmation couvre 25 paires. À 100 knt et 10 %, elle reprend les cinq paires de l’étude d’origine et ajoute les paires 6 à 10. Les trois autres configurations reprennent chacune cinq paires publiées, à 100 knt et 5 %, à 100 knt et 2 %, puis à 1 Mnt et 2 %. Le générateur fixe la graine maîtresse à 2026. Le pilote vérifie les empreintes SHA-256 des entrées publiées et calcule le score de chaque entrée avec Edlib. [Générateur FASTA](../../../scripts/datasets/synthetic.py) · [Pilote de mesure](../../../bench/repeated/divergence/focus.py)

## Temps par paire

Sur le nœud `romeo-a003`, le job Slurm 840653 utilise la tuile de 16 × 128 avec ancres. La confirmation correspond à la reprise 1 du job. Chaque processus chronométré hérite du cœur Grace 72 sur le nœud NUMA 1 du GPU. [Allocation Slurm](confirmation/accounting.tsv) · [Affinité](confirmation/timing/binding.tsv) · [Configuration du noyau](confirmation/configuration.txt) · [Matériel et compilation](confirmation/environment.txt)

Pour chaque paire, huit blocs alternent quatre méthodes selon un plan équilibré. Un processus lit le FASTA, effectue un alignement d’échauffement, puis chronomètre le second appel. Le temps inclut l’indexation effectuée par le solveur, les allocations CUDA et les transferts nécessaires au score. `timing/trials.csv` conserve les deux durées.

Pour une paire $p$ et une méthode de référence $r$, le rapport moyen porte sur les quotients calculés au sein des huit blocs. Les cinq paires ont ensuite le même poids dans la moyenne de la configuration.

$$
R_p=\frac{1}{8}\sum_{b=1}^{8}\frac{T_{p,b,\mathrm{NvPD}}}{T_{p,b,r}},
\qquad
\bar R=\frac{1}{5}\sum_{p=1}^{5} R_p,
\qquad
I_{95}=\bar R\pm t_{0.975,4}\frac{s_R}{\sqrt{5}}.
$$

Une paire compte comme victoire lorsque $R_p<1$. L’écart-type $s_R$ décrit la dispersion entre les cinq paires. Le fichier brut contient les 800 temps mesurés et les 800 temps de chauffe. [Chronométrages](confirmation/timing/trials.csv) · [Calcul des intervalles](../../../bench/repeated/divergence/focus_analysis.py)

À 100 knt et 10 %, le rapport NvPD/WFA2 vaut 0,740 sur les cinq paires initiales (IC à 95 % de 0,733 à 0,747) et 0,742 sur les cinq supplémentaires (de 0,731 à 0,752). NvPD précède WFA2 et Edlib sur chacune des dix paires. À 1 Mnt et 2 %, les temps moyens de NvPD et WFA2 atteignent 0,759 et 0,762 s. Le rapport moyen vaut 0,996 (de 0,986 à 1,007) et deux paires sur cinq favorisent NvPD. [Résultats par paire](confirmation/analysis/pairs.csv) · [Moyennes par configuration](confirmation/analysis/summary.csv) · [Figure des rapports](../../../figures/comparison/nvpd_focused.png)

## Empreinte mémoire et temps de la cohorte

`device_bytes` mesure le maximum des octets alloués aux tampons gérés par NvPD. Les bibliothèques CPU ont une valeur nulle pour cette métrique. Une campagne séparée mesure la mémoire hôte sur les paires à 2 %. Un observateur épinglé sur un autre cœur relève le RSS du processus dans `/proc/PID/smaps_rollup`, avec un intervalle nominal de 1 ms. La plus grande valeur observée couvre la chauffe et le calcul suivant. Le journal donne l’intervalle maximal réellement observé. Les zéros de `getrusage` sur Grace correspondent à des mesures indisponibles. [Collecte du RSS](../../../bench/repeated/divergence/focus_memory.py)

Sur les cinq paires historiques à 100 knt et 10 %, deux modes exécutent le même lot sur un seul GPU. Un processus parcourt les cinq paires avec un espace de travail CUDA. L’autre mode affecte une paire à chacun de cinq processus simultanés. Après la chauffe, les processus attendent une commande. À ce rendez-vous, le pilote additionne les RSS, puis les PSS. Le RSS compte les pages partagées dans chaque processus qui les utilise. Le PSS répartit le coût de chaque page entre les processus qui la partagent. Le chronométrage du lot suivant commence à l’envoi de la commande et s’arrête à la réception du dernier score. Les deux modes alternent sur cinq blocs et utilisent le même ordre de paires dans chaque bloc. L’intervalle du rapport concurrent/séquentiel décrit les répétitions de cette cohorte fixe. [Mesure du lot et des mémoires résidentes](../../../bench/repeated/divergence/cohort_memory.py)

À 1 Mnt et 2 %, le pic moyen des tampons GPU passe de 22,890 à 1,588 Mio. Le maximum RSS hôte observé vaut 159,55 Mio pour le processus NvPD optimisé, 9,63 Mio pour Edlib et 10,07 Mio pour WFA2. Sur le lot de cinq paires à 100 knt et 10 %, un processus traite les cinq scores en 0,3647 s avec 154,92 Mio de PSS. Cinq processus simultanés prennent 0,3890 s avec une somme de PSS de 726,23 Mio. Le rapport apparié des PSS vaut 4,688 (IC à 95 % de 4,687 à 4,689). [Mémoire et temps par méthode](confirmation/analysis/resources.csv) · [Mesures de la cohorte](confirmation/cohort/cohort.csv) · [Intervalles de la cohorte](confirmation/analysis/cohort-summary.csv) · [Figure des mémoires GPU et hôte](../../../figures/comparison/nvpd_focus_memory.png)

## Reproduction sur ROMEO

Depuis un clone du dépôt, préparer le FASTA synthétique et compiler l’exécutable initial sur un nœud Grace Hopper. Le script de confirmation utilise les archives des dépendances du dépôt lorsqu’elles se trouvent dans `deps/`.

```bash
python3 scripts/datasets/synthetic.py
mkdir -p build/focus-repro
git worktree add --detach build/focus-repro/base-src 0c916d5
if [ -d deps ]; then ln -s "$PWD/deps" build/focus-repro/base-src/deps; fi
base_job=$(sbatch --parsable -A "$ROMEO_PROJECT" -p instant -C armgpu -N 1 -c 4 \
    --gpus=1 --mem=8G -t 00:15:00 <<'BATCH'
#!/bin/bash -l
set -eo pipefail
romeo_load_armgpu_env >/dev/null
spack load /phhwuh3
spack load /dvle736
base="$PWD/build/focus-repro/base-src"
cmake -S "$base" -B build/focus-repro/base -DCMAKE_BUILD_TYPE=Release \
    -DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 \
    -DCMAKE_PROJECT_INCLUDE="$base/bench/repeated/targets.cmake"
cmake --build build/focus-repro/base --parallel 4 --target repeat_kernel
BATCH
)
sbatch -A "$ROMEO_PROJECT" --dependency="afterok:${base_job}" \
    bench/repeated/divergence/focus.sbatch \
    build/focus-repro/base/repeat_kernel build/focus-repro/confirmation select
```

Le mode `select` exécute le réglage avant la confirmation. Pour reprendre directement la géométrie publiée, fournir le chemin absolu de `selection/selected-flags.txt` dans `CHPS_NVPD_SELECTION` au moment de la soumission. Les empreintes des exécutables figurent dans `timing/binaries.tsv`. Les cinq combinaisons de longueur, de taux et de graines sont décrites dans `timing/inputs.tsv`. Pour produire les deux figures et les résumés après la fin du job, utiliser Python avec les dépendances d’analyse du projet.

```bash
run=build/focus-repro/confirmation
mkdir -p "$run/analysis"
python3 bench/repeated/divergence/focus_analysis.py "$run/timing" \
    --output "$run/analysis" --smaps "$run/memory/smaps.csv" \
    --cohort "$run/cohort/cohort.csv" --figure "$run/analysis/times.png" \
    --memory-figure "$run/analysis/memory.png"
```
