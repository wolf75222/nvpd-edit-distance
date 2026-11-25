# Accès mémoire entre deux puces Grace

Le lot compare les chromosomes 1 à 4 de CHM13 entre les assemblages v1.0 et v1.1. Il reprend les quatre fichiers de la campagne de passage à l’échelle. Leurs distances valent 176, 16, 659 et 488 éditions. Un lot comporte seize passages sur les quatre paires, soit 64 scores et 28 344 500 032 bases présentées au solveur, répétitions comprises. [Manifeste de référence](../inputs.tsv) · [Empreintes des fichiers](inputs.sha256)

Les cœurs logiques 0 et 72 de `romeo-a002` appartiennent respectivement aux domaines NUMA 0 et 1. Le solveur `nvpd_sparse` conserve l’ordre des alignements de la campagne initiale. Chaque processus démarre avec des pages neuves de 64 Kio. Un fil attaché au cœur voulu y écrit les séquences, puis le pilote relève les pages résidentes de chaque séquence après l’échauffement et à la fin du calcul. La lecture des FASTA et la préparation des copies précèdent le chronométrage. [Topologie](topology.txt) · [Protocole C++](../../../../bench/human/grace-numa/numa.cpp) · [Relevé de résidence sous Linux](https://man7.org/linux/man-pages/man7/numa.7.html)

| Mode | Cœurs de calcul | Copies et domaine de lecture |
| --- | --- | --- |
| `single` | 0 | Le cœur 0 lit une copie sur sa puce. |
| `shared0`, `shared1` | 0 et 72 | Les deux fils partagent une copie sur le domaine 0 ou 1. |
| `twins0`, `twins1` | 0 et 72 | Chaque fil lit sa propre copie. Les deux copies résident sur le domaine 0 ou 1. |
| `local` | 0 et 72 | Chaque fil lit sa propre copie sur sa puce. |
| `remote` | 0 et 72 | Chaque fil lit sa propre copie sur la puce voisine. |

Chaque mode occupe chacune des sept positions deux fois en quatorze blocs. Le temps moyen, l’écart-type et l’intervalle de Student à 95 % portent sur ces blocs. Pour une différence de durée, la soustraction se fait à l’intérieur de chaque bloc. Les copies privées des modes `local` et `remote` occupent le même volume et le même nombre de pages sur chaque puce. Les modes `twins0` et `twins1` permettent de déplacer une copie sans en changer le nombre. [Mesures individuelles](trials.csv) · [Écarts par bloc](paired.csv) · [Pages par séquence](residency.csv)

## Résultats

Les 98 lots totalisent 6 272 scores chronométrés. Les 2 464 relevés par séquence, copie et phase situent chaque page sur le domaine prévu. Avec une copie par lecteur, les temps moyens sur les quatorze blocs valent 0,6831 ± 0,0033 s en lecture locale et 1,4712 ± 0,0057 s en lecture distante. La pénalité appariée vaut 788,10 ms (IC à 95 % de 784,70 à 791,50 ms). Le rapport des temps vaut 2,154 (de 2,147 à 2,161). [Durées et contrastes](summary.csv) · [Résidence des copies](residency-summary.csv) · [Figure](../../../../figures/comparison/grace_numa_cost.png)

À un cœur, le lot prend 1,2574 s. Avec deux cœurs lisant une copie partagée sur le domaine 0, il prend 1,4725 s. Réserver deux copies sur le domaine 0 fait gagner 0,58 ms (IC à 95 % de −3,07 à 4,23 ms). Déplacer ensuite la seconde copie auprès de son lecteur fait gagner 788,81 ms (de 785,38 à 792,23 ms). Les séquences occupent alors 3,30 Gio, contre 1,65 Gio avec une copie partagée. La table d’attribution rapporte ces différences au ralentissement de 215,10 ms observé dans cette campagne. Le gain de placement dépasse ce ralentissement puisque les deux lecteurs locaux terminent avant le témoin à un cœur. [Décomposition appariée](attribution.csv)

## Reproduction

Le script Slurm demande un cœur sur chacune des quatre puces de `romeo-a002`. Il épingle le calcul sur les cœurs 0 et 72. Une première réservation de deux cœurs avait attribué les cœurs 144 et 216. Le contrôle d’affinité a interrompu cette tentative avant toute mesure. Le job 840942 a ensuite obtenu les cœurs 0, 72, 144 et 216. Les alignements ont utilisé les deux premiers. [Comptabilité Slurm](accounting.psv)

Depuis la racine du dépôt sur ROMEO, utiliser le manifeste de la campagne de passage à l’échelle. Il contient les chemins locaux des quatre FASTA, leurs longueurs et les scores de référence. Le script contrôle leurs empreintes SHA-256 avant de lancer le solveur. Les archives d’Edlib et de WFA2-lib se placent dans `deps/`. [Préparation des entrées](../README.md) · [Compilation et dépendances](../../../../README.md)

```bash
sbatch -A "$ROMEO_PROJECT" bench/human/grace-numa/run.sbatch \
    "$PWD/build/interchip-run" "$PWD/build/grace-scaling/inputs.tsv"
```

Pour recalculer les statistiques et la figure à partir des mesures publiées, exécuter les commandes suivantes depuis la racine du dépôt.

```bash
study=results/human/grace-scaling/interchip
python3 bench/human/grace-numa/analyze.py "$study" --manifest results/human/grace-scaling/inputs.tsv
python3 bench/human/grace-numa/plot.py "$study" figures/comparison/grace_numa_cost.png
```
