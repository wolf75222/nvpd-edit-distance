# Amortissement d’une copie par puce Grace

Les cœurs logiques 0 et 72 du nœud `romeo-a002` comparent les mêmes chromosomes 1 à 4 de CHM13 entre les assemblages v1.0 et v1.1. Le cœur 0 ouvre et décode les quatre fichiers FASTA une seule fois par processus. Avec le mode `shared`, les deux cœurs lisent directement les chaînes produites par le décodeur. Avec le mode `local`, le cœur 72 copie les huit séquences dans sa mémoire avant le premier alignement. Le volume supplémentaire des séquences est de 1,65 Gio. [Empreintes des entrées](inputs.sha256) · [Banc C++](../../../../bench/human/grace-numa/amortization/amortization.cpp)

Pour un lot de taille $K$, le premier chronomètre part à l’ouverture du premier FASTA et s’arrête au dernier score vérifié. Le second suit le processus du lancement à la sortie, y compris la libération de la mémoire. Dans chaque comparaison, les deux modes parcourent les chromosomes dans le même ordre. Si $K$ est impair, le cœur 0 calcule un score de plus que le cœur 72.

Les fichiers sont relus sur GPFS après la vérification de leurs empreintes. La sonde Linux `mincore` porte sur les pages du fichier vues par le noyau. GPFS dispose aussi de son propre cache de données. Ce protocole mesure donc les temps obtenus lors de lectures répétées sur ce nœud. [Pages des fichiers](cache.csv) · [Cache de GPFS](https://www.ibm.com/docs/en/linux-on-systems?topic=usage-setting-pagepool-size)

## Durée de la commande et délai du dernier score

La campagne ROMEO 841002 a testé toutes les tailles de 1 à 32, puis 64 et 128 alignements. Après le balayage préliminaire, chaque taille a fait l’objet de seize comparaisons appariées, soit 1 088 processus neufs. À chaque taille, les deux ordres de passage des politiques sont représentés à parts égales. Les intervalles de Student sont simultanés à 95 % sur les 34 tailles pour chaque chronomètre (correction de Bonferroni). [Chronométrages](timings.csv) · [Écarts appariés](paired.csv) · [Moyennes et écarts-types](summary.csv) · [Allocation Slurm](accounting.psv)

Pour la commande entière, le gain moyen apparaît à 25 alignements. Il vaut alors 9,9 ms, avec un intervalle de −38,7 à 58,6 ms. L’intervalle exclut zéro à 28, puis le recouvre à 29 et 30. Il reste positif à partir de 31 sur les tailles mesurées. À 64 alignements, le partage prend 3,356 s et la réplication 2,864 s. Si le délai s’arrête au dernier score, le premier gain moyen apparaît à 24 alignements et la borne inférieure reste positive dès 30. [Seuils par chronomètre](threshold.csv) · [Figure](../../../../figures/comparison/grace_numa_amortization.pdf)

## Reproduction

Depuis la racine du dépôt sur ROMEO, préparer les quatre paires de chromosomes et leur manifeste, puis placer les archives épinglées d’Edlib et WFA2-lib dans `deps/`. L’allocation utilise un cœur de chacune des quatre puces pour obtenir les cœurs physiques de l’étude initiale. Seuls les cœurs 0 et 72 exécutent les alignements. [Préparation du panel](../README.md) · [Configuration du projet](../../../../README.md)

```bash
study=bench/human/grace-numa/amortization
sbatch -A "$ROMEO_PROJECT" "$study/run.sbatch" "$PWD/build/grace-amortization-run" \
    "$PWD/build/grace-scaling/inputs.tsv"
```

Les mesures publiées permettent de recalculer les intervalles et de reproduire la figure depuis la racine du dépôt.

```bash
study=bench/human/grace-numa/amortization
data=results/human/grace-scaling/amortization
python3 "$study/analyze_e2e.py" "$data"
python3 "$study/plot.py" "$data" figures/comparison/grace_numa_amortization.png
python3 "$study/plot.py" "$data" figures/comparison/grace_numa_amortization.pdf
```
