# Parallélisme entre chromosomes sur Grace

Le job ROMEO 840841 compare les chromosomes 1 à 4 entre les assemblages CHM13 v1.0 et v1.1 du même échantillon. Les distances valent respectivement 176, 16, 659 et 488 éditions. Les huit séquences totalisent 1 771 531 252 nt. Chaque cœur exécute le solveur `nvpd_sparse` sur une paire à la fois. [Entrées](inputs.tsv) · [Empreintes du FASTA de chaque paire](inputs.sha256) · [Assemblages et scores de référence](../chromosome-panel/README.md)

À charge fixe, le lot comprend 16 passages sur les quatre paires, soit 64 alignements et 28,345 Gnt de séquences présentées au solveur. À charge par cœur fixe, chaque cœur effectue quatre passages, soit 16 alignements et 7,086 Gnt. Les passages réemploient les mêmes chromosomes. Le débit exprime ce volume logique en Gnt/s.

Chaque bloc démarre un processus sur le nœud `romeo-a002`. Le fil principal charge les séquences, puis WFA2 vérifie les distances. Un échauffement précède les mesures à chaque effectif. Les six blocs couvrent les six ordres possibles pour un, deux et quatre cœurs. Le chronomètre inclut la région OpenMP et les appels au solveur. Les quatre cœurs appartiennent à quatre processeurs Grace du même nœud. [Placement des fils](environment.txt) · [Mesures et échauffements](trials.csv)

Pour chaque bloc, l’accélération à charge fixe vaut $S_p=T_1(W)/T_p(W)$ et l’efficacité vaut $E_p=S_p/p$. À charge par cœur fixe, $E_p=T_1(w)/T_p(pw)$. L’intervalle à 95 % de la moyenne vaut $\bar x\pm t_{0,975;5}s/\sqrt{6}$. Il décrit les répétitions sur ce nœud, avec les quatre mêmes paires. [Rapports par bloc](paired.csv) · [Moyennes, écarts-types et intervalles](summary.csv)

Le lot de 64 alignements prend 1,242 s avec un cœur, 1,467 s avec deux cœurs et 0,746 s avec quatre cœurs. À quatre cœurs, l’accélération atteint 1,665 (écart-type 0,008, IC à 95 % [1,657 ; 1,674]). L’efficacité vaut 41,6 % à charge fixe et 42,0 % à charge par cœur fixe. Sur quatre cœurs d’une même puce, une seconde expérience attribue la lecture de chaque fichier au cœur qui traite la paire. Le temps de lecture passe de 1,826 à 0,619 s. Les pages des séquences résident sur le même nœud NUMA dans les deux modes. [Mesure du placement](placement/README.md) · [Figure](../../../figures/comparison/grace_chromosome_scaling.png)

La comparaison sur les cœurs 0 et 72 reprend ensuite le lot fixe avec sept placements. Le surcoût des lectures distantes atteint 788,1 ms pour 64 alignements. Avec une copie locale auprès de chaque lecteur, le temps à deux cœurs descend sous celui à un cœur. [Comparaison entre puces](interchip/)

Lorsque chaque processus relit les quatre FASTA et prépare ses propres copies, la réplication réduit en moyenne le temps de la commande dès 25 alignements. À 31 alignements, le gain atteint 73,1 ms et sa borne inférieure reste positive pour les tailles supérieures mesurées. [Coût de bout en bout](amortization/)

Depuis la racine du dépôt, la commande suivante recalcule les rapports appariés et la figure.

```sh
python3 bench/repeated/divergence/launch/scaling_analysis.py results/human/grace-scaling \
    --figure figures/comparison/grace_chromosome_scaling.png
```

Pour relancer les mesures sur ROMEO, préparer les quatre FASTA à partir des assemblages publiés, puis soumettre les deux étapes depuis la racine du dépôt. Slurm lance la comparaison des lecteurs lorsque la mesure à un, deux et quatre cœurs a terminé. Le fichier `environment.txt` de chaque étape donne les cœurs effectivement attribués. [Préparation du panel](../chromosome-panel/README.md)

```bash
python3 scripts/datasets/chm13_panel.py --cache build/t2t-cache --out data/human/chm13-panel
campaign=$(mktemp -d "$PWD/build/grace-chromosomes.XXXXXXXX")
launch=bench/repeated/divergence/launch
scaling_job=$(sbatch --parsable -A "$ROMEO_PROJECT" "$launch/scaling.sbatch" "$campaign/scaling" data/human/chm13-panel)
sbatch -A "$ROMEO_PROJECT" --time=00:10:00 --dependency="afterok:${scaling_job%%;*}" "$launch/placement.sbatch" \
    "$campaign/placement" "$campaign/scaling"
```
