# Profilage de NvPD sur ROMEO

Sur ROMEO, Spack fournit MAQAO 2026.1.1 sur EPYC (`x64cpu`, `/inu4pfx`) et sur Grace (`armgpu`, `/ogu2xf5`). Les mesures ci-dessous portent sur un EPYC 9654. [Chronométrages des variantes de NvPD](../dense/README.md)

## Accès aux prédécesseurs

ONEVIEW profile le noyau initial et le noyau à permutation sur deux préfixes homologues de 32 768 bases du chromosome 21. Chacun des 96 alignements donne une distance de 1 479. Sur un cœur, la boucle de mise à jour du noyau initial reçoit 99,18 % des échantillons de cycles. [Profils des boucles dominantes](diagnostics.csv)

Dans la boucle initiale, une instruction rassemble huit scores à des adresses différentes. Le noyau par blocs sélectionne seize scores par permutation en registres. En supposant les données présentes dans le cache L1, le modèle CQA estime le coût d’une itération à 21 cycles pour huit cellules dans le premier cas et à 2,67 cycles pour seize dans le second. [Instructions et diagnostics de boucle](diagnostics.json) · [Inventaire des boucles et des fonctions](diagnostics-objects.csv)

## Avertissements de collecte

À 1 000 Hz, LProf signale que le tampon autorisé par le noyau Linux risque de perdre des échantillons. Cette alerte ne donne pas le nombre d’échantillons perdus. Sur le nœud de calcul, `perf_event_mlock_kb` vaut 516. À 200 Hz, l’alerte disparaît. Charger `hwloc` fournit aussi l’exécutable nécessaire au rapport de topologie. L’échantillonnage de l’espace noyau reste désactivé (`perf_event_paranoid=2`). [Contrôles de collecte](collection.csv)

Le contrôle de durée de ONEVIEW demande au moins dix secondes. À quatre cœurs et 200 Hz, la durée affichée par ONEVIEW passe de 3,82 s pour 96 appels à 15,32 s pour 384 appels au même noyau. ONEVIEW confirme la présence des symboles de débogage et distingue l’absence de pointeur de cadre. Avec une compilation dédiée utilisant `-fno-omit-frame-pointer`, il reconnaît la présence du pointeur de cadre pour reconstruire les piles d’appels. [Contrôles par configuration](collection.csv)

Hors du profileur, huit blocs appariés comparent les deux compilations sur deux paires du chromosome 21, à 32 768 et 100 000 bases, avec un et quatre cœurs. Chaque processus effectue une chauffe. L’ordre des compilations alterne entre les blocs. La variation par bloc est $100(T_{\mathrm{frame}}/T_{\mathrm{ref}}-1)$. Les moyennes par configuration vont de −1,01 % à +1,02 %. Sur la paire homologue de 32 768 bases à quatre cœurs, le surcoût moyen vaut 1,02 % (intervalle de Student à 95 % de 0,13 à 1,92 %). Les sept autres intervalles comprennent zéro. Ces intervalles sont calculés séparément pour chaque configuration. La compilation avec pointeur de cadre sert aux profils de piles d’appels. [128 mesures individuelles](frames-wall.csv) · [Pilote Slurm](../../../../bench/cpu/maqao/frame-cost.sbatch) · [Exécutions sur ROMEO](accounting.tsv)

## Reproduction

Préparer le noyau initial dans `build/maqao/base` et les entrées dans `build/maqao/inputs` selon le protocole de la campagne initiale. Les commandes suivantes se lancent à la racine du dépôt sur ROMEO. [Préparation](../dense/README.md#reproduction)

```sh
romeo_load_x64cpu_env
spack find --long maqao
spack load maqao@2026.1.1
maqao --version

base=build/maqao/base
pair=build/maqao/inputs/near_32768.fa
script=bench/cpu/maqao/versions.sbatch
for config in one four; do
    workers=1
    if [ "$config" = four ]; then workers=4; fi
    sbatch -A "$ROMEO_PROJECT" -p instant -c "$workers" "$script" \
        "$base" . "$pair" "build/maqao/$config" high 2026.1.1
done
sbatch -A "$ROMEO_PROJECT" -p instant -c 4 --export=ALL,MAQAO_HWLOC=/wrpbujm \
    "$script" "$base" . "$pair" build/maqao/control medium 2026.1.1
sbatch -A "$ROMEO_PROJECT" -p instant -c 4 --export=ALL,MAQAO_HWLOC=/wrpbujm,MAQAO_REPEATS=384 \
    "$script" "$base" . "$pair" build/maqao/long medium 2026.1.1
sbatch -A "$ROMEO_PROJECT" -p instant -c 4 \
    --export=ALL,MAQAO_HWLOC=/wrpbujm,MAQAO_REPEATS=384,MAQAO_FRAME_POINTERS=1 \
    "$script" "$base" . "$pair" build/maqao/frames medium 2026.1.1
sbatch -A "$ROMEO_PROJECT" -p instant bench/cpu/maqao/frame-cost.sbatch \
    . build/maqao/inputs build/maqao/frames-wall.csv
```

Chaque archive se décompresse dans un répertoire du même nom, à côté du fichier `.log`. Les noms `one` et `four` identifient les mesures à 1 000 Hz pour un et quatre cœurs. Les répertoires `control`, `long` et `frames` contiennent les contrôles à 200 Hz, avec 96 appels, puis 384, puis 384 avec pointeurs de cadre. [Pilote Slurm](../../../../bench/cpu/maqao/versions.sbatch) · [Extraction des alertes](../../../../bench/cpu/maqao/collection.py)
