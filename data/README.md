# Données et arbres

`python3 scripts/synthetic.py` génère les 60 paires du [manifeste](synthetic/manifest.tsv). Les références contiennent 5 000, 20 000, 100 000 ou 1 000 000 de bases. Chaque longueur comprend cinq paires par taux d’erreur (2, 5 ou 10 %). La graine maîtresse vaut 2026. La mutation choisit des positions sans remise et tire un type d’édition avec probabilité 1/3. Le minimum d’éditions peut être inférieur au nombre d’opérations appliquées.

Le lecteur accepte l’alphabet IUPAC de l’ADN, ignore la casse et compare les symboles littéralement. N/N coûte zéro, N/A coûte un. Les identifiants FASTA sont les premiers mots des en-têtes et doivent être uniques. En mode `matrix`, une seule comparaison par paire alimente la matrice symétrique.

## Cohorte biologique

NCBI Datasets 18.41.1 a fourni les 64 génomes complets du bactériophage φX174. Le manifeste associe une accession versionnée et une empreinte SHA-256 à chaque séquence. Les arbres utilisent les distances d’édition entre les séquences déposées, dans leur orientation et à leur origine publiées. Le volet SARS-CoV-2 contient les métadonnées du panel. [Manifeste φX174](../results/ncbi/phix174.tsv) · [Métadonnées SARS-CoV-2](../results/ncbi/sars_metadata.tsv)

```sh
python3 scripts/ncbi_phix174.py
build/editdist matrix --algorithm wfa2 --threads 8 --input data/phix174/cohort.fa > distances.tsv
python3 scripts/tree.py distances.tsv --method upgma --out upgma.nwk
python3 scripts/tree.py distances.tsv --method nj --out nj.nwk
```

Le téléchargement requiert la commande `datasets`. La construction des arbres utilise la bibliothèque standard de Python. L’option `--groups groupes.tsv` accepte les colonnes `id` et `group` pour vérifier si une branche positive isole chaque groupe. [Documentation NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) · [Matrice et arbres mesurés](../results/romeo/cohort-837271/)
