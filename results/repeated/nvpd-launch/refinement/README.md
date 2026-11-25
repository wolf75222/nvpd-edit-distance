# Confirmation du choix du solveur sur quinze paires indépendantes

À 100 000 bases, deux versions du sélecteur choisissent le moteur qui calcule la distance de Levenshtein. La première retient WFA2 lorsque le coût d’un alignement glouton reste inférieur ou égal à 6 % de la longueur maximale. La seconde essaie aussi un alignement par ancres courtes lorsque la première borne dépasse ce seuil. Si les deux bornes le dépassent, les paires de cette étude passent au noyau CUDA persistant. [Implémentation du sélecteur](../../../../src/hybrid/pair.cpp)

La graine maîtresse 2026 produit cinq nouvelles paires à chacun des trois taux d’éditions introduites, avec les indices 12 à 16. Les deux sélecteurs et les quatre moteurs utilisent les mêmes fichiers. Edlib fournit le score de contrôle. Sur le nœud `romeo-a002`, le pilote épingle le processus au cœur Grace 0, dans le domaine NUMA du GH200. Chaque paire passe dans six blocs. Un processus neuf charge le FASTA, effectue un calcul d’échauffement puis un appel chronométré. La lecture du fichier précède ces appels. Les durées incluent le choix du moteur et les transferts. Le pilote relève le RSS et le PSS entre les deux appels. [Protocole](../../../../bench/repeated/divergence/launch/refine.py) · [Entrées et scores](timing/inputs.tsv) · [Affinité](timing/binding.tsv)

À 5 %, les ancres orientent les cinq paires vers WFA2. Le temps moyen passe de 72,96 à 28,36 ms. La réduction appariée est de 61,12 % (IC à 95 % de 60,86 à 61,39 %). Après échauffement, la mémoire résidente proportionnelle (PSS) passe en moyenne de 151,68 à 3,13 Mio. Sur chacune des cinq paires, WFA2 appelé directement reste plus rapide. Le surcoût moyen du sélecteur atteint 1,18 %. À 2 %, les deux règles choisissent WFA2 pour les cinq paires. La seconde prend 5,34 ms, soit un surcoût de 1,97 % par rapport à WFA2 direct. [Temps et mémoire](analysis/summary.csv)

À 10 %, les cinq paires passent du noyau CUDA ordinaire au noyau persistant. Le sélecteur prend 68,82 ms en moyenne et réduit le temps de 6,31 % par rapport à la première règle (IC à 95 % de 5,95 à 6,67 %). Il devance WFA2 sur les cinq paires de 30,69 % en moyenne (IC à 95 % de 30,07 à 31,31 %) et Edlib de 20,07 %. Sur chaque paire, le rapport des temps est la moyenne des six quotients entre appels d’un même bloc. Les intervalles de Student par taux portent sur les cinq moyennes par paire. [540 mesures et leurs échauffements](timing/trials.csv) · [Rapports par paire](analysis/pairs.csv)

## Reproduction sur ROMEO

Depuis la révision du dépôt associée à cette notice, construire le témoin à partir de `148ecdd` et appliquer le patch de l’allocation initiale. Le fichier `environment.txt` donne les empreintes des sources des deux solveurs. Le script compile les deux versions en mode Release pour `sm_90`, puis vérifie chaque CLI et les tests du sélecteur. Définir `ROMEO_PROJECT` avec le compte Slurm du projet. [Patch du témoin](../job-source.patch) · [Environnement](environment.txt) · [Empreintes des exécutables](timing/binaries.tsv)

```bash
root=$(git rev-parse --show-toplevel)
study=$(mktemp -d "$PWD/refinement.XXXXXXXX")
git worktree add --detach "$study/reference" 148ecdd
git -C "$study/reference" apply "$root/results/repeated/nvpd-launch/job-source.patch"
export CHPS_DEPS="$study/deps"
mkdir "$CHPS_DEPS"
curl -fsSL https://codeload.github.com/Martinsos/edlib/tar.gz/ec2310eda1841ab48c14cd3d866778a4f5eb1491 -o "$CHPS_DEPS/edlib.tar.gz"
curl -fsSL https://codeload.github.com/Marco-Sola-Lab/WFA2-lib/tar.gz/bcf473a6561fa297934a80eaa7c04b4ee412360c -o "$CHPS_DEPS/wfa2.tar.gz"
cd "$root"
sbatch -A "$ROMEO_PROJECT" bench/repeated/divergence/launch/refine.sbatch "$study/reference" "$study/run"
```

Pour recalculer les statistiques à partir des résultats archivés, exécuter depuis la racine du dépôt la commande suivante.

```bash
python3 bench/repeated/divergence/launch/refine.py analyze \
  results/repeated/nvpd-launch/refinement/timing /tmp/nvpd-refinement-analysis
```
