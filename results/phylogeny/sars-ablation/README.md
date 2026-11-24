# Décomposition du calcul par fragments

Le job ROMEO 839296 traite les 2 016 paires du corpus sur 32 cœurs physiques d’un AMD EPYC 9654. Les quatre méthodes utilisent le même binaire et le même ordonnancement OpenMP. [Environnement](environment.txt)

Pour chaque mode de sortie, le bloc zéro échauffe les méthodes. Les cinq blocs suivants chronomètrent le calcul sur les séquences déjà chargées, en faisant varier l’ordre de passage. Les intervalles de Student à 95 % portent sur ces cinq mesures dans une allocation. Les accélérations comparent deux durées du même bloc. [Mesures](trials.csv) · [Moyennes et intervalles](summary.csv)

BiWFA traite la paire entière. `trim-select` enlève les extrémités identiques, puis choisit NvPD pour un fragment dense de 262 144 cases au plus et de côté maximal 20 000. BiWFA traite les fragments plus grands. `split` ajoute les coupes internes. `split-filter` calcule les effectifs de lettres une fois par séquence. Le chronomètre inclut ce comptage.

Soient $a_c$ et $b_c$ les effectifs du symbole $c$. Une édition élimine au plus une unité d’excédent et une unité de déficit. La distance vérifie donc $d\ge L=\max\{\sum_c(a_c-b_c)_+,\sum_c(b_c-a_c)_+\}$.

Si $L>|m-n|$, le coût admissible $U$ dépasse la borne de longueur. Le certificat exige alors $U\le 1\,024$. Chaque coupe doit aussi être distante de ses voisines d’au moins $64(U+1)$ bases dans les deux séquences. Le filtre écarte la découpe lorsque $L>1\,024$ ou lorsque toutes les coupes manquent d’espace à $U=L$.

Pour le score, le filtre ramène le temps moyen de 589,77 à 451,83 ms. BiWFA seule prend 446,61 ms. Les deux variantes conservent 112 coupes sur 86 paires. Le filtre réduit le nombre d’appels de fronts de 6 505 à 3 998, essais de découpe compris.

Pour la reconstruction, le temps moyen passe de 1,220 à 0,929 s. Les 48 matrices coïncident avec la référence edlib. Le pilote rejoue les 48 384 CIGAR obtenus sur les 24 appels avec reconstruction. [Matrice de référence](../sars-edit/matrix.tsv)

Le job 839305 répète le protocole avec les écritures par octets du traceback WFA2. Sans filtre, le score prend 583,51 ms et la reconstruction 1 203,72 ms. Avec le filtre, ces moyennes passent à 448,16 et 918,48 ms. Les accélérations appariées valent respectivement 1,302 et 1,311. Les effectifs de coupes restent identiques entre les deux allocations. [Mesures](839305/trials.csv) · [Intervalles](839305/summary.csv) · [Comparaison du code compilé](../../alternatives/wfa-portability/README.md)

Depuis la racine du dépôt, les commandes suivantes reproduisent la mesure et son résumé.

```sh
sbatch --account="$ROMEO_PROJECT" bench/cpu/surveillance/run.sbatch \
  data/sars/cohort.fa results/phylogeny/sars-edit/matrix.tsv
python3 bench/cpu/surveillance/summarize.py results/phylogeny/sars-ablation/trials.csv \
  > results/phylogeny/sars-ablation/summary.csv
```
