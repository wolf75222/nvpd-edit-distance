# Comparaison des profils du répartiteur

Sur quatre cœurs physiques de l’EPYC 9654 de ROMEO, MAQAO 2025.1.4 et 2026.1.1 ont profilé successivement le même exécutable `editdist-route`. Chaque processus a calculé la matrice des 64 génomes une fois pour la chauffe puis six fois sous la politique `hybrid-counts`. Les deux matrices coïncident octet pour octet avec la référence. Les journaux de décision sont également identiques, avec 610 paires traitées par Edlib et 1 406 par WFA2. [Empreintes des entrées et de l’exécutable](environment.txt) · [Allocation 840307](accounting.tsv)

Myers dans Edlib reçoit 84,90 % puis 84,96 % des échantillons de cycles. L’extension AVX-512 de WFA2 en reçoit 10,97 % puis 11,09 %, et la construction des esquisses 0,143 % puis 0,145 %. Chaque profil fournit une observation du processus entier. Les proportions portent sur les événements classés par LProf dans cette exécution, et la couverture d’une fonction inclut ses boucles internes. La catégorie OpenMP représente 0,293 % puis 0,238 %. [Fonctions dont la couverture atteint 0,1 %](functions.csv) · [Durées des profils et catégories](profiles.csv)

## Diagnostics des boucles

LProf repère 23 boucles communes aux deux profils. Pour les 17 d’entre elles disposant d’une analyse CQA, les fichiers de diagnostic sont identiques octet pour octet, jusque dans le texte des remèdes proposés. Les 12 rapports de fonction communs ont aussi des empreintes identiques. Deux fonctions de WFA2 reçoivent un rapport seulement dans le profil 2026.1.1. Pour le calcul des bornes du front (`wavefront_compute_trim_ends`, 0,0173 % des échantillons), CQA invite à examiner la vectorisation et les conversions. Pour l’allocation des fronts (`wavefront_slab_allocate`, 0,0189 %), il relève aussi les appels de fonctions. Une seule collecte par version ne permet pas d’attribuer cette différence de couverture à une nouvelle règle d’analyse. Les boucles visibles dans un seul profil reçoivent chacune environ 0,003 % des échantillons. [Couvertures, adresses et empreintes des boucles](loops.csv) · [Empreintes des rapports de fonction](function_cqa.csv)

La boucle interne de Myers reçoit 81,17 % des échantillons dans les deux exécutions. Les deux versions de CQA estiment son itération assembleur à 11 cycles et proposent notamment de réorganiser les accès aux structures, d’examiner les dépendances entre itérations et d’essayer le déroulage. Pour la boucle externe d’extension WFA2, le modèle Zen 4 de CQA signale qu’il ne sait pas analyser au moins une instruction. La colonne `model_cycles` reste vide de valeur numérique pour cette boucle. Les indications sur les conversions et les appels de la boucle de construction des esquisses restent les mêmes.

## Diagnostic de la collecte

ONEVIEW 2026.1.1 identifie les symboles de débogage et précise que l’option `-fno-omit-frame-pointer`, utile à la reconstruction des piles d’appels, manque à la compilation. Les deux collectes signalent les restrictions d’accès aux compteurs du noyau et l’absence de `lstopo-no-graphics` dans l’environnement du job. Au démarrage, la version 2026.1.1 signale aussi un risque de perte d’échantillons lié à la limite de mémoire verrouillable du sous-système Linux `perf`. Le journal ne chiffre pas cette perte. [Indicateurs ONEVIEW](oneview.csv)

## Reproduction

Depuis la racine du dépôt sur ROMEO, préparer le FASTA de la cohorte puis lancer la commande suivante dans un répertoire de sortie neuf. Le script compile l’exécutable une seule fois et résout chaque version avec `spack location -i`. Les deux archives contiennent la page d’accueil HTML sous `RESULTS/editdist-route_one_html/index.html` et les diagnostics de chaque boucle sous `static_data/cqa/`. Les colonnes `*_file_sha256` des exports donnent les empreintes SHA-256 des fichiers natifs.

```sh
sbatch -A "$ROMEO_PROJECT" bench/cpu/maqao/route.sbatch \
  build/hybrid-matrix/viral/surveillance.fa results/phylogeny/sars-edit/matrix.tsv \
  build/revue-maqao/route-versions 2025.1.4 2026.1.1
```
