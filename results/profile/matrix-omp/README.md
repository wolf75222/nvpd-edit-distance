# Ordonnancement des matrices NvPD OpenMP

Sur ROMEO, le temps mural médian du lot de 120 paires passe de **227,253 ms** à **43,505 ms** sur huit cœurs EPYC 9654. La médiane des cinq gains appariés est de **5,262×**. Les cinq rapports vont de 5,196× à 5,300×. Le job 838652 exécute les variantes successivement sur les mêmes huit cœurs physiques réservés par Slurm. Chaque fil est lié à un cœur avec OpenMP. Les dix matrices et celle du contrôle local ont la même empreinte edlib. [Mesures EPYC](838652/trials.csv) · [Placement initial](838652/intra-placement.txt) · [Placement inter-paires](838652/inter-placement.txt) · [Compilation](838652/environment.txt) · [Allocation](838652/accounting.tsv)

Pour reproduire l'essai ROMEO depuis le checkout, utiliser `sbatch --account="$ROMEO_PROJECT" --partition=instant bench/cpu/openmp/matrix_omp.sbatch`. Le script utilise GCC 11.5.0 avec `-O3 -DNDEBUG -march=native`. Les deux exécutables partagent les mêmes bibliothèques compilées et ne diffèrent que par leur point d'entrée matriciel. [Script Slurm](../../../bench/cpu/openmp/matrix_omp.sbatch)

La commande `editdist matrix --algorithm nvpd-omp --threads 8` calcule chaque élément du triangle supérieur et écrit la matrice symétrique complète. Si le triangle contient au moins autant de paires que de fils demandés, les fils se partagent des groupes de paires. Une paire est calculée par un seul fil avec la boucle SIMD NvPD. Le réglage du nombre de fils descendants s'applique à la tâche OpenMP courante et empêche la création d'équipes imbriquées actives, y compris avec `OMP_MAX_ACTIVE_LEVELS=4`. Un triangle plus petit conserve une équipe sur les colonnes de chaque paire successive. La commande `pair` suit également ce dernier chemin. [Implémentation](../../../src/cli/main.cpp)

Le banc utilise seize séquences synthétiques de 2 048 à 2 054 bases, soit 120 paires distinctes, et le même exécutable courant pour la matrice de référence edlib. Les essais lancent la CLI en processus neufs. Le chronomètre externe couvre l'ouverture du FASTA, le calcul et la lecture de l'intégralité de la sortie TSV. Chaque variante effectue un premier passage écarté, puis cinq mesures par ordre alterné. Chaque matrice est identique octet par octet à celle d'edlib. Les deux empreintes du CSV identifient l'entrée et le résultat exacts. [Pilote et générateur](../../../bench/cpu/openmp/matrix_omp.py)

## Contrôle local sur poste partagé

Sur le même lot, le temps mural médian passe de **23,191 s** à **0,150 s** sur Apple M1 Pro. Le gain médian des cinq rapports appariés vaut **188,02×**. Les rapports individuels vont de **55,74× à 390,71×**. [Dix mesures locales](trials.csv)

Le M1 Pro possède six cœurs de performance et deux d'efficacité. Les huit fils utilisent libgomp avec placement libre, équipes dynamiques désactivées et imbrication autorisée jusqu'à quatre niveaux. D'autres processus utilisaient la machine pendant le banc. Les rapports décrivent ce lot et cette charge locale. Les deux versions ont été compilées par GCC 16.1.0 avec `-O3 -DNDEBUG -mcpu=native`. La version initiale correspond à `0cc2454b8c4778ddd8599e02d4f6e237c5e02650`. Les deux exécutables utilisent le même `src/cpu/dp.cpp` (SHA-256 `eb8dd38c31d714c0477bc5479a8301b532c6bbfc71cf924738ca9ec57dc70dc1`).

Pour reproduire l'essai, compiler la révision initiale et le checkout courant en mode Release dans deux répertoires, avec le même compilateur et les mêmes options natives, puis exécuter la commande ci-dessous. Le fichier de destination doit être nouveau.

```sh
python3 bench/cpu/openmp/matrix_omp.py build-initial/editdist build-inter/editdist --out /tmp/matrix-omp.csv
```

Un contrôle distinct a comparé à edlib les matrices synthétiques de 1, 2, 4, 5 et 9 séquences, avec des séquences vides, des symboles IUPAC et des longueurs inégales. Un interposeur libgomp en lecture seule a observé `omp_get_active_level()` dans chaque région parallèle des appels de la CLI avec huit fils. Sur deux séquences non vides, il a compté deux régions actives de huit fils, l'index puis les scores. Sur cinq et neuf séquences, il a compté une équipe externe de huit fils et respectivement 12 et 56 régions enfants sérialisées. Le nombre de régions dont le niveau actif dépassait un était nul dans les cinq cas. Ce contrôle a été exécuté séparément du chronométrage.
