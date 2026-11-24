# Vérification du correctif WFA2 sur EPYC

Le job 838814 compare la révision WFA2-lib `bcf473a6561fa297934a80eaa7c04b4ee412360c` avec les deux modifications appliquées par CMake. Le noyau d’extension charge les mots de huit octets avec `memcpy`, qui accepte toute adresse de caractère. L’allocateur maintient son curseur une position après la dernière requête vivante. L’expérience porte sur les exécutables optimisés sans instrumentation. [Transformation des sources](../../../cmake/patches/wfa2-portability.cmake) · [Empreintes des sources et du pilote](838814/sources.sha256)

Les deux builds utilisent le même chemin temporaire, la même archive WFA2 vérifiée par CMake, GCC 11.5 et les options `-O3 -DNDEBUG -march=native`. Pour la variante amont, le script neutralise la transformation dans une copie temporaire du projet. Il extrait chaque section exécutable des deux objets modifiés avec GNU objcopy, puis compare les SHA-256 de ces sections, des fichiers ELF entiers et du désassemblage avec relocalisations. La CLI de production et le pilote chronométré sont aussi comparés. [Empreintes du code machine](838814/machine-code.csv) · [Compilateur et placement](838814/environment.txt)

Pour les cinq paires du job 838762, le pilote vérifie les SHA-256 et les distances globales avec edlib. Huit blocs visitent les cinq entrées en ordre aléatoire. Pour chaque entrée et bloc, deux processus frais sont exécutés successivement sur le même cœur physique. L’ordre amont/corrigé alterne et chaque version passe en premier dans quatre blocs. Chaque processus lit le FASTA, fait un appel d’échauffement, puis chronomètre dix appels WFA2 en mode coût seul avec la même horloge et le même adaptateur. L’unité expérimentale est la moyenne des dix appels. Les intervalles à 95 % sont les intervalles bilatéraux de Student sur huit moyennes de processus, ou sur huit rapports appariés pour la comparaison corrigé/amont. [800 appels individuels](838814/trials.csv) · [Moyennes, écarts-types et intervalles](838814/summary.csv)

Le fichier objet du noyau d’extension est identique octet par octet, y compris sa section machine et les relocalisations du désassemblage. La section machine de l’allocateur diffère, tout comme celle de la CLI et du pilote. Les 800 appels ont donné les distances attendues. Le tableau présente les temps moyens et le rapport calculé au sein de chaque bloc.

| Entrée | Amont (ms) | Corrigé (ms) | Corrigé/amont [IC95] |
| --- | ---: | ---: | ---: |
| 20k-0.02 | 0,2013 | 0,1930 | 0,960 [0,897, 1,023] |
| 20k-0.05 | 1,0480 | 1,0469 | 0,999 [0,989, 1,009] |
| 20k-0.1 | 3,7593 | 3,7479 | 0,997 [0,992, 1,002] |
| chr21-50k | 0,5184 | 0,5190 | 1,001 [0,982, 1,021] |
| chr21-center128 | 0,0399 | 0,0433 | 1,088 [0,735, 1,440] |

Les cinq intervalles de rapport incluent 1. Leur largeur dépend de l’entrée, en particulier pour la paire centrale de distance 72 dont un appel dure environ 30 µs. Ces mesures délimitent la précision de la comparaison sur ce cœur EPYC.

Depuis la racine du dépôt, exécuter `python3 bench/quicked/cases.py data/controls/quicked` après génération des jeux synthétiques et des fenêtres de référence, puis `sbatch -A r260105 -p instant bench/cpu/alternatives/wfa_portability.sbatch` sur ROMEO. `CHPS_WFA_INPUTS` permet d’utiliser un autre emplacement contenant les mêmes FASTA. [Script Slurm](../../../bench/cpu/alternatives/wfa_portability.sbatch) · [Collecteur apparié](../../../bench/cpu/alternatives/wfa_portability.py) · [Manifeste de référence](../romeo/838762/cases.csv)
