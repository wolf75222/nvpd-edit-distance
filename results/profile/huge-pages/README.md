# Taille des pages des entrées chromosomiques

Le travail ROMEO 839316 exécute le score sur un cœur EPYC 9654. Les fichiers FASTA contiennent deux révisions complètes de CHM13 par chromosome. Le [manifeste du panel](../../human/chromosome-panel/manifest.tsv) donne leurs longueurs et leurs empreintes SHA-256. Le solveur par fragments et WFA2 sont ceux de la révision [42b6043](https://github.com/wolf75222/nvpd-edit-distance/tree/42b6043). Le témoin edlib vérifie les trois scores avant les chronométrages.

Le banc réserve deux intervalles anonymes dans un même mappage aligné sur 2 Mio. Une garde isole ce mappage dans `/proc/self/smaps`. Le chronomètre de préparation couvre la réservation, le conseil `madvise`, la copie et le premier accès aux pages. Le lecteur FASTA et les allocations internes des solveurs conservent leur politique Linux. La colonne `huge_bytes` consigne la couverture effective du mappage. Le contrôle à 4 Kio utilise `MADV_NOHUGEPAGE` et le traitement à 2 Mio utilise `MADV_HUGEPAGE`.

Chaque chromosome passe par neuf blocs. Dans chacun, le générateur permute les quatre combinaisons de solveur et de taille de page. Chaque configuration effectue un score d’échauffement puis deux lots de trois scores. Le premier lot compte les cycles et les instructions, les défauts de lecture du cache L1D et les demandes de parcours de tables de pages. Le second compte les lignes reçues par le cœur depuis la DRAM locale ou distante. La source de ces événements peut aussi être une adresse d’entrée-sortie mappée. Le calcul porte ici sur des chaînes en mémoire anonyme. La somme des octets de remplissage par seconde porte sur ce cœur en mode utilisateur. Les écritures vers la DRAM et les accès des autres cœurs relèvent des compteurs du contrôleur mémoire.

Le fichier `trials.csv` conserve les compteurs bruts et les durées **cumulés sur trois appels**. `core_active` et `mem_active` donnent la fraction du temps actif de chaque groupe. Le débit de remplissage se calcule avec le temps du second lot et des lignes de 64 octets. Le fichier `events.txt` donne les descriptions et les masques de la PMU. Le noyau réserve aux utilisateurs autorisés les compteurs de socket UMC et RAPL (voir `socket-counters.txt`). Le champ énergétique de Slurm vaut zéro sur cette allocation.

Sur le chromosome 1, `reuse.csv` mesure le mappage et la résolution répétée en temps mural. Cinq blocs alternent les tailles de pages et les nombres de requêtes. Avant chaque mesure, un appel sur les chaînes d’origine chauffe le solveur. Le premier accès aux nouvelles entrées est inclus. La différence de temps appariée vaut le temps à 4 Kio moins celui à 2 Mio. Les intervalles de Student décrivent les blocs du nœud mesuré.

Le script de lancement prend pour argument le dossier contenant les trois FASTA du panel.

```sh
sbatch -A "$ROMEO_PROJECT" bench/cpu/huge-pages/run.sbatch data/human/chm13-panel
python bench/cpu/huge-pages/summarize.py results/profile/huge-pages/839316
```

La description de la gestion des pages transparentes et du champ `AnonHugePages` figure dans la [documentation du noyau Linux](https://docs.kernel.org/admin-guide/mm/transhuge.html).
