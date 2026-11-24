# Distribution du calcul sur deux nœuds EPYC

Les mesures principales couvrent un et deux nœuds ROMEO à raison de 24 rangs MPI et de huit threads OpenMP par rang. Chaque configuration compte sept lancements. Un lancement exécute deux échauffements sur l’entrée chronométrée puis un calcul. L’ordre des configurations change à chaque bloc. Le programme vérifie les cœurs physiques affectés et compare les distances à edlib. [Protocole Slurm](../../bench/mpi/repeated.sbatch).

Pour le scaling fort, le nombre de cellules reste constant. Pour le scaling faible, le nombre de cellules par cœur reste constant. Les accélérations et les efficacités sont calculées dans chaque bloc avant de prendre la moyenne. Les fichiers de synthèse donnent l’écart-type expérimental et les bornes de l’intervalle de Student bilatéral à 95 %. Les répétitions ont lieu sur la même allocation. [Essais de scaling](scaling/trials.csv) · [Statistiques](scaling/summary.csv) · [Allocation](scaling/accounting.tsv).

La sonde PMPI chronomètre les appels à `MPI_Allgather` sur le rang qui détermine la durée de la paire. Ce temps inclut l’attente à l’arrivée des autres rangs. Les essais sans sonde donnent le temps de référence pour chaque forme de matrice. [Essais avec et sans sonde](collectives/trials.csv) · [Statistiques](collectives/summary.csv) · [Allocation](collectives/accounting.tsv).

```bash
sbatch --account="$ROMEO_PROJECT" --partition=instant bench/mpi/repeated.sbatch scaling
sbatch --account="$ROMEO_PROJECT" --partition=instant bench/mpi/repeated.sbatch collectives
```

Les réservations antérieures regroupaient plusieurs chronométrages dans le même lancement. Leurs résultats restent accessibles par identifiant Slurm. [Scaling fort](838226/trials.csv) · [Scaling faible](weak/) · [Contrôles distribués](838226/mpi-tests.txt) · [Historique des réservations](jobs.tsv).
