# NvPD · Distance d’édition parallèle

Distance de Levenshtein globale en C++20. Le projet compare NvPD sur OpenMP, CUDA, Kokkos et MPI aux anti-diagonales, à la programmation dynamique séquentielle, au calcul en bande, à edlib, à WFA2-lib et à QuickEd. [Article de Sadiq et coll.](https://doi.org/10.1007/s10586-019-02962-w)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python3 scripts/synthetic.py --lengths 5000 --rates 0.05 --pairs 1 --out build/exemple
build/editdist pair --algorithm nvpd-omp --threads 4 --input build/exemple/pairs/*.fa
```

CMake 3.22+ et OpenMP sont requis. Les options `CHPS_CUDA`, `CHPS_KOKKOS`, `CHPS_MPI` et `CHPS_QUICKED` activent les cibles supplémentaires. Pour un GH200, choisir `CMAKE_CUDA_ARCHITECTURES=90`. La commande `pair --algorithm hirschberg --trace` reconstruit le CIGAR avec deux lignes de scores. En mode `matrix`, le programme écrit un TSV directement lisible par `scripts/tree.py`.

Pour reproduire les figures, installer [les dépendances Python](notebooks/requirements.txt) et exécuter le notebook dans l’ordre. La charte des figures se trouve dans [ct_plot.py](python/ct_plot.py). Pour compiler le PDF, lancer `python3 report/phylograms.py` à la racine, puis `tectonic rapport.tex` et `tectonic rapport-court.tex` depuis `report/`.

[Notebook exécuté](notebooks/profil_nvpd.ipynb) · [Rapport PDF](report/rapport.pdf) · [Synthèse courte](report/rapport-court.pdf) · [Données et arbres](data/README.md) · [Benchmarks ROMEO](bench/README.md)
