# NvPD · Distance d’édition parallèle

Distance de Levenshtein globale en C++20. Le projet compare NvPD sur OpenMP, CUDA, Kokkos et MPI aux anti-diagonales, à la programmation dynamique séquentielle, au calcul en bande, à edlib, à WFA2-lib et à QuickEd. [Article de Sadiq et coll.](https://doi.org/10.1007/s10586-019-02962-w)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python3 scripts/synthetic.py --lengths 5000 --rates 0.05 --pairs 1 --out build/exemple
build/editdist pair --algorithm nvpd-omp --threads 4 --input build/exemple/pairs/*.fa
```

CMake 3.22+ et OpenMP sont requis. Les options `CHPS_CUDA`, `CHPS_KOKKOS`, `CHPS_MPI` et `CHPS_QUICKED` activent les cibles supplémentaires. Pour un GH200, choisir `CMAKE_CUDA_ARCHITECTURES=90`. En mode `matrix`, le programme écrit un TSV directement lisible par `scripts/tree.py`.

[Notebook exécuté](notebooks/profil_nvpd.ipynb) · [Rapport PDF](report/rapport.pdf) · [Données et arbres](data/README.md) · [Benchmarks ROMEO](bench/README.md)
