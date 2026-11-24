# NvPD · Distance d’édition parallèle

Distance de Levenshtein exacte en C++20, OpenMP et CUDA. Le programme compare NvPD à la programmation dynamique séquentielle, aux anti-diagonales OpenMP, à edlib et à WFA2-lib. [Article de Sadiq et coll.](https://doi.org/10.1007/s10586-019-02962-w)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python3 scripts/synthetic.py --lengths 5000 --rates 0.05 --pairs 1 --out build/exemple
build/editdist pair --algorithm nvpd-omp --threads 4 --input build/exemple/pairs/*.fa
```

CMake 3.22+ requiert un compilateur C++20 avec OpenMP. Pour un GPU Hopper, ajouter `-DCHPS_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90` à la configuration. Le mode `pair --trace` reconstruit un alignement jusqu’à 20 000 bases avec `seq`, `nvpd-omp` ou `nvpd-cuda`.

[Notebook exécuté](notebooks/profil_nvpd.ipynb) · [Rapport PDF](report/rapport.pdf) · [Données et arbres](data/README.md) · [Benchmarks ROMEO](bench/README.md)
