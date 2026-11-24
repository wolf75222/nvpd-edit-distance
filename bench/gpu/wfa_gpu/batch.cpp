#include "chps/baselines.hpp"
#include "chps/cuda.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) try {
    if (argc != 4) throw std::invalid_argument("batch edlib|wfa2|nvpd-cuda input.seq scores.tsv");
    const std::string backend = argv[1];
    if (backend != "edlib" && backend != "wfa2" && backend != "nvpd-cuda")
        throw std::invalid_argument("unknown backend");
    std::ifstream input(argv[2]);
    std::ofstream output(argv[3]);
    if (!input || !output) throw std::runtime_error("cannot open batch files");
    chps::CudaWorkspace workspace;
    std::string query, target;
    while (std::getline(input, query)) {
        if (!std::getline(input, target) || !query.starts_with('>') || !target.starts_with('<'))
            throw std::invalid_argument("incomplete or malformed .seq pair");
        auto a = std::string_view(query).substr(1), b = std::string_view(target).substr(1);
        auto result = backend == "edlib" ? chps::edlib_distance(a,b)
                    : backend == "wfa2" ? chps::wfa_distance(a,b) : workspace.distance(a,b);
        output << result.distance << '\n';
    }
    output.close();
    if (input.bad() || !output) throw std::runtime_error("batch I/O failed");
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}
