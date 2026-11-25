// Run each loaded pair through a single CUDA workspace on every command.
#include "chps/cuda.hpp"
#include "chps/fasta.hpp"
#include <algorithm>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct Pair {
    std::vector<chps::Record> sequences;
    int expected;
};

int main(int argc, char** argv) {
    try {
        if (argc < 3 || argc % 2 != 1) return 2;
        std::vector<Pair> pairs;
        for (int i = 1; i < argc; i += 2) {
            pairs.push_back({chps::read_fasta(argv[i]), std::stoi(argv[i + 1])});
            if (pairs.back().sequences.size() != 2) throw std::runtime_error("Expected two FASTA records");
        }
        chps::CudaWorkspace workspace;
        std::cout << "ready" << std::endl;
        std::string command;
        while (std::getline(std::cin, command)) {
            if (command != "run") return 2;
            std::size_t peak = 0;
            for (const auto& pair : pairs) {
                const auto result = workspace.distance(pair.sequences[0].bases, pair.sequences[1].bases);
                if (result.distance != pair.expected) throw std::runtime_error("Score differs from Edlib");
                peak = std::max(peak, result.device_bytes);
            }
            std::cout << pairs.size() << ',' << peak << std::endl;
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
