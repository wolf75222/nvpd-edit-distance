#pragma once
#include <string>
#include <vector>

namespace chps {
struct Record {
    std::string id;
    std::string bases;
};
std::vector<Record> read_fasta(const std::string& path);
}
