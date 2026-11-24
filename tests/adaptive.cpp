#include "chps/adaptive.hpp"
#include "chps/baselines.hpp"
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <omp.h>

static void check(const std::string& a, const std::string& b) {
    const int expected = chps::edlib_distance(a,b).distance;
    for (bool reverse : {false, true}) {
        chps::AdaptiveWork work;
        const int got = chps::nvpd_adaptive_score(reverse ? b : a, reverse ? a : b, work);
        if (got != expected) throw std::runtime_error("adaptive mismatch on " + a + " / " + b);
    }
}

int main() {
    omp_set_num_threads(1);
    std::vector<std::string> words{""};
    for (int n=1; n<=6; ++n)
        for (int mask=0; mask<(1<<n); ++mask) {
            std::string word(n,'A');
            for (int j=0; j<n; ++j) if (mask & (1<<j)) word[j]='N';
            words.push_back(word);
        }
    for (const auto& a : words) for (const auto& b : words) check(a,b);
    std::mt19937 random(260127);
    auto sequence = [&](int n) {
        std::string x(n, 'A');
        for (char& c : x) c = "ACGTNRYS"[random()%8];
        return x;
    };
    for (int i=0;i<400;++i) check(sequence(1+random()%550), sequence(1+random()%550));
    std::string core=sequence(1000);
    check("NN"+core+"R", core+"RR"); // translated corridor and absent terminal symbols
    check("N"+std::string(1000,'A')+"C", "T"+std::string(1000,'A')+"G");
    check(std::string(20000,'A'), std::string(20000,'A'));
    check("A",std::string(50000,'N'));
    std::cout << "Adaptive NvPD agrees with edlib in both orientations.\n";
}
