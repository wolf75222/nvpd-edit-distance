// Search deterministic A/C blocks for a counterexample on the measured ISA.
#include "chps/quicked.h"
#include <edlib.h>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    bool found[2]{};
    for (int k : {64, 128, 192, 256, 32, 96, 160})
        for (int width : {k, k + 1, k - 1, 2 * k + 1})
            for (int n : {4*k, 8*k, 12*k+1, 16*k-1, 20*k})
                for (int tail : {0, 1, 63, 64, 2*k})
                    for (int offset : {0, 1, k/2, k-1}) {
                        std::string a, b;
                        for (int i = 0; i < n; ++i) {
                            a += (i+offset) % (width+k) < width ? 'A' : 'C';
                            b += (i+offset) % (width+k) < k ? 'C' : 'A';
                        }
                        a.append(tail, 'G');
                        b.append(tail, 'G');
                        const int len = a.size();
                        auto ref = edlibAlign(a.data(),len,b.data(),len,edlibDefaultAlignConfig());
                        if (ref.status != EDLIB_STATUS_OK || ref.editDistance < 0) return 1;
                        const int d = ref.editDistance;
                        edlibFreeAlignResult(ref);
                        const int got = chps_quicked_default(a.data(),len,b.data(),len);
                        const bool impossible = got < 0 || got > len;
                        if (got == d || found[impossible]) continue;
                        found[impossible] = true;
                        const std::string name = impossible ? "native-impossible" : "native-suboptimal";
                        std::ofstream out(std::string(argv[1]) + "/" + name + ".fa");
                        out << ">a\n" << a << "\n>b\n" << b << "\n";
                        if (!out) return 1;
                        std::cout << name << ',' << len << ',' << d << ',' << got << '\n';
                        if (found[0] && found[1]) return 0;
                    }
}
