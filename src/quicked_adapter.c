#include <quicked.h>

// QuickEd's default exact pipeline constructs a CIGAR internally.
int chps_quicked(const char *a, int m, const char *b, int n) {
    if (!m || !n) return m + n;
    quicked_params_t params = quicked_default_params();
    quicked_aligner_t aligner;
    if (quicked_check_error(quicked_new(&aligner, &params))) return -1;
    const quicked_status_t status = quicked_align(&aligner, a, m, b, n);
    const int score = quicked_check_error(status) ? -1 : aligner.score;
    quicked_free(&aligner);
    return score;
}
