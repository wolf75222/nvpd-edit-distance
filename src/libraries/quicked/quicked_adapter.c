#include "chps/quicked.h"
#include <quicked.h>
#include <edlib.h>
#include <limits.h>
#include <string.h>

static bool valid_input(const char *bases, int length) {
    return length >= 0 && length <= INT_MAX / 100 && (length == 0 || bases != NULL);
}

// Upstream uses M for equality, X for substitution, I for consuming b,
// D for consuming a. Replay the path before accepting its proposed cost.
static bool valid_path(const char *p, const char *a, int m, const char *b, int n, int score) {
    int i = 0, j = 0, edits = 0;
    if (!p) return false;
    while (*p) {
        int k = 0;
        while (*p >= '0' && *p <= '9') {
            const int digit = *p++ - '0';
            if (k > (INT_MAX - digit) / 10) return false;
            k = k * 10 + digit;
        }
        const char op = *p;
        if (!op || !strchr("MXID", op) || !k) return false;
        ++p;
        if ((op != 'I' && k > m - i) || (op != 'D' && k > n - j)) return false;
        if (op == 'M' && memcmp(a + i, b + j, k)) return false;
        if (op == 'X')
            for (int h = 0; h < k; ++h)
                if (a[i + h] == b[j + h]) return false;
        if (op != 'I') i += k;
        if (op != 'D') j += k;
        if (op != 'M') edits += k;
    }
    return i == m && j == n && edits == score;
}

static int run(const char *a, int m, const char *b, int n, quicked_params_t params, bool replay) {
    if (!m || !n) return m + n;
    quicked_aligner_t aligner;
    if (quicked_check_error(quicked_new(&aligner, &params))) return -1;
    const quicked_status_t status = quicked_align(&aligner, a, m, b, n);
    const int score = quicked_check_error(status) ? -1 : aligner.score;
    const bool ok = !replay || (score >= 0 && valid_path(aligner.cigar, a, m, b, n, score));
    quicked_free(&aligner);
    return ok ? score : -1;
}

int chps_quicked_default(const char *a, int m, const char *b, int n) {
    if (!valid_input(a, m) || !valid_input(b, n)) return -1;
    return run(a, m, b, n, quicked_default_params(), false);
}

static int exact_score(const char *a, int m, const char *b, int n) {
    const EdlibAlignConfig cfg = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    EdlibAlignResult ref = edlibAlign(a, m, b, n, cfg);
    const int score = ref.status == EDLIB_STATUS_OK ? ref.editDistance : -1;
    edlibFreeAlignResult(ref);
    return score;
}

chps_quicked_result chps_quicked_certified(const char *a, int m, const char *b, int n) {
    if (!valid_input(a, m) || !valid_input(b, n))
        return (chps_quicked_result){-1, -1, 0, 0, "invalid-input"};
    const int first = run(a, m, b, n, quicked_default_params(), true);
    chps_quicked_result result = {first, first, 0, 0, "quicked-default"};
    if (m == 0 || n == 0) return result;
    const int exact = exact_score(a, m, b, n);
    if (exact < 0) {
        result.distance = -1;
        return result;
    }
    if (first == exact) return result;

    // The default mode chooses its own final cutoff. Select HIRSCHBERG
    // explicitly so the requested, successively larger cutoffs reach traceback.
    quicked_params_t params = quicked_default_params();
    params.algo = HIRSCHBERG;
    const int length = m > n ? m : n;
    int percent = (int)((100LL * exact + length - 1) / length);
    if (percent < 1) percent = 1;
    result.source = "quicked-widened";
    for (;;) {
        params.bandwidth = percent;
        ++result.retries;
        result.band_percent = percent;
        result.distance = run(a, m, b, n, params, true);
        if (result.distance == exact) return result;
        if (percent == 100) break;
        percent = percent < 50 ? 2 * percent : 100;
    }
    // Preserve the exact distance when the upstream path still fails.
    result.source = "edlib-fallback";
    result.distance = exact;
    return result;
}
