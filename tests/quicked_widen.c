// Deliberately start outside the optimal path to exercise production recovery.
// The real CLI uses the unchanged upstream adaptive initial parameters.
#include <quicked.h>
#include <stdio.h>
static quicked_params_t narrow_start(void) {
    quicked_params_t p = quicked_default_params();
    p.algo = HIRSCHBERG;
    p.bandwidth = 1;
    return p;
}
#define quicked_default_params narrow_start
#include "../src/libraries/quicked/quicked_adapter.c"

int main(void) {
    enum { length = 1152 };
    char a[length], b[length];
    for (int i = 0; i < length; ++i) {
        a[i] = i % 384 < 192 ? 'A' : 'C';
        b[i] = i % 384 < 192 ? 'C' : 'A';
    }
    const chps_quicked_result q = chps_quicked_certified(a, length, b, length);
    printf("initial=%d corrected=%d retries=%d percent=%d source=%s\n",
           q.first_score, q.distance, q.retries, q.band_percent, q.source);
    return q.distance != 384 || q.first_score == q.distance || q.retries == 0 ||
           strcmp(q.source, "quicked-widened") != 0;
}
