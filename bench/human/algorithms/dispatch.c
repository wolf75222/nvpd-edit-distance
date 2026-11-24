// Linked only in a separate, untimed GNU ld --wrap probe.
#include "utils/commons.h"
#include "wavefront/wfa.h"

static unsigned entries, bounds;
void __real_wavefront_bialign(wavefront_aligner_t *);
void __wrap_wavefront_bialign(wavefront_aligner_t *aligner) {
    ++entries;
    __real_wavefront_bialign(aligner);
}
void __real_wavefront_bialigner_set_sequences_bounds(wavefront_bialigner_t *, int, int, int, int);
void __wrap_wavefront_bialigner_set_sequences_bounds(wavefront_bialigner_t *aligner,
                                                     int a0, int a1, int b0, int b1) {
    ++bounds;
    __real_wavefront_bialigner_set_sequences_bounds(aligner, a0, a1, b0, b1);
}
__attribute__((destructor)) static void counts(void) {
    fprintf(stderr, "%u,%u\n", entries, bounds);
}
