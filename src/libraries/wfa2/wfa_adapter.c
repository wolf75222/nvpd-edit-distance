#include "utils/commons.h"
#include "wavefront/wfa.h"

int chps_wfa_align(const char *a, int m, const char *b, int n,
                   bool bidirectional, char **path) {
    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    attr.distance_metric = edit;
    attr.alignment_scope = path ? compute_alignment : compute_score;
    attr.alignment_form.span = alignment_end2end;
    attr.heuristic.strategy = wf_heuristic_none;
    attr.memory_mode = bidirectional ? wavefront_memory_ultralow : wavefront_memory_high;
    wavefront_aligner_t *aligner = wavefront_aligner_new(&attr);
    if (aligner == NULL) return -1;
    int distance = -1;
    if (wavefront_align(aligner, a, m, b, n) == WF_STATUS_ALG_COMPLETED) {
        // Both dispatchers store the global unit cost on the CIGAR, even in score-only mode.
        const cigar_t *cigar = aligner->cigar;
        distance = cigar->score;
        if (path) {
            const size_t steps = (size_t)(cigar->end_offset - cigar->begin_offset);
            *path = malloc(2 * steps + 1); // RLE uses at most two characters per operation.
            if (*path == NULL) distance = -1;
            else cigar_sprint(*path, cigar, true);
        }
    }
    wavefront_aligner_delete(aligner);
    return distance;
}
