#include "utils/commons.h"
#include "wavefront/wfa.h"

int chps_wfa_distance(const char *a, int m, const char *b, int n) {
    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    attr.distance_metric = edit;
    attr.alignment_scope = compute_score;
    attr.alignment_form.span = alignment_end2end;
    attr.heuristic.strategy = wf_heuristic_none;
    wavefront_aligner_t *aligner = wavefront_aligner_new(&attr);
    if (aligner == NULL) return -1;
    int status = wavefront_align(aligner, a, m, b, n);
    int distance = status == WF_STATUS_ALG_COMPLETED ? aligner->align_status.score : -1;
    wavefront_aligner_delete(aligner);
    return distance;
}
