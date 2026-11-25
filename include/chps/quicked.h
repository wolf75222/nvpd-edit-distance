#pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int distance;
    int first_score;
    int retries;
    int band_percent;
    const char *source;
} chps_quicked_result;

// Supply ACGT strings. A null pointer denotes an empty sequence only.
// Negative lengths and lengths above INT_MAX / 100 return a distance of -1.
int chps_quicked_default(const char *a, int m, const char *b, int n);
chps_quicked_result chps_quicked_certified(const char *a, int m, const char *b, int n);

#ifdef __cplusplus
}
#endif
