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

// The caller supplies ACGT strings, lengths between zero and INT_MAX / 100.
int chps_quicked_default(const char *a, int m, const char *b, int n);
chps_quicked_result chps_quicked_certified(const char *a, int m, const char *b, int n);

#ifdef __cplusplus
}
#endif
