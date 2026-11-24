# Apply the two boundary checks to the pinned source without modifying the
# extracted vendor tree. The benchmark links the original translation unit.
file(READ "${original_banded}" source)
string(REPLACE [[        const uint64_t mask = 1ULL << (effective_v % BPM_W64_LENGTH);
        const uint64_t mask_r = 1ULL << (effective_v_r % BPM_W64_LENGTH);]]
[[        const int64_t stored_bits = (int64_t)num_words64 * BPM_W64_LENGTH;
        const bool has_v = effective_v >= 0 && effective_v < stored_bits;
        const bool has_v_r = effective_v_r >= 0 && effective_v_r < stored_bits;]] source "${source}")
string(REPLACE [[if (Pv[bdp_idx_r] & mask_r)]]
               [[if (has_v_r && (Pv[bdp_idx_r] & (1ULL << (effective_v_r % BPM_W64_LENGTH))))]]
               source "${source}")
string(REPLACE [[else if (Mv[(bdp_idx)] & mask)]]
               [[else if (has_v && (Mv[bdp_idx] & (1ULL << (effective_v % BPM_W64_LENGTH))))]]
               source "${source}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/quicked_banded.c" "${source}")
