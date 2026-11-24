# The pinned wavefront extension reads words at arbitrary byte offsets.
set(wfa_extend "${wfa2_SOURCE_DIR}/wavefront/wavefront_extend_kernels.c")
set(wfa_alloc "${wfa2_SOURCE_DIR}/system/mm_allocator.c")
set(patched_dir "${CMAKE_CURRENT_BINARY_DIR}/wfa2-portable")
file(READ "${wfa_extend}" source)
string(REPLACE "FORCE_INLINE wf_offset_t wavefront_extend_matches_kernel_blockwise("
  [[static inline uint64_t load_word(const char* p) {
  uint64_t word;
  memcpy(&word,p,sizeof(word));
  return word;
}

FORCE_INLINE wf_offset_t wavefront_extend_matches_kernel_blockwise(]] source "${source}")
foreach(name pattern text)
  string(REPLACE "uint64_t* ${name}_blocks = (uint64_t*)" "const char* ${name}_blocks = " source "${source}")
  string(REPLACE "++${name}_blocks;" "${name}_blocks += sizeof(uint64_t);" source "${source}")
endforeach()
string(REPLACE "*pattern_blocks ^ *text_blocks" "load_word(pattern_blocks) ^ load_word(text_blocks)" source "${source}")
file(GENERATE OUTPUT "${patched_dir}/wavefront_extend_kernels.c" CONTENT "${source}")

# Keep the cursor one past the live requests, including when none remain.
file(READ "${wfa_alloc}" source)
string(REPLACE "vector_get_mem(segment->requests,mm_allocator_request_t) + (num_requests-1);"
               "vector_get_mem(segment->requests,mm_allocator_request_t) + num_requests;" source "${source}")
string(REPLACE "while (num_requests>0 && MM_ALLOCATOR_REQUEST_IS_FREE(request))"
               "while (num_requests>0 && MM_ALLOCATOR_REQUEST_IS_FREE(request-1))" source "${source}")
string(REPLACE "segment->used = request->offset + request->size;"
               "segment->used = request[-1].offset + request[-1].size;" source "${source}")
file(GENERATE OUTPUT "${patched_dir}/mm_allocator.c" CONTENT "${source}")
list(REMOVE_ITEM wfa_sources "${wfa_extend}" "${wfa_alloc}")
list(APPEND wfa_sources "${patched_dir}/wavefront_extend_kernels.c" "${patched_dir}/mm_allocator.c")
