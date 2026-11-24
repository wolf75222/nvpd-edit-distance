# The pinned wavefront extension reads words at arbitrary byte offsets.
set(wfa_alloc "${wfa2_SOURCE_DIR}/system/mm_allocator.c")
set(patched_dir "${CMAKE_CURRENT_BINARY_DIR}/wfa2-portable")
set(word_load [[FORCE_INLINE uint64_t wavefront_load_unaligned_word(const char* p) {
  uint64_t word;
  memcpy(&word,p,sizeof(word));
  return word;
}

]])
foreach(unit wavefront_extend_kernels wavefront_extend_kernels_avx)
  set(input "${wfa2_SOURCE_DIR}/wavefront/${unit}.c")
  file(READ "${input}" source)
  # AVX uses this scalar loop for residual diagonals.
  if(unit STREQUAL "wavefront_extend_kernels")
    set(entry "FORCE_INLINE wf_offset_t wavefront_extend_matches_kernel_blockwise(")
  else()
    set(entry "FORCE_INLINE wf_offset_t wavefront_extend_matches_packed_kernel(")
  endif()
  string(REPLACE "${entry}" "${word_load}${entry}" source "${source}")
  foreach(name pattern text)
    string(REPLACE "uint64_t* ${name}_blocks = (uint64_t*)" "const char* ${name}_blocks = " source "${source}")
    string(REPLACE "++${name}_blocks;" "${name}_blocks += sizeof(uint64_t);" source "${source}")
  endforeach()
  string(REPLACE "*pattern_blocks ^ *text_blocks" "wavefront_load_unaligned_word(pattern_blocks) ^ wavefront_load_unaligned_word(text_blocks)" source "${source}")
  file(GENERATE OUTPUT "${patched_dir}/${unit}.c" CONTENT "${source}")
  list(REMOVE_ITEM wfa_sources "${input}")
  list(APPEND wfa_sources "${patched_dir}/${unit}.c")
endforeach()

# Traceback writes matches into a byte buffer, starting one past its last match.
set(wfa_trace "${wfa2_SOURCE_DIR}/wavefront/wavefront_backtrace.c")
file(READ "${wfa_trace}" source)
string(REPLACE "char* operations = cigar->operations + cigar->begin_offset;"
               "char* operations = cigar->operations + (cigar->begin_offset + 1);" source "${source}")
string(REPLACE "*((uint64_t*)(operations+1)) = matches_lut;"
               "memcpy(operations,&matches_lut,sizeof(matches_lut));" source "${source}")
string(REPLACE "*operations = 'M';\n    --operations;" "*--operations = 'M';" source "${source}")
file(GENERATE OUTPUT "${patched_dir}/wavefront_backtrace.c" CONTENT "${source}")

# Unpack complete 8-byte matches and then the exact tail, even at sequence end.
set(wfa_pcigar "${wfa2_SOURCE_DIR}/wavefront/wavefront_pcigar.c")
file(READ "${wfa_pcigar}" source)
string(FIND "${source}" "  // Fetch pattern/text blocks" begin)
string(FIND "${source}" "  // Return total matches" end)
if(begin LESS 0 OR end LESS_EQUAL begin)
  message(FATAL_ERROR "Cannot locate the pinned WFA2 packed CIGAR extension")
endif()
string(SUBSTRING "${source}" 0 ${begin} prefix)
string(SUBSTRING "${source}" ${end} -1 suffix)
set(pcigar_loop [[  const int available = MIN(pattern_length - v,text_length - h);
  while (available - num_matches >= 8 &&
         memcmp(pattern+v+num_matches,text+h+num_matches,8) == 0) {
    num_matches += 8;
  }
  while (num_matches < available && pattern[v+num_matches] == text[h+num_matches]) {
    ++num_matches;
  }
  memset(cigar_buffer,'M',num_matches);
]])
file(GENERATE OUTPUT "${patched_dir}/wavefront_pcigar.c" CONTENT "${prefix}${pcigar_loop}${suffix}")

# Keep the cursor one past the live requests, including when none remain.
file(READ "${wfa_alloc}" source)
string(REPLACE "vector_get_mem(segment->requests,mm_allocator_request_t) + (num_requests-1);"
               "vector_get_mem(segment->requests,mm_allocator_request_t) + num_requests;" source "${source}")
string(REPLACE "while (num_requests>0 && MM_ALLOCATOR_REQUEST_IS_FREE(request))"
               "while (num_requests>0 && MM_ALLOCATOR_REQUEST_IS_FREE(request-1))" source "${source}")
string(REPLACE "segment->used = request->offset + request->size;"
               "segment->used = request[-1].offset + request[-1].size;" source "${source}")
file(GENERATE OUTPUT "${patched_dir}/mm_allocator.c" CONTENT "${source}")
list(REMOVE_ITEM wfa_sources "${wfa_alloc}" "${wfa_trace}" "${wfa_pcigar}")
list(APPEND wfa_sources "${patched_dir}/mm_allocator.c" "${patched_dir}/wavefront_backtrace.c"
                       "${patched_dir}/wavefront_pcigar.c")
