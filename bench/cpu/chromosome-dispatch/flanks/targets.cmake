if(PROJECT_NAME STREQUAL "chps_edit_distance" AND NOT TARGET flank_trial)
  add_executable(flank_trial "${CMAKE_CURRENT_LIST_DIR}/trial.cpp" "${CMAKE_SOURCE_DIR}/src/io/fasta.cpp"
                 "${CHPS_SPLIT_BASELINE}")
  set_source_files_properties("${CHPS_SPLIT_BASELINE}" PROPERTIES COMPILE_DEFINITIONS
      "SplitComposition=OriginalSplitComposition;nvpd_split_profile=original_split_profile;nvpd_split=original_split;nvpd_sparse=original_sparse")
  set_target_properties(flank_trial PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
  target_include_directories(flank_trial PRIVATE "${CMAKE_SOURCE_DIR}/src/cpu")
  target_link_libraries(flank_trial PRIVATE distances chps_warnings)
endif()
