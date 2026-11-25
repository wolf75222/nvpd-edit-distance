if(PROJECT_NAME STREQUAL "chps_edit_distance" AND NOT TARGET grace_placement)
  add_executable(grace_placement "${CMAKE_CURRENT_LIST_DIR}/placement.cpp"
                 "${CMAKE_SOURCE_DIR}/src/io/fasta.cpp")
  set_target_properties(grace_placement PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  target_link_libraries(grace_placement PRIVATE chps_warnings distances)
endif()
