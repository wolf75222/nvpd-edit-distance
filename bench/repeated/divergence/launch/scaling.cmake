if(PROJECT_NAME STREQUAL "chps_edit_distance" AND NOT TARGET grace_scaling)
  add_executable(grace_scaling "${CMAKE_CURRENT_LIST_DIR}/scaling.cpp"
                 "${CMAKE_SOURCE_DIR}/src/io/fasta.cpp")
  set_target_properties(grace_scaling PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  target_link_libraries(grace_scaling PRIVATE chps_warnings distances)
endif()
