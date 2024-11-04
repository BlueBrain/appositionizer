# This is an automatically generated file.
#
# Please use the script
#   tools/update_cmake_files.py
# to regenerate the list of source files.

target_sources(serial_tests
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/apposition_space.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/arange.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/box.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cell_group.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cell_slicer.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/file_reader.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/index_magic.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/morphology.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/recipe.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/segment_overlap.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/segment_spine.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/segment_storage.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/serialization.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/sorting.cpp
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/spatial_index.cpp
)

add_subdirectory(benchmarks)
