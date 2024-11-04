# This is an automatically generated file.
#
# Please use the script
#   tools/update_cmake_files.py
# to regenerate the list of source files.

target_sources(appositionizer_lib
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/better_mpi.cxx
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cell_slice.cxx
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/cell_slicer.cxx
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/distributed_memory_sorter.cxx
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/memory_tools.cxx
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/partial_collision.cxx
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}/recipe.cxx
)

add_subdirectory(data)
add_subdirectory(io)
add_subdirectory(process)
