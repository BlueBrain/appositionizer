cmake_minimum_required(VERSION 3.19.0 FATAL_ERROR)

project(Appositionizer VERSION 1.0.0)

add_subdirectory(deps/hpc-coding-conventions/cpp)

set(${PROJECT_NAME}_DESCRIPTION
    "Apposition detection software, part of the Simulation tools for the Blue Brain Project")

execute_process(
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tools/update_cmake_files.py" "${CMAKE_CURRENT_SOURCE_DIR}"
  COMMAND_ECHO STDOUT
  RESULT_VARIABLE Appositioner_SELF_MODIFYING_RETURN_CODE
)
if(NOT ${Appositioner_SELF_MODIFYING_RETURN_CODE} EQUAL 0)
  message(AUTHOR_WARNING
    "The file-listing script failed. Therefore, the list of files"
    " might be outdated. If you've added or removed files, you can"
    " run\n"
    "    tools/update_cmake_files.py\n"
    "manually to update the required `CMakeLists.txt`."
  )
endif()

set_property(DIRECTORY APPEND PROPERTY
 CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tools/update_cmake.sh"
)

option(EXTLIB_FROM_SUBMODULES "Use Git submodules for header-only dependencies" OFF)

option(ENABLE_DEBUG "Enable verbose debugging" OFF)
option(ENABLE_CALIPER "Enables profiling with Caliper" OFF)
option(ENABLE_BENCHMARKS "Enables building benchmarks" OFF)
option(ENABLE_TESTS "Enables building and running tests" ON)
option(ENABLE_ASAN "Enables Clang's AddressSanitizer" OFF)
option(ENABLE_UBSAN "Enables Clang's UndefinedBehaviourSanitizer" OFF)

list(APPEND CMAKE_MODULE_PATH
            ${PROJECT_SOURCE_DIR}/CMake
            ${PROJECT_SOURCE_DIR}/CMake/config
)

set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH TRUE)

include("${CMAKE_CURRENT_SOURCE_DIR}/CMake/warnings.cmake")
include(ReleaseDebugAutoFlags)
include(FindPkgConfig)

find_package(MorphIO REQUIRED)
find_package(MPI REQUIRED)
find_package(sonata REQUIRED)
find_package(yaml-cpp REQUIRED)
if(EXTLIB_FROM_SUBMODULES)
  add_subdirectory(deps/eigen)
  add_subdirectory(deps/fmt)
  add_subdirectory(deps/nlohmann-json)
  add_subdirectory(deps/range-v3)
else()
  find_package(Eigen3 REQUIRED NO_MODULE)
  find_package(fmt REQUIRED)
  find_package(nlohmann_json REQUIRED)
  find_package(Range-v3 REQUIRED)
endif()
find_package(TBB REQUIRED)
find_package(ZLIB REQUIRED)
find_package(HighFive REQUIRED)

add_subdirectory(deps/cli11)

find_path(RANDOM123_INCLUDE_PATH NAMES Random123/threefry.h)
add_library(Random123 INTERFACE IMPORTED)
set_target_properties(Random123 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${RANDOM123_INCLUDE_PATH}
)

if(ENABLE_CALIPER)
  find_package(caliper REQUIRED)
  add_definitions(-DAPPOSITIONIZER_CALIPER)
endif()


if(ENABLE_DEBUG)
  add_definitions(-DDEBUG_OUT)
endif()

add_library(${PROJECT_NAME}_deps INTERFACE)
target_link_libraries(${PROJECT_NAME}_deps
  INTERFACE
  CLI11::CLI11
  yaml-cpp
  MorphIO::morphio
  sonata::sonata_shared
  Eigen3::Eigen
  fmt::fmt
  Random123
  range-v3
  TBB::tbb
  nlohmann_json::nlohmann_json
  MPI::MPI_C
  HighFive
)

if(ENABLE_CALIPER)
  target_link_libraries(${PROJECT_NAME}_deps INTERFACE caliper)
endif()


# Configure files -------------------------------------------------------------
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in
               ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/appositionizer/defines.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/appositionizer/defines.h)


# Appositioner properties ---------------------------------------------------------------
function(target_cxx_standard target)
  set_target_properties(${target}
    PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
  )
endfunction()

function (target_sanitizers target)
  if((ENABLE_ASAN OR ENABLE_UBSAN) AND NOT CMAKE_CXX_COMPILER_IS_CLANG)
    message(FATAL_ERROR "Sanitizers require that the compiler is clang.")
  endif()

  if(ENABLE_ASAN)
    target_compile_options(${target} PRIVATE -fsanitize=address)
    target_link_options(${target} PRIVATE -fsanitize=address)
  endif()

  if(ENABLE_UBSAN)
    target_compile_options(${target} PRIVATE -fsanitize=undefined)
    target_link_options(${target} PRIVATE -fsanitize=undefined)
  endif()
endfunction()


# Appositioner library ------------------------------------------------------------------
add_library(appositionizer_lib ${COMPILE_LIBRARY_TYPE})
target_cxx_standard(appositionizer_lib)
target_sanitizers(appositionizer_lib)
set_target_properties(appositionizer_lib PROPERTIES OUTPUT_NAME "appositionizer")
target_include_directories(appositionizer_lib PUBLIC
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/appositionizer>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/appositionizer>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(appositionizer_lib PUBLIC ${PROJECT_NAME}_deps ${PROJECT_NAME}_warnings)
set_target_properties(appositionizer_lib PROPERTIES
  VERSION   ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  CLEAN_DIRECT_OUTPUT 1
)


# Appositioner executable ---------------------------------------------------------------
add_executable(appositionizer_bin
  "${CMAKE_CURRENT_SOURCE_DIR}/appositionizer/appositionizer_main.cxx"
)
target_link_libraries(appositionizer_bin appositionizer_lib)
target_cxx_standard(appositionizer_bin)
target_sanitizers(appositionizer_bin)
set_target_properties(appositionizer_bin
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/appositionizer"
  OUTPUT_NAME "appositionizer"
)


# MM executable
add_executable(morphoi_merge_bin
  "${CMAKE_CURRENT_SOURCE_DIR}/appositionizer/morphoi_merge_main.cxx"
)
target_link_libraries(morphoi_merge_bin
  PRIVATE
  MPI::MPI_C
  CLI11::CLI11
  HighFive
  sonata::sonata_shared
  MorphIO::morphio)
target_cxx_standard(morphoi_merge_bin)
target_sanitizers(morphoi_merge_bin)
set_target_properties(morphoi_merge_bin
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/appositionizer"
  OUTPUT_NAME "morphoi-merge"
)


# Appositioner tests & benchmarks ------------------------------------------------------
if(ENABLE_TESTS)
  # unit tests ---------------------------------------------------------------
  find_package(Catch2 REQUIRED)

  enable_testing()
  add_executable(serial_tests
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/serial_tests_main.cpp"
  )
  target_link_libraries(serial_tests appositionizer_lib Catch2::Catch2)
  target_cxx_standard(serial_tests)
  target_sanitizers(serial_tests)

  set_target_properties(serial_tests
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tests"
  )

  include(Catch)
  catch_discover_tests(serial_tests)


  # benchmarks ----------------------------------------------------------------
  if(ENABLE_BENCHMARKS)
    find_package(benchmark REQUIRED)

    add_executable(overlap
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/benchmarks/overlap_main.cpp"
    )
    target_cxx_standard(overlap)
    target_link_libraries(overlap PRIVATE
      # -lpthread
      benchmark::benchmark
      appositionizer_lib
    )
  endif()


  add_subdirectory(tests)
endif()

# Installation ----------------------------------------------------------------
install(TARGETS appositionizer_lib appositionizer_bin morphoi_merge_bin)

# Recurse into subdirectories.
add_subdirectory(appositionizer)
