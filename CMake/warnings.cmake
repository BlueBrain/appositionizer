include(CompilerFlagsHelpers)

add_library(${PROJECT_NAME}_warnings INTERFACE)


if(CMAKE_CXX_COMPILER_IS_GCC OR CMAKE_CXX_COMPILER_IS_CLANG)
  target_compile_options(${PROJECT_NAME}_warnings
                        INTERFACE
                        -Wall
                        -Wextra
                        -Wcast-align
                        # -Wconversion
                        # -Wdouble-promotion
                        -Wmisleading-indentation
                        -Wnon-virtual-dtor
                        -Wnull-dereference
                        -Wshadow
                        # -Wsign-conversion
                        -Wunused
                        -pedantic
                        )
endif()

if(CMAKE_CXX_COMPILER_IS_GCC)
  target_compile_options(${PROJECT_NAME}_warnings
                         INTERFACE
                         -Wduplicated-cond
                         -Wuseless-cast
                         )
endif()
