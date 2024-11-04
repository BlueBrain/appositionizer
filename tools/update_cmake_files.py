#! /usr/bin/env python3

# This file is a modified version of
# https://github.com/1uc/ZisaMemory/blob/main/bin/update_cmake.py
#
# Copyright (c) 2021 ETH Zurich, Luc Grosheintz-Laval
# Copyright (c) 2022 Blue Brain Project, EPFL
#
# For modifications, all rights are reserved. This file should not
# be distributed without contacting the Blue Brain Project.


import sys
import glob
import os
import errno


def find_files(folder, suffixes):
    files = sum(
        (glob.glob(os.path.join(folder, "*{}".format(s))) for s in suffixes), []
    )
    return list(sorted(files))


def find_source_files(folder):
    suffixes = [".c", ".C", ".cpp", ".cxx", ".c++"]
    return find_files(folder, suffixes)


def find_subdirectories(folder):
    dirs = sorted(glob.glob(os.path.join(folder, "*/")))
    return [d for d in dirs if "CMake" not in d]


def target_sources(target, sources):
    ret = ""
    line_pattern = "  PRIVATE ${{CMAKE_CURRENT_LIST_DIR}}/{:s}\n"

    if sources:
        ret += "".join(
            [
                "target_sources(" + target + "\n",
                "".join(line_pattern.format(os.path.basename(s)) for s in sources),
                ")\n\n",
            ]
        )

    return ret


def add_subdirectory(subfolder):
    line_pattern = "add_subdirectory({:s})\n"
    return line_pattern.format(os.path.basename(subfolder.rstrip("/")))

def self_documenting_comment():
    # It's important that this value doesn't change of unimportant reasons,
    # since that would lead to unnecessary git commits. An example of this
    # problem is
    #   script_name = os.path.relname(__file__)
    # Which isn't desireable, since we have no control over the CWD.

    script_name = "tools/update_cmake_files.py"
    return f"""
# This is an automatically generated file.
#
# Please use the script
#   {script_name}
# to regenerate the list of source files.

""".lstrip()


def remove_file(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise


def append_to_file(filename, text):
    with open(filename, "a") as f:
        f.write(text)


def ignore_directory(directory):
    return "tests/data" in directory


def recurse(base_directory, targets):
    cmake_file = os.path.join(base_directory, "CMakeLists.txt")
    remove_file(cmake_file)
    append_to_file(cmake_file, self_documenting_comment())

    source_files = find_source_files(base_directory)

    for dependency, target in targets.items():
        filtered_sources = list(filter(select_for(dependency), source_files))

        if dependency == "generic":
            append_to_file(cmake_file, target_sources(target, filtered_sources))

        elif dependency == "benchmark":
            append_to_file(cmake_file, "if(benchmark_FOUND)\n\n")
            append_to_file(cmake_file, target_sources(target, filtered_sources))
            append_to_file(cmake_file, "endif()\n")

        else:
            raise Exception("Unknown dependency. [{}]".format(dependency))

    for d in find_subdirectories(base_directory):
        if ignore_directory(d):
            continue

        recurse(d, targets)
        append_to_file(cmake_file, add_subdirectory(d))


def is_executable(path):
    # Anything matching
    #     *_main.{cpp,cxx,...}
    # contains a `main()` by convention.
    return path.rsplit(".", 1)[0].endswith("_main")


def is_benchmark_file(path):
    return "tests/benchmark/" in path


def is_generic_file(path):
    return not any(f(path) for f in [is_benchmark_file, is_executable])


def select_for(dependency):
    if dependency == "generic":
        return is_generic_file

    elif dependency == "benchmark":
        return is_benchmark_file

    else:
        raise Exception(f"Unknown dependency. [{dependency}]")


if __name__ == "__main__":
    appo_dir = "." if len(sys.argv) != 2 else sys.argv[1]

    recurse(os.path.join(appo_dir, "appositionizer"), {"generic": "appositionizer_lib"})
    recurse(os.path.join(appo_dir, "tests"), {"generic": "serial_tests"})
