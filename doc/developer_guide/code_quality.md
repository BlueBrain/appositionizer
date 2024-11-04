# Software Quality Measures {#code_quality}
This page assumes you're using a [modern cluster](@ref modern_cluster).

## Code Formatting {#code_formatting}
Appositionizer uses the HPC Coding Conventions wrapper for `clang-format` to format the
C++ source code. The formatting is asserted during CI.

The HPC Coding Conventions wrapper will check that the correct version of
formatting tools are installed. There's two convenient ways of using this,
either ensure that you've `pip install`ed the correct version, as stated in
`.bbp-project.yaml`; or use the following virtual env:
```
python -m venv .bbp-project-venv
source .bbp-project-venv
```
When using `.bbp-project-venv` the tool will install any required dependencies
automatically.

In principle you can now also simply call `clang-format -i`, e.g., in a
pre-write hook. It's convenient to use the HPC Coding Conventions wrapper:
```
deps/hpc-coding-conventions/bin/format [-n]
```
This will run `clang-format` on the correct source files under version control.
The flag `-n` (short for `--dry-run`) will only assert the formatting but apply
it format.

## Continuous Integration {#continuous_integration}
For continuous integration Appositionizer relies heavily on the [BBP Gitlab
Pipelines][gitlab_pipelines]. The Gitlab Pipelines encapsulate the logic
required to run jobs on BB5 during CI. The project is well-documented.

[gitlab_pipelines]: https://bbpgitlab.epfl.ch/hpc/gitlab-pipelines

### Unit Tests
Appositionizer uses [Catch2][catch2] and keeps all unit test in `tests`. These unit-tests
are run on every commit pushed to the repository. The runtime of unit-tests is to be kept
short/instant.

[catch2]: https://github.com/catchorg/Catch2/tree/v2.x

### Integration Tests
Any tests that can't be performed quickly as part of the unit-tests, will be
considered integration tests. An important part of the integration tests consist
of checking that the output of Appositionizer hasn't changed compared to a ground truth
for selected, small circuits. The corresponding scripts are found in
`.ci/test_circuit_*.sh`.  Further, integration tests are found in the same
folder.

Adding new integration tests is done by creating a script
`.ci/test_${TEST_NAME}.sh` which performs the desired checks. Then register
this test using `${TEST_NAME}` as a job name in `.gitlab-ci.hpc.yml`.

## Sanitizers & Clang-Tidy
The sanitizers ASAN and UBSAN; and static code analysis using `clang-tidy` have
been setup, see [Build Instructions](@ref build_variants) for details on how to
enable them.

## Code Documentation {#code_documentation}
The split between the Developer Guide and the User Guide is that everything
that's user facing. When faced with the decision of where to put a piece of
information, ask yourself: Can a user, i.e., someone that doesn't edit the
code, observe behaviour that requires this piece of information to explain? If
yes, place it in the user guide. If not, it belongs in the Developer Guide.

The User Guide uses the Sphinx documentation system, because this is the
documentation system that is well supported by BBP Software Catalog and
supports BBP themed styling. Since the Developer Guide simply consists of
static HTML pages one can link from the Sphinx documentation into parts of the
Developer Guide with ease and reasonable accuracy.

The Developer Guide uses Doxygen. While there are tools for parsing Doxygen
output to embed it directly into sphinx, in our experience these tools are
exceptionally brittle and hard to configure (2022-05-19). The issue is
sufficiently bad that we prefer the limitations of dealing with two tools. The
vital functionality that Doxygen provides is:
* complete and automatic listing of all functions, classes, members and
  namespaces,
* automatically generated list of classes, etc.
* the ability to embed the source code into the documentation,
* automatic generation of a link to the API documentation for `SomeClass`,
* well-structured API documentation,
* parses Doxygen-style docstrings.

For now, Doxygen is the only feasible choice for documenting C++ codebases. We
use [doxygen-awesome] to generate a modern, unthemed look. In order to have
automatically generated links from regular documentation into the API
documentation, we need to write the developer guide inside the Doxygen
ecosystem and use markdown since this is the format that Doxygen supports.

[doxygen-awesome]: https://jothepro.github.io/doxygen-awesome-css/

### Links to the Developer Guide
It seems that we need to manually maintain the links from the User Guide into
the Developer Guide. In order to keep duplication low the links have been 
collected in a file `doc/source/developer_guide_links.rst`. Please use it when
adding new links.

### Publishing Documentation
You can build the documentation locally using
```
rm -r doc/build; tox -e docs
```
Which will first delete the destination, since otherwise the Doxygen
documentation wont be updated. We use `tox` because it allows us to use the NSE
pipelines to publish the documentation to the BBP Software Catalog. This is
done during CI/CD automatically for every release.
