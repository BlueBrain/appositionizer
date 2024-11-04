# Build Instructions {#build_instructions}
Simply running TD on BB5 does not require building it. Instead simply load the
module `touchdetector`:

    module load touchdetector


## A modern cluster {#modern_cluster}
This documentation assumes you're using modern versions of the standard tools.
You can load them by
```
module load unstable spack git cmake python
```
(or similar).


## Using Spack
As a developer of TD we recommend that you use `spack` to compile TD. If you're
not familiar with `spack` you can find documentation here:
  * BlueBrain Spack [Documentation][bbp_spack_docs] for BB5-centric documentation.
  * Upstream Spack [Documentation][spack_docs] is very extensive.

Using Spack will take care of installing the numerous dependencies
automatically. One convenient way of doing this is by using environments:
```
spack env create TouchDetector
spacktivate -p TouchDetector

spack develop --no-clone -p ${PWD} touchdetector@develop
spack add touchdetector@develop

# Repeat this command:
spack install --jobs $(nproc) [--verbose]
```

### Build Variants {#build_variants}
In order to discover build variants of TD you may read the Spack recipe:
```
spack edit touchdetector
```
which is the definite source of knowledge regarding compiling TD.

The recipe should show
  - `+caliper` to enable the instrumentation of the code with Caliper
  performance measurements,
  - `+asan` to enable Clang's AddressSanitizer,
  - `+ubsan` to enable Clang's UndefinedBehaviorSanitizer,
  - `+clang-tidy` to enable static analysis with `clang-tidy`.

Please note that the sanitizers work best when using Clang as a compiler.

[bbp_spack_docs]: https://github.com/bluebrain/spack
[spack_docs]: https://spack.readthedocs.io/en/latest

### Directly using CMake
Directly using CMake is not recommended due to the number of dependencies.
Nevertheless, TD is CMake based and a build system can be configured using the
usual CMake commands, after figuring out the appropriate value for
`CMAKE_PREFIX_PATH`. Once again `spack` can help with this.
