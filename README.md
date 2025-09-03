![libROM Logo](https://www.librom.net/img/logo-300.png)

# Introduction

[libROM](https://www.librom.net) is a free, lightweight, scalable C++ library for data-driven physical
simulation methods from the intrusive projection-based reduced order models to
non-intrusive black-box approaches.

The best starting point for new users interested in libROM's features 
is to review the [examples](https://www.librom.net/examples.html). 
The [code documentation](https://librom.readthedocs.io/en/latest/index.html) 
provides more details about libROM's classes and functions. 

## Features

- Dynamic data collection
- Dynamic mode decomposition (DMD)
- Data compression
- Physics-informed greedy algorithm
- Projection-based hyper-reduction
- EQP: quadrature-based hyper-reduction sampling algorithm
- [Python interface](https://github.com/LLNL/pylibROM)

## Features to be added

# Installation

To compile libROM with default build settings (Mac and LLNL LC Machines):
```sh
 ./scripts/compile.sh
```

Compilation options:

- -a: Compile a special build for the LLNL codebase: Ardra
- -d: Compile in debug mode.
- -f: Specify additional compiler flags
- -m: Compile with MFEM (required to run the libROM examples)
- -g: Compile MFEM with GSLib (requires -m)
- -r: Compile unit tests (requires Googletest)
- -s: Compile and use a local SCALAPACK
- -t: Use your own cmake/toolchain
- -u: Update all of libROM's dependencies.

## Compiling on LC Machines

libROM provides several CMake toolchains which can be used to compile on LLNL LC machines. 
For more information on installing and using libROM on specific LC machines, 
refer to [this wiki page](https://github.com/LLNL/libROM/wiki/Compiling-on-LC-Machines).

# Installing via Spack

There is a Spack package for libROM; however, the version it installs
is the latest public release. See the [spack
documentation](https://spack.readthedocs.io/en/latest/index.html) for
details on how to use Spack.

To install libROM with default options using spack.

```sh
 spack install librom
```

To install libROM with MFEM using spack.

```sh
 spack install librom +mfem
```

# Compiling and linking with libROM

To compile and link an existing code with libROM, follow these steps:

- Add libROM/lib to the include path
```sh
 -I/path/to/libROM/lib
```
- Add the following to the linker flags (LDFLAGS)
```sh
 -Wl,-rpath,/path/to/libROM/build/lib -L/path/to/libROM/build/lib
```
- Add the following library
```sh
 -lROM
```

For example,
```sh
mpicxx myapp.cpp -I/path/to/libROM/lib -Wl,-rpath,/path/to/libROM/build/lib -L/path/to/libROM/build/lib -lROM -o myapp.out
```

# Using Docker container

Docker container [`librom_env`](https://ghcr.io/llnl/librom/librom_env) provides a containerized environment with all the prerequisites for libROM. For instruction on how to use it, check out [the wiki page](https://github.com/LLNL/libROM/wiki/Using-Docker-container).

# libROM CI

libROM leverages GitHub Actions for CI. The CI currently applies only to commits to pull requests.  Unit tests run for all PR commits. Upon the addition of the `LGTM` label, both the unit tests and regression tests run. While the `LGTM` label is still present, all subsequent commits run both unit tests and regression tests.

To compile and run unit tests locally, build using the `-r` option to `compile.sh` or with `-DENABLE_TESTS=ON`. Building the unit tests will require Googletest to be installed. Unit tests can be run using `ctest` from the root build directory.

# License

libROM is distributed under the terms of both the MIT license and the
Apache License (Version 2.0). Users may choose either license at their
option.

All new contributions must be made under both the MIT and Apache-2.0 licenses.

See
[LICENSE-MIT](https://github.com/LLNL/libROM/blob/master/LICENSE-MIT),
[LICENSE-APACHE](https://github.com/LLNL/libROM/blob/master/LICENSE-APACHE),
[COPYRIGHT](https://github.com/LLNL/libROM/blob/master/COPYRIGHT), and
[NOTICE](https://github.com/LLNL/libROM/blob/master/NOTICE) for
details.

Up to commit 299876e0a0304f25db56f1f9e2eb2c61ef199048, libROM was
previously released under the terms of the BSD-3 license.

SPDX_License-Identifier: (Apache-2.0 OR MIT)

LLNL-CODE-686965 (up to commit 299876e0a0304f25db56f1f9e2eb2c61ef199048)
LLNL-CODE-766763


# Authors
- Robert W. Anderson (LLNL)
- William Anderson (LLNL)
- William Arrighi (LLNL)
- Kyle Chand (LLNL)
- Siu Wun Cheung (LLNL)
- Eric Chin (LLNL)
- Youngsoo Choi (LLNL)
- "Kevin" Seung Whan Chung (LLNL)
- Dylan Copeland (LLNL)
- William Fries (University of Arizona)
- Debojyoti Ghosh (LLNL)
- Xiaolong He (UC San Diego)
- Adrian Humphry (University of Toronto)
- Kevin Huynh (LLNL)
- Coleman Kendrick (LLNL)
- Minji Kim (UNC)
- Tanya Kostova-Vassilevska (LLNL)
- Axel Larsson (Princeton)
- Jessica Lauzon (LLNL)
- Jacob Lotz (TUDelft)
- Sean McBane (UT Austin)
- Geoffrey Oxberry (LLNL)
- Yeonjong Shin (KAIST)
- Seung Won Suh (UIUC)
- Paul Tranquilli (LLNL)
- Chris Vales (Dartmouth)
- Pranav Vempati (LLNL)
- Masayuki Yano (University of Toronto)
