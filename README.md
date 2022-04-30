                        Library for Reduced Order Models
                      _   _   _      ___    ___    __  __ 
                     | | (_) | |__  | _ \  / _ \  |  \/  |
                     | | | | | '_ \ |   / | (_) | | |\/| |
                     |_| |_| |_.__/ |_|_\  \___/  |_|  |_|
                                                                             
                             https://librom.net

# Introduction

libROM is a free, lightweight, scalable C++ library for data-driven physical
simulation methods from the intrusive projection-based reduced order models to
non-intrusive black-box approaches.

## Features

- Dynamic data collection
- Dynamic mode decomposition (DMD)
- Data compression
- Physics-informed Greedy algorithm
- Projection-based hyper-reduction

## Features to be added

- S-OPT: stable sampling algorithm for hyper-reduction
- EQP: quadrature-based hyper-reduction sampling algorithm
- Python interface

# Installation

To compile libROM with default build settings (Mac and LLNL LC Machines):
```sh
 ./scripts/compile.sh
```

Compilation options:

- -a: Compile a special build for the LLNL codebase: Ardra
- -d: Compile in debug mode.
- -m: Compile with MFEM (required to run the libROM examples)
- -t: Use your own cmake/toolchain
- -u: Update all of libROM's dependencies.

# Installing via Spack

There is a Spack package for libROM; however, the version it installs
is the latest public release. See the [spack
documentation](http://spack.readthedocs.io/en/latest/index.html) for
details on how to use Spack.

To install libROM with default options using spack.

```sh
 spack install librom
```

To install libROM with MFEM using spack.

```sh
 spack install librom +mfem
```

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
- William Arrighi (LLNL)
- Kyle Chand (LLNL)
- Siu Wun Cheung (LLNL)
- Eric Chin (LLNL)
- Youngsoo Choi (LLNL)
- Dylan Copeland (LLNL)
- William Fries (University of Arizona)
- Debojyoti Ghosh (LLNL)
- Xiaolong He (UC San Diego)
- Kevin Huynh (LLNL)
- Tanya Kostova-Vassilevska (LLNL)
- Jessica Lauzon (Stanford)
- Sean McBane (UT Austin)
- Geoffrey Oxberry (LLNL)
