# Introduction

libROM is a free, lightweight, scalable C++ library for data-driven physical
simulation methods from the intrusive projection-based reduced order models to
non-intrusive black-box approaches. 


# History of libROM

The original libROM release was written by William Arrighi (retired) mainly to
implement a C++ scalable and incremental singular value decomposition based on
the journal paper, i.e., 

> Geoffrey M. Oxberry, Tanya Kostova-Vassilevska, William Arrighi, and Kyle Chand, [_Limited-memory adaptive snapshot selection for proper orthogonal decomposition_] (https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5283), International Journal of Numerical Methods in Engineering, **109**:198--217.

In addition to the incremental singular value decomposition, many other features
are added:

## Features

- Dynamic data collection
- Data compression
- Greedy algorithm
- Hyper-reduction

# Authors
Robert W. Anderson (LLNL),
Siu Wun Cheung (LLNL),
Youngsoo Choi (LLNL),
Dylan Copeland (LLNL),
Kevin Huynh (LLNL),
Jessica Lauzon (Stanford),
Sean McBane (UT Austin),
Geoffrey Oxberry (LLNL).

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

# Installation

To compile libROM with default build settings (Mac and LLNL LC Machines):
```sh
 ./scripts/compile.sh
```

To compile libROM for Ardra (LLNL LC Machines only):
```sh
./scripts/ardra_compile.sh
```

To compile libROM using a different toolchain within cmake/toolchains (Mac and LLNL LC Machines):
```sh
./scripts/toolchain_compile.sh toolchain.cmake
```

# Installing via Spack

There is a Spack package for libROM; however, the version it installs
is the latest public release. See the [spack
documentation](http://spack.readthedocs.io/en/latest/index.html) for
details on how to use Spack.
