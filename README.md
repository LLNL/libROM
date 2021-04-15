# Introduction

libROM is a library to compute proper orthogonal decomposition-based
reduced order models (POD-based ROMs). It also contains some code to
compute hyperreduced POD-based ROMs using the discrete empirical
interpolation method (DEIM).

# Authors

The original libROM release was written by Bill Arrighi based on
prototype MATLAB code written by Geoffrey Oxberry and Kyle Chand. This
MATLAB code implemented the algorithm presented in

Geoffrey M. Oxberry, Tanya Kostova-Vassilevska, William Arrighi, and
Kyle Chand, [_Limited-memory adaptive snapshot selection for proper
orthogonal
decomposition_](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.5283),
International Journal of Numerical Methods in Engineering,
**109**:198--217.

The key algorithmic idea of computing an SVD incrementally to train
POD-based ROMs on the fly was conceived by Kyle Chand and refined by
Geoffrey Oxberry based on his PhD thesis work, with collaborative
contributions by Bill Arrighi and Tanya Kostova-Vassilevska.

Subsequent commits to libROM have been made by Geoffrey Oxberry,
Robert W. Anderson, Youngsoo Choi, and Dylan Copeland. In particular,
Bill Arrighi and Robert W. Anderson have contributed an implementation
of DEIM to libROM, and Dylan Copeland has contributed an
implementation of GNAT to libROM.

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

To copmile libROM for Laghos (Mac and LLNL LC Machines):
1. ./scripts/laghos_compile.sh

To compile libROM for Ardra (LLNL LC Machines only):
1. ./scripts/ardra_compile.sh

# Installing via Spack

There is a Spack package for libROM; however, the version it installs
is the latest public release. See the [spack
documentation](http://spack.readthedocs.io/en/latest/index.html) for
details on how to use Spack.
