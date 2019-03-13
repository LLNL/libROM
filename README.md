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

See the INSTALL file for details; in brief, it uses the `configure`,
`make`, `make install` sequence of commands, but you probably want to
check out the INSTALL file for some guidance on flags, or the output
of `./configure --help`.

# For developers who need to hack the build system

If you want to add new source files and executables, the best way to
do so is to modify `Makefile.am`, which is a GNU Automake file.

If you want to add new dependencies (e.g., you want to use a new
linear algebra library), the best way to add all of the include/link
flag information is to modify `configure.ac`, which is a GNU Autoconf
file. As part of modifying this file, if you want to add new macros
copied from an m4sh macro library, add them to the `config` folder.

After changing `configure.ac` or `Makefile.am`, run the following
command in the root of this repository:

```console
autoreconf -I config --install
```

A full explanation of GNU Autotools is out of the scope of this
README, but briefly:

* `aclocal` takes `configure.ac` and a bunch of m4sh files -- these
  all have the suffix `.m4` and makes an `aclocal.m4` file containing
  all of the macros needed by `configure.ac`

* `autoconf` takes `configure.ac` and `aclocal.m4` and generates a
  `configure` script

* `automake` takes `Makefile.am`, `configure.ac`, and `aclocal.m4` and
  generates a `Makefile.in` file

* `configure` takes `Makefile.in` and generates a `Makefile`, which
  can be used with GNU Make to build and install libROM

The `autoreconf` command runs these tools -- and others -- in the
correct order to generate a working `configure` script, along with all
of the necessary intermediate files, some of which were mentioned in
the previous list.

The following resources are moderately useful at communicating how to
modify these files, in (roughly) descending order of utility:

* [Wikipedia article on GNU
  Autotools](https://en.wikipedia.org/wiki/GNU_Build_System): Start
  here. Mostly useful for the diagram explaining the inputs and
  outputs of each tool, plus their dependency relationships

* [Autotools Mythbuster](https://autotools.io/index.html): Read after
  the Wikipedia article. Good at communicating the basics, but
  occasionally misses important finer points (e.g., how to add include
  flags and link flags to Automake files).

* [GNU Automake
  Manual](https://www.gnu.org/software/automake/manual/): Ignore all
  of Chapter 1 and Chapter 2, except Section 2.4. Section 2.4
  illustrates the basics covered in Autotools Mythbuster, but more
  verbosely. If you're in a hurry and you're getting errors due to
  `configure.ac`, skip to Section 6.1. If you're in a hurry and you
  want to figure out how to hack `Makefile.am`, skip to Chapter 8 if
  you're trying to add source files or executables. Skip to Chapter 9
  if you're trying to add commands to call scripts, generate headers,
  generate example data, etc. Skip to Chapter 15 if you're trying to
  add test harnesses.

* [GNU Autoconf
  Manual](https://www.gnu.org/software/autoconf/manual/index.html):
  Autotools Mythbuster has a better introduction; use this resource if
  you get stuck.

* [GNU Autoconf
  Archive](https://www.gnu.org/software/autoconf-archive/index.html):
  Collection of prebuilt m4sh macros for autoconf. Useful to search
  for macros that might do what you need (e.g., detect C++
  features). To "install" a macro from here, add a text file with the
  macro to the `config` directory.

* [GNU Manuals](https://www.gnu.org/manual/manual.html): Useful if you
  get really stuck, and the problem is in some other part of the GNU
  Autotools, but try searching for your problem on Google first.

* [Autotools Tutorial](https://www.lrde.epita.fr/~adl/autotools.html):
  Lengthy introduction at a much slower pace than Autotools Mythbuster
  or the GNU Autotools manuals. The format is 500+ LaTeX Beamer slides
  as a big PDF, and is more suitable for a half-day to day-long
  workshop. You're likely to find what you want to know faster by
  searching/reading the previous resources.
