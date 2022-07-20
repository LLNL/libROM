#!/bin/bash

# Check whether Homebrew or wget is installed
if [ "$(uname)" == "Darwin" ]; then
  which -s brew > /dev/null
  if [[ $? != 0 ]] ; then
      # Install Homebrew
      echo "Homebrew installation is required."
      exit 1
  fi

  which -s wget > /dev/null
  # Install wget
  if [[ $? != 0 ]] ; then
      brew install wget
  fi
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
LIB_DIR=$SCRIPT_DIR/../dependencies
mkdir -p $LIB_DIR

# Install SWIG-4.0.2
cd $LIB_DIR
if [ ! -d "swig" ]; then
    wget http://prdownloads.sourceforge.net/swig/swig-4.0.2.tar.gz
    tar -xvf swig-4.0.2.tar.gz
    mv swig-4.0.2 swig
    cd swig
    ./configure --prefix=$LIB_DIR/swig/swig_install
    make
    make install

    # We also need the numpy swig files
    cd swig_install/share/swig/4.0.2
    wget https://raw.githubusercontent.com/numpy/numpy/main/tools/swig/numpy.i
    wget https://raw.githubusercontent.com/numpy/numpy/main/tools/swig/pyfragments.swg

    # We have to delete these two files for some reason or else SWIG doesn't work.
    rm std_except.i
    rm exception.i
fi
