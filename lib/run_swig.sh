#!/bin/bash 

/usr/workspace/huynh24/swig-4.0.2/SWIG_INSTALL/bin/swig -I/usr/workspace/huynh24/swig-4.0.2/SWIG_INSTALL/share/swig/4.0.2 -I/usr/workspace/huynh24/swig-4.0.2/SWIG_INSTALL/share/swig/4.0.2/python -I/usr/workspace/huynh24/swig-4.0.2/SWIG_INSTALL/share/swig/4.0.2/std -c++ -python librom.i
mpic++ -O2 -fPIC -c librom_wrap.cxx -I/usr/include/python2.7 -I/usr/tce/packages/python/python-2.7.16/lib/python2.7/site-packages/mpi4py/include -I/usr/workspace/huynh24/libROM3/lib
mpic++ -shared -Wl,-rpath,/usr/workspace/huynh24/libROM3/build/lib -L/usr/workspace/huynh24/libROM3/build/lib librom_wrap.o -lROM -o _pyROM.so
