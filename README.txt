This is the casarest package, the remainder of the AIPS++
libraries. It consists of the libraries:
   msvis
   calibration
   synthesis
   flagging
   simulators
   ionosphere
The prorgram lwimager (part of synthesis) is the main deliverable.

It can be checked out like:
 svn co svn+ssh://user@lofar9.astron.nl/var/svn/repos/trunk/casarest casarest

It uses cmake as its build system (minimum version is 2.6)
and can be built like:

 cd casarest
 mkdir build
 cd build
 cmake .. -DCASACORE_ROOT_DIR=/home/user
   -DHDF5_ROOT_DIR=/Users/diepen/hdf5-3xx ..
   -DCMAKE_INSTALL_PREFIX=/home/user
   -DLIB_EXTRA_SYNTHESIS=gfortran
   -DBUILD_ALL=1

By default only the msvis, calibration, and synthesis libraries are
built as well as the lwimager program.
If -DBUILD_ALL=1 is given, the flagging and simulators libraries are
built too.
Ionosphere is not built at the moment, because it depends on the PIM package.

The synthesis library is a mix of Fortran and C++ code, hence the
fortran library needs to known when linking synthesis.
Alas cmake versions before 2.8 did not handle it well. Therefore this
library needs to be given using -DLIB_EXTRA_SYNTHESIS. It defaults to
gfortran, which is usually fine. 

HDF5_ROOT_DIR only needs to be given if casacore was built with HDF5
support.

If you don't have your own ~/.casarc, copy the provided casarc file to
~/.casarc.
