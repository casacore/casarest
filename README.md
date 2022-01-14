
# Casarest

This is the casarest package, the remainder of the AIPS++
libraries. It consists of the libraries:
 * msvis
 * calibration
 * synthesis
 * flagging
 * simulators
 * ionosphere

The program lwimager (part of synthesis) is the main deliverable.


# Installation

## Obtaining the source

The casarest source code is maintained on github.

You can obtain it using:

```
$ git clone https://github.com/casacore/casarest
```

## Requirements

To compile casarest you need to meet the following requirements:


* cmake
* g++
* casacore (3.4 or later)
* boost
* wcslib
* cfitsio
* fortran 
* hdf5 (optional)


On Debian / Ubuntu you can install these with:
``` 
$ sudo apt-get install cmake libcasacore2-dev libboost-dev wcslib-dev \
   libcfitsio3-dev libboost-system-dev libboost-thread-dev gfortran \
   g++
```

And the optional dependencies:
``` 
$ sudo apt-get install libhdf5-dev
```

## Compilation

In the casacore source folder run:
```
mkdir build
cd build
cmake ..
make 
make install
```

There are 4 components in the casarest package and you may not need all them. The four targets are: casa_components, casa_msvis, casa_calibration and casa_synthesis. You can make them individuallay by:
```
mkdir build
cd build
cmake ..
make <choice of target>
cmake -DCOMPONENT=<choice of target>  -P cmake_install.cmake 
```
Which will build and install only your choice of target.


## KERN packages

If you run Ubuntuyou can use precompiled binary packages

https://kernsuite.info/


# Problems & bugs

If you have any issues compiling or using casacore, please open an issue on
the issue tracker on github.

