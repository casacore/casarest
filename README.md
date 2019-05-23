
# Casarest

This is the casarest package, the remainder of the AIPS++
libraries. It consists of the libraries:
 * msvis
 * calibration
 * synthesis
 * flagging
 * simulators
 * ionosphere

The prorgram lwimager (part of synthesis) is the main deliverable.


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
* casacore (3.0 or later)
* boost
* wcslib
* cfitsio
* fortran 
* hdf5 (optional)

If you are still at casacore 2.0 then use casarest 1.4.2.


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


## Ubuntu 14.04 packages

If you run Ubuntu 14.04 you can use precompiled binary packages

https://launchpad.net/~radio-astro/+archive/ubuntu/main

installation commands:
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:radio-astro/main
sudo apt-get update
sudo apt-get install casarest
```


# Problems & bugs

If you have any issues compiling or using casacore, please open an issue on
the issue tracker on github.


# travis

[![Build Status](https://travis-ci.org/casacore/casarest.svg?branch=master)](https://travis-ci.org/casacore/casarest)
