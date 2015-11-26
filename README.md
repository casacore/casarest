
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
* casacore (2.0 or later)
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

And the optional libraries:
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
