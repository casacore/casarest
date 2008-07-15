This is the casarest package, the remainder of the AIPS++ libraries

To build and install everything run 'batchbuild.py' like:
   batchbuild.py static install \
                 prefix=/dop131_1/gvandiep \
		 casacoreroot=/dop131_1/gvandiep \
		 casarestroot=/dop131_1/gvandiepen \
		 cfitsioroot=/usr/local cfitsiolibdir=/usr/local/lib64 \
		 wcsroot=/aips++/weekly/code/casa/wcs \
		 lapacklibdir=/usr/lib64 \
		 enable_hdf5 hdf5root=/dop131_1/gvandiep

Only use enable_hdf5 if casacore was built with hdf5.

If you don't have your own ~/.casarc, copy the provided casarc file to
~/.casarc.