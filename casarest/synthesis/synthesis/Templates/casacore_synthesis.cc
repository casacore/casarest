//# casacore_synthesis.cc: casacore template instantiations for synthesis
//# Copyright (C) 2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

// Generate casacore templates if needed.
// undef to auto-instantiate templates used by them.

#ifdef CASACORE_NO_AUTO_TEMPLATES
#undef CASACORE_NO_AUTO_TEMPLATES

// Instantiate all casacore templates needed.

#include <casa/Arrays/Vector.h>
#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/Cube.h>
#include <casa/Arrays/MaskedArray.h>
#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/Containers/SimOrdMap.h>
#include <casa/Containers/Stack.h>
#include <casa/Utilities/GenSort.h>
#include <casa/Utilities/PtrHolder.h>
#include <casa/Utilities/BinarySearch.h>
#include <casa/Quanta/Quantum.h>
#include <scimath/Mathematics/SquareMatrix.h>
#include <scimath/Mathematics/RigidVector.h>
#include <scimath/Mathematics/Gridder.h>
#include <scimath/Mathematics/ConvolveGridder.h>
#include <scimath/Mathematics/Convolver.h>
#include <measures/Measures/MeasBase.h>
#include <measures/Measures/MeasRef.h>
#include <measures/Measures/MeasConvert.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MCEpoch.h>
#include <measures/Measures/MCDoppler.h>
#include <measures/Measures/MRadialVelocity.h>
#include <measures/TableMeasures/ArrayMeasColumn.h>
#include <lattices/Lattices/Lattice.h>
#include <lattices/Lattices/SubLattice.h>
#include <lattices/Lattices/LatticeIterator.h>
#include <lattices/Lattices/LatticeCleaner.h>
#include <lattices/Lattices/LatticeConvolver.h>
#include <lattices/Lattices/LatticeCache.h>
#include <images/Images/PagedImage.h>
#include <images/Images/HDF5Image.h>
#include <images/Images/SubImage.h>
#include <images/Images/TempImage.h>
#include <images/Images/ImageRegrid.h>
#include <components/ComponentModels/Flux.h>

namespace casa { //# NAMESPACE - BEGIN

  template class Array<Complex>;
  template class Array<String>;
  template class Vector<float>;
  template class Matrix<float>;
  template class Matrix<double>;
  template class Cube<bool>;
  template class Cube<int>;
  template class Cube<float>;
  template class Cube<Complex>;
  template class Cube<DComplex>;
  template class MaskedArray<int>;
  template class MaskedArray<float>;
  template class GenSortIndirect<int>;
  template void minMax(float&, float&, IPosition&, IPosition&, const Array<float>&);
  template Array<float> operator*(const float&, const Array<float>&);
  template Array<float> operator*(const Array<float>&, const Array<float>&);
  template Array<Complex> abs(const Array<Complex>&);
  template int min(const Array<int>&);
  template Array<int> min(const Array<int>&, const Array<int>&);
  template unsigned int ntrue(const Array<bool>&);
  template int binarySearch(bool&, const Vector<double>&, const double&, unsigned int, int);
  template class Quantum<double>;
  template class Quantum<Vector<double> >;

  template class SquareMatrix<float,2>;
  template class SquareMatrix<float,4>;
  template class SquareMatrix<Complex,2>;
  template class SquareMatrix<Complex,4>;
  template class RigidVector<Complex,4>;
  template class Vector< SquareMatrix<float,2> >;
  template class Vector< SquareMatrix<Complex,2> >;
  template class Matrix< SquareMatrix<float,2> >;
  template class Matrix< SquareMatrix<Complex,2> >;
  template class Array< SquareMatrix<float,4> >;
  template class Cube< SquareMatrix<float,4> >;
  template class Array< SquareMatrix<Complex,4> >;
  template class Cube< SquareMatrix<Complex,4> >;
  template class Gridder<double, Complex>;
  template class ConvolveGridder<double, Complex>;
  template class Convolver<float>;

  template class MeasBase<MVRadialVelocity, MeasRef<MRadialVelocity> >;
  template class ROArrayMeasColumn<MEpoch>;
  template class MeasConvert<MEpoch>;

  template class Lattice<float>;
  template class Lattice<Complex>;
  template class SubLattice<float>;
  template class SubLattice<Complex>;
  template class LatticeIterator<float>;
  template class RO_LatticeIterator<float>;
  template class RO_LatticeIterator<Complex>;
  template class LatticeCleaner<float>;
  template class LatticeConvolver<float>;
  template class LatticeCache<float>;
  template class LatticeCache<Complex>;

  template class PagedImage<float>;
  template class PagedImage<Complex>;
  template class SubImage<float>;
  template class SubImage<Complex>;
  template class TempImage<float>;
  template class TempImage<Complex>;
  template class ImageRegrid<float>;
  template class PtrHolder<ImageInterface<float> >;
  template class HDF5Image<float>;

  template class Flux<double>;

  template class Matrix<bool>;
  template class SimpleOrderedMap<int, DataType>;
  template class Vector<MDirection>;
  template class Vector<RigidVector<double, 3> >;
  
  template Array<Complex> operator+(Array<Complex> const&, Array<Complex> const&);
  template double max(MaskedArray<double> const&);
  template double min(casa::MaskedArray<double> const&);

  template class Vector<int>;
  template class MeasConvert<MDoppler>;
  template class InterpolateArray1D<float, Complex>;
  template class InterpolateArray1D<float, float>;
  template class Link<void*>;
  template class Stack<void*>;

} //# NAMESPACE - END

#endif
