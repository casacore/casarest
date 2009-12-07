//# casarest_synthesis.cc: casarest template instantiations for synthesis
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

// Generate casarest templates if needed.
// undef to auto-instantiate templates used by them.

#ifdef CASACORE_NO_AUTO_TEMPLATES
#undef CASACORE_NO_AUTO_TEMPLATES
#endif

#ifdef AIPS_NO_TEMPLATE_SRC

// Instantiate all casarest templates needed.

#include <casa/Utilities/CountedPtr.h>
#include <casa/Arrays/Array.h>
#include <lattices/Lattices/Lattice.h>
#include <synthesis/MeasurementComponents/PBMathInterface.h>
#include <synthesis/MeasurementEquations/ArrayModel.tcc>
#include <synthesis/MeasurementEquations/LinearEquation.tcc>
#include <synthesis/MeasurementEquations/LinearModel.tcc>
#include <synthesis/MeasurementEquations/ResidualEquation.tcc>
#include <calibration/CalTables/CalSet.tcc>
#include <calibration/CalTables/SolvableCalSetMCol.tcc>
#include <calibration/CalTables/ROCalMainColumns2.tcc>
#include <calibration/CalTables/CalMainColumns2.tcc>
#include <msvis/MSVis/StokesVector.h>

namespace casa { //# NAMESPACE - BEGIN

  template class CountedPtr<PBMathInterface>;
  template class ArrayModel<Float>;
  template class LinearEquation<Array<Float>, Array<Float> >;
  template class LinearEquation<Lattice<Float>, Lattice<Float> >;
  template class LinearModel<Array<Float> >;
  template class LinearModel<Lattice<Float> >;
  template class ResidualEquation<Array<Float> >;
  template class ResidualEquation<Lattice<Float> >;
  template class Matrix<CStokesVector>;
  template class CalSet<Float>;
  template class CalSet<Complex>;
  template class ROSolvableCalSetMCol<Complex>;
  template class SolvableCalSetMCol<Complex>;
  template class ROCalMainColumns2<float>;
  template class ROCalMainColumns2<Complex>;
  template class CalMainColumns2<float>;
  template class CalMainColumns2<Complex>;

} //# NAMESPACE - END

#endif
