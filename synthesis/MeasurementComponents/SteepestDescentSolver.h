//# ConjugateGradientSolver.h: Definition for an iterative ConjugateGradientSolver
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
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
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#ifndef SYNTHESIS_STEEPESTDESCENT_H
#define SYNTHESIS_STEEPESTDESCENT_H
//----------------------------------------------------------------------------

#include <casacore/casa/aips.h>

#include <casarest/synthesis/MeasurementEquations/VisEquation.h>
#include <casarest/synthesis/MeasurementEquations/Iterate.h>
#include <casarest/synthesis/MeasurementComponents/EPJones.h>
#include <casarest/msvis/MSVis/VisBuffer.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/Timer.h>

#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogMessage.h>

namespace casa {

  class SteepestDescentSolver: public Iterate
  {
  public:
    SteepestDescentSolver(Int nParams, Vector<Int> polMap, 
			  Int nIter=20, Double tol=1E-5);
    ~SteepestDescentSolver(){};

    void setMaxParams(Int nParams){maxParams=nParams;};
    
    Double solve(VisEquation& ve, EPJones& vj, VisBuffer& vb,
		 Int nAnt, Int SlotNo);

    Double solve2(VisEquation& ve, VisIter& vi, EPJones& epj, Int nAnt, Int SlotNo);

    LogIO& logIO() {return logIO_p;};
    
  private:
    Vector<Complex> getVj(const VisBuffer& vb, Int NAnt, Int whichAnt, Int whichPol,
			  Double& sumWt,Int negate=0, Int weighted=1);
    Double getGOF(const VisBuffer& residual,Int& whichPol, Double& sumWt,char *msg="");

    Int maxIter,maxParams;
    
    VisBuffer residual_p,gradient0_p,gradient1_p;
    Matrix<Bool> flags;

    LogIO logIO_p;
    Vector<Int> polMap_p;
  };
};

#endif
