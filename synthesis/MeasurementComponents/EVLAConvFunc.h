// -*- C++ -*-
//# EVLAConvFunc.h: Definition of the EVLAConvFunc class
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
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
//# $Id$
//
#ifndef SYNTHESIS_EVLACONVFUNC_H
#define SYNTHESIS_EVLACONVFUNC_H

#include <casacore/images/Images/ImageInterface.h>
#include <synthesis/MeasurementComponents/Utils.h>
#include <synthesis/MeasurementComponents/BeamCalc.h>
#include <synthesis/MeasurementComponents/CFStore.h>
#include <synthesis/MeasurementComponents/VLACalcIlluminationConvFunc.h>
//#include <synthesis/MeasurementComponents/IlluminationConvFunc.h>
//#include <synthesis/MeasurementComponents/PixelatedConvFunc.h>
#include <synthesis/MeasurementComponents/ConvolutionFunction.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#if defined(casacore)
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#else
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#endif
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogOrigin.h>
#include <casacore/casa/Arrays/ArrayFwd.h>

namespace casacore { //# NAMESPACE CASACORE - BEGIN
  template<class T> class ImageInterface;
  class VisBuffer;
  class EVLAConvFunc : public ConvolutionFunction
  //: public PixelatedConvFunc<Complex>
  {
  public:
    //    EVLAConvFunc(const CountedPtr<IlluminationConvFunc> ATerm):
    //      ConvolutionFunction(),bandID_p(-1), polMap_p(), feedStokes_p(), ATerm_p(ATerm)
    EVLAConvFunc():     
      ConvolutionFunction(),bandID_p(-1), polMap_p(), feedStokes_p()
    {};
    ~EVLAConvFunc() {};
    EVLAConvFunc& operator=(const EVLAConvFunc& other);
    Int getVLABandID(Double& freq,String&telescopeName);
    Bool findSupport(Array<Complex>& func, Float& threshold,Int& origin, Int& R);
    void makeConvFunction(const ImageInterface<Complex>& image,
			  const VisBuffer& vb,
			  const Int wConvSize,
			  const Float pa,
			  CFStore& cfs,
			  CFStore& cfwts);
    int getVisParams(const VisBuffer& vb);
    Int makePBPolnCoords(const VisBuffer&vb,
			 const Vector<Int>& polMap,
			 const Int& convSize,
			 const Int& convSampling,
			 const CoordinateSystem& skyCoord,
			 const Int& skyNx, const Int& skyNy,
			 CoordinateSystem& feedCoord,
			 Vector<Int>& cfStokes);
    //
    // Overloading these functions from ConvolutionFunction class
    //
    void setPolMap(const Vector<Int>& polMap);
    void setFeedStokes(const Vector<Int>& feedStokes);
  private:
    Int bandID_p;
    Float Diameter_p, Nant_p, HPBW, sigma;
    
    LogIO& logIO() {return logIO_p;}
    LogIO logIO_p;
    Vector<Int> polMap_p;
    Vector<Int> feedStokes_p;
    CountedPtr<IlluminationConvFunc> ATerm_p;
  };
};
#endif
