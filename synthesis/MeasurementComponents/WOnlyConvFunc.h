// -*- C++ -*-
//# WOnlyConvFunc.h: Definition of the WOnlyConvFunc class
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
#ifndef SYNTHESIS_WONLYCONVFUNC_H
#define SYNTHESIS_WONLYCONVFUNC_H

#include <synthesis/MeasurementComponents/ConvolutionFunction.h>
#include <synthesis/MeasurementComponents/CFStore.h>
#include <synthesis/MeasurementComponents/WTerm.h>
#include <synthesis/MeasurementComponents/PSTerm.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogOrigin.h>
#include <casacore/casa/Arrays/ArrayFwd.h>

namespace casacore { //# NAMESPACE CASACORE - BEGIN
  template<class T> class ImageInterface;
  class VisBuffer;
  //
  //-------------------------------------------------------------------------------------------
  //
  class WOnlyConvFunc : public ConvolutionFunction
  {
  public:
    WOnlyConvFunc(const CountedPtr<WTerm>& wTerm,
		  const CountedPtr<PSTerm>& psTerm):
      ConvolutionFunction(),wTerm_p(wTerm), psTerm_p(psTerm)
    {};
    ~WOnlyConvFunc() {};
    WOnlyConvFunc& operator=(const WOnlyConvFunc& other);


    void makeConvFunction(const ImageInterface<Complex>& image,
			  const VisBuffer& vb,
			  const Int wConvSize,
			  const Float pa,
			  CFStore& cfs,
			  CFStore& cfwts);

    virtual int getVisParams(const VisBuffer& vb) {return wTerm_p->getVisParams(vb);};
    virtual void setPolMap(const Vector<Int>& polMap) {wTerm_p->setPolMap(polMap);};

    virtual Bool findSupport(Array<Complex>& func, Float& threshold,Int& origin, Int& R);
  protected:
    void setSupport(Array<Complex>& convFunc, CFStore& cfs);

  private:
    CountedPtr<WTerm> wTerm_p;
    CountedPtr<PSTerm> psTerm_p;
  };
  //
  //-------------------------------------------------------------------------------------------
  //
};
#endif
