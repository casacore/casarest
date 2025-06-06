//# tProfileFit1D.cc: test the ProfileFit1D class
//# Copyright (C) 1995,1996,1998,1999,2000,2001,2002,2004
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Containers/Record.h>
#include <components/SpectralComponents/PowerLogPolynomialSpectralElement.h>
#include <components/SpectralComponents/SpectralElementFactory.h>

#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/IO/ArrayIO.h>

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/BasicMath/Functors.h>

#include <casacore/casa/iostream.h>

#include <casacore/casa/namespace.h>

#include <memory>

int main() {
	{
		try {
			cout << "Test constructor" << endl;
			Vector<Double> p(3);
			p[0] = 5.5;
			p[1] = 2.2;
			p[2] = 3.3;
			PowerLogPolynomialSpectralElement plp(p);
			AlwaysAssert(near(plp(2), 123.36739264171082), AipsError);

			cout << plp << endl;
			AlwaysAssert(allTrue(plp.get() == p), AipsError);
			cout << "Test to/from record" << endl;
			Record rec;
			plp.toRecord(rec);
			std::unique_ptr<SpectralElement> el(SpectralElementFactory::fromRecord(rec));
			plp = *dynamic_cast<PowerLogPolynomialSpectralElement *>(el.get());
			AlwaysAssert(allTrue(plp.get() == p), AipsError);
		}
		catch (const AipsError& x) {
			cout << x.getMesg() << endl;
			return 1;
		}
	}

	cout << "ok" << endl;
	return 0;


}
