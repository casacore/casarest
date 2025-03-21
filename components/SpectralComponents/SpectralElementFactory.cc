//# Spectral2Element.cc: Record conversion for SpectralElement
//# Copyright (C) 2001,2004
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
//# $Id: SpectralElementFactory.cc 21465 2014-06-19 05:56:56Z gervandiepen $

//# Includes

#include <casacore/casa/Containers/Record.h>
#include <components/SpectralComponents/CompiledSpectralElement.h>
#include <components/SpectralComponents/GaussianSpectralElement.h>
#include <components/SpectralComponents/GaussianMultipletSpectralElement.h>
#include <components/SpectralComponents/LogTransformedPolynomialSpectralElement.h>
#include <components/SpectralComponents/LorentzianSpectralElement.h>
#include <components/SpectralComponents/PolynomialSpectralElement.h>
#include <components/SpectralComponents/PowerLogPolynomialSpectralElement.h>
#include <components/SpectralComponents/SpectralElementFactory.h>

#include <memory>

namespace casacore { //# NAMESPACE CASACORE - BEGIN

SpectralElement* SpectralElementFactory::fromRecord(
	const RecordInterface &in
) {
	std::unique_ptr<SpectralElement> specEl;
	String origin = "SpectralElementFactory::fromRecord: ";
	if (
		! in.isDefined("type")
		|| in.type(in.idToNumber(RecordFieldId("type"))) != TpString
	) {
		throw AipsError("Record does not represent a SpectralElement");
	}
	String stp;
	SpectralElement::Types tp;
	in.get(RecordFieldId("type"), stp);
	if (!SpectralElement::toType(tp, stp)) {
		throw AipsError("Unknown spectral type in SpectralElement::fromRecord\n");
	}

	Vector<Double> errs;

	// Get the errors if defined in record

	if (in.isDefined("errors")) {
		DataType dataType = in.dataType("errors");
		if (dataType == TpArrayDouble) {
			in.get("errors", errs);
		}
		else if (dataType == TpArrayFloat) {
			Vector<Float> v;
			in.get("errors", v);
			errs.resize(v.nelements());
			convertArray(errs, v);
		}
		else if (dataType == TpArrayInt) {
			Vector<Int> v;
			in.get("errors", v);
			errs.resize(v.nelements());
			convertArray(errs, v);
		}
		else {
			throw AipsError(
				"SpectralElement::fromRecord: errors field "
				"must be double, float or int\n"
			);
		}
	}

	Vector<Double> param;
	if (in.isDefined(("parameters"))) {
		DataType dataType = in.dataType("parameters");
		if (dataType == TpArrayDouble) {
			in.get("parameters", param);
		}
		else if (dataType == TpArrayFloat) {
			Vector<Float> v;
			in.get("parameters", v);
			param.resize(v.nelements());
			convertArray(param, v);
		}
		else if (dataType == TpArrayInt) {
			Vector<Int> v;
			in.get("parameters", v);
			param.resize(v.nelements());
			convertArray(param, v);
		}
		else {
			throw AipsError(
				origin +
				"SpectralElement::fromRecord: parameters field "
				"must be double, float or int\n"
			);
		}
	}

	// Make sizes of errors and params equal
	if (errs.nelements() == 0) {
		errs.resize(param.nelements());
		errs = 0.0;
	}
	if (errs.nelements() != param.nelements()) {
		throw AipsError(
			origin +
			"SpectralElement::fromRecord must have equal lengths "
			"for parameters and errors fields"
		);
	}
	switch (tp) {
	case SpectralElement::GAUSSIAN:
		if (param.nelements() != 3) {
			throw AipsError(
				origin + "Illegal number of parameters for Gaussian element"
			);
		}
		if (param(2) <= 0.0) {
			throw AipsError(
				origin +
				"The width of a Gaussian element must be positive"
			);
		}
		param(2) = GaussianSpectralElement::sigmaFromFWHM (param(2));
		errs(2) = GaussianSpectralElement::sigmaFromFWHM (errs(2));
		specEl.reset(new GaussianSpectralElement(param(0), param(1), param(2)));
		specEl->setError(errs);
		break;
	case SpectralElement::LORENTZIAN:
		if (param.nelements() != 3) {
			throw AipsError(
				"Illegal number of parameters for Lorentzian element"
			);
		}
		if (param(2) <= 0.0) {
			throw AipsError(
				"The width of a Lorentzian element must be positive"
			);
		}
                specEl.reset(new LorentzianSpectralElement(param(0), param(1), param(2)));
                specEl->setError(errs);
                break;
	case SpectralElement::POLYNOMIAL:
		if (param.nelements() == 0) {
			throw AipsError(
				"Polynomial spectral element must have order "
				"of at least zero"
			);
		}
		specEl.reset(new PolynomialSpectralElement(param.nelements() - 1));
		specEl->set(param);
		specEl->setError(errs);
		break;
	case SpectralElement::COMPILED:
		if (
				in.isDefined("compiled")
				&& in.type(in.idToNumber(RecordFieldId("compiled"))) == TpString
		) {
			String function;
			in.get(RecordFieldId("compiled"), function);
			specEl.reset(new CompiledSpectralElement(function, param));
			specEl->setError(errs);
		}
		else {
			throw AipsError(
				"No compiled string in SpectralElement::fromRecord\n"
			);
		}
		break;

	case SpectralElement::GMULTIPLET:
	{
		if (! in.isDefined("gaussians")) {
			throw AipsError("gaussians not defined in record");
		}
		if (! in.isDefined("fixedMatrix")) {
			throw AipsError("fixed matrix not defined in record");
		}
		Record gaussians = in.asRecord("gaussians");
		uInt i = 0;
		std::vector<GaussianSpectralElement> comps(0);
		while(True) {
			String id = "*" + String::toString(i);
			if (gaussians.isDefined(id)) {
                            // MV: I found this code rather untidy w.r.t. pointer management, but leave it as is for now
				std::unique_ptr<SpectralElement> gauss(fromRecord(gaussians.asRecord(id)));
				comps.push_back(
					*dynamic_cast<GaussianSpectralElement*>(
						gauss.get()
					)
				);
				i++;
			}
			else {
				break;
			}
		}
		Matrix<Double> fixedMatrix = in.asArrayDouble("fixedMatrix");
		fixedMatrix.reform(IPosition(2, comps.size()-1, 3));
		specEl.reset(new GaussianMultipletSpectralElement(comps, fixedMatrix));
	}
	break;

    case SpectralElement::POWERLOGPOLY: {
		specEl.reset(new PowerLogPolynomialSpectralElement(param));
		specEl->set(param);
		specEl->setError(errs);
	}
	break;

    case SpectralElement::LOGTRANSPOLY: {
                specEl.reset(new LogTransformedPolynomialSpectralElement(param));
                specEl->set(param);
    		specEl->setError(errs);
    	}
    	break;

	default:
		throw AipsError(
			"Unhandled or illegal spectral element record in "
			"SpectralElementFactory::fromRecord\n"
		);
	}

	if (in.isDefined("fixed")) {
		specEl->fix(in.asArrayBool("fixed"));
	}
    // ready to return, fish out the pointer and return it without deleting it
    SpectralElement *sp = specEl.release();
    return sp;
}


} //# NAMESPACE CASACORE - END

