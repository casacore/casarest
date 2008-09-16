//# lwimager.cc: Program to create and/or clean an image
//# Copyright (C) 2008
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
//# $Id: lwimager.cc 20245 2008-02-14 11:05:10Z gervandiepen $

//# Includes
#include <casa/aips.h>
#include <synthesis/MeasurementEquations/Imager.h>
#include <images/Images/PagedImage.h>
#ifdef HAVE_HDF5
# include <images/Images/HDF5Image.h>
#endif
#include <images/Images/ImageFITSConverter.h>
#include <casa/Inputs.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Utilities/Regex.h>
#include <casa/Utilities/Assert.h>
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <casa/sstream.h>

using namespace casa;

IPosition handlePos (const IPosition& pos, const IPosition& def)
{
  if (pos.nelements() == 0) {
    return def;
  }
  if (pos.nelements() != 2) {
    throw AipsError("Give 0 or 2 values in maskblc and masktrc");
  }
  IPosition npos(def);
  int n = npos.nelements();
  npos[n-1] = pos[1];
  npos[n-2]= pos[0];
  return npos;
}

IPosition readIPosition (const String& in)
{
  if (in.empty()) {
    return IPosition();
  }
  Vector<String> strs = stringToVector (in);
  IPosition ipos(strs.nelements());
  for (uInt i=0; i<strs.nelements(); ++i) {
    istringstream iss(strs[i]);
    iss >> ipos[i];
  }
  return ipos;
}

Quantity readQuantity (const String& in)
{
  Quantity res;
  if (!Quantity::read(res, in)) {
    throw AipsError (in + " is an illegal quantity");
  }
  return res;
}

MDirection readDirection (const String& in)
{
  Vector<String> vals = stringToVector(in);
  if (vals.size() > 3) {
    throw AipsError ("MDirection value " + in + " is invalid;"
		     " up to 3 values can be given");
  }
  MDirection::Types tp;
  if (! MDirection::getType (tp, vals[0])) {
    throw AipsError(vals[0] + " is an invalid MDirection type");
  }
  Quantity v0(0, "deg");
  Quantity v1(90, "deg");     // same default as in measures.g
  if (vals.size() > 1  &&  !vals[1].empty()) {
    v0 = readQuantity(vals[1]);
  }
  if (vals.size() > 2  &&  !vals[2].empty()) {
    v1 = readQuantity(vals[2]);
  }
  return MDirection(v0, v1, tp);
}

int main (Int argc, char** argv)
{
  try {
    Input inputs(1);
    // define the input structure
    inputs.version("20080818-GvD");
    inputs.create ("ms", "",
		   "Name of input MeasurementSet",
		   "string");
    inputs.create ("image", "",
		   "Name of output image file (default is <msname-stokes-mode-nchan>.img)",
		   "string");
    inputs.create ("fits", "no",
		   "Name of output image fits file ('no' means no fits file) empty is <imagename>.fits",
		   "string");
    inputs.create ("hdf5", "no",
    		   "Name of output image HDF5 file ('no' means no HDF5 file) empty is <imagename>.hdf5",
    		   "string");
    inputs.create ("model", "",
		   "Name of model image file (default is <imagename>.model",
		   "string");
    inputs.create ("restored", "",
		   "Name of restored image file (default is <imagename>.restored",
		   "string");
    inputs.create ("residual", "",
		   "Name of residual image file (default is <imagename>.residual",
		   "string");
    inputs.create ("data", "DATA",
		   "Name of DATA column to use",
		   "string");
    inputs.create ("mode", "mfs",
		   "Imaging mode (mfs, channel, or velocity)",
		   "string");
    inputs.create ("weight", "briggs",
		   "Weighting scheme (uniform, superuniform, natural, briggs (robust), or radial",
		   "string");
    inputs.create ("wprojplanes", "0",
		   "if >0 specifies nr of convolution functions to use in W-projection",
		   "int");
    inputs.create ("padding", "1.0",
		   "padding factor in image plane (>=1.0)",
		   "float");
    inputs.create ("cachesize", "512",
		   "maximum size of gridding cache (in MBytes)",
		   "int");
    inputs.create ("stokes", "I",
		   "Stokes parameters to image (e.g. IQUV)",
		   "string");
    inputs.create ("npix", "256",
		   "number of image pixels in x and y direction",
		   "int");
    inputs.create ("cellsize", "1arcsec",
		   "pixel width in x and y direction",
		   "quantity string");
    inputs.create ("phasecenter", "",
		   "phase center to be used (e.g. 'j2000, 05h30m, -30.2deg')",
		   "direction string");
    inputs.create ("field", "0",
		   "field id to be used",
		   "int");
    inputs.create ("spwid", "0",
		   "spectral window id to be used",
		   "int");
    inputs.create ("chanmode", "channel",
		   "frequency channel mode",
		   "string");
    inputs.create ("nchan", "1",
		   "number of frequency channels in MS",
		   "int");
    inputs.create ("chanstart", "0",
		   "first frequency channel in MS (0-relative)",
		   "int");
    inputs.create ("chanstep", "1",
		   "frequency channel step in MS",
		   "int");
    inputs.create ("img_nchan", "1",
		   "number of frequency channels in image",
		   "int");
    inputs.create ("img_chanstart", "0",
		   "first frequency channel in image (0-relative)",
		   "int");
    inputs.create ("img_chanstep", "1",
		   "frequency channel step in image",
		   "int");
    inputs.create ("select", "",
		   "TaQL selection string for MS",
		   "string");
    inputs.create ("operation", "image",
		   "Operation (image,clark,hogbom)",
		   "string");
    inputs.create ("niter", "1000",
		   "Number of clean iterations",
		   "int");
    inputs.create ("gain", "0.1",
		   "Loop gain for cleaning",
		   "float");
    inputs.create ("threshold", "0Jy",
		   "Flux level at which to stop cleaning",
		   "quantity string");
    inputs.create ("fixed", "False",
		   "Keep clean model fixed",
		   "bool");
    inputs.create ("mask", "",
		   "Name of the mask to use in cleaning",
		   "string");
    inputs.create ("maskblc", "0,0",
		   "bottom-left corner of mask region",
		   "int vector");
    inputs.create ("masktrc", "image shape",
		   "top-right corner of mask region",
		   "int vector");
    inputs.create ("maskvalue", "-1.0",
		   "Value to store in mask region; if given, mask is created; if mask not exists, defaults to 1.0",
		   "float");

    // Fill the input structure from the command line.
    inputs.readArguments (argc, argv);

    // Get the input specification.
    Bool fixed       = inputs.getBool("fixed");
    Long cachesize   = inputs.getInt("cachesize");
    Int fieldid      = inputs.getInt("field");
    Int spwid        = inputs.getInt("spwid");
    Int npix         = inputs.getInt("npix");
    Int nchan        = inputs.getInt("nchan");
    Int chanstart    = inputs.getInt("chanstart");
    Int chanstep     = inputs.getInt("chanstep");
    Int img_nchan    = inputs.getInt("img_nchan");
    Int img_start    = inputs.getInt("img_chanstart");
    Int img_step     = inputs.getInt("img_chanstep");
    Int wplanes      = inputs.getInt("wprojplanes");
    Int niter        = inputs.getInt("niter");
    Double padding   = inputs.getDouble("padding");
    Double gain      = inputs.getDouble("gain");
    Double maskValue = inputs.getDouble("maskvalue");
    String mode      = inputs.getString("mode");
    String operation = inputs.getString("operation");
    String weight    = inputs.getString("weight");
    String stokes    = inputs.getString("stokes");
    String chanmode  = inputs.getString("chanmode");
    String cellsize  = inputs.getString("cellsize");
    String phasectr  = inputs.getString("phasecenter");
    String threshStr = inputs.getString("threshold");
    String msName    = inputs.getString("ms");
    String imgName   = inputs.getString("image");
    String fitsName  = inputs.getString("fits");
    String hdf5Name  = inputs.getString("hdf5");
    String modelName = inputs.getString("model");
    String restoName = inputs.getString("restored");
    String residName = inputs.getString("residual");
    String imageType = inputs.getString("data");
    String select    = inputs.getString("select");
    String maskName  = inputs.getString("mask");
    String mstrBlc   = inputs.getString("maskblc");
    String mstrTrc   = inputs.getString("masktrc");

    // Check and interpret input values.
    Quantity qcellsize = readQuantity (cellsize);
    if (msName.empty()) {
      throw AipsError("An MS name must be given like ms=test.ms");
    }
    imageType.downcase();
    if (imageType == "data") {
      imageType = "observed";
    } else if (imageType == "corrected_data") {
      imageType = "corrected";
    } else if (imageType == "model_data") {
      imageType = "model";
    } else if (imageType == "residual_data") {
      imageType = "residual";
    }
    if (select.empty()) {
      select = "ANTENNA1 != ANTENNA2";
    } else {
      select = '(' + select + ") && ANTENNA1 != ANTENNA2";
    }
    if (imgName.empty()) {
      imgName = msName;
      imgName.gsub (Regex("\\..*"), "");
      imgName.gsub (Regex(".*/"), "");
      imgName += '-' + stokes + '-' + mode + String::toString(img_nchan)
	+ ".img";
    }
    if (fitsName == "no") {
      fitsName = String();
    } else if (fitsName.empty()) {
      fitsName = imgName + ".fits";
    }
    if (hdf5Name == "no") {
      hdf5Name = String();
#ifdef HAVE_HDF5
    } else if hdf5Name.empty()) {
      hdf5Name = imgName + ".hdf5";
#else
    } else {
      std::cerr << "Cannot create image in HDF5 format"
		<< " (HDF5 support is not compiled in)" << std::endl;
      hdf5Name = String();
#endif
    }
    if (modelName.empty()) {
      modelName = imgName + ".model";
    }
    if (restoName.empty()) {
      restoName = imgName + ".restored";
    }
    if (residName.empty()) {
      residName = imgName + ".residual";
    }
    if (weight == "robust") {
      weight = "briggs";
    }
    bool doShift = False;
    MDirection phaseCenter;
    if (! phasectr.empty()) {
      doShift = True;
      phaseCenter = readDirection (phasectr);
    }
    operation.downcase();
    AlwaysAssertExit (operation=="image" || operation=="hogbom" || operation=="clark");
    IPosition maskBlc, maskTrc;
    Quantity  threshold;
    if (operation != "image") {
      maskBlc = readIPosition (mstrBlc);
      maskTrc = readIPosition (mstrTrc);
      threshold = readQuantity (threshStr);
    }

    // Set the various imager variables.
    // The non-parameterized values used are the defaults in imager.g.
    MeasurementSet ms(msName, Table::Update);
    Imager imager(ms);
    imager.setdata (chanmode,                       // mode
		    Vector<Int>(1, nchan),          // nchan
		    Vector<Int>(1, chanstart),      // start
		    Vector<Int>(1, chanstep),       // step
		    MRadialVelocity(),              // mStart
		    MRadialVelocity(),              // mStep
		    Vector<Int>(1, spwid),          // spectralwindowsids
		    Vector<Int>(1, fieldid),        // fieldids
		    select);                        // msSelect
    imager.setimage (npix,                          // nx
		     npix,                          // ny
		     qcellsize,                     // cellx
		     qcellsize,                     // celly
		     stokes,                        // stokes
		     doShift,                       // doShift
		     phaseCenter,                   // phaseCenter
		     Quantity(0, "arcsec"),         // shiftx
		     Quantity(0, "arcsec"),         // shifty
		     mode,                          // mode
		     img_nchan,                     // nchan
		     img_start,                     // start
		     img_step,                      // step
		     MRadialVelocity(),             // mStart
		     MRadialVelocity(),             // mStep
		     Vector<Int>(1, spwid),         // spectralwindowids
		     fieldid,                       // fieldid
		     1,                             // facets
		     Quantity(0, "m"),              // distance
		     5.0,                           // paStep
		     5e-2);                         // pbLimit
    if (weight != "default") {
      imager.weight (weight,                        // type
		     "none",                        // rmode
		     Quantity(0, "Jy"),             // noise
		     0.0,                           // robust
		     Quantity(0, "rad"),            // fieldofview
		     0);                            // npixels
    }
    String ftmachine("ft");
    if (wplanes > 0) {
      ftmachine = "wproject";
    }
    imager.setoptions(ftmachine,                    // ftmachine
		      cachesize*1024*(1024/8),      // cache
		      16,                           // tile
		      "SF",                         // gridfunction
		      MPosition(),                  // mLocation
		      padding,                      // padding
		      True,                         // usemodelcol
		      wplanes,                      // wprojplanes
		      "",                           // epJTableName
		      True,                         // applyPointingOffsets
		      True,                         // doPointingCorrection
		      "");                          // cfCacheDirName

    // Do the imaging.
    if (operation == "image") {
      imager.makeimage (imageType, imgName, "");

      // Convert to HDF5 if needed.
#ifdef HAVE_HDF5
      if (! hdf5Name.empty()) {
	PagedImage<float> pimg(imgName);
	HDF5Image<float>  himg(pimg.shape(), pimg.coordinates(), hdf5Name);
	himg.copyData (pimg);
	himg.setUnits     (pimg.units());
	himg.setImageInfo (pimg.imageInfo());
	himg.setMiscInfo  (pimg.miscInfo());
      }
#endif
      // Convert result to fits if needed.
      if (! fitsName.empty()) {
	String error;
	PagedImage<float> img(imgName);
	if (! ImageFITSConverter::ImageToFITS (error, img, fitsName)) {
	  throw AipsError(error);
	}
      }

//      // Delete PagedImage if HDF5 was used.
//      if (! hdf5Name.empty()) {
//        Table::deleteTable (imgName);
//      }

    } else {
    // Do the cleaning.
      if (! maskName.empty()) {
	if (maskValue >= 0) {
	  PagedImage<float> pimg(imgName);
	  maskBlc = handlePos (maskBlc, IPosition(pimg.ndim(), 0));
	  maskTrc = handlePos (maskTrc, pimg.shape() - 1);
	  imager.boxmask (maskName,
			  maskBlc.asVector(),
			  maskTrc.asVector(),
			  maskValue);
	}
      }
      imager.clean (operation,                      // algorithm,
		    niter,                          // niter
		    gain,                           // gain
		    threshold,                      // threshold
		    False,                          // displayProgress
		    Vector<String>(1, modelName),   // model
		    Vector<Bool>(1, fixed),         // fixed
		    "",                             // complist
		    Vector<String>(1, maskName),    // mask
		    Vector<String>(1, restoName),   // restored
		    Vector<String>(1, residName));  // residual
      // Convert result to fits if needed.
       
      if (! fitsName.empty()) {
	String error;
	PagedImage<float> img(restoName);
	if (! ImageFITSConverter::ImageToFITS (error, img, fitsName)) {
	  throw AipsError(error);
	}
      }
    }

  } catch (AipsError x) {
    cout << x.getMesg() << endl;
    return 1;
  } 
  cout << "lwimager normally ended" << endl;
  return 0;
}
