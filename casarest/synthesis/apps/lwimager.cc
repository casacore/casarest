//# lwimager.cc: Program to convert an image to FITS format
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
#include <images/Images/HDF5Image.h>
#endif
#include <images/Images/ImageFITSConverter.h>
#include <casa/Inputs.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Utilities/Regex.h>
#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>

using namespace casa;

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
    inputs.version("20080516-GvD");
    inputs.create ("ms", "",
		   "Name of input MeasurementSet",
		   "string");
    inputs.create ("image", "",
		   "Name of output image file (default is <msname-stokes-mode-nchan>.img)",
		   "string");
    inputs.create ("fits", "",
		   "Name of output fits file ('no' means no fits file) default is <imagename>.fits",
		   "string");
    inputs.create ("data", "DATA",
		   "Name of DATA column to use",
		   "string");
    inputs.create ("hdf5", "false",
    		   "Create the image in HDF5 format?"
    		   "bool");
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
		   "string");
    inputs.create ("phasecenter", "",
		   "phase center to be used (e.g. 'j2000, 05h30m, -30.2deg')",
		   "string");
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

    // Fill the input structure from the command line.
    inputs.readArguments (argc, argv);

    // Get the input specification.
    Bool hdf5        = inputs.getBool("hdf5");
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
    Double padding   = inputs.getDouble("padding");
    String mode      = inputs.getString("mode");
    String weight    = inputs.getString("weight");
    String stokes    = inputs.getString("stokes");
    String chanmode  = inputs.getString("chanmode");
    String cellsize  = inputs.getString("cellsize");
    String phasectr  = inputs.getString("phasecenter");
    String msName    = inputs.getString("ms");
    String imgName   = inputs.getString("image");
    String fitsName  = inputs.getString("fits");
    String imageType = inputs.getString("data");
    String select    = inputs.getString("select");

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
    if (fitsName.empty()) {
      fitsName = imgName;
      fitsName.gsub (Regex("\\..*"), ".fits");
    }
    String hdf5Name;
    if (hdf5) {
#ifndef HAVE_HDF5
      std::cerr << "Cannot create image in HDF5 format"
		<< " (HDF5 support is not compiled in)" << std::endl;
      hdf5 = False;
#else
      hdf5Name = imgName;
      imgName += "_tmp";
#endif
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

    std::cerr << phaseCenter<< std::endl;

    // Do the imaging.
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
    imager.makeimage (imageType, imgName, "");

    // Convert to HDF5 if needed.
#ifdef HAVE_HDF5
    if (hdf5) {
      PagedImage<float> pimg(imgName);
      HDF5Image<float>  himg(pimg.shape(), pimg.coordinates(), hdf5Name);
      himg.copyData (pimg);
      himg.setUnits     (pimg.units());
      himg.setImageInfo (pimg.imageInfo());
      himg.setMiscInfo  (pimg.miscInfo());
    }
#endif

    // Convert to fits if needed.
    if (fitsName != "no") {
      String error;
      PagedImage<float> img(imgName);
      if (! ImageFITSConverter::ImageToFITS (error, img, fitsName)) {
	throw AipsError(error);
      }
    }

    // Delete PagedImage if HDF5 was used.
    if (hdf5) {
      Table::deleteTable (imgName);
    }

  } catch (AipsError x) {
    cout << x.getMesg() << endl;
    return 1;
  } 
  cout << "lwimager normally ended" << endl;
  return 0;
}
