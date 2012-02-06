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
#include <images/Images/HDF5Image.h>
#include <images/Images/ImageFITSConverter.h>
#include <casa/Inputs.h>
#include <casa/Arrays/ArrayUtil.h>
#include <casa/Arrays/ArrayMath.h>
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

void readFilter (const String& filter,
                 Quantity& bmajor, Quantity& bminor, Quantity& bpa)
{
  if (filter.empty()) {
    return;
  }
  Vector<String> strs = stringToVector(filter);
  if (strs.size() != 3) {
    throw AipsError("Specify gaussian tapering filter as bmajor,bminor,bpa");
  }
  if (! strs[0].empty()) {
    bmajor = readQuantity (strs[0]);
  }
  if (! strs[1].empty()) {
    bminor = readQuantity (strs[1]);
  }
  if (! strs[2].empty()) {
    bpa = readQuantity (strs[2]);
  }
}

void readCleanBeam (const String& clean_beam,
                 Quantity& c_bmaj, Quantity& c_bmin, Quantity& c_bpa)
{
  if (clean_beam.empty()) {
    return;
  }
  Vector<String> strs = stringToVector(clean_beam);
  if (strs.size() != 3) {
    throw AipsError("Specify clean restoring beam bmajor,bminor,bpa");
  }
  if (! strs[0].empty()) {
    c_bmaj = readQuantity (strs[0]);
  }
  if (! strs[1].empty()) {
    c_bmin = readQuantity (strs[1]);
  }
  if (! strs[2].empty()) {
    c_bpa = readQuantity (strs[2]);
  }
// c_bmaj.setUnit("arcsec");
// c_bmin.setUnit("arcsec");
// c_bpa.setUnit("deg");
}

void makeEmpty (Imager& imager, const String& imgName, Int fieldid)
{
  CoordinateSystem coords;
  AlwaysAssert (imager.imagecoordinates(coords), AipsError);
  String name(imgName);
  imager.makeEmptyImage(coords, name, fieldid);
  imager.unlock();
}

int main (Int argc, char** argv)
{
  try {
    Input inputs(1);
    // define the input structure
    inputs.version("20110628-GvD");
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
    inputs.create ("prior", "",
		   "Name of prior image file (default is <imagename>.prior)",
		   "string");
    inputs.create ("model", "",
		   "Name of model image file (default is <imagename>.model)",
		   "string");
    inputs.create ("restored", "",
		   "Name of restored image file (default is <imagename>.restored)",
		   "string");
    inputs.create ("residual", "",
		   "Name of residual image file (default is <imagename>.residual)",
		   "string");
    inputs.create ("data", "DATA",
		   "Name of DATA column to use",
		   "string");
    inputs.create ("mode", "mfs",
		   "Imaging mode (mfs, channel, or velocity)",
		   "string");
    inputs.create ("filter", "",
                   "Apply gaussian tapering filter; specify as major,minor,pa",
                   "string");
    inputs.create ("clean_beam", "",
                   "Specify clean restoring beam; specify as major axis,minor axis,position angle with units",
                   "string");
    inputs.create ("weight", "briggs",
		   "Weighting scheme (uniform, superuniform, natural, briggs (robust), briggsabs, or radial)",
		   "string");
    inputs.create ("noise", "1.0",
		   "Noise (in Jy) for briggsabs weighting"
		   "float");
    inputs.create ("robust", "0.0",
		   "Robust parameter",
		   "float");
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
    inputs.create ("nfacets", "1",
                   "number of facets in x or y",
                   "int");
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
		   "spectral window id(s) to be used",
		   "int vector");
    inputs.create ("chanmode", "channel",
		   "frequency channel mode",
		   "string");
    inputs.create ("nchan", "1",
		   "number of frequency channels to select from each spectral window (one number per spw)",
		   "int vector");
    inputs.create ("chanstart", "0",
		   "first frequency channel per each spw (0-relative)",
		   "int vector");
    inputs.create ("chanstep", "1",
		   "frequency channel step per each spw",
		   "int vector");
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
		   "Operation (empty,image,clark,hogbom,csclean,multiscale,entropy)",
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
    inputs.create ("targetflux", "1.0Jy",
		   "Target flux for maximum entropy",
		   "quantity string");
    inputs.create ("sigma", "0.001Jy",
		   "deviation for maximum entropy",
		   "quantity string");
    inputs.create ("fixed", "False",
		   "Keep clean model fixed",
		   "bool");
    inputs.create ("constrainflux", "False",
		   "Constrain image to match target flux? For max entropy",
		   "bool");
    inputs.create ("prefervelocity", "True",
                   "Should FITS image spectral axis be velocity or frequency",
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
    inputs.create ("nscales", "5",
                   "Scales for MultiScale Clean",
                   "int");
    inputs.create ("uservector", "0",
		   "user-defined scales for MultiScale clean",
		   "float vector");
    inputs.create ("maskvalue", "-1.0",
		   "Value to store in mask region; if given, mask is created; if mask not exists, defaults to 1.0",
		   "float");
    inputs.create ("cyclefactor", "1.5",
                   "multi-field deconvolution parameter; see casapy's imager.setmfcontrol",
                   "float");
    inputs.create ("cyclespeedup", "-1",
                   "multi-field deconvolution parameter; see casapy's imager.setmfcontrol",
                   "float");
    inputs.create ("cyclemaxpsffraction", "0.8",
                   "multi-field deconvolution parameter; see casapy's imager.setmfcontrol",
                   "float");
    inputs.create ("stoplargenegatives", "2",
                   "multi-field deconvolution parameter; see casapy's imager.setmfcontrol",
                   "int");
    inputs.create ("stoppointmode", "-1",
                   "multi-field deconvolution parameter; see casapy's imager.setmfcontrol",
                   "setmfcontrol parameter; see casapy",
                   "int");

    // Fill the input structure from the command line.
    inputs.readArguments (argc, argv);

    // Get the input specification.
    Bool fixed       = inputs.getBool("fixed");
    Bool constrainFlux  = inputs.getBool("constrainflux");
    Bool preferVelocity = inputs.getBool("prefervelocity");
    Long cachesize   = inputs.getInt("cachesize");
    Int fieldid      = inputs.getInt("field");
    Vector<Int> spwid(inputs.getIntArray("spwid"));
    Int npix         = inputs.getInt("npix");
    Int nfacet       = inputs.getInt("nfacets");
    Vector<Int> nchan(inputs.getIntArray("nchan"));
    Vector<Int> chanstart(inputs.getIntArray("chanstart"));
    Vector<Int> chanstep(inputs.getIntArray("chanstep"));
    Int img_nchan    = inputs.getInt("img_nchan");
    Int img_start    = inputs.getInt("img_chanstart");
    Int img_step     = inputs.getInt("img_chanstep");
    Int wplanes      = inputs.getInt("wprojplanes");
    Int niter        = inputs.getInt("niter");
    Int nscales      = inputs.getInt("nscales");
    Vector<Double> userScaleSizes(inputs.getDoubleArray("uservector"));
    Double padding   = inputs.getDouble("padding");
    Double gain      = inputs.getDouble("gain");
    Double maskValue = inputs.getDouble("maskvalue");
    String mode      = inputs.getString("mode");
    String operation = inputs.getString("operation");
    String weight    = inputs.getString("weight");
    Double noise     = inputs.getDouble("noise");
    Double robust    = inputs.getDouble("robust");
    String filter    = inputs.getString("filter");
    String clean_beam = inputs.getString("clean_beam");
    String stokes    = inputs.getString("stokes");
    String chanmode  = inputs.getString("chanmode");
    String cellsize  = inputs.getString("cellsize");
    String phasectr  = inputs.getString("phasecenter");
    String sigmaStr  = inputs.getString("sigma");
    String targetStr = inputs.getString("targetflux");
    String threshStr = inputs.getString("threshold");
    String msName    = inputs.getString("ms");
    String imgName   = inputs.getString("image");
    String fitsName  = inputs.getString("fits");
    String hdf5Name  = inputs.getString("hdf5");
    String modelName = inputs.getString("model");
    String priorName = inputs.getString("prior");
    String restoName = inputs.getString("restored");
    String residName = inputs.getString("residual");
    String imageType = inputs.getString("data");
    String select    = inputs.getString("select");
    String maskName  = inputs.getString("mask");
    String mstrBlc   = inputs.getString("maskblc");
    String mstrTrc   = inputs.getString("masktrc");
    Double cycleFactor   = inputs.getDouble ("cyclefactor");
    Double cycleSpeedup  = inputs.getDouble ("cyclespeedup");
    Double cycleMaxPsfFr = inputs.getDouble ("cyclemaxpsffraction");
    Int    stopLargeNeg  = inputs.getInt    ("stoplargenegatives");
    Int    stopPointMode = inputs.getInt    ("stoppointmode");

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
    }
    if (priorName.empty()) {
      priorName = imgName + ".prior";
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
    } else if (weight == "robustabs") {
      weight = "briggsabs";
    }
    string rmode = "norm";
    if (weight == "briggsabs") {
      weight = "briggs";
      rmode  = "abs";
    }
    bool doShift = False;
    MDirection phaseCenter;
    if (! phasectr.empty()) {
      doShift = True;
      phaseCenter = readDirection (phasectr);
    }
    operation.downcase();
    AlwaysAssertExit (operation=="empty" || operation=="image" || operation=="hogbom" || operation=="clark" || operation=="csclean" || operation=="multiscale" || operation =="entropy");
    IPosition maskBlc, maskTrc;
    Quantity threshold;
    Quantity sigma;
    Quantity targetFlux;
    Bool doClean = (operation != "empty"  &&  operation != "image");
    if (doClean) {
      maskBlc = readIPosition (mstrBlc);
      maskTrc = readIPosition (mstrTrc);
      threshold = readQuantity (threshStr);
      sigma = readQuantity (sigmaStr);
      targetFlux = readQuantity (targetStr);
    }
    // Get axis specification from filter.
    Quantity bmajor, bminor, bpa;
    readFilter (filter, bmajor, bminor, bpa);

    // Get restoring beam specification from clean_beam.
    Quantity c_bmaj, c_bmin, c_bpa;
    readCleanBeam (clean_beam, c_bmaj, c_bmin, c_bpa);

    // Set the various imager variables.
    // The non-parameterized values used are the defaults in imager.g.
    MeasurementSet ms(msName, Table::Update);
    Imager imager(ms);
    // Use channel mode if only one data channel per image-channel.
    if (nchan.size() == 1  &&  nchan[0] == img_nchan) {
      mode = "channel";
    }
    imager.setdata (chanmode,                       // mode
		    nchan,
		    chanstart,
                    chanstep,
		    MRadialVelocity(),              // mStart
		    MRadialVelocity(),              // mStep
		    spwid,
		    Vector<Int>(1,fieldid),
		    select,                         // msSelect
                    String(),                       // timerng
                    String(),                       // fieldnames
                    Vector<Int>(),                  // antIndex
                    String(),                       // antnames
                    String(),                       // spwstring
                    String(),                       // uvdist
                    String(),                       // scan
                    True);                          // useModelCol

    imager.defineImage (npix,                       // nx
                        npix,                       // ny
                        qcellsize,                  // cellx
                        qcellsize,                  // celly
                        stokes,                     // stokes
                        phaseCenter,                // phaseCenter
                        doShift  ?  -1 : fieldid,   // fieldid
                        mode,                       // mode
                        img_nchan,                  // nchan
                        img_start,                  // start
                        img_step,                   // step
                        MFrequency(),               // mFreqstart
                        MRadialVelocity(),          // mStart
                        Quantity(1,"km/s"),         // qstep, Def=1 km/s
                        spwid,                      // spectralwindowids
                        nfacet);                    // facets

    // Create empty image?
    if (operation == "empty" ) {
      makeEmpty (imager, imgName, fieldid);
    } else {

      // Define weighting.
      if (weight != "default") {
        imager.weight (weight,                      // type
                       rmode,                       // rmode
                       Quantity(noise, "Jy"),       // briggsabs noise
                       robust,                      // robust
                       Quantity(0, "rad"),          // fieldofview
                       0);                          // npixels
      }

      // If multiscale, set its parameters.
      if (operation == "multiscale") {
        String scaleMethod;
        Vector<Float> userVector(userScaleSizes.shape());
        convertArray (userVector, userScaleSizes);
        if (userScaleSizes.size() > 1) {
          scaleMethod = "uservector";
        } else {
          scaleMethod = "nscales";
        }
        imager.setscales(scaleMethod, nscales, userVector);
      }
      if (! filter.empty()) {
        imager.filter ("gaussian", bmajor, bminor, bpa);
      }
      if (! clean_beam.empty()) {
        imager.setbeam (c_bmaj, c_bmin, c_bpa);
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
                        wplanes);                     // wprojplanes
      // Do the imaging.
      if (operation == "image") {
        imager.makeimage (imageType, imgName);

        // Convert result to fits if needed.
        if (! fitsName.empty()) {
          String error;
          PagedImage<float> img(imgName);
          if (! ImageFITSConverter::ImageToFITS (error,
                                                 img,
                                                 fitsName,
                                                 64,         // memoryInMB
                                                 preferVelocity)) {
            throw AipsError(error);
          }
        }

        // Convert to HDF5 if needed.
        if (! hdf5Name.empty()) {
          PagedImage<float> pimg(imgName);
          HDF5Image<float>  himg(pimg.shape(), pimg.coordinates(), hdf5Name);
          himg.copyData (pimg);
          himg.setUnits     (pimg.units());
          himg.setImageInfo (pimg.imageInfo());
          himg.setMiscInfo  (pimg.miscInfo());
          // Delete PagedImage if HDF5 is used.
          Table::deleteTable (imgName);
        }

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
        imager.setmfcontrol (cycleFactor,
                             cycleSpeedup,
                             cycleMaxPsfFr,
                             stopLargeNeg, 
                             stopPointMode,
                             "SAULT",                 // scaleType
                             0.1,                     // minPB
                             0.4,                     // constPB
                             Vector<String>(),        // fluxscale
                             True);                   // flatnoise
        if (operation == "entropy") {
          imager.mem(operation,                       // algorithm
                     niter,                           // niter
                     sigma,                           // sigma
                     targetFlux,                      // targetflux
                     constrainFlux,                   // constrainflux
                     False,                           // displayProgress
                     Vector<String>(1, modelName),    // model
                     Vector<Bool>(1, fixed),          // fixed
                     "",                              // complist
                     Vector<String>(1, priorName),    // prior
                     Vector<String>(1, maskName),     // mask
                     Vector<String>(1, restoName),    // restored
                     Vector<String>(1, residName));   // residual

        } else {
          imager.clean(operation,                     // algorithm,
                       niter,                         // niter
                       gain,                          // gain
                       threshold,                     // threshold
                       False,                         // displayProgress
                       Vector<String>(1, modelName),  // model
                       Vector<Bool>(1, fixed),        // fixed
                       "",                            // complist
                       Vector<String>(1, maskName),   // mask
                       Vector<String>(1, restoName),  // restored
                       Vector<String>(1, residName)); // residual
        }
        // Convert result to fits if needed.
        if (! fitsName.empty()) {
          String error;
          PagedImage<float> img(restoName);
          if (! ImageFITSConverter::ImageToFITS (error,
                                                 img,
                                                 fitsName,
                                                 64,         // memoryInMB
                                                 preferVelocity)) {
            throw AipsError(error);
          }
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
