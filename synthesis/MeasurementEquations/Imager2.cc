//All helper functions of imager moved here for readability
//# Imager.cc: Implementation of Imager.h
//# Copyright (C) 1997-2008
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
//# $Id$


#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/iostream.h>
#include <synthesis/MeasurementEquations/Imager.h>
#include <synthesis/MeasurementComponents/EPJones.h>

#include <casacore/ms/MeasurementSets/MSHistoryHandler.h>

#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <casacore/casa/Logging.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogMessage.h>

#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/OS/Path.h>

#include <casacore/casa/OS/HostInfo.h>
#include <casacore/tables/Tables/RefRows.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/TaQL/TableParse.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/TableLock.h>
#include <casacore/tables/TaQL/ExprNode.h>

#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Utilities/Fallible.h>
#include <casacore/casa/Utilities/CompositeNumber.h>

#include <casacore/casa/BasicSL/Constants.h>

#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogMessage.h>

#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Slice.h>
//#include <casacore/images/Images/ImageAnalysis.h>
#include <casacore/images/Images/ImageExpr.h>
//#include <casacore/images/Images/ImagePolarimetry.h>
#include <synthesis/MeasurementEquations/ClarkCleanProgress.h>
#if defined(casacore)
#include <casacore/lattices/LatticeMath/LatticeCleanProgress.h>
#else
#include <casacore/lattices/LatticeMath/LatticeCleanProgress.h>
#endif
#include <msvis/MSVis/VisSet.h>
#include <msvis/MSVis/VisSetUtil.h>
#include <msvis/MSVis/VisImagingWeight.h>
#include <msvis/MSVis/SubMS.h>
// Disabling Imager::correct() (gmoellen 06Nov20)
//#include <synthesis/MeasurementComponents/TimeVarVisJones.h>

#include <casacore/measures/Measures/Stokes.h>
#include <casacore/casa/Quanta/UnitMap.h>
#include <casacore/casa/Quanta/UnitVal.h>
#include <casacore/casa/Quanta/MVAngle.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MeasTable.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
#include <casacore/ms/MeasurementSets/MSDopplerUtil.h>
#if defined(casacore)
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSDataDescIndex.h>
#include <casacore/ms/MSSel/MSSourceIndex.h>
#include <casacore/ms/MSOper/MSSummary.h>
#else
#include <casacore/ms/MSSel/MSSelection.h>
#include <casacore/ms/MSSel/MSDataDescIndex.h>
#include <casacore/ms/MSSel/MSSourceIndex.h>
#include <casacore/ms/MSOper/MSSummary.h>
#endif
#include <synthesis/MeasurementEquations/CubeSkyEquation.h>
// Disabling Imager::correct() (gmoellen 06Nov20)
//#include <synthesis/MeasurementEquations/VisEquation.h>
#include <synthesis/MeasurementComponents/ImageSkyModel.h>
#include <synthesis/MeasurementComponents/CEMemImageSkyModel.h>
#include <synthesis/MeasurementComponents/MFCEMemImageSkyModel.h>
#include <synthesis/MeasurementComponents/MFCleanImageSkyModel.h>
#include <synthesis/MeasurementComponents/CSCleanImageSkyModel.h>
#include <synthesis/MeasurementComponents/MFMSCleanImageSkyModel.h>
#include <synthesis/MeasurementComponents/HogbomCleanImageSkyModel.h>
#include <synthesis/MeasurementComponents/MSCleanImageSkyModel.h>
#include <synthesis/MeasurementComponents/NNLSImageSkyModel.h>
#include <synthesis/MeasurementComponents/WBCleanImageSkyModel.h>
#include <synthesis/MeasurementComponents/MultiThreadedVisResampler.h>
#include <synthesis/MeasurementComponents/GridBoth.h>
#include <synthesis/MeasurementComponents/rGridFT.h>
#include <synthesis/MeasurementComponents/MosaicFT.h>
#include <synthesis/MeasurementComponents/WProjectFT.h>
#include <synthesis/MeasurementComponents/nPBWProjectFT.h>
#include <synthesis/MeasurementComponents/PBMosaicFT.h>
#include <synthesis/MeasurementComponents/PBMath.h>
#include <synthesis/MeasurementComponents/SimpleComponentFTMachine.h>
#include <synthesis/MeasurementComponents/SimpCompGridMachine.h>
#include <synthesis/MeasurementComponents/VPSkyJones.h>
#include <synthesis/MeasurementComponents/SynthesisError.h>
#include <synthesis/MeasurementComponents/HetArrayConvFunc.h>
#include <synthesis/MeasurementComponents/VisibilityResamplerBase.h>

#include <synthesis/DataSampling/SynDataSampling.h>
#include <synthesis/DataSampling/SDDataSampling.h>
#include <synthesis/DataSampling/ImageDataSampling.h>
#include <synthesis/DataSampling/PixonProcessor.h>

#include <synthesis/MeasurementEquations/StokesImageUtil.h>
#include <casacore/lattices/Lattices/TiledLineStepper.h> 
#include <casacore/lattices/Lattices/LatticeIterator.h> 
#if defined(casacore)
#include <casacore/lattices/LEL/LatticeExpr.h> 
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#include <casacore/lattices/LRegions/LattRegionHolder.h>
#include <casacore/lattices/LRegions/LCEllipsoid.h>
#include <casacore/lattices/LRegions/LCRegion.h>
#include <casacore/lattices/LRegions/LCBox.h>
#include <casacore/lattices/LRegions/LCIntersection.h>
#include <casacore/lattices/LRegions/LCUnion.h>
#include <casacore/lattices/LRegions/LCExtension.h>
#else
#include <casacore/lattices/LEL/LatticeExpr.h> 
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#include <casacore/lattices/LRegions/LattRegionHolder.h>
#include <casacore/lattices/LRegions/LCEllipsoid.h>
#include <casacore/lattices/LRegions/LCRegion.h>
#include <casacore/lattices/LRegions/LCBox.h>
#include <casacore/lattices/LRegions/LCIntersection.h>
#include <casacore/lattices/LRegions/LCUnion.h>
#include <casacore/lattices/LRegions/LCExtension.h>
#endif

#include <casacore/images/Images/ImageRegrid.h>
#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/images/Regions/RegionManager.h>
#include <casacore/images/Regions/WCBox.h>
#include <casacore/images/Regions/WCUnion.h>
#include <casacore/images/Regions/WCIntersection.h>
#include <synthesis/MeasurementComponents/PBMath.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/ImageInfo.h>
#include <casacore/images/Images/SubImage.h>
///#include <casacore/images/Images/ImageMetaData.h>
#include <casacore/images/Images/ImageUtilities.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/coordinates/Coordinates/Projection.h>
#include <casacore/coordinates/Coordinates/ObsInfo.h>

#include <components/ComponentModels/ComponentList.h>
#include <components/ComponentModels/ConstantSpectrum.h>
#include <components/ComponentModels/TabularSpectrum.h>
#include <components/ComponentModels/Flux.h>
#include <components/ComponentModels/FluxStandard.h>
#include <components/ComponentModels/PointShape.h>
#include <components/ComponentModels/DiskShape.h>

#include <casacore/casa/OS/HostInfo.h>

#include <components/ComponentModels/ComponentList.h>

#include <casacore/measures/Measures/UVWMachine.h>

#include <casacore/casa/sstream.h>

#include <sys/types.h>
#include <unistd.h>
#include <limits>

#include <synthesis/MeasurementComponents/AWProjectFT.h>
#include <synthesis/MeasurementComponents/AWProjectWBFT.h>
#include <synthesis/MeasurementComponents/MultiTermFT.h>
#include <synthesis/MeasurementComponents/AWConvFunc.h>

using namespace std;


#ifdef PABLO_IO
#include "PabloTrace.h"
#endif


namespace casacore { //# NAMESPACE CASACORE - BEGIN

String Imager::imageName()
{
  LogIO os(LogOrigin("imager", "imageName()", WHERE));
  try {
    lock();
    String name(msname_p);
    ROMSColumns msc(*ms_p);
    if(datafieldids_p.shape() !=0) {
      name=msc.field().name()(datafieldids_p(0));
    }
    else if(fieldid_p > -1) {
       name=msc.field().name()(fieldid_p);
    }
    unlock();
    return name;
  } catch (AipsError x) {
    unlock();
    os << LogIO::SEVERE << "Caught Exception: "<< x.getMesg() << LogIO::EXCEPTION; 
    return "";
  } 
  return String("imagerImage");
}

// imagecoordinates2 (use subMS method to get freq vectors)
// Make standard choices for coordinates
Bool Imager::imagecoordinates2(CoordinateSystem& coordInfo, const Bool verbose) 
{  
  Bool rstate;
  if(!valid()) return False;
  if(!assertDefinedImageParameters()) return False;
  LogIO os(LogOrigin("Imager", "imagecoordinates()", WHERE));
  
  //===Adjust for some setimage defaults if imageNchan=-1 and spwids=-1
  if(imageNchan_p <=0){
    imageNchan_p=1;
  }
  if(spectralwindowids_p.nelements()==1){
    if(spectralwindowids_p[0]<0){
      spectralwindowids_p.resize();
      if(dataspectralwindowids_p.nelements()==0){
         Int nspwinms=ms_p->spectralWindow().nrow();
         dataspectralwindowids_p.resize(nspwinms);
         indgen(dataspectralwindowids_p);
      }
      spectralwindowids_p=dataspectralwindowids_p;
    }
    
  }
  if(fieldid_p < 0){
     if(datafieldids_p.shape() !=0) {
       fieldid_p=datafieldids_p[0];
     }
     else{
       fieldid_p=0; //best default if nothing is specified
     }
  }
  //===end of default
  Vector<Double> deltas(2);
  deltas(0)=-mcellx_p.get("rad").getValue();
  deltas(1)=mcelly_p.get("rad").getValue();
  
  ROMSColumns msc(*ms_p);
  MFrequency::Types obsFreqRef=MFrequency::DEFAULT;
  ScalarColumn<Int> measFreqRef(ms_p->spectralWindow(),
				  MSSpectralWindow::columnName(MSSpectralWindow::MEAS_FREQ_REF));
  //using the first frame of reference; TO DO should do the right thing 
  //for different frames selected. 
  //Int eh = spectralwindowids_p(0);
  if(spectralwindowids_p.size() && measFreqRef(spectralwindowids_p(0)) >=0) {
     obsFreqRef=(MFrequency::Types)measFreqRef(spectralwindowids_p(0));
  }
			    

  MVDirection mvPhaseCenter(phaseCenter_p.getAngle());
  // Normalize correctly
  MVAngle ra=mvPhaseCenter.get()(0);
  ra(0.0);
  MVAngle dec=mvPhaseCenter.get()(1);
  Vector<Double> refCoord(2);
  refCoord(0)=ra.get().getValue();    
  refCoord(1)=dec;    
  
  Vector<Double> refPixel(2); 
  refPixel(0) = Double(nx_p / 2);
  refPixel(1) = Double(ny_p / 2);
  
  //defining observatory...needed for position on earth
  String telescop = msc.observation().telescopeName()(0);

  // defining epoch as begining time from timerange in OBSERVATION subtable
  // Using first observation for now
  //MEpoch obsEpoch = msc.observation().timeRangeMeas()(0)(IPosition(1,0));
  // modified to use main table's TIME column for better match with what
  // VisIter does.
  MEpoch obsEpoch = msc.timeMeas()(0);

  //Now finding the position of the telescope on Earth...needed for proper
  //frequency conversions

  MPosition obsPosition;
  if(!(MeasTable::Observatory(obsPosition, telescop))){
    os << LogIO::WARN << "Did not get the position of " << telescop 
       << " from data repository" << LogIO::POST;
    os << LogIO::WARN 
       << "Please contact CASA to add it to the repository."
       << LogIO::POST;
    os << LogIO::WARN << "Frequency conversion will not work " << LogIO::POST;
    freqFrameValid_p = False;
  }
  else{
    mLocation_p = obsPosition;
    freqFrameValid_p = True;
  }
  // Now find the projection to use: could probably also use
  // max(abs(w))=0.0 as a criterion
  Projection projection(Projection::SIN);
  if(telescop == "ATCASCP" || telescop == "WSRT" || telescop == "DRAO") {
    os << LogIO::NORMAL // Loglevel NORMAL
       << "Using SIN image projection adjusted for "
       << (telescop == "ATCASCP" ? 'S' : 'N') << "CP" 
       << LogIO::POST;
    Vector<Double> projectionParameters(2);
    projectionParameters(0) = 0.0;
    if(sin(dec) != 0.0){
      projectionParameters(1) = cos(dec)/sin(dec);
      projection=Projection(Projection::SIN, projectionParameters);
    }
    else {
      os << LogIO::WARN
         << "Singular projection for " << telescop << ": using plain SIN"
         << LogIO::POST;
      projection=Projection(Projection::SIN);
    }
  }
  else {
    os << LogIO::DEBUGGING << "Using SIN image projection" << LogIO::POST;
  }
  os << LogIO::NORMAL;
  
  Matrix<Double> xform(2,2);
  xform=0.0;xform.diagonal()=1.0;
  DirectionCoordinate
    myRaDec(MDirection::Types(phaseCenter_p.getRefPtr()->getType()),
	    projection,
	    refCoord(0), refCoord(1),
	    deltas(0), deltas(1),
	    xform,
	    refPixel(0), refPixel(1));
  
  // Now set up spectral coordinate
  SpectralCoordinate* mySpectral=0;
  Double refChan=0.0;
  
  // Spectral synthesis
  // For mfs band we set the window to include all spectral windows
  Int nspw=spectralwindowids_p.nelements();
  if (imageMode_p=="MFS") {
    Double fmin=C::dbl_max;
    Double fmax=-(C::dbl_max);
    Double fmean=0.0;
    for (Int i=0;i<nspw;++i) {
      Int spw=spectralwindowids_p(i);
      Vector<Double> chanFreq=msc.spectralWindow().chanFreq()(spw); 
      Vector<Double> freqResolution = msc.spectralWindow().chanWidth()(spw); 
      
      if(dataMode_p=="none"){
      
	if(i==0) {
	  fmin=min(chanFreq-abs(0.5*freqResolution));
	  fmax=max(chanFreq+abs(0.5*freqResolution));
	}
	else {
	  fmin=min(fmin,min(chanFreq-abs(0.5*freqResolution)));
	  fmax=max(fmax,max(chanFreq+abs(0.5*freqResolution)));
	}
      }
      else if(dataMode_p=="channel"){
	// This needs some careful thought about roundoff - it is likely 
	// still adding an extra half-channel at top and bottom but 
	// if the freqResolution is nonlinear, there are subtleties
	Int elnchan=chanFreq.nelements();
	Int firstchan=0;
        Int elstep=1;
	for (uInt jj=0; jj < dataspectralwindowids_p.nelements(); ++jj){
	  if(dataspectralwindowids_p[jj]==spw){
	    firstchan=dataStart_p[jj];
	    elnchan=dataNchan_p[jj];
	    elstep=dataStep_p[jj];
	  }	
	}
	Int lastchan=firstchan+ elnchan*elstep;
        for(Int k=firstchan ; k < lastchan ;  k+=elstep){
	  fmin=min(fmin,chanFreq[k]-abs(freqResolution[k]*(elstep-0.5)));
	  fmax=max(fmax,chanFreq[k]+abs(freqResolution[k]*(elstep-0.5)));
        }
      }
      else{
	this->unlock();
	os << LogIO::SEVERE 
	   << "setdata has to be in 'channel' or 'none' mode for 'mfs' imaging to work"
	   << LogIO::EXCEPTION;
      return False;
      }
 
    }

    fmean=(fmax+fmin)/2.0;
    Vector<Double> restFreqArray;
    Double restFreq=fmean;
    if(getRestFreq(restFreqArray, spectralwindowids_p(0))){
      restFreq=restFreqArray[0];
    }
    imageNchan_p=1;
    Double finc=(fmax-fmin); 
    mySpectral = new SpectralCoordinate(freqFrame_p,  fmean//-finc/2.0
					, finc,
      					refChan, restFreq);
    os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3) // Loglevel INFO
       << "Center frequency = "
       << MFrequency(Quantity(fmean, "Hz")).get("GHz").getValue()
       << " GHz, synthesized continuum bandwidth = "
       << MFrequency(Quantity(finc, "Hz")).get("GHz").getValue()
       << " GHz" << LogIO::POST;

    if(ntaylor_p>1 && reffreq_p==0.0) 
    {
	    reffreq_p = fmean;
	    os << "Setting center frequency as MS-MFS reference frequency" << LogIO::POST;
    }
  }
  
  else if(imageMode_p.contains("FREQ")) {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      Double restFreq=mfImageStart_p.get("Hz").getValue();
      Vector<Double> restFreqVec;
      if(getRestFreq(restFreqVec, spectralwindowids_p(0))){
	restFreq=restFreqVec[0];
      }
      MFrequency::Types mfreqref=(obsFreqRef==(MFrequency::REST)) ? MFrequency::REST : MFrequency::castType(mfImageStart_p.getRef().getType()) ; 

      Vector<Double> imgridfreqs;
      Vector<Double> imfreqres;
      //rstate=calcImFreqs(imgridfreqs, imfreqres, mfreqref, obsEpoch, obsPosition,restFreq);
      // should use obsFreqRef
      rstate=calcImFreqs(imgridfreqs, imfreqres, obsFreqRef, obsEpoch, obsPosition,restFreq);
      //cerr<<"imfreqres(0)="<<imfreqres(0)<<endl;


      if (imageNchan_p==1) {
        mySpectral = new SpectralCoordinate(mfreqref,
      					  mfImageStart_p.get("Hz").getValue(),
      					  mfImageStep_p.get("Hz").getValue(),
      					  refChan, restFreq);
      }
      else {
        Double finc= imgridfreqs(1)-imgridfreqs(0); 
        mySpectral = new SpectralCoordinate(mfreqref, imgridfreqs(0), finc, refChan, restFreq);

        //debug TT
        //Double wrld,pixl;
        //pixl=0.0;
        //mySpectral->toWorld(wrld,pixl);
        //cerr<<"world="<<wrld<<" pixel="<<pixl;
      }
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Start frequency = " // Loglevel INFO
	 << mfImageStart_p.get("GHz").getValue()
	 << ", channel increment = "
	 << mfImageStep_p.get("GHz").getValue() 
	 << "GHz, frequency frame = "
         << MFrequency::showType(mfreqref)
         << endl;
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Rest frequency is "  // Loglevel INFO
	 << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	 << "GHz" << LogIO::POST;
      
  }
  

  else {
    Vector<Double> chanFreq;
    Vector<Double> freqResolution;
    //starting with a default rest frequency to be ref 
    //in case none is defined
    Double restFreq=
      msc.spectralWindow().refFrequency()(spectralwindowids_p(0));

    for (Int spwIndex=0; spwIndex < nspw; ++spwIndex){
 
      Int spw=spectralwindowids_p(spwIndex);
      
      Vector<Double> restFreqArray;
      if(getRestFreq(restFreqArray, spw)){
	if(spwIndex==0){
	  restFreq = restFreqArray(0);
	}
	else{
	  if(restFreq != restFreqArray(0)){
	    os << LogIO::WARN << "Rest frequencies are different for  spectralwindows selected " 
	       << LogIO::POST;
	    os << LogIO::WARN 
	       <<"Will be using the restFreq defined in spectralwindow "
	       << spectralwindowids_p(0) << LogIO::POST;
	  }
	  
	}	
      }
    }
  
    // use data frame here (for channel mode)
    rstate=calcImFreqs(chanFreq, freqResolution, obsFreqRef, obsEpoch, obsPosition,restFreq);

    if(imageMode_p=="CHANNEL") {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      if(imageStep_p==0)
	imageStep_p=1;
//	TT: commented these out otherwise the case for multiple MSes would not work
//	Int nsubchans=
//	(chanFreq.shape()(0) - Int(imageStart_p)+1)/Int(imageStep_p)+1;
//      if((nsubchans >0) && (imageNchan_p>nsubchans)) imageNchan_p=nsubchans;

      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Image spectral coordinate: "<< imageNchan_p
         << " channels, starting at visibility channel "
	 << imageStart_p << " stepped by " << imageStep_p << LogIO::POST;
     
      Double finc;
      if(imageNchan_p > 1){
	  finc=chanFreq(1)-chanFreq(0);
      }
      else if(imageNchan_p==1) {
	finc=freqResolution(IPosition(1,0));
      }

      // Now use outframe (instead of data frame) as the rest of
      // the modes do
      //
      mySpectral = new SpectralCoordinate(freqFrame_p,
      //mySpectral = new SpectralCoordinate(obsFreqRef,
      					  chanFreq(0),
                                          finc,  
      					  refChan, restFreq);
        
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Frequency = " // Loglevel INFO
	 << MFrequency(Quantity(chanFreq(0), "Hz")).get("GHz").getValue()
	 << ", channel increment = "
	 << MFrequency(Quantity(finc, "Hz")).get("GHz").getValue() 
	 << "GHz" << endl;
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Rest frequency is "  // Loglevel INFO
	 << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	 << "GHz" << LogIO::POST;
      
    }
    // Spectral channels resampled at equal increments in optical velocity
    // Here we compute just the first two channels and use increments for
    // the others
    else if (imageMode_p=="VELOCITY" || imageMode_p.contains("RADIO")) {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      {
	ostringstream oos;
	oos << "Image spectral coordinate:"<< imageNchan_p 
	    << " channels, starting at radio velocity " << mImageStart_p
	    << "  stepped by " << mImageStep_p << endl;
	os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
           << String(oos); // Loglevel INFO
      }
     
      MFrequency::Types mfreqref=MFrequency::LSRK;
      //Can't convert to frame in mImageStart
      if(!MFrequency::getType(mfreqref, (MRadialVelocity::showType(mImageStart_p.getRef().getType()))))
	mfreqref=freqFrame_p;
      mfreqref=(obsFreqRef==(MFrequency::REST)) ? MFrequency::REST : mfreqref; 
      //rstate=calcImFreqs(chanFreq, freqResolution, mfreqref, obsEpoch, obsPosition,restFreq);
      rstate=calcImFreqs(chanFreq, freqResolution, obsFreqRef, obsEpoch, obsPosition,restFreq);


      Double finc;
      if(imageNchan_p ==1) {
        finc = freqResolution(0);
      }
      else {
        finc = chanFreq(1)-chanFreq(0);
      }
      //mySpectral = new SpectralCoordinate(obsFreqRef,
      mySpectral = new SpectralCoordinate(mfreqref,
      					  chanFreq(0),
                                          finc,  
      					  refChan, restFreq);
     
      {
	ostringstream oos;
	oos << "Reference Frequency = "
	    << MFrequency(Quantity(chanFreq(0), "Hz")).get("GHz")
	    << ", spectral increment = "
	    << MFrequency(Quantity(finc, "Hz")).get("GHz") 
	    << ", frequency frame = "
            << MFrequency::showType(mfreqref)
            << endl; 
	oos << "Rest frequency is " 
	    << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	    << " GHz" << endl;
	os << LogIO::NORMAL << String(oos) << LogIO::POST; // Loglevel INFO
      }
      
    }
    // Since optical/relativistic velocity is non-linear in frequency, we have to
    // pass in all the frequencies. For radio velocity we can use 
    // a linear axis.
    else if (imageMode_p=="OPTICALVELOCITY" || imageMode_p.contains("OPTICAL") || imageMode_p.contains("TRUE") 
	     || imageMode_p.contains("BETA") ||  imageMode_p.contains("RELATI") ) {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      {
	ostringstream oos;
	oos << "Image spectral coordinate: "<< imageNchan_p 
	    << " channels, starting at optical velocity " << mImageStart_p
	    << "  stepped by " << mImageStep_p << endl;
	os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
           << String(oos); // Loglevel INFO
      }
      
      // Use this next line when non-linear is working
      // when selecting in velocity its specfied freqframe or REST 
      MFrequency::Types imfreqref=(obsFreqRef==MFrequency::REST) ? MFrequency::REST : freqFrame_p;
      //rstate=calcImFreqs(chanFreq, freqResolution, imfreqref, obsEpoch, obsPosition,restFreq);
      rstate=calcImFreqs(chanFreq, freqResolution, obsFreqRef, obsEpoch, obsPosition,restFreq);
  
      if (imageNchan_p==1) {
	mySpectral = new SpectralCoordinate(imfreqref,
					    chanFreq(0),
					    freqResolution(0),
					    refChan, restFreq);
      }
      else {
        mySpectral = new SpectralCoordinate(imfreqref, chanFreq, restFreq);
      }

      {
	ostringstream oos;
	oos << "Reference Frequency = "
	    //<< MFrequency(Quantity(freqs(0), "Hz")).get("GHz")
	    << MFrequency(Quantity(chanFreq(0), "Hz")).get("GHz")
	    << " Ghz, " 
            <<" frequency frame= "<<MFrequency::showType(imfreqref)<<endl;
	os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
           << String(oos) << LogIO::POST; // Loglevel INFO
      }
    }
    else {
      this->unlock();
      os << LogIO::SEVERE << "Unknown mode " << imageMode_p
	 << LogIO::EXCEPTION;
      return False;
    }
        
    
  }
 
    //In FTMachine lsrk is used for channel matching with data channel 
    //hence we make sure that
    // we convert to lsrk when dealing with the channels
  freqFrameValid_p=freqFrameValid_p && (obsFreqRef !=MFrequency::REST);
  if(freqFrameValid_p){
      mySpectral->setReferenceConversion(MFrequency::LSRK, obsEpoch, 
					 obsPosition,
					 phaseCenter_p);
  }

  // Polarization
  Vector<String> polType=msc.feed().polarizationType()(0);
  if (polType(0)!="X" && polType(0)!="Y" &&
      polType(0)!="R" && polType(0)!="L") {
    os << LogIO::WARN << "Unknown stokes types in feed table: ["
       << polType(0) << ", " << polType(1) << "]" << endl
       << "Results open to question!" << LogIO::POST;
  }
  
  if (polType(0)=="X" || polType(0)=="Y") {
    polRep_p=SkyModel::LINEAR;
    os << LogIO::DEBUG1 
       << "Preferred polarization representation is linear" << LogIO::POST;
  }
  else {
    polRep_p=SkyModel::CIRCULAR;
    os << LogIO::DEBUG1
       << "Preferred polarization representation is circular" << LogIO::POST;
  }

  // Compare user input with whatever is allowed by the data. 
  // If possible, allow.

  Vector<Int> whichStokes = decideNPolPlanes(True);
  if(whichStokes.nelements()==0 || (whichStokes.nelements()==1 && whichStokes[0]==0) ) 
    {
      if(polRep_p==SkyModel::CIRCULAR) 
	os << LogIO::SEVERE << "Stokes selection of " << stokes_p << " is not valid for Circular feeds." << LogIO::EXCEPTION;
      else 
	os << LogIO::SEVERE << "Stokes selection of " << stokes_p << " is not valid for Linear feeds." << LogIO::EXCEPTION;
      return False;
    }
 
  StokesCoordinate myStokes(whichStokes);
  
  //  os << LogIO::DEBUG1 << "imagecoordinate : " << (coordInfo).stokesCoordinate((coordInfo).findCoordinate(Coordinate::STOKES)).stokes() << LogIO::POST;

  //Set Observatory info
  ObsInfo myobsinfo;
  myobsinfo.setTelescope(telescop);
  myobsinfo.setPointingCenter(mvPhaseCenter);
  myobsinfo.setObsDate(obsEpoch);
  myobsinfo.setObserver(msc.observation().observer()(0));
  this->setObsInfo(myobsinfo);

  //Adding everything to the coordsystem
  coordInfo.addCoordinate(myRaDec);
  coordInfo.addCoordinate(myStokes);
  coordInfo.addCoordinate(*mySpectral);
  coordInfo.setObsInfo(myobsinfo);

  if(mySpectral) delete mySpectral;

  return True;
}


// Make standard choices for coordinates
Bool Imager::imagecoordinates(CoordinateSystem& coordInfo, const Bool verbose) 
{  
  if(!valid()) return False;
  if(!assertDefinedImageParameters()) return False;
  LogIO os(LogOrigin("Imager", "imagecoordinates()", WHERE));
  
  //===Adjust for some setimage defaults if imageNchan=-1 and spwids=-1
  if(imageNchan_p <=0){
    imageNchan_p=1;
  }
  if(spectralwindowids_p.nelements()==1){
    if(spectralwindowids_p[0]<0){
      spectralwindowids_p.resize();
      if(dataspectralwindowids_p.nelements()==0){
         Int nspwinms=ms_p->spectralWindow().nrow();
         dataspectralwindowids_p.resize(nspwinms);
         indgen(dataspectralwindowids_p);
      }
      spectralwindowids_p=dataspectralwindowids_p;
    }
    
  }
  if(fieldid_p < 0){
     if(datafieldids_p.shape() !=0) {
       fieldid_p=datafieldids_p[0];
     }
     else{
       fieldid_p=0; //best default if nothing is specified
     }
  }
  //===end of default
  Vector<Double> deltas(2);
  deltas(0)=-mcellx_p.get("rad").getValue();
  deltas(1)=mcelly_p.get("rad").getValue();
  
  ROMSColumns msc(*ms_p);
  MFrequency::Types obsFreqRef=MFrequency::DEFAULT;
  ScalarColumn<Int> measFreqRef(ms_p->spectralWindow(),
				  MSSpectralWindow::columnName(MSSpectralWindow::MEAS_FREQ_REF));
  //using the first frame of reference; TO DO should do the right thing 
  //for different frames selected. 
  //Int eh = spectralwindowids_p(0);
  if(spectralwindowids_p.size() && measFreqRef(spectralwindowids_p(0)) >=0) {
     obsFreqRef=(MFrequency::Types)measFreqRef(spectralwindowids_p(0));
  }
			    

  MVDirection mvPhaseCenter(phaseCenter_p.getAngle());
  // Normalize correctly
  MVAngle ra=mvPhaseCenter.get()(0);
  ra(0.0);
  MVAngle dec=mvPhaseCenter.get()(1);
  Vector<Double> refCoord(2);
  refCoord(0)=ra.get().getValue();    
  refCoord(1)=dec;    
  
  Vector<Double> refPixel(2); 
  refPixel(0) = Double(nx_p / 2);
  refPixel(1) = Double(ny_p / 2);
  
  //defining observatory...needed for position on earth
  String telescop = msc.observation().telescopeName()(0);

  // defining epoch as begining time from timerange in OBSERVATION subtable
  // Using first observation for now
  //MEpoch obsEpoch = msc.observation().timeRangeMeas()(0)(IPosition(1,0));
  // modified to use main table's TIME column for better match with what
  // VisIter does.
  MEpoch obsEpoch = msc.timeMeas()(0);

  //Now finding the position of the telescope on Earth...needed for proper
  //frequency conversions

  MPosition obsPosition;
  if(!(MeasTable::Observatory(obsPosition, telescop))){
    os << LogIO::WARN << "Did not get the position of " << telescop 
       << " from data repository" << LogIO::POST;
    os << LogIO::WARN 
       << "Please contact CASA to add it to the repository."
       << LogIO::POST;
    os << LogIO::WARN << "Frequency conversion will not work " << LogIO::POST;
    freqFrameValid_p = False;
  }
  else{
    mLocation_p = obsPosition;
    freqFrameValid_p = True;
  }
  // Now find the projection to use: could probably also use
  // max(abs(w))=0.0 as a criterion
  Projection projection(Projection::SIN);
  if(telescop == "ATCASCP" || telescop == "WSRT" || telescop == "DRAO") {
    os << LogIO::NORMAL // Loglevel NORMAL
       << "Using SIN image projection adjusted for "
       << (telescop == "ATCASCP" ? 'S' : 'N') << "CP" 
       << LogIO::POST;
    Vector<Double> projectionParameters(2);
    projectionParameters(0) = 0.0;
    if(sin(dec) != 0.0){
      projectionParameters(1) = cos(dec)/sin(dec);
      projection=Projection(Projection::SIN, projectionParameters);
    }
    else {
      os << LogIO::WARN
         << "Singular projection for " << telescop << ": using plain SIN"
         << LogIO::POST;
      projection=Projection(Projection::SIN);
    }
  }
  else {
    os << LogIO::DEBUGGING << "Using SIN image projection" << LogIO::POST;
  }
  os << LogIO::NORMAL;
  
  Matrix<Double> xform(2,2);
  xform=0.0;xform.diagonal()=1.0;
  DirectionCoordinate
    myRaDec(MDirection::Types(phaseCenter_p.getRefPtr()->getType()),
	    projection,
	    refCoord(0), refCoord(1),
	    deltas(0), deltas(1),
	    xform,
	    refPixel(0), refPixel(1));
  
  // Now set up spectral coordinate
  SpectralCoordinate* mySpectral=0;
  Double refChan=0.0;
  
  // Spectral synthesis
  // For mfs band we set the window to include all spectral windows
  Int nspw=spectralwindowids_p.nelements();
  if (imageMode_p=="MFS") {
    Double fmin=C::dbl_max;
    Double fmax=-(C::dbl_max);
    Double fmean=0.0;
    for (Int i=0;i<nspw;++i) {
      Int spw=spectralwindowids_p(i);
      Vector<Double> chanFreq=msc.spectralWindow().chanFreq()(spw); 
      Vector<Double> freqResolution = msc.spectralWindow().chanWidth()(spw); 
      
      if(dataMode_p=="none"){
      
	if(i==0) {
	  fmin=min(chanFreq-abs(0.5*freqResolution));
	  fmax=max(chanFreq+abs(0.5*freqResolution));
	}
	else {
	  fmin=min(fmin,min(chanFreq-abs(0.5*freqResolution)));
	  fmax=max(fmax,max(chanFreq+abs(0.5*freqResolution)));
	}
      }
      else if(dataMode_p=="channel"){
	// This needs some careful thought about roundoff - it is likely 
	// still adding an extra half-channel at top and bottom but 
	// if the freqResolution is nonlinear, there are subtleties
	Int elnchan=chanFreq.nelements();
	Int firstchan=0;
        Int elstep=1;
	for (uInt jj=0; jj < dataspectralwindowids_p.nelements(); ++jj){
	  if(dataspectralwindowids_p[jj]==spw){
	    firstchan=dataStart_p[jj];
	    elnchan=dataNchan_p[jj];
	    elstep=dataStep_p[jj];
	  }	
	}
	Int lastchan=firstchan+ elnchan*elstep;
        for(Int k=firstchan ; k < lastchan ;  k+=elstep){
	  fmin=min(fmin,chanFreq[k]-abs(freqResolution[k]*(elstep-0.5)));
	  fmax=max(fmax,chanFreq[k]+abs(freqResolution[k]*(elstep-0.5)));
        }
      }
      else{
	this->unlock();
	os << LogIO::SEVERE 
	   << "setdata has to be in 'channel' or 'none' mode for 'mfs' imaging to work"
	   << LogIO::EXCEPTION;
      return False;
      }
 
    }

    fmean=(fmax+fmin)/2.0;
    Vector<Double> restFreqArray;
    Double restFreq=fmean;
    if(getRestFreq(restFreqArray, spectralwindowids_p(0))){
      restFreq=restFreqArray[0];
    }
    imageNchan_p=1;
    Double finc=(fmax-fmin); 
    mySpectral = new SpectralCoordinate(freqFrame_p,  fmean//-finc/2.0
					, finc,
      					refChan, restFreq);
    os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3) // Loglevel INFO
       << "Center frequency = "
       << MFrequency(Quantity(fmean, "Hz")).get("GHz").getValue()
       << " GHz, synthesized continuum bandwidth = "
       << MFrequency(Quantity(finc, "Hz")).get("GHz").getValue()
       << " GHz" << LogIO::POST;

    if(ntaylor_p>1 && reffreq_p==0.0) 
    {
	    reffreq_p = fmean;
	    os << "Setting center frequency as MS-MFS reference frequency" << LogIO::POST;
    }
  }
  
  else if(imageMode_p.contains("FREQ")) {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      Double restFreq=mfImageStart_p.get("Hz").getValue();
      Vector<Double> restFreqVec;
      if(getRestFreq(restFreqVec, spectralwindowids_p(0))){
	restFreq=restFreqVec[0];
      }
      MFrequency::Types mfreqref=(obsFreqRef==(MFrequency::REST)) ? MFrequency::REST : MFrequency::castType(mfImageStart_p.getRef().getType()) ; 
      //  mySpectral = new SpectralCoordinate(mfreqref,
      //					  mfImageStart_p.get("Hz").getValue()+
      //mfImageStep_p.get("Hz").getValue()/2.0,
      //					  mfImageStep_p.get("Hz").getValue(),
      //					  refChan, restFreq);
      mySpectral = new SpectralCoordinate(mfreqref,
      					  mfImageStart_p.get("Hz").getValue(),
      					  mfImageStep_p.get("Hz").getValue(),
      					  refChan, restFreq);
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Start frequency = " // Loglevel INFO
	 << mfImageStart_p.get("GHz").getValue()
	 << ", channel increment = "
	 << mfImageStep_p.get("GHz").getValue() 
	 << "GHz, frequency frame = "
         << MFrequency::showType(mfreqref)
         << endl;
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Rest frequency is "  // Loglevel INFO
	 << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	 << "GHz" << LogIO::POST;
      
  }
  

  else {
    //    if(nspw>1) {
    //      os << LogIO::SEVERE << "Line modes allow only one spectral window"
    //	 << LogIO::POST;
    //      return False;
    //    }
    Vector<Double> chanFreq;
    Vector<Double> freqResolution;
    //starting with a default rest frequency to be ref 
    //in case none is defined
    Double restFreq=
      msc.spectralWindow().refFrequency()(spectralwindowids_p(0));

    for (Int spwIndex=0; spwIndex < nspw; ++spwIndex){
 
      Int spw=spectralwindowids_p(spwIndex);
      Int origsize=chanFreq.shape()(0);
      Int newsize=origsize+msc.spectralWindow().chanFreq()(spw).shape()(0);
      chanFreq.resize(newsize, True);
      chanFreq(Slice(origsize, newsize-origsize))=msc.spectralWindow().chanFreq()(spw);
      freqResolution.resize(newsize, True);
     freqResolution(Slice(origsize, newsize-origsize))=
	msc.spectralWindow().chanWidth()(spw); 
      

     
      Vector<Double> restFreqArray;
      if(getRestFreq(restFreqArray, spw)){
	if(spwIndex==0){
	  restFreq = restFreqArray(0);
	}
	else{
	  if(restFreq != restFreqArray(0)){
	    os << LogIO::WARN << "Rest frequencies are different for  spectralwindows selected " 
	       << LogIO::POST;
	    os << LogIO::WARN 
	       <<"Will be using the restFreq defined in spectralwindow "
	       << spectralwindowids_p(0) << LogIO::POST;
	  }
	  
	}	
      }
    }
  

    if(imageMode_p=="CHANNEL") {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      Vector<Double> freqs;
      if(imageStep_p==0)
	imageStep_p=1;
//	TT: commented these out otherwise the case for multiple MSes would not work
//	Int nsubchans=
//	(chanFreq.shape()(0) - Int(imageStart_p)+1)/Int(imageStep_p)+1;
//      if((nsubchans >0) && (imageNchan_p>nsubchans)) imageNchan_p=nsubchans;

      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Image spectral coordinate: "<< imageNchan_p
         << " channels, starting at visibility channel "
	 << imageStart_p << " stepped by " << imageStep_p << LogIO::POST;
      freqs.resize(imageNchan_p);
      for (Int chan=0;chan<imageNchan_p;chan++) {
	freqs(chan)=chanFreq(Int(imageStart_p)+Int(Float(chan)*Float(imageStep_p)));
      }
      // Use this next line when non-linear working
      //    mySpectral = new SpectralCoordinate(obsFreqRef, freqs,
      //					restFreq);
      // Since we are taking the frequencies as is, the frame must be
      // what is specified in the SPECTRAL_WINDOW table
      //      Double finc=(freqs(imageNchan_p-1)-freqs(0))/(imageNchan_p-1);
      Double finc;
      if(imageNchan_p > 1){
	finc=freqs(1)-freqs(0);
      }
      else if(imageNchan_p==1) {
	finc=freqResolution(IPosition(1,0))*imageStep_p;
      }


	  //in order to outframe to work need to set here original freq frame
      //mySpectral = new SpectralCoordinate(freqFrame_p, freqs(0)-finc/2.0, finc,
      mySpectral = new SpectralCoordinate(obsFreqRef, freqs(0)//-finc/2.0
					  , finc,
					  refChan, restFreq);
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Frequency = " // Loglevel INFO
	 << MFrequency(Quantity(freqs(0), "Hz")).get("GHz").getValue()
	 << ", channel increment = "
	 << MFrequency(Quantity(finc, "Hz")).get("GHz").getValue() 
	 << "GHz" << endl;
      os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
         << "Rest frequency is "  // Loglevel INFO
	 << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	 << "GHz" << LogIO::POST;
      
    }
    // Spectral channels resampled at equal increments in optical velocity
    // Here we compute just the first two channels and use increments for
    // the others
    else if (imageMode_p=="VELOCITY" || imageMode_p.contains("RADIO")) {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      {
	ostringstream oos;
	oos << "Image spectral coordinate:"<< imageNchan_p 
	    << " channels, starting at radio velocity " << mImageStart_p
	    << "  stepped by " << mImageStep_p << endl;
	os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
           << String(oos); // Loglevel INFO
      }
      Vector<Double> freqs(2);
      freqs=0.0;
      if(Double(mImageStep_p.getValue())!=0.0) {
	MRadialVelocity mRadVel=mImageStart_p;
	for (Int chan=0;chan<2;chan++) {
	  MDoppler mdoppler(mRadVel.getValue().get(), MDoppler::RADIO);
	  freqs(chan)=
	    MFrequency::fromDoppler(mdoppler, 
				    restFreq).getValue().getValue();
	  Quantity vel=mRadVel.get("m/s");
	  Quantity inc=mImageStep_p.get("m/s");
	  vel+=inc;
	  mRadVel=MRadialVelocity(vel, MRadialVelocity::LSRK);
	}
      }
      else {
	for (Int chan=0;chan<2;++chan) {
	  freqs(chan)=chanFreq(chan);
	}
      }

      MFrequency::Types mfreqref=MFrequency::LSRK;
      //Can't convert to frame in mImageStart
      if(!MFrequency::getType(mfreqref, (MRadialVelocity::showType(mImageStart_p.getRef().getType()))))
	mfreqref=freqFrame_p;
      mfreqref=(obsFreqRef==(MFrequency::REST)) ? MFrequency::REST : mfreqref; 
      mySpectral = new SpectralCoordinate(mfreqref, freqs(0)//-(freqs(1)-freqs(0))/2.0
					  , freqs(1)-freqs(0), refChan,
					  restFreq);
     
      {
	ostringstream oos;
	oos << "Reference Frequency = "
	    << MFrequency(Quantity(freqs(0), "Hz")).get("GHz")
	    << ", spectral increment = "
	    << MFrequency(Quantity(freqs(1)-freqs(0), "Hz")).get("GHz") 
	    << ", frequency frame = "
            << MFrequency::showType(mfreqref)
            << endl; 
	oos << "Rest frequency is " 
	    << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	    << " GHz" << endl;
	os << LogIO::NORMAL << String(oos) << LogIO::POST; // Loglevel INFO
      }
      
    }
    // Since optical/relativistic velocity is non-linear in frequency, we have to
    // pass in all the frequencies. For radio velocity we can use 
    // a linear axis.
    else if (imageMode_p=="OPTICALVELOCITY" || imageMode_p.contains("OPTICAL") || imageMode_p.contains("TRUE") 
	     || imageMode_p.contains("BETA") ||  imageMode_p.contains("RELATI") ) {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::EXCEPTION;
	return False;
      }
      {
	ostringstream oos;
	oos << "Image spectral coordinate: "<< imageNchan_p 
	    << " channels, starting at optical velocity " << mImageStart_p
	    << "  stepped by " << mImageStep_p << endl;
	os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
           << String(oos); // Loglevel INFO
      }
      Vector<Double> freqs(imageNchan_p);
      freqs=0.0;
      Double chanVelResolution=0.0;
      if(Double(mImageStep_p.getValue())!=0.0) {
	MRadialVelocity mRadVel=mImageStart_p;
	MDoppler mdoppler;
	for (Int chan=0;chan<imageNchan_p;++chan) {
	  if(imageMode_p.contains("OPTICAL")){
	    mdoppler=MDoppler(mRadVel.getValue().get(), MDoppler::OPTICAL);
	  } 
	  else{
	    mdoppler=MDoppler(mRadVel.getValue().get(), MDoppler::BETA);
	  }
	  freqs(chan)=
	    MFrequency::fromDoppler(mdoppler, restFreq).getValue().getValue();
	  mRadVel.set(mRadVel.getValue()+mImageStep_p.getValue());
	  if(imageNchan_p==1)
	    chanVelResolution=MFrequency::fromDoppler(mdoppler, restFreq).getValue().getValue()-freqs(0);
	}
      }
      else {
	for (Int chan=0;chan<imageNchan_p;++chan) {
	    freqs(chan)=chanFreq(chan);
	}
      }
      // Use this next line when non-linear is working
      // when selecting in velocity its specfied freqframe or REST 
      MFrequency::Types imfreqref=(obsFreqRef==MFrequency::REST) ? MFrequency::REST : freqFrame_p;
  
      if(imageNchan_p==1){
	if(chanVelResolution==0.0)
	  chanVelResolution=freqResolution(0);
	mySpectral = new SpectralCoordinate(imfreqref,
					    freqs(0)//-chanVelResolution/2.0,
					    ,
					    chanVelResolution,
					    refChan, restFreq);
      }
      else{
	mySpectral = new SpectralCoordinate(imfreqref, freqs,
					    restFreq);
      }
      // mySpectral = new SpectralCoordinate(MFrequency::DEFAULT, freqs(0),
      //				       freqs(1)-freqs(0), refChan,
      //				        restFreq);
      {
	ostringstream oos;
	oos << "Reference Frequency = "
	    << MFrequency(Quantity(freqs(0), "Hz")).get("GHz")
	    << " Ghz, " 
            <<" frequency frame= "<<MFrequency::showType(imfreqref)<<endl;
	os << (verbose ? LogIO::NORMAL : LogIO::NORMAL3)
           << String(oos) << LogIO::POST; // Loglevel INFO
      }
    }
    else {
      this->unlock();
      os << LogIO::SEVERE << "Unknown mode " << imageMode_p
	 << LogIO::EXCEPTION;
      return False;
    }
        
    
  }
 
 
    //In FTMachine lsrk is used for channel matching with data channel 
    //hence we make sure that
    // we convert to lsrk when dealing with the channels
  freqFrameValid_p=freqFrameValid_p && (obsFreqRef !=MFrequency::REST);
  if(freqFrameValid_p){
      mySpectral->setReferenceConversion(MFrequency::LSRK, obsEpoch, 
					 obsPosition,
					 phaseCenter_p);
  }

  // Polarization
  Vector<String> polType=msc.feed().polarizationType()(0);
  if (polType(0)!="X" && polType(0)!="Y" &&
      polType(0)!="R" && polType(0)!="L") {
    os << LogIO::WARN << "Unknown stokes types in feed table: ["
       << polType(0) << ", " << polType(1) << "]" << endl
       << "Results open to question!" << LogIO::POST;
  }
  
  if (polType(0)=="X" || polType(0)=="Y") {
    polRep_p=SkyModel::LINEAR;
    os << LogIO::DEBUG1 
       << "Preferred polarization representation is linear" << LogIO::POST;
  }
  else {
    polRep_p=SkyModel::CIRCULAR;
    os << LogIO::DEBUG1
       << "Preferred polarization representation is circular" << LogIO::POST;
  }

  // Compare user input with whatever is allowed by the data. 
  // If possible, allow.

  Vector<Int> whichStokes = decideNPolPlanes(True);
  if(whichStokes.nelements()==0 || (whichStokes.nelements()==1 && whichStokes[0]==0) ) 
    {
      if(polRep_p==SkyModel::CIRCULAR) 
	os << LogIO::SEVERE << "Stokes selection of " << stokes_p << " is not valid for Circular feeds." << LogIO::EXCEPTION;
      else 
	os << LogIO::SEVERE << "Stokes selection of " << stokes_p << " is not valid for Linear feeds." << LogIO::EXCEPTION;
      return False;
    }
  /* STOKESDBG */  //cout << "Imager2::imagecoordinate : Image 'whichStokes' : " << whichStokes << endl;

  StokesCoordinate myStokes(whichStokes);
  
  //Set Observatory info
  ObsInfo myobsinfo;
  myobsinfo.setTelescope(telescop);
  myobsinfo.setPointingCenter(mvPhaseCenter);
  myobsinfo.setObsDate(obsEpoch);
  myobsinfo.setObserver(msc.observation().observer()(0));
  this->setObsInfo(myobsinfo);

  //Adding everything to the coordsystem
  coordInfo.addCoordinate(myRaDec);
  coordInfo.addCoordinate(myStokes);
  coordInfo.addCoordinate(*mySpectral);
  coordInfo.setObsInfo(myobsinfo);

  if(mySpectral) delete mySpectral;

  return True;
}

/* Match the user input stokes_p string to all available options.
   If it matches none, return False - exceptions are thrown.
   If it matches an option, set npol_p : number of desired image planes.

  If checkwithMS = True, check with the PolRep and valid combinations of what
                                     is asked for and what the data allows.

  If checkwithMS = False, check only that the inputs are among the valid list
                                     and give npol=1,2,3 or 4.

  There are two levels of checks because when called from defineimage, it does
  not know what PolRep is and cannot check with the data. This used to be in
  two different functions earlier. It is now easier to keep things consistent.

  PLEASE NOTE :  If any stokes options are to be allowed or disallowed,
                           this is the place to do it. StokesImageUtil handles all combinations,
*/
Vector<Int> Imager::decideNPolPlanes(Bool checkwithMS)
{
    Vector<Int> whichStokes(0);

    // Check with options for npol=1
    if(stokes_p=="I" || stokes_p=="Q" || stokes_p=="U" || stokes_p=="V" || 
       stokes_p=="RR" ||stokes_p=="LL" || //stokes_p=="RL" ||stokes_p=="LR" || 
       stokes_p=="XX" || stokes_p=="YY" ) //|| stokes_p=="XY" ||stokes_p=="YX" ) 
      {
        npol_p=1;

        if(checkwithMS)
	{
 
           // Fill in the stokes vector for the output image
           whichStokes.resize(npol_p);
           // The first 8 depend on circular vs linear
           if(polRep_p==SkyModel::CIRCULAR && stokes_p=="RR")  whichStokes(0)=Stokes::RR;
           else if(polRep_p==SkyModel::CIRCULAR && stokes_p=="LL")  whichStokes(0)=Stokes::LL;
           ////else if(polRep_p==SkyModel::CIRCULAR && stokes_p=="RL")  whichStokes(0)=Stokes::RL;
           ////else if(polRep_p==SkyModel::CIRCULAR && stokes_p=="LR")  whichStokes(0)=Stokes::LR;
           else if(polRep_p==SkyModel::LINEAR && stokes_p=="XX") whichStokes(0)=Stokes::XX;
           else if(polRep_p==SkyModel::LINEAR && stokes_p=="YY")  whichStokes(0)=Stokes::YY;
           ////else if(polRep_p==SkyModel::LINEAR && stokes_p=="XY")  whichStokes(0)=Stokes::XY;
           ////else if(polRep_p==SkyModel::LINEAR && stokes_p=="YX")  whichStokes(0)=Stokes::YX;
           // these next 4 don't depend on circular vs linear
           else if(stokes_p=="I") whichStokes(0)=Stokes::I;
           else if(stokes_p=="Q") whichStokes(0)=Stokes::Q;
           else if(stokes_p=="U")  whichStokes(0)=Stokes::U;
           else if(stokes_p=="V")  whichStokes(0)=Stokes::V;
           else {whichStokes.resize(0);}

	}// end of checkwithMS

    }// end of npol=1
    // Check for options with npol=2
    else if(stokes_p=="IV" || stokes_p=="IQ" || 
              stokes_p=="RRLL" || stokes_p=="XXYY" ||
	    ////stokes_p=="RLLR" || stokes_p=="XYYX" ||
	      stokes_p=="QU" || stokes_p=="UV") 
      {
        npol_p=2;

        if(checkwithMS)
	{
            // Check with polrep and fill in the stokes vector for the output image
           whichStokes.resize(npol_p);
           // The first 4 depend on circular vs linear
           if(polRep_p==SkyModel::CIRCULAR && stokes_p=="RRLL")  
                 {whichStokes(0)=Stokes::RR; whichStokes(1)=Stokes::LL;}
           ////else if(polRep_p==SkyModel::CIRCULAR && stokes_p=="RLLR") 
           ////      {whichStokes(0)=Stokes::RL; whichStokes(1)=Stokes::LR;}
           else if(polRep_p==SkyModel::LINEAR && stokes_p=="XXYY") 
                  {whichStokes(0)=Stokes::XX; whichStokes(1)=Stokes::YY;}
           ////else if(polRep_p==SkyModel::LINEAR && stokes_p=="XYYX") 
	   ////      {whichStokes(0)=Stokes::XY; whichStokes(1)=Stokes::YX;}
           // These 4 don't care about circular vs linear
           else if(stokes_p=="IV") 
                  {whichStokes(0)=Stokes::I; whichStokes(1)=Stokes::V;}
           else if(stokes_p=="QU") 
                  {whichStokes(0)=Stokes::Q; whichStokes(1)=Stokes::U;}
           else if(stokes_p=="IQ") 
                  {whichStokes(0)=Stokes::I; whichStokes(1)=Stokes::Q;}
           else if(stokes_p=="UV") 
                  {whichStokes(0)=Stokes::U; whichStokes(1)=Stokes::V;}
           else 
                  {whichStokes.resize(0);}

	}// end of checkwithMS

    }
    // Check for options with npol=3
    else if(stokes_p=="IQU" || stokes_p=="IUV") 
    {
        npol_p=3;
 
        if(checkwithMS)
	{
          whichStokes.resize(npol_p);
           // Check with polRep_p
	  if(stokes_p=="IUV")
	    { whichStokes(0)=Stokes::I; whichStokes(1)=Stokes::U; whichStokes(2)=Stokes::V; }
          else if(stokes_p=="IQU")
	    { whichStokes(0)=Stokes::I; whichStokes(1)=Stokes::Q; whichStokes(2)=Stokes::U; }
          else 
            { whichStokes.resize(0);}

	}// end of checkwithMS

    }//end of npol=3
    // Check with options for npol=4
    else if(stokes_p=="IQUV") // NOT allowing RRLLRLLR and XXYYXYYX !!
    {
        npol_p=4;
 
       if(checkwithMS)
	{
	  whichStokes.resize(npol_p);
          whichStokes(0)=Stokes::I;
          whichStokes(1)=Stokes::Q;
          whichStokes(2)=Stokes::U;
          whichStokes(3)=Stokes::V;
	}// end of checkwithMS

    }// end of npol=4
    else // If none of the options are valid (checked in setimage/defineimage)
    {
      whichStokes.resize(1); whichStokes[0]=0; // undefined.
    }

    // If no match has been found, or an incompatible match is found, this returns Vector<Int>(0);
    return whichStokes;
}// end of decideNPolPlanes




String Imager::state() 
{
  ostringstream os;
  
  try {
    this->lock();
    os << "General: " << endl;
    os << "  MeasurementSet is " << ms_p->tableName() << endl;
    if(beamValid_p) {
      os << "  Beam fit: " << bmaj_p.get("arcsec").getValue() << " by "
	 << bmin_p.get("arcsec").getValue() << " (arcsec) at pa " 
	 << bpa_p.get("deg").getValue() << " (deg) " << endl;
    }
    else {
      os << "  Beam fit is not valid" << endl;
    }
    
    ROMSColumns msc(*ms_p);
    MDirection mDesiredCenter;
    if(setimaged_p) {
      os << "Image definition settings: "
	"(use setimage in Function Group <setup> to change)" << endl;
      os << "  nx=" << nx_p << ", ny=" << ny_p
	 << ", cellx=" << mcellx_p << ", celly=" << mcelly_p
	 << ", Stokes axes : " << stokes_p << endl;
      Int widthRA=20;
      Int widthDec=20;
      if(doShift_p) {
	os << "  doshift is True: Image phase center will be as specified "
	   << endl;
      }
      else {
	os << "  doshift is False: Image phase center will be that of field "
	   << fieldid_p << " :" << endl;
      }
      
      if(shiftx_p.get().getValue()!=0.0||shifty_p.get().getValue()!=0.0) {
	os << "  plus the shift: longitude: " << shiftx_p
	   << " / cos(latitude) : latitude: " << shifty_p << endl;
      }
      
      MVAngle mvRa=phaseCenter_p.getAngle().getValue()(0);
      MVAngle mvDec=phaseCenter_p.getAngle().getValue()(1);
      os << "     ";
      os.setf(ios::left, ios::adjustfield);
      os.width(widthRA);  os << mvRa(0.0).string(MVAngle::TIME,8);
      os.width(widthDec); os << mvDec.string(MVAngle::DIG2,8);
      os << "     " << MDirection::showType(phaseCenter_p.getRefPtr()->getType())
	 << endl;
      
      if(distance_p.get().getValue()!=0.0) {
	os << "  Refocusing to distance " << distance_p << endl;
      }
      
      if(imageMode_p=="MFS") {
	os << "  Image mode is mfs: Image will be frequency synthesised from spectral windows : ";
	for (uInt i=0;i<spectralwindowids_p.nelements();++i) {
	  os << spectralwindowids_p(i) << " ";
	}
	os << endl;
      }
      
      else {
	os << "  Image mode is " << imageMode_p
	   << "  Image number of spectral channels ="
	   << imageNchan_p << endl;
      }
      
    }
    else {
      os << "Image: parameters not yet set (use setimage "
	"in Function Group <setup> )" << endl;
    }
    
    os << "Data selection settings: (use setdata in Function Group <setup> "
      "to change)" << endl;
    if(dataMode_p=="none") {
      if(mssel_p->nrow() == ms_p->nrow()){
	os << "  All data selected" << endl;
      }
      else{
        os << " Number of rows of data selected= " << mssel_p->nrow() << endl;
      }
    }
    else {
      os << "  Data selection mode is " << dataMode_p << ": have selected "
	 << dataNchan_p << " channels";
      if(dataspectralwindowids_p.nelements()>0) {
	os << " spectral windows : ";
	for (uInt i=0;i<dataspectralwindowids_p.nelements();++i) {
	  os << dataspectralwindowids_p(i) << " ";
	}
      }
      if(datafieldids_p.nelements()>0) {
	os << "  Data selected includes fields : ";
	for (uInt i=0;i<datafieldids_p.nelements();++i) {
	  os << datafieldids_p(i) << " ";
	}
      }
      os << endl;
    }
    os << "Options settings: (use setoptions in Function Group <setup> "
      "to change) " << endl;
    os << "  Gridding cache has " << cache_p << " complex pixels, in tiles of "
       << tile_p << " pixels on a side" << endl;
    os << "  Gridding convolution function is ";
    
    if(gridfunction_p=="SF") {
      os << "Spheroidal wave function";
    }
    else if(gridfunction_p=="BOX") {
      os << "Box car convolution";
    }
    else if(gridfunction_p=="PB") {
      os << "Using primary beam for convolution";
    }
    else {
      os << "Unknown type : " << gridfunction_p;
    }
    os << endl;
    
    if(doVP_p) {
      os << "  Primary beam correction is enabled" << endl;
      //       Table vpTable( vpTableStr_p );   could fish out info and summarize
    }
    else {
      os << "  No primary beam correction will be made " << endl;
    }
    os << "  Image plane padding : " << padding_p << endl;
    
    this->unlock();
  } catch (AipsError x) {
    this->unlock();
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << LogIO::EXCEPTION;

  } 
  return String(os);
}




// Apply a primary beam or voltage pattern to an image
Bool Imager::pb(const String& inimage, 
		const String& outimage,
		const String& incomps,
		const String& /*outcomps*/,
		const String& operation, 
		const MDirection& pointingCenter,
		const Quantity& pa,
		const String& pborvp)

{
  if(!valid()) return False;
  
  LogIO os(LogOrigin("imager", "pb()", WHERE));
  
  PagedImage<Float> * inImage_pointer = 0;
  PagedImage<Float> * outImage_pointer = 0;
  ComponentList * inComps_pointer = 0;
  ComponentList * outComps_pointer = 0;
  PBMath * myPBp = 0;
  try {

    if ( ! doVP_p ) {
      this->unlock();
      os << LogIO::SEVERE << 
	"Must invoke setvp() first in order to apply the primary beam" << LogIO::EXCEPTION;
      return False;
    }
    
    if (pborvp == "vp") {
      this->unlock();
      os << LogIO::SEVERE << "VP application is not yet implemented in DOimager" << LogIO::EXCEPTION;
      return False;
    }

    if (operation == "apply") {
      os << LogIO::DEBUG1
         << "function pb will apply " << pborvp << LogIO::POST;
    } else if (operation=="correct") {
      os << LogIO::DEBUG1
         << "function pb will correct for " << pborvp << LogIO::POST;
    } else {
      this->unlock();
      os << LogIO::SEVERE << "Unknown pb operation " << operation 
	 << LogIO::EXCEPTION;
      return False;
    }
    
    // Get initial image and/or SkyComponents

    if (incomps!="") {
      if(!Table::isReadable(incomps)) {
	this->unlock();
	os << LogIO::SEVERE << "ComponentList " << incomps
	   << " not readable" << LogIO::EXCEPTION;
	return False;
      }
      inComps_pointer = new ComponentList(incomps);
      outComps_pointer = new ComponentList( inComps_pointer->copy() );
    }
    if (inimage !="") {
      if(!Table::isReadable(inimage)) {
	this->unlock();
	os << LogIO::SEVERE << "Image " << inimage << " not readable" 
	   << LogIO::EXCEPTION;
	return False;
      }
      inImage_pointer = new PagedImage<Float>( inimage );
      if (outimage != "") {
	outImage_pointer = new PagedImage<Float>( TiledShape(inImage_pointer->shape(), inImage_pointer->niceCursorShape()), 
						  inImage_pointer->coordinates(), outimage);
      }
    }
    // create the PBMath object, needed to apply PB 
    // to make high res Fourier weight image
    Quantity qFreq;
    if (doDefaultVP_p && inImage_pointer==0) {
      this->unlock();
      os << LogIO::SEVERE << 
	"There is no default telescope associated with a componentlist" 
	 << LogIO::POST;
      os << LogIO::SEVERE << 
	"Either specify the PB/VP via a vptable or supply an image as well" 
	 << LogIO::EXCEPTION;
	return False;
    } else if (doDefaultVP_p && inImage_pointer!=0) {
      // look up the telescope in ObsInfo
      ObsInfo oi = inImage_pointer->coordinates().obsInfo();
      String myTelescope = oi.telescope();
      if (myTelescope == "") {
	this->unlock();
	os << LogIO::SEVERE << "No telescope imbedded in image" 
	   << LogIO::EXCEPTION;
	return False;
      }
      {
	Int spectralIndex=inImage_pointer->coordinates().findCoordinate(Coordinate::SPECTRAL);
	AlwaysAssert(spectralIndex>=0, AipsError);
	SpectralCoordinate
	  spectralCoord=inImage_pointer->coordinates().spectralCoordinate(spectralIndex);
	Vector<String> units(1); units = "Hz";
	spectralCoord.setWorldAxisUnits(units);	
	Vector<Double> spectralWorld(1);
	Vector<Double> spectralPixel(1);
	spectralPixel(0) = 0;
	spectralCoord.toWorld(spectralWorld, spectralPixel);  
	Double freq  = spectralWorld(0);
	qFreq = Quantity( freq, "Hz" );
      }
      String band;
      PBMath::CommonPB whichPB;
      String pbName;
      // get freq from coordinates
      PBMath::whichCommonPBtoUse (myTelescope, qFreq, band, whichPB, pbName);
      if (whichPB  == PBMath::UNKNOWN) {
	this->unlock();
	os << LogIO::SEVERE << "Unknown telescope for PB type: " 
	   << myTelescope << LogIO::EXCEPTION;
	return False;
      }
      myPBp = new PBMath(whichPB);
    } else {
      // get the PB from the vpTable
      Table vpTable( vpTableStr_p );
      ScalarColumn<TableRecord> recCol(vpTable, (String)"pbdescription");
      myPBp = new PBMath(recCol(0));
    }
    AlwaysAssert((myPBp != 0), AipsError);


    // Do images (if indeed we have any)
    if (outImage_pointer!=0) {
      Vector<Int> whichStokes;
      CoordinateSystem cCoords;
      cCoords=StokesImageUtil::CStokesCoord(//inImage_pointer->shape(),
					    inImage_pointer->coordinates(),
					    whichStokes,
					    polRep_p);
      TempImage<Complex> cIn(inImage_pointer->shape(),
			     cCoords);
      StokesImageUtil::From(cIn, *inImage_pointer);
      if (operation=="apply") {
	myPBp->applyPB(cIn, cIn, pointingCenter, 
		       pa, squintType_p, False);
      } else {
	myPBp->applyPB(cIn, cIn, pointingCenter, 
		       pa, squintType_p, True);
      }
      StokesImageUtil::To(*outImage_pointer, cIn);
    }
    // Do components (if indeed we have any)
    if (inComps_pointer!=0) {
      if (inImage_pointer==0) {
	this->unlock();
	os << LogIO::SEVERE << 
	  "No input image was given for the componentList to get the frequency from" 
	   << LogIO::EXCEPTION;
	return False;
      }
      Int ncomponents = inComps_pointer->nelements();
      for (Int icomp=0;icomp<ncomponents;++icomp) {
	SkyComponent component=outComps_pointer->component(icomp);
	if (operation=="apply") {
	  myPBp->applyPB(component, component, pointingCenter, 
			 qFreq, pa, squintType_p, False);
	} else {
	  myPBp->applyPB(component, component, pointingCenter, 
			 qFreq, pa, squintType_p, True);
	}
      }
    }
    if (myPBp) delete myPBp;
    if (inImage_pointer) delete inImage_pointer;
    if (outImage_pointer) delete outImage_pointer; 
    if (inComps_pointer) delete inComps_pointer; 
    if (outComps_pointer) delete outComps_pointer; 
    return True;
  } catch (AipsError x) {
    if (myPBp) delete myPBp;
    if (inImage_pointer) delete inImage_pointer;
    if (outImage_pointer) delete outImage_pointer; 
    if (inComps_pointer) delete inComps_pointer; 
    if (outComps_pointer) delete outComps_pointer; 
    this->unlock();
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << LogIO::EXCEPTION;
    return False;
  }
  return True;
}


// Apply a primary beam or voltage pattern to an image
Bool Imager::pbguts(ImageInterface<Float>& inImage, 
		    ImageInterface<Float>& outImage,
		    const MDirection& pointingDirection,
		    const Quantity& pa)
{
  if(!valid()) return False;
  
  LogIO os(LogOrigin("imager", "pbguts()", WHERE));
  
  try {
    if ( ! doVP_p ) {
      os << LogIO::SEVERE << 
	"Must invoke setvp() first in order to apply the primary beam" << LogIO::POST;
      return False;
    }
    String operation = "apply";  // could have as input in the future!

    // create the PBMath object, needed to apply PB 
    // to make high res Fourier weight image
    Quantity qFreq;
    PBMath * myPBp = 0;

    if (doDefaultVP_p) {
      // look up the telescope in ObsInfo
      ObsInfo oi = inImage.coordinates().obsInfo();
      String myTelescope = oi.telescope();
      if (myTelescope == "") {
	os << LogIO::SEVERE << "No telescope embedded in image" << LogIO::POST;
	return False;
      }
      {
	Int spectralIndex=inImage.coordinates().findCoordinate(Coordinate::SPECTRAL);
	AlwaysAssert(spectralIndex>=0, AipsError);
	SpectralCoordinate
	  spectralCoord=inImage.coordinates().spectralCoordinate(spectralIndex);
	Vector<String> units(1); units = "Hz";
	spectralCoord.setWorldAxisUnits(units);	
	Vector<Double> spectralWorld(1);
	Vector<Double> spectralPixel(1);
	spectralPixel(0) = 0;
	spectralCoord.toWorld(spectralWorld, spectralPixel);  
	Double freq  = spectralWorld(0);
	qFreq = Quantity( freq, "Hz" );
      }
      String band;
      PBMath::CommonPB whichPB;
      String pbName;
      // get freq from coordinates
      PBMath::whichCommonPBtoUse (myTelescope, qFreq, band, whichPB, pbName);
      if (whichPB  != PBMath::UNKNOWN) {
      
	myPBp = new PBMath(whichPB);
      }
      else{
	ROMSAntennaColumns ac(ms_p->antenna());
	Double dishDiam=ac.dishDiameter()(0);
	myPBp= new PBMath(dishDiam, True, qFreq);

      }
    } else {
      // get the PB from the vpTable
      Table vpTable( vpTableStr_p );
      ScalarColumn<TableRecord> recCol(vpTable, (String)"pbdescription");
      myPBp = new PBMath(recCol(0));
    }
    AlwaysAssert((myPBp != 0), AipsError);

    Vector<Int> whichStokes;
    CoordinateSystem cCoords;
    cCoords=StokesImageUtil::CStokesCoord(//inImage.shape(),
					  inImage.coordinates(),
					  whichStokes,
					  polRep_p);
    TempImage<Complex> cIn(inImage.shape(),
			   cCoords);
    StokesImageUtil::From(cIn, inImage);
    if (operation=="apply") {
      myPBp->applyPB(cIn, cIn, pointingDirection, 
		     pa, squintType_p, False);
    } else {
      myPBp->applyPB(cIn, cIn, pointingDirection, 
		     pa, squintType_p, True);
    }
    StokesImageUtil::To(outImage, cIn);
    delete myPBp;    
    return True;
  } catch (AipsError x) {
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << LogIO::POST;
    return False;
  } 
  
  return True;
}


Bool Imager::pixon(const String& algorithm,
		   const Quantity& sigma, 
		   const String& model)
{
  if(!valid()) {
    return False;
  }
  LogIO os(LogOrigin("imager", "pixon()", WHERE));
  
  this->lock();

  try {

    if(algorithm=="singledish") {

      if(!assertDefinedImageParameters()) {
	return False;
      }
      String modelName=model;
      if(modelName=="") modelName=imageName()+".pixon";
      make(modelName);
      
      PagedImage<Float> modelImage(modelName);
      
      os << LogIO::NORMAL << "Single dish pixon processing" << LogIO::POST; // Loglevel PROGRESS
      os << LogIO::NORMAL // Loglevel INFO
         << "Using defaults for primary beams in pixon processing" << LogIO::POST;
      gvp_p=new VPSkyJones(*mssel_p, True, parAngleInc_p, squintType_p,
                           skyPosThreshold_p);
      os << LogIO::NORMAL << "Calculating data sampling, etc." << LogIO::POST; // Loglevel PROGRESS
      SDDataSampling ds(*mssel_p, *gvp_p, modelImage.coordinates(),
			modelImage.shape(), sigma);
      
      os << LogIO::NORMAL << "Finding pixon solution" << LogIO::POST; // Loglevel PROGRESS
      PixonProcessor pp;

      IPosition zero(4, 0, 0, 0, 0);
      Array<Float> result;
      if(pp.calculate(ds, result)) {
	os << LogIO::NORMAL << "Pixon solution succeeded" << LogIO::POST; // Loglevel INFO
	modelImage.putSlice(result, zero);
      }
      else {
	os << LogIO::WARN << "Pixon solution failed" << LogIO::POST;
      }
    }
    else if(algorithm=="synthesis") {

      if(!assertDefinedImageParameters()) {
	return False;
      }
      String modelName=model;
      if(modelName=="") modelName=imageName()+".pixon";
      make(modelName);
      
      PagedImage<Float> modelImage(modelName);
      
      os << LogIO::NORMAL << "Synthesis pixon processing" << LogIO::POST; // Loglevel INFO
      os << LogIO::NORMAL << "Calculating data sampling, etc." << LogIO::POST; // Loglevel PROGRESS
      SynDataSampling ds(*mssel_p, modelImage.coordinates(),
			 modelImage.shape(), sigma);
      
      os << LogIO::NORMAL << "Finding pixon solution" << LogIO::POST; // Loglevel PROGRESS
      PixonProcessor pp;
      
      IPosition zero(4, 0, 0, 0, 0);
      Array<Float> result;
      if(pp.calculate(ds, result)) {
	os << LogIO::NORMAL << "Pixon solution succeeded" << LogIO::POST; // Loglevel INFO
	modelImage.putSlice(result, zero);
      }
      else {
	os << LogIO::WARN << "Pixon solution failed" << LogIO::POST;
      }
    }
    else if(algorithm=="synthesis-image") {

      if(!assertDefinedImageParameters()) {
	return False;
      }
      String modelName=model;
      if(modelName=="") modelName=imageName()+".pixon";
      make(modelName);
      
      PagedImage<Float> modelImage(modelName);
      
      os << LogIO::NORMAL << "Synthesis image pixon processing" << LogIO::POST; // Loglevel PROGRESS
      String dirtyName=modelName+".dirty";
      Imager::makeimage("corrected", dirtyName);
      String psfName=modelName+".psf";
      Imager::makeimage("psf", psfName);
      PagedImage<Float> dirty(dirtyName);
      PagedImage<Float> psf(psfName);
      os << "Calculating data sampling, etc." << LogIO::POST;
      ImageDataSampling imds(dirty, psf, sigma.getValue());
      os << "Finding pixon solution" << LogIO::POST;
      PixonProcessor pp;
      IPosition zero(4, 0, 0, 0, 0);
      Array<Float> result;
      if(pp.calculate(imds, result)) {
	os << "Pixon solution succeeded" << LogIO::POST;
	modelImage.putSlice(result, zero);
      }
      else {
	os << LogIO::WARN << "Pixon solution failed" << LogIO::POST;
      }
    }

    else if(algorithm=="test") {

      os << LogIO::NORMAL << "Pixon standard test" << LogIO::POST; // Loglevel INFO
      PixonProcessor pp;

      return pp.standardTest();
      
    } else {
      this->unlock();
      os << LogIO::SEVERE << "Unknown algorithm: " << algorithm << LogIO::POST;
      return False;

    }

    this->unlock();

    return True;
  } catch (AipsError x) {
    this->unlock();
    os << LogIO::SEVERE << "Exception: " << x.getMesg() << LogIO::POST;

    return False;
  } 
  this->unlock();
  return True;

}

void Imager::printbeam(CleanImageSkyModel *sm_p, LogIO &os, const Bool firstrun)
{
  //Use predefined beam for restoring or find one by fitting
  Bool printBeam = false;
  if(!beamValid_p){
    Vector<Float> beam(3);
    beam=sm_p->beam(0);
    if(beam[0] > 0.0){
      bmaj_p=Quantity(abs(beam(0)), "arcsec"); 
      bmin_p=Quantity(abs(beam(1)), "arcsec");
      bpa_p=Quantity(beam(2), "deg");
      beamValid_p=True;
      printBeam = true;
      os << LogIO::NORMAL << "Fitted beam used in restoration: " ;	 // Loglevel INFO
    }
  }
  else if(firstrun){
    printBeam = true;
    os << LogIO::NORMAL << "Beam used in restoration: "; // Loglevel INFO
  }
  if(printBeam)
    os << LogIO::NORMAL << bmaj_p.get("arcsec").getValue() << " by " // Loglevel INFO
       << bmin_p.get("arcsec").getValue() << " (arcsec) at pa " 
       << bpa_p.get("deg").getValue() << " (deg) " << LogIO::POST;
}

Bool Imager::restoreImages(const Vector<String>& restoredNames)
{

  LogIO os(LogOrigin("imager", "restoreImages()", WHERE));
  try{
    // It's important that we use the congruent images in both
    // cases. This means that we must use the residual image as
    // passed to the SkyModel and not the one returned. This 
    // distinction is important only (currently) for WFCleanImageSkyModel
    // which has a different representation for the image internally.
    Vector<String> residualNames(images_p.nelements());
    Vector<String> modelNames(images_p.nelements());
    for(uInt k=0; k < modelNames.nelements() ; ++k){
      residualNames[k]=residuals_p[k]->name();
      modelNames[k]=images_p[k]->name();
    }

    /*
    if((nx_p*ny_p > 7000*7000) && (ft_p->name() != "MosaicFT")){
      // very large for convolution ...destroy Skyequations to release memory
      // need to fix the convolution to be leaner
      destroySkyEquation();

    }
    */
    Bool dorestore=False;
    if( (bmaj_p.getValue() >0.0) && (bmin_p.getValue() > 0.0))
      dorestore=True;
    
    if(restoredNames.nelements()>0) {
      for (Int thismodel=0;thismodel<Int(restoredNames.nelements()); 
	   ++thismodel) {
	if(restoredNames(thismodel)!="") {
	  PagedImage<Float> modelIm(modelNames[thismodel]);
	  PagedImage<Float> residIm(residualNames[thismodel]);
	  TempImage<Float> restored;
	  if(dorestore){
	    restored=TempImage<Float>(modelIm.shape(),
				modelIm.coordinates());
	    restored.copyData(modelIm);
	   
	    StokesImageUtil::Convolve(restored, bmaj_p, bmin_p, bpa_p);
	  }
	  // We can work only if the residual image was defined.
	  if(residIm.name() != "") {
	    if(dorestore){
	      LatticeExpr<Float> le(restored+(residIm)); 
	      restored.copyData(le);
	    }
	    if(freqFrameValid_p){
	      CoordinateSystem cs=residIm.coordinates();
	      String errorMsg;
	      if (CoordinateUtil::setSpectralConversion (errorMsg, cs,MFrequency::showType(freqFrame_p))) {
		residIm.setCoordinateInfo(cs);
		if(dorestore)
		  restored.setCoordinateInfo(cs);
	      }
	    } 
	    
	    //should be able to do that only on testing dofluxscale
	    // ftmachines or sm_p should tell us that
	    
	    // this should go away...e
	    //special casing like this gets hard to maintain
	    // need to redo how interactive is done so that it is outside 
	    //of this code 


	    //
	    // Using minPB_p^2 below to make it consistent with the normalization in SkyEquation.
	    //
	    Float cutoffval=minPB_p;
	    if(ft_p->name()=="MosaicFT")	      
	      cutoffval=minPB_p*minPB_p;
	    
	    if (sm_p->doFluxScale(thismodel)) {
	      TempImage<Float> cover(modelIm.shape(),modelIm.coordinates());
	      if(ft_p->name()=="MosaicFT")
		se_p->getCoverageImage(thismodel, cover);
              else
                  cover.copyData(sm_p->fluxScale(thismodel));
	      if(scaleType_p=="NONE"){
		if(dorestore){
		  LatticeExpr<Float> le(iif(cover < minPB_p, 
					    0.0,(restored/(sm_p->fluxScale(thismodel)))));
		  restored.copyData(le);
		}
		LatticeExpr<Float> le1(iif(cover < minPB_p, 
					   0,(residIm/(sm_p->fluxScale(thismodel)))));
		residIm.copyData(le1);
		
	      }
	      
	      //Setting the bit-mask for mosaic image
	      LatticeExpr<Bool> lemask(iif((cover < cutoffval) , 
					   False, True));
	      if(dorestore){
		ImageRegion outreg=restored.makeMask("mask0", False, True);
		LCRegion& outmask=outreg.asMask();
		outmask.copyData(lemask);
		restored.defineRegion("mask0", outreg, RegionHandler::Masks, True);
		restored.setDefaultMask("mask0");
	  
	      }
	      
	    }
	    if(dorestore){
	      PagedImage<Float> diskrestore(restored.shape(), restored.coordinates(), 
					    restoredNames(thismodel));
	      diskrestore.copyData(restored);
	      ImageInfo ii = modelIm.imageInfo();
	      ii.setRestoringBeam(bmaj_p, bmin_p, bpa_p); 
	      diskrestore.setImageInfo(ii);
	      diskrestore.setUnits(Unit("Jy/beam"));
	      copyMask(diskrestore, restored, "mask0");
		 
	    }
	    
	  }
	  else {
	    os << LogIO::SEVERE << "No residual image for model "
	       << thismodel << ", cannot restore image" << LogIO::POST;
	  }
	  
	  if(residuals_p[thismodel]->ok()){
	    residuals_p[thismodel]->table().relinquishAutoLocks(True);
	    residuals_p[thismodel]->table().unlock();
	    //detaching residual so that other processes can use it
	    residuals_p[thismodel]=0;
	  }
	}
	
      }

    }
  }
catch (exception &x) { 
    os << LogIO::SEVERE << "Exception: " << x.what() << LogIO::POST;
    return False;
  }
  return True;
}

Bool Imager::copyMask(ImageInterface<Float>& out, const ImageInterface<Float>& in, String maskname, Bool setdef){
  try{
    if(in.hasRegion(maskname) && !out.hasRegion(maskname)){
      ImageRegion outreg=out.makeMask(maskname, False, True);
      LCRegion& outmask=outreg.asMask();
      outmask.copyData(in.getRegion("mask0").asLCRegion());
      LatticeExpr<Float> myexpr(iif(outmask, out, 0.0) );
      out.copyData(myexpr);
      out.defineRegion(maskname, outreg, RegionHandler::Masks, True);
      if(setdef)
	out.setDefaultMask(maskname);
    }
    else
      return False;
    
  }
  catch(exception &x){
    throw(AipsError(x.what()));
    return False;
  }


  return True;
}

Bool Imager::writeFluxScales(const Vector<String>& fluxScaleNames)
{
  LogIO os(LogOrigin("imager", "writeFluxScales()", WHERE));
  Bool answer = False;
  ImageInterface<Float> *cover;
  if(fluxScaleNames.nelements()>0) {
    for (Int thismodel=0;thismodel<Int(fluxScaleNames.nelements());++thismodel) {
      if(fluxScaleNames(thismodel)!="") {
        PagedImage<Float> fluxScale(images_p[thismodel]->shape(),
                                    images_p[thismodel]->coordinates(),
                                    fluxScaleNames(thismodel));
	PagedImage<Float> coverimage(images_p[thismodel]->shape(),
				     images_p[thismodel]->coordinates(),
				     fluxScaleNames(thismodel)+".pbcoverage");
        coverimage.table().markForDelete();
	if(freqFrameValid_p){
	  CoordinateSystem cs=fluxScale.coordinates();
	  String errorMsg;
	  if (CoordinateUtil::setSpectralConversion (errorMsg, cs,MFrequency::showType(freqFrame_p))) {
	    fluxScale.setCoordinateInfo(cs);
            coverimage.setCoordinateInfo(cs);
	  }
        }
        if (sm_p->doFluxScale(thismodel)) {
	  cover=&(sm_p->fluxScale(thismodel));
	  answer = True;
          fluxScale.copyData(sm_p->fluxScale(thismodel));
	  Float cutoffval=minPB_p;
	  if(ft_p->name()=="MosaicFT"){
	    cutoffval=minPB_p*minPB_p;
	    se_p->getCoverageImage(thismodel, coverimage);
            cover=&(coverimage);
	    //Do the sqrt
	    coverimage.copyData(( LatticeExpr<Float> )(iif(coverimage > 0.0, sqrt(coverimage), 0.0)));
            coverimage.table().unmarkForDelete();
	    LatticeExpr<Bool> lemask(iif((*cover) < sqrt(cutoffval), 
				       False, True));
	    ImageRegion outreg=coverimage.makeMask("mask0", False, True);
	    LCRegion& outmask=outreg.asMask();
	    outmask.copyData(lemask);
	    coverimage.defineRegion("mask0", outreg,RegionHandler::Masks, True); 
	    coverimage.setDefaultMask("mask0");
	  }
	  LatticeExpr<Bool> lemask(iif((*cover) < minPB_p, 
				       False, True));
	  ImageRegion outreg=fluxScale.makeMask("mask0", False, True);
	  LCRegion& outmask=outreg.asMask();
	  outmask.copyData(lemask);
	  fluxScale.defineRegion("mask0", outreg,RegionHandler::Masks, True); 
	  fluxScale.setDefaultMask("mask0");


        } else {
	  answer = False;
          os << LogIO::NORMAL // Loglevel INFO
             << "No flux scale available (or required) for model " << thismodel
             << LogIO::POST;
          os << LogIO::NORMAL // Loglevel INFO
             << "(This is only pertinent to mosaiced images)" << LogIO::POST;
          os << LogIO::NORMAL // Loglevel INFO
             << "Writing out image of constant 1.0" << LogIO::POST;
          fluxScale.set(1.0);
        }
      }
    }
  }
  return answer;
}
// Supports the "[] or -1 => everything" convention using the rule:
// If v is empty or only has 1 element, and it is < 0, 
//     replace v with 0, 1, ..., nelem - 1.
// Returns whether or not it modified v.
//   If so, v is modified in place.
Bool Imager::expand_blank_sel(Vector<Int>& v, const uInt nelem)
{
  if((v.nelements() == 1 && v[0] < 0) || v.nelements() == 0){
    v.resize(nelem);
    indgen(v);
    return true;
  }
  return false;
}


// Correct the data using a plain VisEquation.
// This just moves data from observed to corrected.
// Eventually we should pass in a calibrater
// object to do the work.
Bool Imager::correct(const Bool /*doparallactic*/, const Quantity& /*t*/)
{
  if(!valid()) return False;

  LogIO os(LogOrigin("imager", "correct()", WHERE));
  
  // (gmoellen 06Nov20)
  // VisEquation/VisJones have been extensively revised, so
  //  so we are disabling this method.  Will probably soon
  /// delete it entirely, since this is a calibrater responsibility....

  this->lock();
  try {

    // Warn users we are disabled.
    throw(AipsError("Imager::correct() has been disabled (gmoellen 06Nov20)."));

 /* Commenting out old VisEquation/VisJones usage
    os << "Correcting data: CORRECTED_DATA column will be replaced"
       << LogIO::POST;
    
    
    if(doparallactic) {
      os<<"Correcting parallactic angle variation"<<LogIO::POST;
      VisEquation ve(*vs_p);
      Float ts=t.get("s").getValue();
      Float dt=ts/10.0;
      PJones pj(*vs_p, ts, dt);
      ve.setVisJones(pj);
      ve.correct();
    }
    else {
      VisEquation ve(*vs_p);
      ve.correct();
    }
    this->unlock();
    return True;
 */


  } catch (AipsError x) {
    this->unlock();
    os << LogIO::SEVERE << "Exception: " << x.getMesg() << LogIO::POST;
    return False;
  } 
  this->unlock();
  
  return True;
}

uInt Imager::count_visibilities(ROVisibilityIterator *rvi,
                                const Bool unflagged_only,
                                const Bool must_have_imwt)
{
  if(!valid()) return 0;
  LogIO os(LogOrigin("imager", "count_visibilities()", WHERE));
  
  this->lock();
  ROVisIter& vi(*rvi);
  VisBuffer vb(vi);
    
  uInt nVis = 0;
  for(vi.originChunks(); vi.moreChunks(); vi.nextChunk()){
    for(vi.origin(); vi.more(); ++vi){
      uInt nRow = vb.nRow();
      uInt nChan = vb.nChannel();
      for(uInt row = 0; row < nRow; ++row)
        for(uInt chn = 0; chn < nChan; ++chn)
          if((!unflagged_only || !vb.flag()(chn,row)) &&
             (!must_have_imwt || vb.imagingWeight()(chn,row) > 0.0))
            ++nVis;
    }
  }
  this->unlock();
  return nVis;
}
// Create the FTMachine as late as possible
Bool Imager::createFTMachine()
{
  

  if(ft_p) {delete ft_p; ft_p=0;}
  if(gvp_p) {delete gvp_p; gvp_p=0;}
  if(cft_p) {delete cft_p; cft_p=0;}

  //For ftmachines that can use double precision gridders 
  Bool useDoublePrecGrid=False;
  //few channels use Double precision
  //till we find a better algorithm to determine when to use Double prec gridding
  if((imageNchan_p < 5) && !(singlePrec_p))
    useDoublePrecGrid=True;

  LogIO os(LogOrigin("imager", "createFTMachine()", WHERE));
  
  // This next line is only a guess
  Int numberAnt=((MeasurementSet&)*ms_p).antenna().nrow();
  
  Int imageVolume=imageshape().product()/square(facets_p);
  if(imageVolume) {
    //    Int minimumCachesize=
    //      tile_p*tile_p*npol_p*imageNchan_p*numberAnt*(numberAnt-1)/2;
    if(imageVolume>cache_p/2) {
      os<<"Cache too small to fit entire image, will use tiled gridding"
	<< endl;  
      os<<"MS has at most "<<numberAnt<<" antennas"<<endl;
      os<<"Optimum cache size to do all baselines is "<< imageVolume
	<< " (complex words)" << LogIO::POST;
    }
  }
  

  Float padding;
  padding=1.0;
  if(doMultiFields_p||(facets_p>1)) {
    padding = padding_p;
    os << LogIO::NORMAL // Loglevel INFO
       << "Multiple fields or facets: transforms will be padded by a factor "
       << padding << LogIO::POST;
  }

  if(ftmachine_p=="sd") {
    os << LogIO::NORMAL // Loglevel INFO
       << "Performing single dish gridding..." << LogIO::POST;
    os << LogIO::NORMAL1 // gridfunction_p is too cryptic for most users.
       << "with convolution function " << gridfunction_p << LogIO::POST;

    // Now make the Single Dish Gridding
    os << LogIO::NORMAL // Loglevel INFO
       << "Gridding will use specified common tangent point:" << LogIO::POST;
    os << LogIO::NORMAL << tangentPoint() << LogIO::POST; // Loglevel INFO
    if(gridfunction_p=="pb") {
      if(!gvp_p) {
	if (doDefaultVP_p) {
	  os << LogIO::NORMAL // Loglevel INFO
             << "Using defaults for primary beams used in gridding" << LogIO::POST;
	  gvp_p=new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p,
                               skyPosThreshold_p);
	} else {
	  os << LogIO::NORMAL // Loglevel INFO
             << "Using VP as defined in " << vpTableStr_p <<  LogIO::POST;
	  Table vpTable( vpTableStr_p ); 
	  gvp_p=new VPSkyJones(*ms_p, vpTable, parAngleInc_p, squintType_p,
                               skyPosThreshold_p);
	}
      } 
      ft_p = new SDGrid(mLocation_p, *gvp_p, cache_p/2, tile_p, gridfunction_p,
                        sdConvSupport_p);
    }
    else {
      ft_p = new SDGrid(mLocation_p, cache_p/2, tile_p, gridfunction_p,
                        sdConvSupport_p);
    }
    ft_p->setPointingDirColumn(pointingDirCol_p);

    ROVisIter& vi(*rvi_p);
    // Get bigger chunks o'data: this should be tuned some time
    // since it may be wrong for e.g. spectral line
    vi.setRowBlocking(100);
    
    AlwaysAssert(ft_p, AipsError);
    
    cft_p = new SimpleComponentGridMachine();
    AlwaysAssert(cft_p, AipsError);
  }
  else if(ftmachine_p=="mosaic") {
    os << LogIO::NORMAL << "Performing mosaic gridding" << LogIO::POST; // Loglevel PROGRESS
   
    setMosaicFTMachine(useDoublePrecGrid);

    // VisIter& vi(vs_p->iter());
    //   vi.setRowBlocking(100);
    
    AlwaysAssert(ft_p, AipsError);
    
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);
  }
  //
  // Make WProject FT machine (for non co-planar imaging)
  //
  else if (ftmachine_p == "wproject"){
    os << LogIO::NORMAL << "Performing w-plane projection" // Loglevel PROGRESS
       << LogIO::POST;
    if(wprojPlanes_p<64) {
      os << LogIO::WARN
	 << "No of WProjection planes set too low for W projection - recommend at least 128"
	 << LogIO::POST;
    }
    
    ft_p = new WProjectFT(wprojPlanes_p,  mLocation_p,
			  cache_p/2, tile_p, True, padding_p, useDoublePrecGrid);
    AlwaysAssert(ft_p, AipsError);
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);
  }
  //
  // Make PBWProject FT machine (for non co-planar imaging with
  // antenna based PB corrections)
  //
  else if (ftmachine_p == "pbwproject"){
    if (wprojPlanes_p<=1)
      {
	os << LogIO::NORMAL
	   << "You are using wprojplanes=1. Doing co-planar imaging (no w-projection needed)" 
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing pb-projection" << LogIO::POST; // Loglevel PROGRESS
      }
    if((wprojPlanes_p>1)&&(wprojPlanes_p<64)) 
      {
	os << LogIO::WARN
	   << "No. of w-planes set too low for W projection - recommend at least 128"
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing pb + w-plane projection" // Loglevel PROGRESS
	   << LogIO::POST;
      }


    if(!gvp_p) 
      {
	os << LogIO::NORMAL // Loglevel INFO
           << "Using defaults for primary beams used in gridding" << LogIO::POST;
	gvp_p = new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p,
                               skyPosThreshold_p);
      }

    ft_p = new nPBWProjectFT(//*ms_p, 
			     wprojPlanes_p, cache_p/2,
                            cfCacheDirName_p, doPointing, doPBCorr,
                             tile_p, paStep_p, pbLimit_p, True);
    //
    // Explicit type casting of ft_p does not look good.  It does not
    // pick up the setPAIncrement() method of PBWProjectFT without
    // this
    //
    ((nPBWProjectFT *)ft_p)->setPAIncrement(parAngleInc_p);
    if (doPointing)
      {
        try
          {
	    //            epJ = new EPJones(*vs_p, *ms_p);
	    VisSet elVS(*rvi_p);
            epJ = new EPJones(elVS);
            RecordDesc applyRecDesc;
            applyRecDesc.addField("table", TpString);
            applyRecDesc.addField("interp",TpString);
            Record applyRec(applyRecDesc);
            applyRec.define("table",epJTableName_p);
            applyRec.define("interp", "nearest");
            epJ->setApply(applyRec);
            ((nPBWProjectFT *)ft_p)->setEPJones(epJ);
          }
        catch(AipsError& x)
          {
            //
            // Add some more useful info. to the message and translate
            // the generic AipsError exception object to a more specific
            // SynthesisError object.
            //
            String mesg = x.getMesg();
            mesg += ". Error in loading pointing offset table.";
            SynthesisError err(mesg);
            throw(err);
          }
      }

    AlwaysAssert(ft_p, AipsError);
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);
  }
  else if (ftmachine_p == "pbmosaic"){

    if (wprojPlanes_p<=1)
      {
	os << LogIO::NORMAL
	   << "You are using wprojplanes=1. Doing co-planar imaging (no w-projection needed)" 
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing pb-mosaic" << LogIO::POST; // Loglevel PROGRESS
      }
    if((wprojPlanes_p>1)&&(wprojPlanes_p<64)) 
      {
	os << LogIO::WARN
	   << "No. of w-planes set too low for W projection - recommend at least 128"
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing pb + w-plane projection" << LogIO::POST; // Loglevel PROGRESS
      }

    if(!gvp_p) 
      {
	os << LogIO::NORMAL // Loglevel INFO
           << "Using defaults for primary beams used in gridding" << LogIO::POST;
	gvp_p = new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p,
                               skyPosThreshold_p);
      }
    ft_p = new PBMosaicFT(*ms_p, wprojPlanes_p, cache_p/2, 
			  cfCacheDirName_p, /*True */doPointing, doPBCorr, 
			  tile_p, paStep_p, pbLimit_p, True);
    ((PBMosaicFT *)ft_p)->setObservatoryLocation(mLocation_p);
    //
    // Explicit type casting of ft_p does not look good.  It does not
    // pick up the setPAIncrement() method of PBWProjectFT without
    // this
    //
    os << LogIO::NORMAL << "Setting PA increment to " << parAngleInc_p.getValue("deg") << " deg" << endl;
    ((nPBWProjectFT *)ft_p)->setPAIncrement(parAngleInc_p);

    if (doPointing) 
      {
	try
	  {
	    // Warn users we are have temporarily disabled pointing cal
	    //	    throw(AipsError("Pointing calibration temporarily disabled (gmoellen 06Nov20)."));
	    //  TBD: Bring this up-to-date with new VisCal mechanisms
	    VisSet elVS(*rvi_p);
	    epJ = new EPJones(elVS, *ms_p);
	    RecordDesc applyRecDesc;
	    applyRecDesc.addField("table", TpString);
	    applyRecDesc.addField("interp",TpString);
	    Record applyRec(applyRecDesc);
	    applyRec.define("table",epJTableName_p);
	    applyRec.define("interp", "nearest");
	    epJ->setApply(applyRec);
	    ((nPBWProjectFT *)ft_p)->setEPJones(epJ);
	  }
	catch(AipsError& x)
	  {
	    //
	    // Add some more useful info. to the message and translate
	    // the generic AipsError exception object to a more specific
	    // SynthesisError object.
	    //
	    String mesg = x.getMesg();
	    mesg += ". Error in loading pointing offset table.";
	    SynthesisError err(mesg);
	    throw(err);
	  }
      }
    AlwaysAssert(ft_p, AipsError);
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);

  }
  else if(ftmachine_p=="both") {
      
    os << LogIO::NORMAL // Loglevel INFO
       << "Performing single dish gridding with convolution function "
       << gridfunction_p << LogIO::POST;
    os << LogIO::NORMAL // Loglevel INFO
       << "and interferometric gridding with the prolate spheroidal convolution function"
       << LogIO::POST;
    
    // Now make the Single Dish Gridding
    os << LogIO::NORMAL // Loglevel INFO
       << "Gridding will use specified common tangent point:" << LogIO::POST;
    os << LogIO::NORMAL << tangentPoint() << LogIO::POST; // Loglevel INFO
    if(!gvp_p) {
      os << LogIO::NORMAL // Loglevel INFO
         << "Using defaults for primary beams used in gridding" << LogIO::POST;
      gvp_p = new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p,
                             skyPosThreshold_p);
    }
    if(sdScale_p != 1.0)
      os << LogIO::NORMAL // Loglevel INFO
         << "Multiplying single dish data by factor " << sdScale_p << LogIO::POST;
    if(sdWeight_p != 1.0)
      os << LogIO::NORMAL // Loglevel INFO
         << "Multiplying single dish weights by factor " << sdWeight_p
         << LogIO::POST;
    ft_p = new GridBoth(*gvp_p, cache_p/2, tile_p,
			mLocation_p, phaseCenter_p,
			gridfunction_p, "SF", padding,
			sdScale_p, sdWeight_p);
    
    //VisIter& vi(vs_p->iter());
    // Get bigger chunks o'data: this should be tuned some time
    // since it may be wrong for e.g. spectral line
    rvi_p->setRowBlocking(100);
    
    AlwaysAssert(ft_p, AipsError);
    
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);
    
  }  
  else if(ftmachine_p=="nift") {
    os << LogIO::NORMAL // Loglevel INFO
       << "Using FTMachine " << ftmachine_p << LogIO::POST
       << "Performing interferometric gridding..."
       << LogIO::POST;
    os << LogIO::NORMAL1 // gridfunction_p is too cryptic for most users.
       << "...with convolution function " << gridfunction_p << LogIO::POST;

    // Make the re-gridder components.  Here, make the basic
    // re-sampler.
    CountedPtr<VisibilityResamplerBase> visResamplerCtor = new VisibilityResampler();
    // Make the multi-threaded re-sampler and supply the basic
    // re-sampler used in the worklet threads.
    CountedPtr<VisibilityResamplerBase> mthVisResampler = new MultiThreadedVisibilityResampler(useDoublePrecGrid,
											       visResamplerCtor);
    // Now make the FTMachine
    if(facets_p>1) {
      os << LogIO::NORMAL // Loglevel INFO
         << "Multi-facet Fourier transforms will use specified common tangent point:"
	 << LogIO::POST;
      os << LogIO::NORMAL << tangentPoint() << LogIO::POST; // Loglevel INFO
      ft_p = new rGridFT(cache_p / 2, tile_p, mthVisResampler, gridfunction_p, mLocation_p,
			 phaseCenter_p, padding, False, useDoublePrecGrid);
      
    }
    else {
      os << LogIO::DEBUG1
         << "Single facet Fourier transforms will use image center as tangent points"
	 << LogIO::POST;
      ft_p = new rGridFT(cache_p/2, tile_p, mthVisResampler, gridfunction_p, mLocation_p,
			padding, False, useDoublePrecGrid);

    }
    AlwaysAssert(ft_p, AipsError);
    
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);
    
  }
  //===============================================================
  // A-Projection FTMachine code start here
  //===============================================================
  else if (ftmachine_p == "wbawp"){

    if (wprojPlanes_p<=1)
      {
	os << LogIO::NORMAL
	   << "You are using wprojplanes=1. Doing co-planar imaging (no w-projection needed)" 
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing WBA-Projection" << LogIO::POST; // Loglevel PROGRESS
      }
    if((wprojPlanes_p>1)&&(wprojPlanes_p<64)) 
      {
	os << LogIO::WARN
	   << "No. of w-planes set too low for W projection - recommend at least 128"
	   << LogIO::POST;
	os << LogIO::NORMAL << "Performing WBAW-Projection" << LogIO::POST; // Loglevel PROGRESS
      }

    // if(!gvp_p) 
    //   {
    // 	os << LogIO::NORMAL // Loglevel INFO
    //        << "Using defaults for primary beams used in gridding" << LogIO::POST;
    // 	gvp_p = new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p,
    //                            skyPosThreshold_p);
    //   }
    useDoublePrecGrid=False;
    CountedPtr<ATerm> apertureFunction = createTelescopeATerm(*ms_p);
    CountedPtr<ConvolutionFunction> awConvFunc = new AWConvFunc(apertureFunction);
    CountedPtr<VisibilityResamplerBase> visResampler = new AWVisResampler();
    CountedPtr<VisibilityResamplerBase> mthVisResampler = new MultiThreadedVisibilityResampler(useDoublePrecGrid,
     											       visResampler);
    CountedPtr<CFCache> cfcache = new CFCache();
    cfcache->setCacheDir(cfCacheDirName_p.data());
    cfcache->initCache();
    ft_p = new AWProjectWBFT(wprojPlanes_p, cache_p/2, 
			     cfcache, awConvFunc, mthVisResampler,//visResampler,
			     /*True */doPointing, doPBCorr, 
			     tile_p, paStep_p, pbLimit_p, True);
    ((AWProjectWBFT *)ft_p)->setObservatoryLocation(mLocation_p);
    //
    // Explicit type casting of ft_p does not look good.  It does not
    // pick up the setPAIncrement() method of PBWProjectFT without
    // this
    //
    // os << LogIO::NORMAL << "Setting PA increment to " << parAngleInc_p.getValue("deg") << " deg" << endl;
    ((AWProjectFT *)ft_p)->setPAIncrement(parAngleInc_p);

    if (doPointing) 
      {
	try
	  {
	    // Warn users we are have temporarily disabled pointing cal
	    //	    throw(AipsError("Pointing calibration temporarily disabled (gmoellen 06Nov20)."));
	    //  TBD: Bring this up-to-date with new VisCal mechanisms
	    VisSet elVS(*rvi_p);
	    epJ = new EPJones(elVS, *ms_p);
	    RecordDesc applyRecDesc;
	    applyRecDesc.addField("table", TpString);
	    applyRecDesc.addField("interp",TpString);
	    Record applyRec(applyRecDesc);
	    applyRec.define("table",epJTableName_p);
	    applyRec.define("interp", "nearest");
	    epJ->setApply(applyRec);
	    ((AWProjectFT *)ft_p)->setEPJones(epJ);
	  }
	catch(AipsError& x)
	  {
	    //
	    // Add some more useful info. to the message and translate
	    // the generic AipsError exception object to a more specific
	    // SynthesisError object.
	    //
	    String mesg = x.getMesg();
	    mesg += ". Error in loading pointing offset table.";
	    SynthesisError err(mesg);
	    throw(err);
	  }
      }
    AlwaysAssert(ft_p, AipsError);
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);

  }
  else if (ftmachine_p == "awp")
    {
      if (wprojPlanes_p<=1)
	{
	  os << LogIO::NORMAL
	     << "You are using wprojplanes=1. Doing co-planar imaging (no w-projection needed)" 
	     << LogIO::POST;
	  os << LogIO::NORMAL << "Performing A-Projection" << LogIO::POST; // Loglevel PROGRESS
	}
      if((wprojPlanes_p>1)&&(wprojPlanes_p<64)) 
	{
	  os << LogIO::WARN
	     << "No. of w-planes set too low for W projection - recommend at least 128"
	     << LogIO::POST;
	  os << LogIO::NORMAL << "Performing AW-Projection"
	     << LogIO::POST; // Loglevel PROGRESS
	}
      // if(!gvp_p) 
      // 	{
      // 	  os << LogIO::NORMAL // Loglevel INFO
      // 	     << "Using defaults for primary beams used in gridding" << LogIO::POST;
      // 	  gvp_p = new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p,
      // 				 skyPosThreshold_p);
      // 	}
      //      CountedPtr<ATerm> evlaAperture = new EVLAAperture();
      useDoublePrecGrid=False;
      CountedPtr<ATerm> apertureFunction = createTelescopeATerm(*ms_p);
      CountedPtr<ConvolutionFunction> awConvFunc=new AWConvFunc(apertureFunction);
      CountedPtr<VisibilityResamplerBase> visResampler = new AWVisResampler();
      CountedPtr<VisibilityResamplerBase> mthVisResampler = new MultiThreadedVisibilityResampler(useDoublePrecGrid,
												 visResampler);
      CountedPtr<CFCache> cfcache=new CFCache();
      cfcache->setCacheDir(cfCacheDirName_p.data());
      cfcache->initCache();
      ft_p = new AWProjectFT(wprojPlanes_p, cache_p/2,
			     cfcache, awConvFunc, mthVisResampler,
			     doPointing, doPBCorr,
			     tile_p, pbLimit_p, True);
      ((AWProjectFT *)ft_p)->setObservatoryLocation(mLocation_p);
      //
      // Explicit type casting of ft_p does not look good.  It does not
      // pick up the setPAIncrement() method of PBWProjectFT without
      // this
      //
      Quantity paInc(paStep_p,"deg");
      // os << LogIO::NORMAL << "Setting PA increment to " 
      // 	 << paInc.getValue("deg") << " deg" << endl;
      ((AWProjectFT *)ft_p)->setPAIncrement(paInc);

      if (doPointing) 
	{
	  try
	    {
	      VisSet elVS(*rvi_p);
	      epJ = new EPJones(elVS, *ms_p);
	      RecordDesc applyRecDesc;
	      applyRecDesc.addField("table", TpString);
	      applyRecDesc.addField("interp",TpString);
	      Record applyRec(applyRecDesc);
	      applyRec.define("table",epJTableName_p);
	      applyRec.define("interp", "nearest");
	      epJ->setApply(applyRec);
	      ((AWProjectFT *)ft_p)->setEPJones(epJ);
	  }
	  catch(AipsError& x)
	    {
	      //
	      // Add some more useful info. to the message and translate
	      // the generic AipsError exception object to a more specific
	      // SynthesisError object.
	      //
	      String mesg = x.getMesg();
	      mesg += ". Error in loading pointing offset table.";
	      SynthesisError err(mesg);
	      throw(err);
	    }
	}
      AlwaysAssert(ft_p, AipsError);
      cft_p = new SimpleComponentFTMachine();
      AlwaysAssert(cft_p, AipsError);
    }
  //===============================================================
  // A-Projection FTMachine code end here
  //===============================================================
  else {
    os << LogIO::NORMAL // Loglevel INFO
       << "Performing interferometric gridding..."
       << LogIO::POST;
    os << LogIO::NORMAL1 // gridfunction_p is too cryptic for most users.
       << "...with convolution function " << gridfunction_p << LogIO::POST;
    // Now make the FTMachine
    if(facets_p>1) {
      os << LogIO::NORMAL // Loglevel INFO
         << "Multi-facet Fourier transforms will use specified common tangent point:"
	 << LogIO::POST;
      os << LogIO::NORMAL << tangentPoint() << LogIO::POST; // Loglevel INFO
      ft_p = new GridFT(cache_p / 2, tile_p, gridfunction_p, mLocation_p,
                        phaseCenter_p, padding, False, useDoublePrecGrid);
      
    }
    else {
      os << LogIO::DEBUG1
         << "Single facet Fourier transforms will use image center as tangent points"
	 << LogIO::POST;
      ft_p = new GridFT(cache_p/2, tile_p, gridfunction_p, mLocation_p,
			padding, False, useDoublePrecGrid);

    }
    AlwaysAssert(ft_p, AipsError);
    
    cft_p = new SimpleComponentFTMachine();
    AlwaysAssert(cft_p, AipsError);
    
  }

  /******* Start MTFT code ********/
  // MultiTermFT is a container for an FTMachine of any type.
  //    It will apply Taylor-polynomial weights during gridding and degridding
  //    and will do multi-term grid-correction (normalizations).
  //    (1) ft_p already holds an FT of the correct type
  //    (2) If nterms>1, create a new MultiTermFT using ft_p, and reassign ft_p. 
  // Currently, Multi-Term applies only to wideband imaging.
  if( ntaylor_p > 1 )
  { 
    //cout << "Creating a Multi-Term FT machine containing " << ftmachine_p << endl;
     FTMachine *tempftm = new MultiTermFT(ft_p, ftmachine_p, ntaylor_p, reffreq_p);
     ft_p = tempftm;
  }
  /******* End MTFT code ********/

  ft_p->setSpw(dataspectralwindowids_p, freqFrameValid_p);
  ft_p->setFreqInterpolation(freqInterpMethod_p);
  if(doTrackSource_p){
    ft_p->setMovingSource(trackDir_p);
  }
  ft_p->setSpwChanSelection(spwchansels_p);

  return True;
}

String Imager::tangentPoint()
{
  MVAngle mvRa=phaseCenter_p.getAngle().getValue()(0);
  MVAngle mvDec=phaseCenter_p.getAngle().getValue()(1);
  ostringstream oos;
  oos << "     ";
  Int widthRA=20;
  Int widthDec=20;
  oos.setf(ios::left, ios::adjustfield);
  oos.width(widthRA);  oos << mvRa(0.0).string(MVAngle::TIME,8);
  oos.width(widthDec); oos << mvDec.string(MVAngle::DIG2,8);
  oos << "     "
      << MDirection::showType(phaseCenter_p.getRefPtr()->getType());
  return String(oos);
}

Bool Imager::removeTable(const String& tablename) {
  
  LogIO os(LogOrigin("imager", "removeTable()", WHERE));
  
  if(Table::isReadable(tablename)) {
    if (! Table::isWritable(tablename)) {
      os << LogIO::SEVERE << "Table " << tablename
	 << " is not writable!: cannot alter it" << LogIO::POST;
      return False;
    }
    else {
      if (Table::isOpened(tablename)) {
	os << LogIO::SEVERE << "Table " << tablename
	   << " is already open in the process. It needs to be closed first"
	   << LogIO::POST;
	  return False;
      } else {
	Table table(tablename, Table::Update);
	if (table.isMultiUsed()) {
	  os << LogIO::SEVERE << "Table " << tablename
	     << " is already open in another process. It needs to be closed first"
	     << LogIO::POST;
	    return False;
	} else {
	  Table table(tablename, Table::Delete);
	}
      }
    }
  }
  return True;
}

Bool Imager::updateSkyModel(const Vector<String>& model,
			    const String complist) {
  LogIO os(LogOrigin("imager", "updateSkyModel()", WHERE));
  if(redoSkyModel_p)
    throw(AipsError("Programming error: update skymodel is called without a valid skymodel"));
  Bool coordMatch=True; 
  for (Int thismodel=0;thismodel<Int(model.nelements());++thismodel) {
    CoordinateSystem cs=(sm_p->image(thismodel)).coordinates();
    coordMatch= coordMatch || checkCoord(cs, model(thismodel));
    ///return False if any fails anyways
    if(!coordMatch)
      return False;
    if(model(thismodel)=="") {
      os << LogIO::SEVERE << "Need a name for model "
	 << model << LogIO::POST;
      return False;
    }
    else {
      if(!Table::isReadable(model(thismodel))) {
	os << LogIO::SEVERE << model(thismodel) << "is unreadable"
	   << model << LogIO::POST;
	return False;
      }
    }
    images_p[thismodel]=0;
    images_p[thismodel]=new PagedImage<Float>(model(thismodel));
    AlwaysAssert(!images_p[thismodel].null(), AipsError);
    sm_p->updatemodel(thismodel, *images_p[thismodel]);
  } 
  if((complist !="") && Table::isReadable(complist)){
      ComponentList cl(Path(complist), True);
      sm_p->updatemodel(cl);
    }
  return True;
}

Bool Imager::createSkyEquation(const String complist) 
{
  Vector<String> image;
  Vector<String> mask;
  Vector<String> fluxMask;
  Vector<Bool> fixed;
  return createSkyEquation(image, fixed, mask, fluxMask, complist);
}

Bool Imager::createSkyEquation(const Vector<String>& image,
			       const String complist) 
{
  Vector<Bool> fixed(image.nelements()); fixed=False;
  Vector<String> mask(image.nelements()); mask="";
  Vector<String> fluxMask(image.nelements()); fluxMask="";
  return createSkyEquation(image, fixed, mask, fluxMask, complist);
}

Bool Imager::createSkyEquation(const Vector<String>& image, const Vector<Bool>& fixed,
			       const String complist) 
{
  Vector<String> mask(image.nelements()); mask="";
  Vector<String> fluxMask(image.nelements()); fluxMask="";
  return createSkyEquation(image, fixed, mask, fluxMask, complist);
}

Bool Imager::createSkyEquation(const Vector<String>& image, const Vector<Bool>& fixed,
			       const Vector<String>& mask,
			       const String complist) 
{
  Vector<String> fluxMask(image.nelements()); fluxMask="";
  return createSkyEquation(image, fixed, mask, fluxMask, complist);
}

Bool Imager::createSkyEquation(const Vector<String>& image,
			       const Vector<Bool>& fixed,
			       const Vector<String>& mask,
			       const Vector<String>& fluxMask,
			       const String complist)
{
 
  if(!valid()) return False;

  LogIO os(LogOrigin("imager", "createSkyEquation()", WHERE));
  
  // If there is no sky model, we'll make one:

  if(sm_p==0) {
    if((facets_p >1)){
      // Support serial and parallel specializations
      setWFCleanImageSkyModel();
    }
    else if(ntaylor_p>1){
      // Init the msmfs sky-model so that Taylor weights are triggered in CubeSkyEqn
      sm_p = new WBCleanImageSkyModel(ntaylor_p, 1 ,reffreq_p);
    }
    else {
      sm_p = new CleanImageSkyModel();
    }
  }
  AlwaysAssert(sm_p, AipsError);

  // Read data pol rep from the MS. This is also done in imagecoordinates().
  // This code segment is repeated here because imagecoordinates
  // is not always called before SkyModel is made (im.ft).
  // Note : This should go into a function.
  {
    ROMSColumns msc(*ms_p);
    Vector<String> polType=msc.feed().polarizationType()(0);
    if (polType(0)!="X" && polType(0)!="Y" &&
        polType(0)!="R" && polType(0)!="L") {
       os << LogIO::WARN << "Unknown stokes types in feed table: ["
       << polType(0) << ", " << polType(1) << "]" << endl
       << "Results open to question!" << LogIO::POST;
    }
  
    if (polType(0)=="X" || polType(0)=="Y") {
      polRep_p=SkyModel::LINEAR;
      os << LogIO::DEBUG1 
          << "Preferred polarization representation is linear" << LogIO::POST;
    }
    else {
      polRep_p=SkyModel::CIRCULAR;
      os << LogIO::DEBUG1
         << "Preferred polarization representation is circular" << LogIO::POST;
    }
  }// end of reading-in polRep_p from the MS

  // Send the data correlation type to the SkyModel
  // This needs to be done before the first call to 'ImageSkyModel::cImage()'
   os << LogIO::DEBUG1 << "Data PolRep in Imager2.cc::createSkyEquation : " << polRep_p << LogIO::POST;
   sm_p->setDataPolFrame(polRep_p);


  // Add the componentlist
  if(complist!="") {
    if(!Table::isReadable(complist)) {
      os << LogIO::SEVERE << "ComponentList " << complist
	 << " not readable" << LogIO::POST;
      return False;
    }
    componentList_p=new ComponentList(complist, True);
    if(componentList_p==0) {
      os << LogIO::SEVERE << "Cannot create ComponentList from " << complist
	 << LogIO::POST;
      return False;
    }
    if(!sm_p->add(*componentList_p)) {
      os << LogIO::SEVERE << "Cannot add ComponentList " << complist
	 << " to SkyModel" << LogIO::POST;
      return False;
    }
    os << LogIO::NORMAL // Loglevel INFO
       << "Processing after subtracting componentlist " << complist << LogIO::POST;
  }
  else {
    componentList_p=0;
  }
 
  // Make image with the required shape and coordinates only if
  // they don't exist yet
  nmodels_p=image.nelements();

  // Remove degenerate case (due to user interface?)
  if((nmodels_p==1)&&(image(0)=="")) {
    nmodels_p=0;
  }
  if(nmodels_p>0) {
    images_p.resize(nmodels_p); 
    masks_p.resize(nmodels_p);  
    fluxMasks_p.resize(nmodels_p); 
    residuals_p.resize(nmodels_p); 
    for (Int model=0;model<Int(nmodels_p);++model) {
      if(image(model)=="") {
	os << LogIO::SEVERE << "Need a name for model "
	   << model << LogIO::POST;
	return False;
      }
      else {
	if(!Table::isReadable(image(model))) {
	  if(!assertDefinedImageParameters()) return False;
	  make(image(model));
	}
      }
      images_p[model]=0;
      images_p[model]=new PagedImage<Float>(image(model));
      AlwaysAssert(!images_p[model].null(), AipsError);

      //Determining the number of XFR
      Int numOfXFR=nmodels_p+1;
      if(datafieldids_p.nelements() >0)
	numOfXFR=datafieldids_p.nelements()*nmodels_p + 1;
      if(squintType_p != BeamSquint::NONE){
	if(parAngleInc_p.getValue("deg") >0 ){
	  numOfXFR= numOfXFR* Int(360/parAngleInc_p.getValue("deg"));
	}	
	else{
	numOfXFR= numOfXFR*10;
	}
      }
      if((sm_p->add(*images_p[model], numOfXFR))!=model) {
	os << LogIO::SEVERE << "Error adding model " << model << LogIO::POST;
	return False;
      }
 
      fluxMasks_p[model]=0;
      if(fluxMask(model)!=""&&Table::isReadable(fluxMask(model))) {
	fluxMasks_p[model]=new PagedImage<Float>(fluxMask(model));
	AlwaysAssert(!fluxMasks_p[model].null(), AipsError);
        if(!sm_p->addFluxMask(model, *fluxMasks_p[model])) {
	  os << LogIO::SEVERE << "Error adding flux mask " << model
	     << " : " << fluxMask(model) << LogIO::POST;
	  return False;
	}
      }
      residuals_p[model]=0;
    }
    addMasksToSkyEquation(mask, fixed);
  }
  
  // Always need a VisSet and an FTMachine
  if (!ft_p)
    createFTMachine();
  
  // Now set up the SkyEquation
  AlwaysAssert(sm_p, AipsError);
  // AlwaysAssert(vs_p, AipsError);
  AlwaysAssert(rvi_p, AipsError);
  AlwaysAssert(ft_p, AipsError);
  AlwaysAssert(cft_p, AipsError);

  // This block determines which SkyEquation is to be used.
  // We are using a mf* algorithm and there is more than one image
  if (doMultiFields_p && multiFields_p) {
    // Mosaicing case
    if(doVP_p){
      //bypassing the minimum size FT stuff as its broken till its fixed
      //se_p=new MosaicSkyEquation(*sm_p, *vs_p, *ft_p, *cft_p);
      
      setSkyEquation();
      if(ft_p->name() != "MosaicFT") 
	sm_p->mandateFluxScale(0);
      os << LogIO::NORMAL // Loglevel INFO
         << "Mosaicing multiple fields with simple sky equation" << LogIO::POST;
    }
    // mosaicing with no vp correction
    else{
      setSkyEquation();
      os << LogIO::NORMAL // Loglevel INFO
         << "Processing multiple fields with simple sky equation" << LogIO::POST;
      os << LogIO::WARN
         << "Voltage Pattern is not set: will not correct for primary beam"
	 << LogIO::POST;
      doMultiFields_p=False;
    }
  }
  // We are not using an mf* algorithm or there is only one image
  else {
    // Support serial and parallel specializations
    if((facets_p >1)){
	setSkyEquation();
	//se_p=new SkyEquation(*sm_p, *vs_p, *ft_p, *cft_p);
	os << LogIO::NORMAL // Loglevel INFO
           << "Processing multiple facets with simple sky equation" << LogIO::POST;
    }
    // Mosaicing
    else if(doVP_p) {
      //Bypassing the mosaicskyequation to the slow version for now.
      //      se_p=new MosaicSkyEquation(*sm_p, *vs_p, *ft_p, *cft_p);
      
      setSkyEquation();
      if(ft_p->name() != "MosaicFT") 
	sm_p->mandateFluxScale(0);
      os << LogIO::NORMAL // Loglevel PROGRESS
         << "Mosaicing single field with simple sky equation" << LogIO::POST;      
    }
    // Default
    else {
      setSkyEquation();
      os << LogIO::DEBUG1
         << "Processing single field with simple sky equation" << LogIO::POST;    
    } 
  }
  //os.localSink().flush();
  //For now force none sault weighting with mosaic ft machine
  
  String scaleType=scaleType_p;
  if(ft_p->name()=="MosaicFT")
    scaleType="NONE";
  
  se_p->setImagePlaneWeighting(scaleType, minPB_p, constPB_p);
  se_p->doFlatNoise(flatnoise_p);
  
  AlwaysAssert(se_p, AipsError);

  // Now add any SkyJones that are needed
  if(doVP_p && (ft_p->name()!="MosaicFT")) {
    if (doDefaultVP_p) {
      vp_p=new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p, skyPosThreshold_p);
    } else { 
      Table vpTable( vpTableStr_p ); 
      vp_p=new VPSkyJones(*ms_p, vpTable, parAngleInc_p, squintType_p, skyPosThreshold_p);
    }
    se_p->setSkyJones(*vp_p);
  }
  else {
    vp_p=0;
  }
  return True;  
}

Bool Imager::addResiduals(const Vector<String>& imageNames) {
  Bool retval=True;
  residuals_p.resize(imageNames.nelements(), True, False);
  for (Int thismodel=0;thismodel<Int(imageNames.nelements());++thismodel) {
    if(imageNames(thismodel)!="") {
      residuals_p[thismodel]=0;
      if(Table::isWritable(imageNames(thismodel))){
	residuals_p[thismodel]=new PagedImage<Float>(imageNames(thismodel));
	if(!(residuals_p[thismodel]->shape()).isEqual(images_p[thismodel]->shape())){
	  residuals_p[thismodel]=0;
	  removeTable(imageNames(thismodel));
	}
      }
      if(residuals_p[thismodel].null()){
	if(Table::isReadable(imageNames(thismodel)))
	   removeTable(imageNames(thismodel));
	residuals_p[thismodel]=
	  new PagedImage<Float> (TiledShape(images_p[thismodel]->shape(), 
					    images_p[thismodel]->niceCursorShape()),
			       images_p[thismodel]->coordinates(),
				 imageNames(thismodel));
	AlwaysAssert(!residuals_p[thismodel].null(), AipsError);
	residuals_p[thismodel]->setUnits(Unit("Jy/beam"));
      }
      if(residuals_p[thismodel].null()) 
	retval=False;
    }
    else{
      retval=False;
    }
  }
  return retval;
} 
// Tell the sky model to use the specified images as the residuals    
Bool Imager::addResidualsToSkyEquation(const Vector<String>& imageNames) {
  
  addResiduals(imageNames);
  for (Int thismodel=0;thismodel<Int(imageNames.nelements());++thismodel) {
    if(imageNames(thismodel)!="") 
      sm_p->addResidual(thismodel, *residuals_p[thismodel]);   
  }
  return True;
} 

void Imager::destroySkyEquation() 
{
  if(se_p) delete se_p; se_p=0;
  if(sm_p) delete sm_p; sm_p=0;
  if(vp_p) delete vp_p; vp_p=0;
  if(gvp_p) delete gvp_p; gvp_p=0;
  
  if(componentList_p) delete componentList_p; componentList_p=0;
 
  
  for (Int model=0;model<Int(nmodels_p); ++model) {
    //As these are CountedPtrs....just assigning them to NULL 
    //get the objects out of context
    if(!images_p[model].null())  
      images_p[model]=0;
    
    if(!masks_p[model].null()) 
      masks_p[model]=0;
 
   if(!fluxMasks_p[model].null()) 
     fluxMasks_p[model]=0;
    
   if(!residuals_p[model].null()) 
     residuals_p[model]=0;
  }
  
  redoSkyModel_p=True;
}

Bool Imager::assertDefinedImageParameters() const
{
  LogIO os(LogOrigin("imager", "if(!assertDefinedImageParameters()", WHERE));
  if(!setimaged_p) { 
    os << LogIO::SEVERE << "Image parameters not yet set: use im.defineimage."
       << LogIO::POST;
    return False;
  }
  return True;
}

Bool Imager::valid() const {
  LogIO os(LogOrigin("imager", "if(!valid()) return False", WHERE));
  if(ms_p.null()) {
    os << LogIO::SEVERE << "Program logic error: MeasurementSet pointer ms_p not yet set"
       << LogIO::POST;
    return False;
  }
  if(mssel_p.null()) {
    os << LogIO::SEVERE << "Program logic error: MeasurementSet pointer mssel_p not yet set"
       << LogIO::POST;
    return False;
  }
  if(!rvi_p) {
    os << LogIO::SEVERE << "Program logic error: VisibilityIterator not yet set"
       << LogIO::POST;
    return False;
  }
  return True;
}


Bool Imager::addMasksToSkyEquation(const Vector<String>& mask, const Vector<Bool>& fixed){
  LogIO os(LogOrigin("imager", "addMasksToSkyEquation()", WHERE));

  for(Int model=0 ;model < nmodels_p; ++model){

    
    if((Int(fixed.nelements())>model) && fixed(model)) {
      os << LogIO::NORMAL // Loglevel INFO
         << "Model " << model << " will be held fixed" << LogIO::POST;
      sm_p->fix(model);
    }     
    /*
    if(!(masks_p[model].null())) delete masks_p[model];
    masks_p[model]=0;
    */
      if(mask(model)!=""&&Table::isReadable(mask(model))) {
	masks_p[model]=new PagedImage<Float>(mask(model));
	AlwaysAssert(!masks_p[model].null(), AipsError);
        if(!sm_p->addMask(model, *masks_p[model])) {
	  os << LogIO::SEVERE << "Error adding mask " << model
	     << " : " << mask(model) << LogIO::POST;
	  return False;
	}
      }
  }
  return True;
}


Bool Imager::openSubTables(){

 antab_p=Table(ms_p->antennaTableName(),
	       TableLock(TableLock::UserNoReadLocking));
 datadesctab_p=Table(ms_p->dataDescriptionTableName(),
	       TableLock(TableLock::UserNoReadLocking));
 feedtab_p=Table(ms_p->feedTableName(),
		 TableLock(TableLock::UserNoReadLocking));
 fieldtab_p=Table(ms_p->fieldTableName(),
		  TableLock(TableLock::UserNoReadLocking));
 obstab_p=Table(ms_p->observationTableName(),
		TableLock(TableLock::UserNoReadLocking));
 poltab_p=Table(ms_p->polarizationTableName(),
		TableLock(TableLock::UserNoReadLocking));
 proctab_p=Table(ms_p->processorTableName(),
		TableLock(TableLock::UserNoReadLocking));
 spwtab_p=Table(ms_p->spectralWindowTableName(),
		TableLock(TableLock::UserNoReadLocking));
 statetab_p=Table(ms_p->stateTableName(),
		TableLock(TableLock::UserNoReadLocking));

 if(Table::isReadable(ms_p->dopplerTableName()))
   dopplertab_p=Table(ms_p->dopplerTableName(),
		      TableLock(TableLock::UserNoReadLocking));

 if(Table::isReadable(ms_p->flagCmdTableName()))
   flagcmdtab_p=Table(ms_p->flagCmdTableName(),
		      TableLock(TableLock::UserNoReadLocking));
 if(Table::isReadable(ms_p->freqOffsetTableName()))
   freqoffsettab_p=Table(ms_p->freqOffsetTableName(),
			 TableLock(TableLock::UserNoReadLocking));

 if(ms_p->isWritable()){
   if(!(Table::isReadable(ms_p->historyTableName()))){
     // setup a new table in case its not there
     TableRecord &kws = ms_p->rwKeywordSet();
     SetupNewTable historySetup(ms_p->historyTableName(),
				MSHistory::requiredTableDesc(),Table::New);
     kws.defineTable(MS::keywordName(MS::HISTORY), Table(historySetup));
     
   }
   historytab_p=Table(ms_p->historyTableName(),
		      TableLock(TableLock::UserNoReadLocking), Table::Update);
 }
 if(Table::isReadable(ms_p->pointingTableName()))
   pointingtab_p=Table(ms_p->pointingTableName(), 
		       TableLock(TableLock::UserNoReadLocking));
 
 if(Table::isReadable(ms_p->sourceTableName()))
   sourcetab_p=Table(ms_p->sourceTableName(),
		     TableLock(TableLock::UserNoReadLocking));

 if(Table::isReadable(ms_p->sysCalTableName()))
 syscaltab_p=Table(ms_p->sysCalTableName(),
		   TableLock(TableLock::UserNoReadLocking));
 if(Table::isReadable(ms_p->weatherTableName()))
   weathertab_p=Table(ms_p->weatherTableName(),
		      TableLock(TableLock::UserNoReadLocking));
 if(ms_p->isWritable()){
   hist_p= new MSHistoryHandler(*ms_p, "imager");
 }
return True;

}

Bool Imager::lock(){

  Bool ok; 
  ok=True;
  if(lockCounter_p == 0){

    ok= ok && (ms_p->lock());
    ok= ok && antab_p.lock(False);
    ok= ok && datadesctab_p.lock(False);
    ok= ok && feedtab_p.lock(False);
    ok= ok && fieldtab_p.lock(False);
    ok= ok && obstab_p.lock(False);
    ok= ok && poltab_p.lock(False);
    ok= ok && proctab_p.lock(False);
    ok= ok && spwtab_p.lock(False);
    ok= ok && statetab_p.lock(False);
    if(!dopplertab_p.isNull())
      ok= ok && dopplertab_p.lock(False);
    if(!flagcmdtab_p.isNull())
      ok= ok && flagcmdtab_p.lock(False);
    if(!freqoffsettab_p.isNull())
      ok= ok && freqoffsettab_p.lock(False);
    if(!historytab_p.isNull())
      ok= ok && historytab_p.lock(False);
    if(!pointingtab_p.isNull())
      ok= ok && pointingtab_p.lock(False);
    if(!sourcetab_p.isNull())
      ok= ok && sourcetab_p.lock(False);
    if(!syscaltab_p.isNull())
      ok= ok && syscaltab_p.lock(False);
    if(!weathertab_p.isNull())
      ok= ok && weathertab_p.lock(False);
 
  }
  ++lockCounter_p;

  return ok ; 
}

Bool Imager::unlock(){

  if(lockCounter_p==1){
    ms_p->unlock();
    antab_p.unlock();
    datadesctab_p.unlock();
    feedtab_p.unlock();
    fieldtab_p.unlock();
    obstab_p.unlock();
    poltab_p.unlock();
    proctab_p.unlock();
    spwtab_p.unlock();
    statetab_p.unlock();
    if(!dopplertab_p.isNull())
      dopplertab_p.unlock();
    if(!flagcmdtab_p.isNull())
      flagcmdtab_p.unlock();
    if(!freqoffsettab_p.isNull())
    freqoffsettab_p.unlock();
    if(!historytab_p.isNull())
      historytab_p.unlock();
    if(!pointingtab_p.isNull())
      pointingtab_p.unlock();
    if(!sourcetab_p.isNull())
      sourcetab_p.unlock();
    if(!syscaltab_p.isNull())
      syscaltab_p.unlock();
    if(!weathertab_p.isNull())
      weathertab_p.unlock();
  }
  for (Int thismodel=0;thismodel<Int(images_p.nelements());++thismodel) {
    if ((images_p.nelements() > uInt(thismodel)) && (!images_p[thismodel].null())) {
      images_p[thismodel]->table().relinquishAutoLocks(True);
      images_p[thismodel]->table().unlock();
    }
    if ((residuals_p.nelements()> uInt(thismodel)) && (!residuals_p[thismodel].null())) {
      residuals_p[thismodel]->table().relinquishAutoLocks(True);
      residuals_p[thismodel]->table().unlock();
    }
  }
  if(lockCounter_p > 0 )
    --lockCounter_p;
  return True ; 
}

Bool Imager::selectDataChannel(Vector<Int>& spectralwindowids, 
			       String& dataMode, 
			       Vector<Int>& dataNchan, 
			       Vector<Int>& dataStart, Vector<Int>& dataStep,
			       MRadialVelocity& mDataStart, 
			       MRadialVelocity& mDataStep){



  LogIO os(LogOrigin("Imager", "selectDataChannel()", WHERE));



  if(dataMode=="channel") {
      if (dataNchan.nelements() != spectralwindowids.nelements()){
	if(dataNchan.nelements()==1){
	  dataNchan.resize(spectralwindowids.nelements(), True);
	  for(uInt k=1; k < spectralwindowids.nelements(); ++k){
	    dataNchan[k]=dataNchan[0];
	  }
	}
	else{
	  os << LogIO::SEVERE 
	     << "Vector of nchan has to be of size 1 or be of the same shape as spw " 
	     << LogIO::POST;
	  return False; 
	}
      }
      if (dataStart.nelements() != spectralwindowids.nelements()){
	if(dataStart.nelements()==1){
	  dataStart.resize(spectralwindowids.nelements(), True);
	  for(uInt k=1; k < spectralwindowids.nelements(); ++k){
	    dataStart[k]=dataStart[0];
	  }
	}
	else{
	  os << LogIO::SEVERE 
	     << "Vector of start has to be of size 1 or be of the same shape as spw " 
	     << LogIO::POST;
	  return False; 
	}
      }
      if (dataStep.nelements() != spectralwindowids.nelements()){
	if(dataStep.nelements()==1){
	  dataStep.resize(spectralwindowids.nelements(), True);
	  for(uInt k=1; k < spectralwindowids.nelements(); ++k){
	    dataStep[k]=dataStep[0];
	  }
	}
	else{
	  os << LogIO::SEVERE 
	     << "Vector of step has to be of size 1 or be of the same shape as spw " 
	     << LogIO::POST;
	  return False; 
	}
      }

      if(spectralwindowids.nelements()>0) {
	Int nch=0;
	for(uInt i=0;i<spectralwindowids.nelements();++i) {
	  Int spwid=spectralwindowids(i);
	  Int numberChan=rvi_p->msColumns().spectralWindow().numChan()(spwid);
	  if(dataStart[i]<0) {
	    os << LogIO::SEVERE << "Illegal start pixel = " 
	       << dataStart[i]  << " for spw " << spwid
	       << LogIO::POST;
	    return False;
	  }
	 
	  if(dataNchan[i]<=0){ 
	    if(dataStep[i] <= 0)
	      dataStep[i]=1;
	    nch=(numberChan-dataStart[i])/Int(dataStep[i])+1;
	  }
	  else nch = dataNchan[i];
	  while((nch*dataStep[i]+dataStart[i]) > numberChan){
	    --nch;
	  }
	  Int end = Int(dataStart[i]) + Int(nch-1) * Int(dataStep[i]);
	  if(end < 0 || end > (numberChan)-1) {
	    os << LogIO::SEVERE << "Illegal step pixel = " << dataStep[i]
	       << " for spw " << spwid
	       << "\n end channel " << end 
	       << " is out of range " << dataStart[i] << " to " 
	       << (numberChan-1)
	       << LogIO::POST;
	    return False;
	  }

          os << LogIO::NORMAL // Too contentious for DEBUG1
             << "Selecting ";
          if(nch > 1)
            os << nch << " channels, starting at "
               << dataStart[i]  << ", stepped by " << dataStep[i] << ",";
          else
            os << "channel " << dataStart[i];
          os << " for spw " << spwid << LogIO::POST;
	  rvi_p->selectChannel(1, Int(dataStart[i]), Int(nch),
				     Int(dataStep[i]), spwid);
	  dataNchan[i]=nch;
	}
      }	else {
        VisBufferAutoPtr vb (rvi_p);
        rvi_p->originChunks ();
        Int numberChan=vb->msColumns().spectralWindow().numChan()(0);

	if(dataNchan[0]<=0){
	  if(dataStep[0] <=0)
	    dataStep[0]=1;
	  dataNchan[0]=(numberChan-dataStart[0])/Int(dataStep[0])+1;
	  
	}
	while((dataNchan[0]*dataStep[0]+dataStart[0]) > numberChan)
	  --dataNchan[0];

	Int end = Int(dataStart[0]) + Int(dataNchan[0]-1) 
	  * Int(dataStep[0]);
	if(end < 0 || end > (numberChan)-1) {
	  os << LogIO::SEVERE << "Illegal step pixel = " << dataStep[0]
	     << "\n end channel " << end << " is out of range 1 to " 
	     << (numberChan-1)
	     << LogIO::POST;
	  return False;
	}
	os << LogIO::NORMAL << "Selecting "<< dataNchan[0] // Loglevel INFO
	   << " channels, starting at visibility channel "
	 << dataStart[0]  << " stepped by "
	   << dataStep[0] << LogIO::POST;
      }
  }
  else if (dataMode=="velocity") {
      MVRadialVelocity mvStart(mDataStart.get("m/s"));
      MVRadialVelocity mvStep(mDataStep.get("m/s"));
      MRadialVelocity::Types
	vType((MRadialVelocity::Types)mDataStart.getRefPtr()->getType());
      os << LogIO::NORMAL << "Selecting "<< dataNchan[0] // Loglevel INFO
	 << " channels, starting at radio velocity " << mvStart
	 << " stepped by " << mvStep << ", reference frame is "
	 << MRadialVelocity::showType(vType) << LogIO::POST;
      rvi_p->selectVelocity(Int(dataNchan[0]), mvStart, mvStep,
				  vType, MDoppler::RADIO);
  }
  else if (dataMode=="opticalvelocity") {
      MVRadialVelocity mvStart(mDataStart.get("m/s"));
      MVRadialVelocity mvStep(mDataStep.get("m/s"));
      MRadialVelocity::Types
	vType((MRadialVelocity::Types)mDataStart.getRefPtr()->getType());
      os << LogIO::NORMAL << "Selecting "<< dataNchan[0] // Loglevel INFO
	 << " channels, starting at optical velocity " << mvStart
	 << " stepped by " << mvStep << ", reference frame is "
	 << MRadialVelocity::showType(vType) << LogIO::POST;
      rvi_p->selectVelocity(Int(dataNchan[0]), mvStart, mvStep,
				  vType, MDoppler::OPTICAL);
  }

  return True;

}


Bool Imager::checkCoord(const CoordinateSystem& coordsys,  
			const String& imageName){ 

  PagedImage<Float> image(imageName);
  CoordinateSystem imageCoord= image.coordinates();
  Vector<Int> imageShape= image.shape().asVector();

  if(imageShape.nelements() > 3){
    if(imageShape(3) != imageNchan_p)
      return False;
  }
  else{
    if(imageNchan_p >1)
      return False;
  }

  if(imageShape.nelements() > 2){
    if(imageShape(2) != npol_p)
      return False;
  } 
  else{
    if(npol_p > 1)
      return False;
  }
  if(imageShape(0) != nx_p)
    return False;
  if(imageShape(1) != ny_p)
    return False;


 
  if(!imageCoord.near(coordsys)){
    return False;
  }
  
  /*
  DirectionCoordinate dir1(coordsys.directionCoordinate(0));
  DirectionCoordinate dir2(imageCoord.directionCoordinate(0));
  if(dir1.increment()(0) != dir2.increment()(0))
    return False;
  if(dir1.increment()(1) != dir2.increment()(1))
    return False;
  SpectralCoordinate sp1(coordsys.spectralCoordinate(2));
  SpectralCoordinate sp2(imageCoord.spectralCoordinate(2));
  if(sp1.increment()(0) != sp2.increment()(0))
    return False;
  */
  return True;
}

void Imager::setImageParam(Int& nx, Int& ny, Int& npol, Int& nchan){

  nx_p=nx;
  ny_p=ny;
  npol_p=npol;
  nchan_p=nchan;

}

void Imager::makeVisSet(MeasurementSet& ms, 
			Bool compress, Bool mosaicOrder){

  if(rvi_p) {
    delete rvi_p;
    rvi_p=0;
    wvi_p=0;
  }

  Block<Int> sort(0);
  if(mosaicOrder){
    sort.resize(4);
    sort[0] = MS::FIELD_ID;
    sort[1] = MS::ARRAY_ID;
    sort[2] = MS::DATA_DESC_ID;
    sort[3] = MS::TIME;
 
  }
  //else use default sort order
  else{
    sort.resize(4);
    sort[0] = MS::ARRAY_ID;
    sort[1] = MS::FIELD_ID;
    sort[2] = MS::DATA_DESC_ID;
    sort[3] = MS::TIME;
  }
  Matrix<Int> noselection;
  Double timeInterval=0;
  //if you want to use scratch col...make sure they are there
  if(useModelCol_p)
    VisSet(ms,sort,noselection,useModelCol_p,timeInterval,compress);
  
  if(imwgt_p.getType()=="none"){
      imwgt_p=VisImagingWeight("natural");
  }

  if(!useModelCol_p){
    rvi_p=new ROVisibilityIterator(ms, sort, timeInterval);
  }
  else{
    wvi_p=new VisibilityIterator(ms, sort, timeInterval);
    rvi_p=wvi_p;    
  }
  rvi_p->useImagingWeight(imwgt_p);

}
/*
void Imager::makeVisSet(MeasurementSet& ms, 
			Bool compress, Bool mosaicOrder){


  Block<Int> sort(0);
  if(mosaicOrder){
    sort.resize(4);
    sort[0] = MS::FIELD_ID;
    sort[1] = MS::ARRAY_ID;
    sort[2] = MS::DATA_DESC_ID;
    sort[3] = MS::TIME;
  }
  //else use default sort order
  else{
  
    sort.resize(4);
    sort[0] = MS::ARRAY_ID;
    sort[1] = MS::FIELD_ID;
    sort[2] = MS::DATA_DESC_ID;
    sort[3] = MS::TIME;
  }
  Matrix<Int> noselection;
  Double timeInterval=0;

  VisSet vs(ms,sort,noselection,useModelCol_p,timeInterval,compress);

}
*/

void Imager::writeHistory(LogIO& os){

  LogIO oslocal(LogOrigin("Imager", "writeHistory"));
  try{

    os.postLocally();
    if(!hist_p.null())
      hist_p->addMessage(os);
  }catch (AipsError x) {
    oslocal << LogIO::SEVERE << "Caught exception: " << x.getMesg()
	    << LogIO::POST;
  } 
}

void Imager::writeCommand(LogIO& os){

  LogIO oslocal(LogOrigin("Imager", "writeCommand"));
  try{
    os.postLocally();
    if(!hist_p.null())
      hist_p->cliCommand(os);
  }catch (AipsError x) {
    oslocal << LogIO::SEVERE << "Caught exception: " << x.getMesg()
	    << LogIO::POST;
  } 
}

Bool Imager::makePBImage(ImageInterface<Float>& pbImage, 
			 Bool useSymmetricBeam){

  LogIO os(LogOrigin("Imager", "makePBImage()", WHERE));
  CoordinateSystem imageCoord=pbImage.coordinates();
   Int spectralIndex=imageCoord.findCoordinate(Coordinate::SPECTRAL);
  SpectralCoordinate
    spectralCoord=imageCoord.spectralCoordinate(spectralIndex);
  Vector<String> units(1); units = "Hz";
  spectralCoord.setWorldAxisUnits(units);	
  Vector<Double> spectralWorld(1);
  Vector<Double> spectralPixel(1);
  spectralPixel(0) = 0;
  spectralCoord.toWorld(spectralWorld, spectralPixel);  
  Double freq  = spectralWorld(0);
  Quantity qFreq( freq, "Hz" );
  String telName=imageCoord.obsInfo().telescope();
  if(telName=="UNKNOWN"){
    os << LogIO::SEVERE << "Telescope encoded in image in not known " 
       << LogIO::POST;
	  return False;
  }

    
  PBMath myPB(telName, useSymmetricBeam, qFreq);
  return makePBImage(myPB, pbImage);

}

Bool Imager::makePBImage(const CoordinateSystem& imageCoord, 
			 const String& telescopeName, 
			 const String& diskPBName, 
			 Bool useSymmetricBeam, Double diam){

  LogIO os(LogOrigin("Imager", "makePBImage()", WHERE));
  Int spectralIndex=imageCoord.findCoordinate(Coordinate::SPECTRAL);
  SpectralCoordinate
    spectralCoord=imageCoord.spectralCoordinate(spectralIndex);
  Vector<String> units(1); units = "Hz";
  spectralCoord.setWorldAxisUnits(units);	
  Vector<Double> spectralWorld(1);
  Vector<Double> spectralPixel(1);
  spectralPixel(0) = 0;
  spectralCoord.toWorld(spectralWorld, spectralPixel);  
  Double freq  = spectralWorld(0);
  Quantity qFreq( freq, "Hz" );
  String telName=telescopeName;
  PBMath::CommonPB whichPB;
  PBMath::enumerateCommonPB(telName, whichPB);  
  PBMath myPB;
  if(whichPB!=PBMath::UNKNOWN && whichPB!=PBMath::NONE){
    
    myPB=PBMath(telName, useSymmetricBeam, qFreq);
  }
  else if(diam > 0.0){
    myPB=PBMath(diam,useSymmetricBeam, qFreq);
  }
  else{
    os << LogIO::WARN << "Telescope " << telName << " is not known\n "
       << "Not making the PB  image" 
       << LogIO::POST;
    return False; 
  }
  return makePBImage(imageCoord, myPB, diskPBName);

}

Bool Imager::makePBImage(const CoordinateSystem& imageCoord, 
			 const Table& vpTable, const String& diskPBName){
  ScalarColumn<TableRecord> recCol(vpTable, (String)"pbdescription");
  PBMath myPB(recCol(0));
  return makePBImage(imageCoord, myPB, diskPBName);

}


Bool Imager::makePBImage(const Table& vpTable, ImageInterface<Float>& pbImage){
  ScalarColumn<TableRecord> recCol(vpTable, (String)"pbdescription");
  PBMath myPB(recCol(0));
  return makePBImage(myPB, pbImage);

}

Bool Imager::makePBImage(const CoordinateSystem& /*imageCoord*/, PBMath& pbMath,
			 const String& diskPBName){

  make(diskPBName);
  PagedImage<Float> pbImage(diskPBName);
  return makePBImage(pbMath, pbImage);
}


Bool Imager::makePBImage(PBMath& pbMath, ImageInterface<Float>& pbImage){

  CoordinateSystem imageCoord=pbImage.coordinates();
  pbImage.set(0.0);
  MDirection wcenter;  
  Int directionIndex=imageCoord.findCoordinate(Coordinate::DIRECTION);
  DirectionCoordinate
    directionCoord=imageCoord.directionCoordinate(directionIndex);

  IPosition imShape=pbImage.shape();
  //Vector<Double> pcenter(2);
  // pcenter(0) = imShape(0)/2;
  // pcenter(1) = imShape(1)/2;    
  //directionCoord.toWorld( wcenter, pcenter );
  VisBuffer vb(*rvi_p);
  Int fieldCounter=0;
  Vector<Int> fieldsDone;
  
  for(rvi_p->originChunks(); rvi_p->moreChunks(); rvi_p->nextChunk()){
    Bool fieldDone=False;
    for (uInt k=0;  k < fieldsDone.nelements(); ++k)
      fieldDone=fieldDone || (vb.fieldId()==fieldsDone(k));
    if(!fieldDone){
      ++fieldCounter;
      fieldsDone.resize(fieldCounter, True);
      fieldsDone(fieldCounter-1)=vb.fieldId();
      wcenter=vb.msColumns().field().phaseDirMeas(vb.fieldId());
      TempImage<Float> pbTemp(imShape, imageCoord);
      TempImage<Complex> ctemp(imShape, imageCoord);
      ctemp.set(1.0);
      pbMath.applyPB(ctemp, ctemp, wcenter, Quantity(0.0, "deg"), BeamSquint::NONE);
      StokesImageUtil::To(pbTemp, ctemp);
      pbImage.copyData(  (LatticeExpr<Float>)(pbImage+pbTemp) );
    }
  }
  LatticeExprNode elmax= max( pbImage );
  Float fmax = abs(elmax.getFloat());
  //If there are multiple overlap of beam such that the peak is larger than 1 then normalize
  //otherwise leave as is
  if(fmax>1.0)
    pbImage.copyData((LatticeExpr<Float>)(pbImage/fmax));

  Float cutoffval=minPB_p;
  LatticeExpr<Bool> lemask(iif(pbImage < cutoffval, 
			       False, True));
  ImageRegion outreg=pbImage.makeMask("mask0", False, True);
  LCRegion& outmask=outreg.asMask();
  outmask.copyData(lemask);
  pbImage.defineRegion("mask0", outreg,RegionHandler::Masks, True); 
  pbImage.setDefaultMask("mask0");



  return True;
}
void Imager::setObsInfo(ObsInfo& obsinfo){

  latestObsInfo_p=obsinfo;
}

ObsInfo& Imager::latestObsInfo(){
  return latestObsInfo_p;
}

Bool Imager::makeEmptyImage(CoordinateSystem& coords, String& name, Int fieldID){

  Int tilex=32;
  if(imageTileVol_p >0){
    tilex=static_cast<Int>(ceil(sqrt(imageTileVol_p/min(4, npol_p)/min(32, imageNchan_p))));
    if(tilex >0){
      if(tilex > min(nx_p, ny_p))
	tilex=min(nx_p, ny_p);
      else
	tilex=nx_p/Int(nx_p/tilex);
    }    
    //Not too small in x-y tile
    if(tilex < 10)
      tilex=10;
  }
  IPosition tileShape(4, min(tilex, nx_p), min(tilex, ny_p),
		     min(4, npol_p), min(32, imageNchan_p));
  IPosition imageShape(4, nx_p, ny_p, npol_p, imageNchan_p);
  PagedImage<Float> modelImage(TiledShape(imageShape, tileShape), coords, name);
  modelImage.set(0.0);
  modelImage.table().markForDelete();
    
  // Fill in miscellaneous information needed by FITS
  ROMSColumns msc(*ms_p);
  Record info;
  String object=msc.field().name()(fieldID)
;  //defining object name
  String objectName=msc.field().name()(fieldID);
  ImageInfo iinfo=modelImage.imageInfo();
  iinfo.setObjectName(objectName);
  modelImage.setImageInfo(iinfo);
  String telescop=msc.observation().telescopeName()(0);
  info.define("OBJECT", object);
  info.define("TELESCOP", telescop);
  info.define("INSTRUME", telescop);
  info.define("distance", 0.0);
  modelImage.setMiscInfo(info);
  modelImage.table().tableInfo().setSubType("GENERIC");
  modelImage.setUnits(Unit("Jy/beam"));
  modelImage.table().unmarkForDelete();
  modelImage.table().relinquishAutoLocks(True);
  modelImage.table().unlock();

  return True;
  
}

String Imager::frmtTime(const Double time) {
  MVTime mvtime(Quantity(time, "s"));
  Time t=mvtime.getTime();
  ostringstream os;
  os << t;
  return os.str();
}

Bool Imager::getRestFreq(Vector<Double>& restFreq, const Int& spw){
  // MS Doppler tracking utility
  MSDopplerUtil msdoppler(*ms_p);
  restFreq.resize();
  if(restFreq_p.getValue() > 0){// User defined restfrequency
    restFreq.resize(1);
    restFreq[0]=restFreq_p.getValue("Hz");
  }
  else{
    // Look up first rest frequency found (for now)

    Int fieldid = (datafieldids_p.nelements()>0 ? datafieldids_p(0) : 
		   fieldid_p);
    try{
      msdoppler.dopplerInfo(restFreq ,spw,fieldid);
    }
    catch(...){
      restFreq.resize();
    }
  }
  if(restFreq.nelements() >0) 
    return True;
  return False;
}


void Imager::setSkyEquation(){
  /*  if(sm_p->nmodels() >0){
    Long npix=0;
    for (Int model=0; model < sm_p->nmodels(); ++model){
      Long pixmod=sm_p->image(model).product();
      npix=max(pixmod, npix);
    }
    Long pixInMem=(HostInfo::memoryTotal()/8)*1024;
    if(npix > (pixInMem/2)){
      se_p = new CubeSkyEquation(*sm_p, *vs_p, *ft_p, *cft_p, !useModelCol_p); 
      return;
    }
    //else lets make the base SkyEquation for now
  }
  
  se_p = new SkyEquation(*sm_p, *vs_p, *ft_p, *cft_p, !useModelCol_p);

  */
//  if (ft_p->name() == "PBWProjectFT")
//    {
//      logSink_p.clearLocally();
//      LogIO os(LogOrigin("imager", "setSkyEquation()"), logSink_p);
//      os << "Creating SkyEquation for PBWProjectFT" << LogIO::POST;
//     se_p = new SkyEquation(*sm_p, *vs_p, *ft_p, *cft_p, !useModelCol_p);
//    }
//  else
  se_p = new CubeSkyEquation(*sm_p, *rvi_p, *ft_p, *cft_p, !useModelCol_p);
  return;
}

void Imager::savePSF(const Vector<String>& psf){

  if( (psf.nelements() > 0) && (Int(psf.nelements()) <= sm_p->numberOfModels())){

    for (Int thismodel=0;thismodel<Int(psf.nelements());++thismodel) {
      if(removeTable(psf(thismodel))) {
	Int whichmodel=thismodel;
	if(facets_p >1 && thismodel > 0)
	  whichmodel=facets_p*facets_p-1+thismodel;
	IPosition shape=images_p[thismodel]->shape();
	PagedImage<Float> psfimage(shape,
				   images_p[thismodel]->coordinates(),
				   psf(thismodel));
	if(freqFrameValid_p){
	  CoordinateSystem cs=psfimage.coordinates();
	  String errorMsg;
	  if (CoordinateUtil::setSpectralConversion (errorMsg, cs,MFrequency::showType(freqFrame_p))) {
	    psfimage.setCoordinateInfo(cs);
	  }
        }
	psfimage.set(0.0);
	if((shape[0]*shape[1]) > ((sm_p->PSF(whichmodel)).shape()[0]*(sm_p->PSF(whichmodel)).shape()[1])){
	  IPosition blc(4, 0, 0, 0, 0);
	  IPosition trc=shape-1;
	  blc[0]=(shape[0]-(sm_p->PSF(whichmodel)).shape()[0])/2;
	  trc[0]=(sm_p->PSF(whichmodel)).shape()[0]+blc[0]-1;
	  blc[1]=(shape[1]-(sm_p->PSF(whichmodel)).shape()[1])/2;
	  trc[1]=(sm_p->PSF(whichmodel)).shape()[1]+blc[1]-1;
	  Slicer sl(blc, trc, Slicer::endIsLast);
	  SubImage<Float> sub(psfimage, sl, True);
	  sub.copyData(sm_p->PSF(whichmodel));	  
	}
	else{
	  psfimage.copyData(sm_p->PSF(whichmodel));
	}
      }
    }



  }  

}

void Imager::setMosaicFTMachine(Bool useDoublePrec){
  LogIO os(LogOrigin("Imager", "setmosaicftmachine", WHERE));
  ROMSColumns msc(*ms_p);
  String telescop=msc.observation().telescopeName()(0);
  PBMath::CommonPB kpb;
  PBMath::enumerateCommonPB(telescop, kpb);
  if(!((kpb == PBMath::UNKNOWN) || 
       (kpb==PBMath::OVRO) || (kpb==PBMath::ALMA) || (kpb==PBMath::ACA))){
    
    if(!gvp_p) {
      if (doDefaultVP_p) {
	os << LogIO::NORMAL // Loglevel INFO
           << "Using defaults for primary beams used in gridding" << LogIO::POST;
	gvp_p=new VPSkyJones(*ms_p, True, parAngleInc_p, squintType_p,
                             skyPosThreshold_p);
      } else {
	os << LogIO::NORMAL // Loglevel INFO
           << "Using VP as defined in " << vpTableStr_p <<  LogIO::POST;
	Table vpTable( vpTableStr_p ); 
	gvp_p=new VPSkyJones(*ms_p, vpTable, parAngleInc_p, squintType_p,
                             skyPosThreshold_p);
      }
    } 
    gvp_p->setThreshold(minPB_p);
  }
  
  ft_p = new MosaicFT(gvp_p, mLocation_p, stokes_p, cache_p/2, tile_p, True, 
		      useDoublePrec);

  if((kpb == PBMath::UNKNOWN) || (kpb==PBMath::OVRO) || (kpb==PBMath::ACA)
     || (kpb==PBMath::ALMA)){
    os << LogIO::NORMAL // Loglevel INFO
       << "Using antenna diameters for determining beams for gridding"
       << LogIO::POST;
    CountedPtr<SimplePBConvFunc> mospb=new HetArrayConvFunc();
    static_cast<MosaicFT &>(*ft_p).setConvFunc(mospb);
  }
  
}
ATerm* Imager::createTelescopeATerm(MeasurementSet& ms)
{
  LogIO log_l(LogOrigin("Imager", "createTelescopeATerm"));
  ROMSObservationColumns msoc(ms.observation());
  String ObsName=msoc.telescopeName()(0);
  if ((ObsName == "EVLA") || (ObsName == "VLA"))
    return new EVLAAperture();
  else
    log_l << "Telescope name ('"+
      ObsName+"') in the MS not recognized to create the telescope specific ATerm"
	  << LogIO::EXCEPTION;

  return NULL;
}

//use SubMS::calcChanFreqs to calculate spectral gridding 
//call from imagecoodinates2
Bool Imager::calcImFreqs(Vector<Double>& imgridfreqs,
                          Vector<Double>& imfreqres,
			  const MFrequency::Types& oldRefFrame,
			  const MEpoch& obsEpoch,
			  const MPosition& obsPosition,
			  const Double& restFreq 
                          )
{ 

  logSink_p.clearLocally();
  LogIO os(LogOrigin("imager", "setGridFreqs()"), logSink_p);

  ROMSColumns msc(*ms_p);
  Vector<Double> oldChanFreqs;
  Vector<Double> oldFreqResolution;
  String veltype;
  String mode;
  String restfreq;
  String start;
  String width;
  String outframe;
  Bool reversevec(False);
  Bool descendfreq(False);
  
  if (imageMode_p.contains("RADIO")) {
    veltype="radio";
    mode="velocity";
    start=dQuantitytoString(mImageStart_p.get("m/s"));
    width=dQuantitytoString(mImageStep_p.get("m/s"));
    if (!width.contains(casacore::Regex("^-"))) {
      //positive vel. width (descending frequencies) 
      //reversevec=True;
      descendfreq=True;
    }
  }
  else if (imageMode_p.contains("OPTICAL")) {
    veltype="optical";
    mode="velocity";
    start=dQuantitytoString(mImageStart_p.get("m/s"));
    width=dQuantitytoString(mImageStep_p.get("m/s"));
    //cerr<<"optical vel width USED="<<width<<endl;
    if (!width.contains(casacore::Regex("^-"))) {
      //positive vel. width (descending frequencies)
      //reversevec=True;
      descendfreq=True;
    }
  }
  else if (imageMode_p.contains("FREQ")) {
    veltype="radio";
    mode="frequency";
    start=dQuantitytoString(mfImageStart_p.get("Hz"));
    width=dQuantitytoString(mfImageStep_p.get("Hz"));
    if (width.contains(casacore::Regex("^-"))) {
      //reversevec=True;
      descendfreq=True;
    }
  }
  else if (imageMode_p.contains("CHANNEL")) {
    veltype="radio";
    mode="channel";
    start=String::toString(imageStart_p);
    width=String::toString(imageStep_p);
    if (width.contains(casacore::Regex("^-"))) {
      // means here going to lower chan index
      descendfreq=True;
    }
  }
  else if (imageMode_p.contains("MFS")) {
    veltype="radio";
    mode="mfs";
    start=String::toString(imageStart_p);
  }
 
  restfreq = dQuantitytoString(Quantity(restFreq,"Hz"));
  //MFrequency::getType(freqFrame_p, outframe);
  if (freqFrame_p!=oldRefFrame) {
    outframe=MFrequency::showType(freqFrame_p);
  }
  else {
    outframe="";
  }

 
  try {
    if(spectralwindowids_p.nelements()==1){
      if(spectralwindowids_p[0]<0){
        spectralwindowids_p.resize();
        if(dataspectralwindowids_p.nelements()==0){
          Int nspwinms=ms_p->spectralWindow().nrow();
          dataspectralwindowids_p.resize(nspwinms);
          indgen(dataspectralwindowids_p);
        }
        spectralwindowids_p=dataspectralwindowids_p;
      }
    }

    if(spectralwindowids_p.nelements()==1) {
      oldChanFreqs=msc.spectralWindow().chanFreq()(spectralwindowids_p[0]);  
      oldFreqResolution=msc.spectralWindow().chanWidth()(spectralwindowids_p[0]);
    }
    else {
      SubMS thems(*ms_p);
      if(!thems.combineSpws(spectralwindowids_p,True,oldChanFreqs,oldFreqResolution)){
        os << LogIO::SEVERE << "Error combining SpWs" << LogIO::POST;
      }
    }
    Bool isDescendingData=False;
    // Some descending order data has positive channel widths...so check chan freqs
    // first...
    //if (oldFreqResolution(0) < 0) isDescendingData=True;
    if (oldChanFreqs.nelements()>1) {
      if ((oldChanFreqs[1] - oldChanFreqs[0])<0) isDescendingData=True;
    }
    else if (oldFreqResolution(0) < 0) {
      isDescendingData=True;
    }
    
    // need theOldRefFrame,theObsTime,mObsPos,mode,nchan,start,width,restfreq,
    // outframe,veltype
    //

    Bool rst=SubMS::calcChanFreqs(os,
			   imgridfreqs, 
			   imfreqres,
			   oldChanFreqs, 
			   oldFreqResolution,
			   phaseCenter_p,
			   oldRefFrame,
			   obsEpoch,
			   obsPosition,
			   mode, 
			   imageNchan_p, 
			   start, 
			   width,
			   restfreq, 
			   outframe,
			   veltype
			   );
    
    if (!rst) {
      os << LogIO::SEVERE << "calcChanFreqs failed, check input start and width parameters"
         << LogIO::EXCEPTION;
      return False;
    }

    //cout<<"=================calcChanFreqs arguments ======================"<<endl;
    //cout.precision(10);
    //cout<<"imgridfreqs(0)="<<imgridfreqs(0)<<endl;
    //cout<<"imgridfreqs("<<imgridfreqs.nelements()-1<<")="<<imgridfreqs(imgridfreqs.nelements()-1)<<endl;
    //cout<<"oldChanFreqs(0)="<<oldChanFreqs(0)<<endl;
    //cout<<"oldChanFreqs("<<oldChanFreqs.nelements()-1<<")="<<oldChanFreqs(oldChanFreqs.nelements()-1)<<endl;
    //cout<<"phaseCenter_p="<<phaseCenter_p<<endl;
    //cout<<"oldRefFrame="<<oldRefFrame<<endl;
    //cout<<"outframe="<<outframe<<endl;
    //cout<<"obsEpoch="<<obsEpoch<<endl;
    //cout<<"obsPosition="<<obsPosition<<endl;
    //cout<<"start="<<start<<" width="<<width<<endl;
    //cout<<"restfreq="<<restfreq<<endl;
    //cout<<"veltype="<<veltype<<endl;
    //cout<<"=================calcChanFreqs arguments end==================="<<endl;
    Bool isDescendingNewData=False;
    if (imgridfreqs(0)-imgridfreqs(1)>0) isDescendingNewData=True;
    //reverse frequency vector? 
    //evaluate reversing condition differ for chan mode from other modes
    if(mode.contains("channel")) {
      if ((descendfreq && !isDescendingNewData && !isDescendingData) |
          (descendfreq && isDescendingNewData && isDescendingData) |
          (!descendfreq && !isDescendingNewData && isDescendingData)) {
        reversevec = True;
      }
    }
    else {
      if ((descendfreq && !isDescendingNewData) | 
          (!descendfreq && isDescendingNewData)) { 
        reversevec = True;
      }
    }
    if (reversevec) {
    //if(reversevec && isAscendingData ) {
      //Int ndata=imgridfreqs.nelements();
      //tempimgridfreqs.resize(ndata);
      /**
      for (Int i=0;i<ndata;i++) {
        tempimgridfreqs[i] = imgridfreqs[ndata - 1 - i];
      }    
      for (Int i=0;i<ndata;i++) {
        std::swap(imgridfreqs[ndata-1-i],tempimgridfreqs[i]);
      }  
      **/
      std::vector<double>  stlimgridfreqs;
      imgridfreqs.tovector(stlimgridfreqs);
      std::reverse(stlimgridfreqs.begin(),stlimgridfreqs.end());  
      imgridfreqs=Vector<double>(stlimgridfreqs);
    }
    //cerr<<"Final imgridfreqs(0)="<<imgridfreqs(0)<<endl;
    
  } catch (AipsError x) {
    this->unlock();
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << LogIO::EXCEPTION;
    return False;
  } 
  return True;
}//end of calcImFreqs

// convert a double precision quanity to a String
String Imager::dQuantitytoString(const Quantity& dq) {
  std::ostringstream ss;
  ss.precision(std::numeric_limits<double>::digits10);
  ss << dq;
  return ss.str();
} 

} //# NAMESPACE CASACORE - END





