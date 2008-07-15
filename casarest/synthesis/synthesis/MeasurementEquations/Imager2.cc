//# Imager2.cc: Imager.cc is split in 3 parts
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003
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
//# $Id: Imager.cc,v 19.43 2006/10/05 20:37:13 rurvashi Exp $



#include <casa/Exceptions/Error.h>
#include <casa/iostream.h>
#include <synthesis/MeasurementEquations/Imager.h>

#include <ms/MeasurementSets/MSHistoryHandler.h>

#include <casa/Arrays/Matrix.h>
#include <casa/Arrays/ArrayMath.h>

#include <casa/Logging.h>
#include <casa/Logging/LogIO.h>
#include <casa/Logging/LogMessage.h>

#include <casa/OS/File.h>
#include <casa/OS/HostInfo.h>
#include <casa/Containers/Record.h>


#include <tables/Tables/Table.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/TableParse.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/TableLock.h>
#include <tables/Tables/ExprNode.h>

#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/Utilities/Fallible.h>
#include <casa/Utilities/CompositeNumber.h>

#include <casa/BasicSL/Constants.h>

#include <casa/Logging/LogSink.h>
#include <casa/Logging/LogMessage.h>

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/Slice.h>
#include <synthesis/MeasurementEquations/ClarkCleanProgress.h>
#include <lattices/Lattices/LatticeCleanProgress.h>
#include <msvis/MSVis/VisSet.h>
#include <msvis/MSVis/VisSetUtil.h>
#include <synthesis/MeasurementComponents/TimeVarVisJones.h>

#include <measures/Measures/Stokes.h>
#include <casa/Quanta/UnitMap.h>
#include <casa/Quanta/UnitVal.h>
#include <casa/Quanta/MVAngle.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MPosition.h>
#include <casa/Quanta/MVEpoch.h>
#include <casa/Quanta/MVTime.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasTable.h>

#include <ms/MeasurementSets/MeasurementSet.h>
#include <ms/MeasurementSets/MSColumns.h>

#include <ms/MeasurementSets/MSDopplerUtil.h>
#include <ms/MeasurementSets/MSSourceIndex.h>
#include <ms/MeasurementSets/MSSummary.h>
#include <synthesis/MeasurementEquations/MosaicSkyEquation.h>
#include <synthesis/MeasurementEquations/WFSkyEquation.h>
#include <synthesis/MeasurementEquations/WBSkyEquation.h>
#include <synthesis/MeasurementEquations/VisEquation.h>
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
#include <synthesis/MeasurementComponents/GridBoth.h>
#include <synthesis/MeasurementComponents/WFGridFT.h>
#include <synthesis/MeasurementComponents/MosaicFT.h>
#include <synthesis/MeasurementComponents/WProjectFT.h>
#include <synthesis/MeasurementComponents/PBWProjectFT.h>
#include <synthesis/MeasurementComponents/WideBandFT.h>
#include <synthesis/MeasurementComponents/SimpleComponentFTMachine.h>
#include <synthesis/MeasurementComponents/SimpCompGridMachine.h>
#include <synthesis/MeasurementComponents/VPSkyJones.h>
#include <synthesis/MeasurementComponents/SynthesisError.h>

#include <synthesis/DataSampling/SynDataSampling.h>
#include <synthesis/DataSampling/SDDataSampling.h>
#include <synthesis/DataSampling/ImageDataSampling.h>
#include <synthesis/DataSampling/PixonProcessor.h>

#include <synthesis/MeasurementEquations/StokesImageUtil.h>
#include <lattices/Lattices/TiledLineStepper.h> 
#include <lattices/Lattices/LatticeIterator.h> 
#include <lattices/Lattices/LatticeExpr.h> 
#include <lattices/Lattices/LCBox.h> 
#include <lattices/Lattices/LatticeFFT.h>
#include <images/Images/ImageRegrid.h>
#include <synthesis/MeasurementComponents/PBMath.h>


#include <images/Images/PagedImage.h>
#include <images/Images/ImageInfo.h>
#include <images/Images/SubImage.h>
#include <images/Images/ImageUtilities.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <coordinates/Coordinates/DirectionCoordinate.h>
#include <coordinates/Coordinates/SpectralCoordinate.h>
#include <coordinates/Coordinates/StokesCoordinate.h>
#include <coordinates/Coordinates/Projection.h>
#include <coordinates/Coordinates/ObsInfo.h>

#include <components/ComponentModels/ComponentList.h>
#include <components/ComponentModels/ConstantSpectrum.h>
#include <components/ComponentModels/Flux.h>
#include <components/ComponentModels/PointShape.h>
#include <components/ComponentModels/FluxStandard.h>


#include <casa/OS/HostInfo.h>
#include <casa/System/PGPlotter.h>

#include <components/ComponentModels/ComponentList.h>

#include <measures/Measures/UVWMachine.h>

#include <casa/sstream.h>

#ifdef PABLO_IO
#include "PabloTrace.h"
#endif

namespace casa { //# NAMESPACE CASA - BEGIN

// Make standard choices for coordinates
Bool Imager::imagecoordinates(CoordinateSystem& coordInfo) 
{  
  if(!valid()) return False;
  if(!assertDefinedImageParameters()) return False;
  LogIO os(LogOrigin("Imager", "imagecoordinates()", WHERE));
  
  Vector<Double> deltas(2);
  deltas(0)=-mcellx_p.get("rad").getValue();
  deltas(1)=mcelly_p.get("rad").getValue();
  
  MSColumns msc(*ms_p);
  MFrequency::Types obsFreqRef=MFrequency::DEFAULT;
  ROScalarColumn<Int> measFreqRef(ms_p->spectralWindow(),
				  MSSpectralWindow::columnName(MSSpectralWindow::MEAS_FREQ_REF));
  //using the first frame of reference; TO DO should do the right thing 
  //for different frames selected. 
  if(measFreqRef(spectralwindowids_p(0)) >=0) 
     obsFreqRef=(MFrequency::Types)measFreqRef(spectralwindowids_p(0));
			    

  // MS Doppler tracking utility
  MSDopplerUtil msdoppler(*ms_p);

  MVDirection mvPhaseCenter(phaseCenter_p.getAngle());
  // Normalize correctly
  MVAngle ra=mvPhaseCenter.get()(0);
  ra(0.0);
  MVAngle dec=mvPhaseCenter.get()(1);
  Vector<Double> refCoord(2);
  refCoord(0)=ra.get().getValue();    
  refCoord(1)=dec;    
  
  Vector<Double> refPixel(2); 
  refPixel(0)=Double(nx_p/2);
  refPixel(1)=Double(ny_p/2);
  
  //defining observatory...needed for position on earth
  String telescop=msc.observation().telescopeName()(0);

  // defining epoch as begining time from timerange in OBSERVATION subtable
  // Using first observation for now
  MEpoch obsEpoch=msc.observation().timeRangeMeas()(0)(IPosition(1,0));

  //Now finding the position of the telescope on Earth...needed for proper
  //frequency conversions

  MPosition obsPosition;
  if(! (MeasTable::Observatory(obsPosition, telescop))){
    os << LogIO::WARN << "Did not get the position of " << telescop 
       << " from data repository" << LogIO::POST ;
    os << LogIO::WARN 
       << "Please do inform aips++  to put in the repository "
       << LogIO::POST;
    os << LogIO::WARN << "Frequency conversion will not work " << LogIO::POST;
    freqFrameValid_p=False;
  }
  else{
    mLocation_p=obsPosition;
    freqFrameValid_p=True;
  }
  // Now find the projection to use: could probably also use
  // max(abs(w))=0.0 as a criterion
  Projection projection(Projection::SIN);
  if(telescop=="ATCASCP") {
    os << LogIO::NORMAL << "Using SIN image projection adjusted for SCP" 
       << LogIO::POST;
    Vector<Double> projectionParameters(2);
    projectionParameters(0)=0.0;
    if(sin(dec)!=0.0) {
      projectionParameters(1)=cos(dec)/sin(dec);
      projection=Projection(Projection::SIN, projectionParameters);
    }
    else {
      os << LogIO::WARN << "Singular projection for ATCA: using plain SIN" << LogIO::POST;
      projection=Projection(Projection::SIN);
    }
  }
  else if(telescop=="WSRT") {
    os << LogIO::NORMAL << "Using SIN image projection adjusted for NCP" 
       << LogIO::POST;
    Vector<Double> projectionParameters(2);
    projectionParameters(0)=0.0;
    if(sin(dec)!=0.0) {
      projectionParameters(1)=cos(dec)/sin(dec);
      projection=Projection(Projection::SIN, projectionParameters);
    }
    else {
      os << LogIO::WARN << "Singular projection for WSRT: using plain SIN" 
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
  if (imageMode_p=="mfs") {
    Double fmin=C::dbl_max;
    Double fmax=-(C::dbl_max);
    Double fmean=0.0;
    for (Int i=0;i<nspw;++i) {
      Int spw=spectralwindowids_p(i);
      Vector<Double> chanFreq=msc.spectralWindow().chanFreq()(spw); 
      Vector<Double> freqResolution=msc.spectralWindow().resolution()(spw); 
      
      if(dataMode_p=="none"){
      

	if(i==0) {
	  fmin=min(chanFreq-abs(freqResolution));
	  fmax=max(chanFreq+abs(freqResolution));
	}
	else {
	  fmin=min(fmin,min(chanFreq-abs(freqResolution)));
	  fmax=max(fmax,max(chanFreq+abs(freqResolution)));
	}
      }
      else if(dataMode_p=="channel"){
	Int lastchan=dataStart_p[i]+ dataNchan_p[i]*dataStep_p[i];
        for(Int k=dataStart_p[i] ; k < lastchan ;  k+=dataStep_p[i]){
	  fmin=min(fmin,chanFreq[k]-abs(freqResolution[k]*dataStep_p[i]));
	  fmax=max(fmax,chanFreq[k]+abs(freqResolution[k]*dataStep_p[i]));
        }
      }
      else{
	this->unlock();
	os << LogIO::SEVERE 
	   << "setdata has to be in 'channel' or 'none' mode for 'mfs' imaging to work"
	   << LogIO::POST;
      return False;
      }
 
    }

    fmean=(fmax+fmin)/2.0;
    // Look up first rest frequency found (for now)
    Vector<Double> restFreqArray;
    Double restFreq=fmean;
    Int fieldid = (datafieldids_p.nelements()>0 ? datafieldids_p(0) : fieldid_p);
    if (msdoppler.dopplerInfo(restFreqArray,spectralwindowids_p(0),fieldid)) {
      restFreq = restFreqArray(0);    
    } 
    imageNchan_p=1;
    Double finc=(fmax-fmin); 
    mySpectral = new SpectralCoordinate(obsFreqRef,  fmean, finc,
      					refChan, restFreq);
    os <<  "Frequency = "
       << MFrequency(Quantity(fmean, "Hz")).get("GHz").getValue()
       << " GHz, synthesized continuum bandwidth = "
       << MFrequency(Quantity(finc, "Hz")).get("GHz").getValue()
       << " GHz" << LogIO::POST;
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
	msc.spectralWindow().resolution()(spw); 
      
      // Look up first rest frequency found (for now)
     
      Vector<Double> restFreqArray;
      Int fieldid = (datafieldids_p.nelements()>0 ? datafieldids_p(0) : 
		     fieldid_p);
      if (msdoppler.dopplerInfo(restFreqArray,spw,fieldid)) {
	if(spwIndex==0){
	  restFreq = restFreqArray(0);
	}
	else{
	  if(restFreq != restFreqArray(0)){
	    os << LogIO::WARN << "Rest frequencies are different for  spectralwindows selected " 
	       << LogIO::POST;
	    os << LogIO::WARN 
	       <<"Will be using the restFreq defined in spectralwindow "
	       << spectralwindowids_p(0)+1 << LogIO::POST;
	  }
	  
	}	
      }
    }

    if(imageMode_p=="channel") {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::POST;
	return False;
      }
      Vector<Double> freqs;
      Int nsubchans=
	(chanFreq.shape()(0) - Int(imageStart_p)+1)/Int(imageStep_p);
      if(imageNchan_p>nsubchans) imageNchan_p=nsubchans;
      os << "Image spectral coordinate: "<< imageNchan_p
	   << " channels, starting at visibility channel "
	 << imageStart_p+1 << " stepped by "
	 << imageStep_p << endl;
      freqs.resize(imageNchan_p);
      for (Int chan=0;chan<imageNchan_p;chan++) {
	freqs(chan)=chanFreq(Int(imageStart_p)+Int(Float(chan+0.5)*Float(imageStep_p)-0.5));
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
      mySpectral = new SpectralCoordinate(obsFreqRef, freqs(0), finc,
					  refChan, restFreq);
      os <<  "Frequency = "
	 << MFrequency(Quantity(freqs(0), "Hz")).get("GHz").getValue()
	 << ", channel increment = "
	 << MFrequency(Quantity(finc, "Hz")).get("GHz").getValue() 
	 << "GHz" << endl;
      os << LogIO::NORMAL << "Rest frequency is " 
	 << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	 << "GHz" << LogIO::POST;
      
    }
    // Spectral channels resampled at equal increments in optical velocity
    // Here we compute just the first two channels and use increments for
    // the others
    else if (imageMode_p=="velocity") {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::POST;
	return False;
      }
      {
	ostringstream oos;
	oos << "Image spectral coordinate:"<< imageNchan_p 
	    << " channels, starting at radio velocity " << mImageStart_p
	    << "  stepped by " << mImageStep_p << endl;
	os << String(oos);
      }
      Vector<Double> freqs(2);
      freqs=0.0;
      if(Double(mImageStep_p.getValue())!=0.0) {
	MRadialVelocity mRadVel=mImageStart_p;
	for (Int chan=0;chan<2;chan++) {
	  MDoppler mdoppler=mRadVel.toDoppler();
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

      // when setting in velocity its in LSRK
      mySpectral = new SpectralCoordinate(MFrequency::LSRK, freqs(0),
					  freqs(1)-freqs(0), refChan,
					  restFreq);
      {
	ostringstream oos;
	oos << "Reference Frequency = "
	    << MFrequency(Quantity(freqs(0), "Hz")).get("GHz")
	    << ", spectral increment = "
	    << MFrequency(Quantity(freqs(1)-freqs(0), "Hz")).get("GHz") 
	    << endl;
	oos << "Rest frequency is " 
	    << MFrequency(Quantity(restFreq, "Hz")).get("GHz").getValue()
	    << " GHz" << endl;
	os << String(oos) << LogIO::POST;
      }
      
    }
    // Since optical velocity is non-linear in frequency, we have to
    // pass in all the frequencies. For radio velocity we can use 
    // a linear axis.
    else if (imageMode_p=="opticalvelocity") {
      if(imageNchan_p==0) {
	this->unlock();
	os << LogIO::SEVERE << "Must specify number of channels" 
	   << LogIO::POST;
	return False;
      }
      {
	ostringstream oos;
	oos << "Image spectral coordinate: "<< imageNchan_p 
	    << " channels, starting at optical velocity " << mImageStart_p
	    << "  stepped by " << mImageStep_p << endl;
	os << String(oos);
      }
      Vector<Double> freqs(imageNchan_p);
      freqs=0.0;
      if(Double(mImageStep_p.getValue())!=0.0) {
	MRadialVelocity mRadVel=mImageStart_p;
	for (Int chan=0;chan<imageNchan_p;++chan) {
	  MDoppler mdoppler=mRadVel.toDoppler();
	  freqs(chan)=
	    MFrequency::fromDoppler(mdoppler, restFreq).getValue().getValue();
	  mRadVel.set(mRadVel.getValue()+mImageStep_p.getValue());
	}
      }
      else {
	for (Int chan=0;chan<imageNchan_p;++chan) {
	    freqs(chan)=chanFreq(chan);
	}
      }
      // Use this next line when non-linear is working
      // when selecting in velocity its LSRK
      mySpectral = new SpectralCoordinate(MFrequency::LSRK, freqs,
					  restFreq);
      // mySpectral = new SpectralCoordinate(MFrequency::DEFAULT, freqs(0),
      //				       freqs(1)-freqs(0), refChan,
      //				        restFreq);
      {
	ostringstream oos;
	oos << "Reference Frequency = "
	    << MFrequency(Quantity(freqs(0), "Hz")).get("GHz")
	    << " Ghz" << endl;
	os << String(oos) << LogIO::POST;
      }
    }
    else {
      this->unlock();
      os << LogIO::SEVERE << "Unknown mode " << imageMode_p
	 << LogIO::POST;
      return False;
    }
        
    
  }
 
 
    //In FTMachine lsrk is used for channel matching with data channel 
    //hence we make sure that
    // we convert to lsrk when dealing with the channels
    
  if(freqFrameValid_p){
      mySpectral->setReferenceConversion(MFrequency::LSRK, obsEpoch, 
					 obsPosition,
					 phaseCenter_p);
  }


 
  // Polarization
  Vector<String> polType=msc.feed().polarizationType()(0);
  if (polType(0)!="X" && polType(0)!="Y" &&
      polType(0)!="R" && polType(0)!="L") {
    os << "Warning: Unknown stokes types in feed table: ["
       << polType(0) << ", " << polType(1) << "]" << endl
       << "Results open to question!" << LogIO::POST;
  }
  
  if (polType(0)=="X" || polType(0)=="Y") {
    polRep_p=SkyModel::LINEAR;
    os << "Preferred polarization representation is linear" << LogIO::POST;
  }
  else {
    polRep_p=SkyModel::CIRCULAR;
    os << "Preferred polarization representation is circular" << LogIO::POST;
  }

  Vector<Int> whichStokes(npol_p);
  switch(npol_p) {
  case 1:
    whichStokes.resize(1);
    whichStokes(0)=Stokes::I;
    os <<  "Image polarization = Stokes I" << LogIO::POST;
    break;
  case 2:
    whichStokes.resize(2);
    whichStokes(0)=Stokes::I;
    if (polRep_p==SkyModel::LINEAR) {
      whichStokes(1)=Stokes::Q;
      os <<  "Image polarization = Stokes I,Q" << LogIO::POST;
    }
    else {
      whichStokes(1)=Stokes::V;
      os <<  "Image polarization = Stokes I,V" << LogIO::POST;
    }
    break;
  case 3:
    whichStokes.resize(3);
    whichStokes(0)=Stokes::I;
    whichStokes(1)=Stokes::Q;
    whichStokes(1)=Stokes::U;
    os <<  "Image polarization = Stokes I,Q,U" << LogIO::POST;
    break;
  case 4:
    whichStokes.resize(4);
    whichStokes(0)=Stokes::I;
    whichStokes(1)=Stokes::Q;
    whichStokes(2)=Stokes::U;
    whichStokes(3)=Stokes::V;
    os <<  "Image polarization = Stokes I,Q,U,V" << LogIO::POST;
    break;
  default:
    this->unlock();
    os << LogIO::SEVERE << "Illegal number of Stokes parameters: " << npol_p
       << LogIO::POST;
    return False;
  };
  
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

IPosition Imager::imageshape() const
{
  return IPosition(4, nx_p, ny_p, npol_p, imageNchan_p);
}

Bool Imager::summary() 
{
  if(!valid()) return False;
  LogOrigin OR("imager", "Imager::summary()", WHERE);
  
  LogIO los(OR);
  
  los << "Logging summary" << LogIO::POST;
  try {
    
    this->lock();
    MSSummary mss(*ms_p);
    mss.list(los, True);
    
    los << endl << state() << LogIO::POST;
    this->unlock();
    return True;
  } catch (AipsError x) {
    los << LogIO::SEVERE << "Caught Exception: " << x.getMesg()
	<< LogIO::POST;
    this->unlock();
    return False;
  } 
  
  return True;
}


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
    
    MSColumns msc(*ms_p);
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
	   << fieldid_p+1 << " :" << endl;
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
      
      if(imageMode_p=="mfs") {
	os << "  Image mode is mfs: Image will be frequency synthesised from spectral windows : ";
	for (uInt i=0;i<spectralwindowids_p.nelements();++i) {
	  os << spectralwindowids_p(i)+1 << " ";
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
	  os << dataspectralwindowids_p(i)+1 << " ";
	}
      }
      if(datafieldids_p.nelements()>0) {
	os << "  Data selected includes fields : ";
	for (uInt i=0;i<datafieldids_p.nelements();++i) {
	  os << datafieldids_p(i)+1 << " ";
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
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << endl;
    this->unlock();
  } 
  return String(os);
}

Bool Imager::setimage(const Int nx, const Int ny,
		      const Quantity& cellx, const Quantity& celly,
		      const String& stokes,
		      Bool doShift,
		      const MDirection& phaseCenter, 
		      const Quantity& shiftx, const Quantity& shifty,
		      const String& mode, const Int nchan,
		      const Int start, const Int step,
		      const MRadialVelocity& mStart, const MRadialVelocity& mStep,
		      const Vector<Int>& spectralwindowids,
		      const Int fieldid,
		      const Int facets,
		      const Quantity& distance,
		      const Float &paStep, const Float &pbLimit)
{



#ifdef PABLO_IO
  traceEvent(1,"Entering Imager::setimage",26);
#endif

  if(!valid())
    {

#ifdef PABLO_IO
      traceEvent(1,"Exiting Imager::setimage",25);
#endif

      return False;
    }

  //Clear the sink 
  logSink_p.clearLocally();
  LogIO os(LogOrigin("imager", "setimage()"), logSink_p);

  os << "nx=" << nx << " ny=" << ny
     << " cellx='" << cellx.getValue() << cellx.getUnit()
     << "' celly='" << celly.getValue() << celly.getUnit()
     << "' stokes=" << stokes << " doShift=" << doShift
     << " shiftx='" << shiftx.getValue() << shiftx.getUnit()
     << "' shifty='" << shifty.getValue() << shifty.getUnit()
     << "' mode=" << mode << " nchan=" << nchan
     << " start=" << start << " step=" << step
     << " spwids=" << spectralwindowids
     << " fieldid=" <<   fieldid << " facets=" << facets
     << " distance='" << distance.getValue() << distance.getUnit() <<"'";
  ostringstream clicom;
  clicom << " phaseCenter='" << phaseCenter;
  clicom << "' mStart='" << mStart << "' mStep='" << mStep << "'";
  os << String(clicom);
  
  try {
    
    this->lock();
    this->writeCommand(os);

    os << "Defining image properties" << LogIO::POST;
  
    /**** this check is not really needed here especially for SD imaging
    if(2*Int(nx/2)!=nx) {
      this->unlock();
      os << LogIO::SEVERE << "nx must be even" << LogIO::POST;
      return False;
    }
    if(2*Int(ny/2)!=ny) {
      this->unlock();
      os << LogIO::SEVERE << "ny must be even" << LogIO::POST;
      return False;
    }

    */
    {
      CompositeNumber cn(nx);
      if (! cn.isComposite(nx)) {
	Int nxc = (Int)cn.nextLargerEven(nx);
	Int nnxc = (Int)cn.nearestEven(nx);
	if (nxc == nnxc) {
	  os << LogIO::WARN << "nx = " << nx << " is not composite; nx = " 
	     << nxc << " will be more efficient" << LogIO::POST;
	} else {
	  os <<  LogIO::WARN << "nx = " << nx << " is not composite; nx = " 
	     << nxc <<  " or " << nnxc << " will be more efficient" << LogIO::POST;
	}
      }
      if (! cn.isComposite(ny)) {
	Int nyc = (Int)cn.nextLargerEven(ny);
	Int nnyc = (Int)cn.nearestEven(ny);
	if (nyc == nnyc) {
	  os <<  LogIO::WARN << "ny = " << ny << " is not composite; ny = " 
	     << nyc << " will be more efficient" << LogIO::POST;
	} else {
	  os <<  LogIO::WARN << "ny = " << ny << " is not composite; ny = " << nyc << 
	      " or " << nnyc << " will be more efficient" << LogIO::POST;
	}
      }
      os << LogIO::WARN 
	 << "You may safely ignore this message for single dish imaging" 
	 << LogIO::POST;

    }

    paStep_p = paStep;
    pbLimit_p = pbLimit;
    nx_p=nx;
    ny_p=ny;
    mcellx_p=cellx;
    mcelly_p=celly;
    distance_p=distance;
    stokes_p=stokes;
    imageMode_p=mode;
    imageNchan_p=nchan;
    imageStart_p=start;
    imageStep_p=step;
    mImageStart_p=mStart;
    mImageStep_p=mStep;
    spectralwindowids_p.resize(spectralwindowids.nelements());
    spectralwindowids_p=spectralwindowids;
    fieldid_p=fieldid;
    facets_p=facets;
    redoSkyModel_p=True;
    destroySkyEquation();    

    // Now make the derived quantities 
    if(stokes_p=="I") {
      npol_p=1;
    }
    else if(stokes_p=="IQ") {
      npol_p=2;
    }
    else if(stokes_p=="IV") {
      npol_p=2;
    }
    else if(stokes_p=="IQU") {
      npol_p=3;
    }
    else if(stokes_p=="IQUV") {
      npol_p=4;
    }
    else {
      this->unlock();
      os << LogIO::SEVERE << "Illegal Stokes string " << stokes_p
	 << LogIO::POST;

      

#ifdef PABLO_IO
      traceEvent(1,"Exiting Imager::setimage",25);
#endif

      return False;
    };
    

    //THIS NEEDS TO GO
    this->setImageParam(nx_p, ny_p, npol_p, imageNchan_p);

    // Now do the shifts
    //    MSColumns msc(*ms_p);

    doShift_p=doShift;
    if(doShift_p) {
      phaseCenter_p=phaseCenter;
    }
    else {

      ROMSFieldColumns msfield(ms_p->field());
      phaseCenter_p=msfield.phaseDirMeas(fieldid_p);
      //    phaseCenter_p=msc.field().phaseDirMeas(fieldid_p);
    }
    
    // Now add the optional shifts
    shiftx_p=shiftx;
    shifty_p=shifty;
    if(shiftx_p.get().getValue()!=0.0||shifty_p.get().getValue()!=0.0) {
      Vector<Double> vPhaseCenter(phaseCenter_p.getAngle().getValue());
      if(cos(vPhaseCenter(1))!=0.0) {
	vPhaseCenter(0)+=shiftx_p.get().getValue()/cos(vPhaseCenter(1));
      }
      vPhaseCenter(1)+=shifty_p.get().getValue();
      phaseCenter_p.set(MVDirection(vPhaseCenter));
    }
    
    // Now we have set the image parameters
    setimaged_p=True;
    beamValid_p=False;
    
    this->unlock();

#ifdef PABLO_IO
    traceEvent(1,"Exiting Imager::setimage",25);
#endif

    return True;
  } catch (AipsError x) {
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << LogIO::POST;
    this->unlock();

#ifdef PABLO_IO
    traceEvent(1,"Exiting Imager::setimage",25);
#endif

    return False;
  } 

#ifdef PABLO_IO
  traceEvent(1,"Exiting Imager::setimage",25);
#endif

  return True;
}

Bool Imager::advise(const Bool takeAdvice, const Float amplitudeLoss,
		    const Quantity& fieldOfView, Quantity& cell,
		    Int& pixels, Int& facets, MDirection& phaseCenter)
{
  if(!valid()) return False;
  
  LogIO os(LogOrigin("imager", "advise()", WHERE));
  
  try {
    
    os << "Advising image properties" << LogIO::POST;
    
    Float maxAbsUV=0.0;
    Float maxWtAbsUV=0.0;
    // To determine the number of facets, we need to fit w to
    // a.u + b.v. The misfit from this (i.e. the dispersion 
    // will determine the error beam due to the non-coplanar
    // baselines. We'll do both cases: where the position
    // errors are important and where they are not. We'll use
    // the latter.
    Double sumWt = 0.0;

    Double sumUU=0.0;
    Double sumUV=0.0;
    Double sumUW=0.0;
    Double sumVV=0.0;
    Double sumVW=0.0;
    Double sumWW=0.0;

    Double sumWtUU=0.0;
    Double sumWtUV=0.0;
    Double sumWtUW=0.0;
    Double sumWtVV=0.0;
    Double sumWtVW=0.0;
    Double sumWtWW=0.0;

    Double sum = 0.0;

    this->lock();
    VisIter& vi(vs_p->iter());
    VisBuffer vb(vi);
    
    for (vi.originChunks();vi.moreChunks();vi.nextChunk()) {
      for (vi.origin();vi.more();vi++) {
	Int nRow=vb.nRow();
	Int nChan=vb.nChannel();
	for (Int row=0; row<nRow; ++row) {
	  for (Int chn=0; chn<nChan; ++chn) {
	    if(!vb.flag()(chn,row)) {
	      Float f=vb.frequency()(chn)/C::c;
	      Float u=vb.uvw()(row)(0)*f;
	      Float v=vb.uvw()(row)(1)*f;
	      Float w=vb.uvw()(row)(2)*f;
              Double wt=vb.imagingWeight()(chn,row);
	      if(wt>0.0) {
		if(abs(u)>maxWtAbsUV) maxWtAbsUV=abs(u);
		if(abs(v)>maxWtAbsUV) maxWtAbsUV=abs(v);
		sumWt += wt;
                sumWtUU += wt * u * u;
                sumWtUV += wt * u * v;
                sumWtUW += wt * u * w;
                sumWtVV += wt * v * v;
                sumWtVW += wt * v * w;
                sumWtWW += wt * w * w;
	      }
	      sum += 1;
	      if(abs(u)>maxAbsUV) maxAbsUV=abs(u);
	      if(abs(v)>maxAbsUV) maxAbsUV=abs(v);
	      sumUU += u * u;
	      sumUV += u * v;
	      sumUW += u * w;
	      sumVV += v * v;
	      sumVW += v * w;
	      sumWW += w * w;
	    }
	  }
	}
      }
    }
    
    if(sumWt==0.0) {
      os << LogIO::WARN << "Visibility data are not yet weighted: using unweighted values" << LogIO::POST;
      sumWt = sum;
    }
    else {
      sumUU = sumWtUU;
      sumUV = sumWtUV;
      sumUW = sumWtUW;
      sumVV = sumWtVV;
      sumVW = sumWtVW;
      sumWW = sumWtWW;
      maxAbsUV = maxWtAbsUV;
    }

    // First find the cell size
    if(maxAbsUV==0.0) {
      this->unlock();
      os << LogIO::SEVERE << "Maximum uv distance is zero" << LogIO::POST;
      return False;
    }
    else {
      cell=Quantity(0.5/maxAbsUV, "rad").get("arcsec");
      os << "Maximum uv distance = " << maxAbsUV << " wavelengths" << endl;
      os << "Recommended cell size < " << cell.get("arcsec").getValue()
	 << " arcsec" << LogIO::POST;
    }

    // Now we can find the number of pixels for the specified field of view
    pixels = 2*Int((fieldOfView.get("rad").getValue()/cell.get("rad").getValue())/2.0);
    {
      CompositeNumber cn(pixels);
      pixels = (Int) (cn.nextLargerEven(pixels));
    }
    if(pixels < 64) pixels = 64;
    os << "Recommended number of pixels = " << pixels << endl;

      // Rough rule for number of facets:
      // For the specified facet size, the loss in amplitude
      // due to the peeling of facets from the sphere should 
      // be equal to the amplitude error.
      Int worstCaseFacets=1;
      if(sumWt<=0.0||sumUU<=0.0||(sumUU+sumVV)<=0.0) {
	this->unlock();
	os << LogIO::SEVERE << "Sum of imaging weights is zero" << LogIO::POST;
	return False;
      }
      else {
	Double rmsUV  = sqrt((sumUU + sumVV)/sumWt);
	Double rmsW = sqrt(sumWW/sumWt);
	os << "Dispersion in uv, w distance = " << rmsUV << ", "<< rmsW
	   << " wavelengths" << endl;
	if(rmsW>0.0&&rmsUV>0.0&&amplitudeLoss>0.0) {
	  worstCaseFacets =
	    Int (pixels * (abs(cell.get("rad").getValue())*
				  sqrt(C::pi*rmsW/(sqrt(32.0*amplitudeLoss)))));
	}
	else {
	  os << LogIO::WARN << "Cannot calculate number of facets: using 1"
	     << LogIO::POST;
	  worstCaseFacets = 1;
	}
	// Solve for the parameters:
	Double Determinant = sumUU * sumVV - square(sumUV);
	Double rmsFittedW = rmsW;
	if(Determinant > 0.0) {
	  Double a = ( sumVV * sumUW - sumUV * sumVW)/Determinant;
	  Double b = (-sumUV * sumUW + sumUU * sumVW)/Determinant;
	  os << "Best fitting plane is w = " << a << " * u + "
	     << b << " * v" << endl;
	  Double FittedWW =
	    sumWW  + square(a) * sumUU + square(b) * sumVV +
	    + 2.0 * a * b * sumUV - 2.0 * (a * sumUW + b * sumVW);
	  rmsFittedW  = sqrt(FittedWW/sumWt);
	  os << "Dispersion in fitted w = " << rmsFittedW
	     << " wavelengths" << endl;
	  facets = Int (pixels * (abs(cell.get("rad").getValue())*
				  sqrt(C::pi*rmsFittedW/(sqrt(32.0*amplitudeLoss)))));
          if (facets<1) facets = 1;
	}
	else {
	  os << "Error in fitting plane to uvw data" << LogIO::POST;
	}
	if(worstCaseFacets<1) worstCaseFacets=1;
	if(worstCaseFacets>1) {
	  os << "imager recommends that you use the wide field clean" << endl
	     << "For accurate positions, use " << worstCaseFacets
	     << " facets on each axis" << endl
	     << "For accurate removal of sources, you only need "
	     << facets << " facets on each axis" << LogIO::POST;
	}
	else {
	  os << "Wide field cleaning is not necessary"
	     << LogIO::POST;
	}
      }

    MSColumns msc(*mssel_p);
    if(datafieldids_p.shape()!=0){
      //If setdata has been used prior to this
    phaseCenter=msc.field().phaseDirMeas(datafieldids_p(0));
    }
    else{
    phaseCenter=msc.field().phaseDirMeas(fieldid_p);   
    }

    
    // Now we have set the image parameters
    if(takeAdvice) {
      os << "Using advised image properties" << LogIO::POST;
      mcellx_p=cell;
      mcelly_p=cell;
      phaseCenter_p=phaseCenter;
      setimaged_p=True;
      beamValid_p=False;
      facets_p=facets;
      nx_p=ny_p=pixels;
    }
    
    this->unlock();
    return True;
  } catch (AipsError x) {
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << LogIO::POST;
    this->unlock();
    return False;
  } 
  
  return True;
}



Bool Imager::setDataPerMS(const String& msname, const String& mode, 
			  const Vector<Int>& nchan, 
			  const Vector<Int>& start,
			  const Vector<Int>& step,
			  const Vector<Int>& spectralwindowids,
			  const Vector<Int>& fieldids,
			  const String& msSelect)
{
  LogIO os(LogOrigin("imager", "setdata()"), logSink_p);
  if(msname != ""){

    LogIO os(LogOrigin("imager", "setdata()"), logSink_p);
    os << LogIO::WARN
       << "Ignoring that ms" << msname << "specified here"
       << LogIO::POST;
    os << LogIO::WARN
       << "Imager was constructed with an ms "
       << LogIO::POST;
    os << LogIO::WARN
       << "if multi-ms are to be used please construct imager without parameters and use setdata to specify the ms's and selection"
       << LogIO::POST;

  }
  MRadialVelocity dummy;
  //Calling the old setdata
  return   setdata(mode, nchan, start, step, dummy, dummy, spectralwindowids, 
		   fieldids, msSelect);

}


Bool Imager::setdata(const String& mode, const Vector<Int>& nchan,
		     const Vector<Int>& start, const Vector<Int>& step,
		     const MRadialVelocity& mStart,
		     const MRadialVelocity& mStep,
		     const Vector<Int>& spectralwindowids,
		     const Vector<Int>& fieldids,
		     const String& msSelect)
  
{
  logSink_p.clearLocally();
  LogIO os(LogOrigin("imager", "setdata()"), logSink_p);

  if(!ms_p) {
    os << LogIO::SEVERE << "Program logic error: MeasurementSet pointer ms_p not yet set"
       << LogIO::POST;
    return False;
  }

  os << "mode=" << mode << " nchan=" << nchan 
     <<  " start=" << start << " step=" << step;
  ostringstream clicom;
  clicom <<  " mstart='" << mStart << "' mstep='" << mStep;
  os << String(clicom) ;
  os <<  "' spectralwindowids=" << spectralwindowids
     << " fieldids=" << fieldids << " msselect=" << msSelect;

  try {
    
    this->lock();
    this->writeCommand(os);

    os << "Selecting data" << LogIO::POST;
    nullSelect_p=False;
    dataMode_p=mode;
    dataNchan_p.resize();
    dataStart_p.resize();
    dataStep_p.resize();
    dataNchan_p=nchan;
    dataStart_p=start;
    dataStep_p=step;
    mDataStart_p=mStart;
    mDataStep_p=mStep;
    dataspectralwindowids_p.resize(spectralwindowids.nelements());
    dataspectralwindowids_p=spectralwindowids;
    datafieldids_p.resize(fieldids.nelements());
    datafieldids_p=fieldids;
    
    if (fieldids.nelements() > 1) {
      os << "Multiple fields specified via fieldids" << LogIO::POST;
      multiFields_p = True;
    }

   // Map the selected spectral window ids to data description ids
    MSDataDescColumns dataDescCol(ms_p->dataDescription());
    Vector<Int> ddSpwIds=dataDescCol.spectralWindowId().getColumn();

    datadescids_p.resize(0);
    for (uInt row=0; row<ddSpwIds.nelements(); ++row) {
      Bool found=False;
      for (uInt j=0; j<dataspectralwindowids_p.nelements(); ++j) {
	if (ddSpwIds(row)==dataspectralwindowids_p(j)) found=True;
      };
      if (found) {
	datadescids_p.resize(datadescids_p.nelements()+1,True);
	datadescids_p(datadescids_p.nelements()-1)=row;
      };
    };

    // If a selection has been made then close the current MS
    // and attach to a new selected MS. We do this on the original
    // MS. 
    if(datafieldids_p.nelements()>0||datadescids_p.nelements()>0) {
      os << "Performing selection on MeasurementSet" << LogIO::POST;
      if(vs_p) delete vs_p; vs_p=0;
      if(mssel_p) delete mssel_p; mssel_p=0;
      
      // check that sorted table exists (it should), if not, make it now.
      this->makeVisSet(*ms_p);
      
      Table sorted=ms_p->keywordSet().asTable("SORTED_TABLE");
      
      
      // Now we make a condition to do the old FIELD_ID, SPECTRAL_WINDOW_ID
      // selection
      TableExprNode condition;
      String colf=MS::columnName(MS::FIELD_ID);
      String cols=MS::columnName(MS::DATA_DESC_ID);
      if(datafieldids_p.nelements()>0&&datadescids_p.nelements()>0){
	condition=sorted.col(colf).in(datafieldids_p)&&
	  sorted.col(cols).in(datadescids_p);
        os << "Selecting on field and spectral window ids" << LogIO::POST;
      }
      else if(datadescids_p.nelements()>0) {
	condition=sorted.col(cols).in(datadescids_p);
        os << "Selecting on spectral window id" << LogIO::POST;
      }
      else if(datafieldids_p.nelements()>0) {
	condition=sorted.col(colf).in(datafieldids_p);
        os << "Selecting on field id" << LogIO::POST;
      }
      
      // Now remake the selected ms
      mssel_p = new MeasurementSet(sorted(condition));

      AlwaysAssert(mssel_p, AipsError);
      mssel_p->rename(msname_p+"/SELECTED_TABLE", Table::Scratch);
      if(mssel_p->nrow()==0) {
	delete mssel_p; mssel_p=0;
	os << LogIO::WARN
	   << "Selection is empty: reverting to sorted MeasurementSet"
	   << LogIO::POST;
	mssel_p=new MeasurementSet(sorted);
	nullSelect_p=True;
      }
      else {
	mssel_p->flush();
	nullSelect_p=False;
      }
      if (nullSelect_p) {
	Table mytab(msname_p+"/FIELD", Table::Old);
	if (mytab.nrow() > 1) {
	  os << "Multiple fields selected" << LogIO::POST;
	   multiFields_p = True;
	} else {
	  os << "Single field selected" << LogIO::POST;
	   multiFields_p = False;
	}
      }

      Int len = msSelect.length();
      Int nspace = msSelect.freq (' ');
      Bool nullSelect=(msSelect.empty() || nspace==len);
      if (!nullSelect) {
	MeasurementSet* mssel_p2;
	// Apply the TAQL selection string, to remake the selected MS
	String parseString="select from $1 where " + msSelect;
	mssel_p2=new MeasurementSet(tableCommand(parseString,*mssel_p));
	AlwaysAssert(mssel_p2, AipsError);
	// Rename the selected MS as */SELECTED_TABLE2
	mssel_p2->rename(msname_p+"/SELECTED_TABLE2", Table::Scratch); 
	if (mssel_p2->nrow()==0) {
	  os << LogIO::WARN
	     << "Selection string results in empty MS: "
	     << "reverting to sorted MeasurementSet"
	     << LogIO::POST;
	  delete mssel_p2;
	} else {
	  if (mssel_p) {
	    delete mssel_p; 
	    mssel_p=mssel_p2;
	    mssel_p->flush();
	  }
	}
      } else {
	os << "No selection string given" << LogIO::POST;
      }

      if(mssel_p->nrow()!=ms_p->nrow()) {
	os << "By selection " << ms_p->nrow() << " rows are reduced to "
	   << mssel_p->nrow() << LogIO::POST;
      }
      else {
	os << "Selection did not drop any rows" << LogIO::POST;
      }
    }
    
    // Now create the VisSet
    this->makeVisSet(vs_p, *mssel_p); 
    AlwaysAssert(vs_p, AipsError);
    
    // Now we do a selection to cut down the amount of information
    // passed around.

    this->selectDataChannel(*vs_p, dataspectralwindowids_p, dataMode_p,
			    dataNchan_p, dataStart_p, dataStep_p,
			    mDataStart_p, mDataStep_p);

    // Guess that the beam is no longer valid
    beamValid_p=False;
    destroySkyEquation();
    if(!valid()){ 
      this->unlock();
      os << LogIO::SEVERE << "Check your data selection or Measurement set " << LogIO::POST;
      return False;
    }
    this->unlock();
    return True;
  } catch (AipsError x) {
    os << LogIO::SEVERE << "Caught exception: " << x.getMesg()
       << LogIO::POST;
    this->unlock();
    return False;
  } 
  return True;
}


Bool Imager::setmfcontrol(const Float cyclefactor,
			  const Float cyclespeedup,
			  const Int stoplargenegatives, 
			  const Int stoppointmode,
			  const String& scaleType,
			  const Float minPB,
			  const Float constPB,
			  const Vector<String>& fluxscale)
{  
  cyclefactor_p = cyclefactor;
  cyclespeedup_p =  cyclespeedup;
  stoplargenegatives_p = stoplargenegatives;
  stoppointmode_p = stoppointmode;
  fluxscale_p.resize( fluxscale.nelements() );
  fluxscale_p = fluxscale;
  scaleType_p = scaleType;
  minPB_p = minPB;
  constPB_p = constPB;
  return True;
}  


Bool Imager::setvp(const Bool dovp,
		   const Bool doDefaultVPs,
		   const String& vpTable,
		   const Bool doSquint,
		   const Quantity &parAngleInc,
		   const Quantity &skyPosThreshold,
		   String defaultTel)
{

#ifdef PABLO_IO
  traceEvent(1,"Entering Imager::setvp",23);
#endif

  //  if(!valid())
  //    {

  //#ifdef PABLO_IO
  //      traceEvent(1,"Exiting Imager::setvp",22);
  //#endif

  //     return False;
  //    }
  LogIO os(LogOrigin("Imager", "setvp()", WHERE));
  
  os << "Setting voltage pattern parameters" << LogIO::POST;
  
  doVP_p=dovp;
  doDefaultVP_p = doDefaultVPs;
  vpTableStr_p = vpTable;
  telescope_p= defaultTel;
  if (doSquint) {
    squintType_p = BeamSquint::GOFIGURE;
  } else {
    squintType_p = BeamSquint::NONE;
  }

  parAngleInc_p = parAngleInc;

  skyPosThreshold_p = skyPosThreshold;
  os<<"Sky position tolerance is "<<skyPosThreshold_p.getValue("deg")<<
      " degrees" << LogIO::POST;

  if (doDefaultVP_p) {
    os << "Using system default voltage patterns for each telescope"  << LogIO::POST;
  } else {
    os << "Using user defined voltage patterns in Table "<<  vpTableStr_p << LogIO::POST;
  }
  if (doSquint) {
    os << "Beam Squint will be included in the VP model" <<  LogIO::POST;
    os << "and the Parallactic Angle increment is " 
       << parAngleInc_p.getValue("deg") << " degrees"  << LogIO::POST;
  }

#ifdef PABLO_IO
  traceEvent(1,"Exiting Imager::setvp",22);
#endif

  // muddled with the state of SkyEquation..so redo it
  destroySkyEquation();
  return True;

}

Bool Imager::setoptions(const String& ftmachine, const Long cache, const Int tile,
			const String& gridfunction, const MPosition& mLocation,
			const Float padding, const Bool usemodelcol, 
			const Int wprojplanes,
			const String& epJTableName,
			const Bool applyPointingOffsets,
			const Bool doPointingCorrection,
			const String& cfCacheDirName)
{

#ifdef PABLO_IO
  traceEvent(1,"Entering Imager::setoptions",28);
#endif

  if(!valid()) 
    {

#ifdef PABLO_IO
      traceEvent(1,"Exiting Imager::setoptions",27);
#endif

      return False;
    }
  if(!assertDefinedImageParameters())
    {

#ifdef PABLO_IO
      traceEvent(1,"Exiting Imager::setoptions",27);
#endif

      return False;
    }
  LogIO os(LogOrigin("imager", "setoptions()", WHERE));
  
  os << "Setting processing options" << LogIO::POST;
  useModelCol_p=usemodelcol;

  ftmachine_p=downcase(ftmachine);
  if(ftmachine_p=="gridft") {
    os << "FT machine gridft is now called ft - please use the new name in future" << endl;
    ftmachine_p="ft";
  }

  if(ftmachine_p=="wfmemoryft"){
    wfGridding_p=True;
    ftmachine_p="ft";
  }

  wprojPlanes_p=wprojplanes;
  epJTableName_p = epJTableName;
  cfCacheDirName_p = cfCacheDirName;

  if(cache>0) cache_p=cache;
  if(tile>0) tile_p=tile;
  gridfunction_p=downcase(gridfunction);
  mLocation_p=mLocation;
  if(padding>=1.0) {
    padding_p=padding;
  }
  // Destroy the FTMachine
  if(ft_p) {delete ft_p; ft_p=0;}
  if(gvp_p) {delete gvp_p; gvp_p=0;}
  if(cft_p) {delete cft_p; cft_p=0;}

#ifdef PABLO_IO
  traceEvent(1,"Exiting Imager::setoptions",27);
#endif

  doPointing = applyPointingOffsets;
  doPBCorr = doPointingCorrection;

  return True;
}

Bool Imager::setsdoptions(const Float scale, const Float weight, 
			  const Int convsupport)
{

#ifdef PABLO_IO
  traceEvent(1,"Entering Imager::setsdoptions",28);
#endif

  if(!valid()) 
    {

#ifdef PABLO_IO
      traceEvent(1,"Exiting Imager::setsdoptions",27);
#endif

      return False;
    }

  LogIO os(LogOrigin("imager", "setsdoptions()", WHERE));
  
  os << "Setting single dish processing options" << LogIO::POST;
  
  sdScale_p=scale;
  sdWeight_p=weight;
  sdConvSupport_p=convsupport;

  // Destroy the FTMachine
  if(ft_p) {delete ft_p; ft_p=0;}
  if(gvp_p) {delete gvp_p; gvp_p=0;}
  if(cft_p) {delete cft_p; cft_p=0;}

#ifdef PABLO_IO
  traceEvent(1,"Exiting Imager::setsdoptions",27);
#endif

  return True;
}


} //# NAMESPACE CASA - END

