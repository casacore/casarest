//# Imager.cc: Implementation of Imager.h
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

Imager::Imager() 
  :  ms_p(0),msname_p(""), mssel_p(0), vs_p(0), ft_p(0), cft_p(0), se_p(0),
    sm_p(0), vp_p(0), gvp_p(0), setimaged_p(False), nullSelect_p(False), pgplotter_p(0)
{
  defaults();
};


void Imager::defaults() 
{

#ifdef PABLO_IO
traceEvent(1,"Entering imager::defaults",25);
#endif

  setimaged_p=False;
  nullSelect_p=False;
  nx_p=128; ny_p=128; facets_p=1;
  wprojPlanes_p=1;
  mcellx_p=Quantity(1, "arcsec"); mcelly_p=Quantity(1, "arcsec");
  shiftx_p=Quantity(0.0, "arcsec"); shifty_p=Quantity(0.0, "arcsec");
  distance_p=Quantity(0.0, "m");
  stokes_p="I"; npol_p=1;
  nscales_p=5;
  ntaylor_p=2;
  scaleMethod_p="nscales";  
  scaleInfoValid_p=False;
  dataMode_p="none";
  imageMode_p="mfs";
  dataNchan_p=0;
  imageNchan_p=0;
  doVP_p=False;
  doDefaultVP_p = True;
  parAngleInc_p=Quantity(360.,"deg");
  skyPosThreshold_p=Quantity(180.,"deg");
  telescope_p="";
  gridfunction_p="SF";
  doMultiFields_p=False;
  doWideBand_p=False;
  multiFields_p=False;
  // Use half the machine memory as cache. The user can override
  // this via the setoptions function().
  cache_p=(HostInfo::memoryTotal()/8)*1024;
  //On 32 bit machines with more than 2G of mem this can become negative
  // overriding it to 2 Gb.
  if(cache_p <=0 )
    cache_p=2000000000/8;
  tile_p=16;
  ftmachine_p="ft";
  wfGridding_p=False;
  padding_p=1.2;
  sdScale_p=1.0;
  sdWeight_p=1.0;
  sdConvSupport_p=-1;

  doShift_p=False;
  spectralwindowids_p.resize(1); 
  spectralwindowids_p=0;
  fieldid_p=0;
  dataspectralwindowids_p.resize(0); 
  datadescids_p.resize(0);
  datafieldids_p.resize(0);
  mImageStart_p=MRadialVelocity(Quantity(0.0, "km/s"), MRadialVelocity::LSRK);
  mImageStep_p=MRadialVelocity(Quantity(0.0, "km/s"), MRadialVelocity::LSRK);
  mDataStart_p=MRadialVelocity(Quantity(0.0, "km/s"), MRadialVelocity::LSRK);
  mDataStep_p=MRadialVelocity(Quantity(0.0, "km/s"), MRadialVelocity::LSRK);
  beamValid_p=False;
  bmaj_p=Quantity(0, "arcsec");
  bmin_p=Quantity(0, "arcsec");
  bpa_p=Quantity(0, "deg");
  images_p.resize(0);
  masks_p.resize(0);
  fluxMasks_p.resize(0);
  residuals_p.resize(0);
  componentList_p=0;

  cyclefactor_p = 1.5;
  cyclespeedup_p =  -1;
  stoplargenegatives_p = 2;
  stoppointmode_p = -1;
  fluxscale_p.resize(0);
  scaleType_p = "NONE";
  minPB_p = 0.1;
  constPB_p = 0.4;
  redoSkyModel_p=True;
  nmodels_p=0;
  useModelCol_p=True;  
  freqFrameValid_p=False;
  logSink_p=LogSink(LogMessage::NORMAL, False);
#ifdef PABLO_IO
  traceEvent(1,"Exiting imager::defaults",24);
#endif

}


Imager::Imager(MeasurementSet& theMS,  Bool compress)
: ms_p(0), msname_p(""), mssel_p(0), vs_p(0), ft_p(0), cft_p(0), se_p(0),
    sm_p(0), vp_p(0), gvp_p(0), setimaged_p(False), nullSelect_p(False), pgplotter_p(0)
{
  lockCounter_p=0;
  LogIO os(LogOrigin("Imager", "Imager(MeasurementSet &theMS)", WHERE));
  if(!open(theMS, compress)) {
    os << LogIO::SEVERE << "Open of MeasurementSet failed" << LogIO::EXCEPTION;
  };

  defaults();
  latestObsInfo_p=ObsInfo();
}



Imager::Imager(MeasurementSet& theMS, PGPlotter& thePlotter, Bool compress)
: ms_p(0), msname_p(""), mssel_p(0), vs_p(0), ft_p(0), cft_p(0), se_p(0),
    sm_p(0), vp_p(0), gvp_p(0), setimaged_p(False), nullSelect_p(False), pgplotter_p(0)
{
  lockCounter_p=0;
  LogIO os(LogOrigin("Imager", "Imager(MeasurementSet &theMS)", WHERE));
  if(!open(theMS, compress)) {
    os << LogIO::SEVERE << "Open of MeasurementSet failed" << LogIO::EXCEPTION;
  };

  defaults();

  pgplotter_p=&thePlotter;
  latestObsInfo_p=ObsInfo();
}

Imager::Imager(const Imager & other)
  :  ms_p(0),msname_p(""), mssel_p(0), vs_p(0), ft_p(0), cft_p(0), se_p(0),
    sm_p(0), vp_p(0), gvp_p(0), setimaged_p(False), nullSelect_p(False), pgplotter_p(0)
{
  operator=(other);
}

Imager &Imager::operator=(const Imager & other)
{
  if (ms_p && this != &other) {
    *ms_p = *(other.ms_p);
  }
  //Equating the table and ms parameters
  antab_p=other.antab_p;
  datadesctab_p=other.datadesctab_p;
  feedtab_p=other.feedtab_p;
  fieldtab_p=other.fieldtab_p;
  obstab_p=other.obstab_p;
  pointingtab_p=other.pointingtab_p;
  poltab_p=other.poltab_p;
  proctab_p=other.proctab_p;
  spwtab_p=other.spwtab_p;
  statetab_p=other.statetab_p;
  latestObsInfo_p=other.latestObsInfo_p;
  parAngleInc_p=other.parAngleInc_p;
  skyPosThreshold_p=other.skyPosThreshold_p;
  if (mssel_p && this != &other) {
    *mssel_p = *(other.mssel_p);
  }
  if (vs_p && this != &other) {
    *vs_p = *(other.vs_p);
  }
  if (ft_p && this != &other) {
    *ft_p = *(other.ft_p);
  }
  if (cft_p && this != &other) {
    *cft_p = *(other.cft_p);
  }
  if (se_p && this != &other) {
    *se_p = *(other.se_p);
  }
  if (sm_p && this != &other) {
    *sm_p = *(other.sm_p);
  }
  if (vp_p && this != &other) {
    *vp_p = *(other.vp_p);
  }
  if (gvp_p && this != &other) {
    *gvp_p = *(other.gvp_p);
  }
  if (pgplotter_p && this != &other) {
    *pgplotter_p = *(other.pgplotter_p);
  }

  return *this;
}

Imager::~Imager()
{

  destroySkyEquation();
  this->unlock(); //unlock things if they are in a locked state

  if (mssel_p) {
    delete mssel_p;
  }
  mssel_p = 0;
  if (ms_p) {
    delete ms_p;
  }
  ms_p = 0;
  if (vs_p) {
    delete vs_p;
  }
  vs_p = 0;
  if (ft_p) {
    delete ft_p;
  }
  ft_p = 0;
  if (cft_p) {
    delete cft_p;
  }
  cft_p = 0;

  //Note we don't deal with pgplotter here.
  

}


Bool Imager::open(MeasurementSet& theMs, Bool compress)
{

#ifdef PABLO_IO
  traceEvent(1,"Entering Imager::open",21);
#endif

  LogIO os(LogOrigin("Imager", "open()", WHERE));
  
  if (ms_p) {
    *ms_p = theMs;
  } else {
    ms_p = new MeasurementSet(theMs);
    AlwaysAssert(ms_p, AipsError);
  }
  

  try {
    this->openSubTables();
    this->lock();
    msname_p = ms_p->tableName();
    
    os << "Opening MeasurementSet " << msname_p << LogIO::POST;

    // Check for DATA or FLOAT_DATA column
    if(!ms_p->tableDesc().isColumn("DATA") && 
       !ms_p->tableDesc().isColumn("FLOAT_DATA")) {
      os << LogIO::SEVERE
	 << "Missing DATA or FLOAT_DATA column: imager cannot be run"
	 << LogIO::POST;
      ms_p->unlock();
      delete ms_p; ms_p=0;
      return False;
    }
    
    Bool initialize=(!ms_p->tableDesc().isColumn("CORRECTED_DATA"));
    
    if(vs_p) {
      delete vs_p; vs_p=0;
    }
    
    // Now open the selected MeasurementSet to be initially the
    // same as the original MeasurementSet

    mssel_p=new MeasurementSet(*ms_p);
    
    
    // Now create the VisSet
    this->makeVisSet(vs_p, *mssel_p);
    AlwaysAssert(vs_p, AipsError);
    
    // Polarization
    MSColumns msc(*mssel_p);
    Vector<String> polType=msc.feed().polarizationType()(0);
    if (polType(0)!="X" && polType(0)!="Y" &&
	polType(0)!="R" && polType(0)!="L") {
      this->unlock();
      os << LogIO::SEVERE << "Warning: Unknown stokes types in feed table: "
	 << polType(0) << endl
	 << "Results open to question!" << LogIO::POST;
    }
    
    
    // Initialize the weights if the IMAGING_WEIGHT column
    // was just created
    if(initialize) {
      os << LogIO::NORMAL
	 << "Initializing natural weights"
	 << LogIO::POST;
      Double sumwt=0.0;
      VisSetUtil::WeightNatural(*vs_p, sumwt);
    }
    this->unlock();

#ifdef PABLO_IO
    traceEvent(1,"Exiting Imager::open",21);
#endif

    return True;
  } catch (AipsError x) {
    this->unlock();
    os << LogIO::SEVERE << "Caught Exception: "<< x.getMesg() << LogIO::POST;

#ifdef PABLO_IO
    traceEvent(1,"Exiting Imager::open",21);
#endif

    return False;
  } 

#ifdef PABLO_IO
  traceEvent(1,"Exiting Imager::open",21);
#endif

  return True;
}

Bool Imager::close()
{
  if(!valid()) return False;
  if (detached()) return True;
  LogIO os(LogOrigin("imager", "close()", WHERE));
  os << "Closing MeasurementSet and detaching from imager"
     << LogIO::POST;
  this->unlock();
  if(ft_p) delete ft_p; ft_p = 0;
  if(cft_p) delete cft_p; cft_p = 0;
  if(vs_p) delete vs_p; vs_p = 0;
  if(mssel_p) delete mssel_p; mssel_p = 0;
  if(ms_p) delete ms_p; ms_p = 0;

  if(se_p) delete se_p; se_p = 0;

  if(vp_p) delete vp_p; vp_p = 0;
  if(gvp_p) delete gvp_p; gvp_p = 0;

  // if(pgplotter_p) delete pgplotter_p; pgplotter_p = 0;
  destroySkyEquation();

  return True;
}

String Imager::name() const
{
  if (detached()) {
    return "none";
  }
  return msname_p;
}

String Imager::imageName()
{
  LogIO os(LogOrigin("imager", "imageName()", WHERE));
  try {
    lock();
    String name(msname_p);
    MSColumns msc(*ms_p);
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
    os << LogIO::SEVERE << "Caught Exception: "<< x.getMesg() << LogIO::POST; return "";
  } 
  return String("imagerImage");
}


} //# NAMESPACE CASA - END

