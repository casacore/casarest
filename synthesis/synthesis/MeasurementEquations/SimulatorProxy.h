//# SimulatorProxy.h: Based on DOnewsimulator.h
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002
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
//#
//# $Id: SimulatorProxy.h,v 1.1 2007/02/12 00:05:58 gvandiep Exp $

#ifndef SYNTHESIS_SIMULATORPROXY_H
#define SYNTHESIS_SIMULATORPROXY_H

#include <casa/aips.h>
#include <casa/Quanta/Quantum.h>
#include <measures/Measures/MPosition.h>
#include <synthesis/MeasurementComponents/BeamSquint.h>
#include <synthesis/MeasurementComponents/VPSkyJones.h>
#include <synthesis/MeasurementComponents/EPJones.h>

#include <casa/namespace.h>
namespace casa { //# NAMESPACE CASA - BEGIN
class MeasurementSet;
class VisSet;
class VisJones;
class ACoh;
class SkyEquation;
class ComponentList;
class CleanImageSkyModel;
class FTMachine;
class ComponentFTMachine;
class MEpoch;
class NewMSSimulator;


template<class T> class PagedImage;



// <summary>Simulates MeasurementSets from SkyModel and SkyEquation</summary>


// <use visibility=export>

// <reviewed reviewer="" date="yyyy/mm/dd" tests="" demos="">
// </reviewed>

// <prerequisite>
//   <li> <linkto class="MeasurementSet">MeasurementSet</linkto>
//   <li> <linkto class="SkyEquation">SkyEquation</linkto>
//   <li> <linkto class="SkyModel">SkyModel</linkto>
// </prerequisite>
//
// <synopsis>
// </synopsis>
//
// <motivation> 
// This class was written to make the simulation capability useful from python
// </motivation>
//
// <thrown>
// <li> AipsError - if an error occurs
// </thrown>
//
// <todo asof="1999/07/22">
//   <li> everything
// </todo>

class SimulatorProxy
{

public:
  SimulatorProxy() {};

  // Construct from string
  SimulatorProxy(const String& msname);

  // "SimulatorProxy" ctor
  SimulatorProxy(MeasurementSet &thems);
  
  SimulatorProxy(const SimulatorProxy &other);

  SimulatorProxy &operator=(const SimulatorProxy &other);
 
  ~SimulatorProxy();
  
  // Select the data to be predicted or corrupted
  Bool setdata(const Vector<Int>& spectralwindowids,
	       const Vector<Int>& fieldids,
	       const String& msSelect);
  
  Bool settimes(const Record& integrationTime, 
		const Bool      useHourAngle,
		const Record&   refTime);

  // set the configuration; NOTE: the telname used
  // here will determine the Voltage Pattern to be used in the
  // simulation
  Bool setconfig(const String& telname,
		 const Vector<Double>& x, 
		 const Vector<Double>& y, 
		 const Vector<Double>& z,
		 const Vector<Double>& dishDiameter,
		 const Vector<Double>& offset,
		 const Vector<String>& mount,
		 const Vector<String>& antName,
		 const String& coordsystem,
		 //const MPosition& referenceLocation);
		 const Record& recRefLocation);
  
  // set the observed fields for the simulation
  Bool setfield(const String& sourceName,           
		const Record& sourceDirection,  
		const String& calCode,
		const Record& distance);
  
  // set one or more spectral windows and their characteristics
  Bool setspwindow(const String& spwName,
		   const Record& freq,
		   const Record& deltafreq,
		   const Record& freqresolution,
		   const Int nchannels,
		   const String& stokes);
  
  // Set the simulated feed characteristics
  Bool setfeed(const String& mode,
	       const Vector<Double>& x,
	       const Vector<Double>& y,
	       const Vector<String>& pol);		

  // Set limits
  Bool setlimits(const Double shadowFraction,
		 const Record& elevationLimit);

  // Set the voltage pattern
  Bool setvp(const Bool dovp,
             const Bool defaultVP,
             const String& vpTable,
             const Bool doSquint,
             const Record &parAngleInc,
	     const Record &skyPosThreshold,
	     const Float &pbLimit);

  
  // Set the random number generator seed for the addition of errors
  Bool setseed(const Int seed);

  // Apply antenna-based gain errors
  Bool setgain(const String& mode, const String& table,
	       const Record& interval, const Vector<Double>& amplitude);

  /*
  // Apply antenna pointing and squint errors
  Bool setpointingerror(const String& epJTableName,
			const Bool applyPointingOffsets,
			const Bool doPBCorrection);
  */

  // Apply polarization leakage errors
  Bool setleakage(const String& mode, const String& table,
		  const Record& interval, const Double amplitude);

  // Apply bandpass errors
  Bool setbandpass(const String& mode, const String& table,
		   const Record& interval, const Vector<Double>& amplitude);

  /*
  // Simulate the parallactic angle phase effect
  Bool setpa(const String& mode, const String& table,
	     const Quantity& interval);
  */

  // Simulate quasi-realistic thermal noise, which can depend upon
  // elevation, bandwidth, antenna diameter, as expected
  Bool setnoise(const String& mode, 
		const Record& simplenoise,
		const String& table,
		const Float antefficiency,
		const Float correfficiency,
		const Float spillefficiency,
		const Float tau,
		const Float trx,
		const Float tatmos, 
		const Float tcmb);
		// const Quantity& trx,
		// const Quantity& tatmos, 
		// const Quantity& tcmb);

  // calculate errors and apply them to the data in our MS
  Bool corrupt();

  /*
  // Set autocorrelation weight
  Bool setauto(const Float autocorrwt);
  */

  // add new visibilities as described by the set methods to an existing
  // or just created measurement set
  Bool observe(const String& sourcename, const String& spwname,
	       const Record& startTime, 
	       const Record& stopTime);

  
  // Given a model image, predict the visibilities onto the (u,v) coordinates
  // of our MS
  Bool predict(const Vector<String>& modelImage, 
	       const String& compList, 
	       const Bool incremental);

  /*
  // Set the processing options
  Bool setoptions(const String& ftmachine, const Int cache, const Int tile,
		  const String& gridfunction, const MPosition& mLocation,
		  const Float padding, const Int facets,
		  const Double maxData,const Int wprojPlanes);

  Bool setauto(const Double autocorrwt);
  */

private:
  // Set up some defaults
  void defaults();

  // Make a VisSet if needed
  void makeVisSet();

  // Print out some help about create()
  // Format direction nicely
  String formatDirection(const MDirection& direction);

  // Format time nicely
  String formatTime(const Double time);
  
  // SkyEquation management
  // <group>
  Bool createSkyEquation( const Vector<String>& image, const String complist);
  void destroySkyEquation();
  // </group>

  String msname_p;
  MeasurementSet* ms_p;
  MeasurementSet* mssel_p;
  VisSet* vs_p;

  Int seed_p;

  VisJones *gj_p;
  VisJones *pj_p;
  VisJones *dj_p;
  VisJones *bj_p;
  ACoh     *ac_p;

  SkyEquation* se_p;
  CleanImageSkyModel* sm_p;
  FTMachine *ft_p;
  ComponentFTMachine *cft_p;

  Int nmodels_p;
  PtrBlock<PagedImage<Float>* > images_p;
  ComponentList *componentList_p;

  String ftmachine_p, gridfunction_p;
  Int cache_p, tile_p;
  MPosition mLocation_p;
  Float padding_p;
  Bool MSMayBeOK;
  Int facets_p;
  Int wprojPlanes_p;
  Long maxData_p;

  // Info for coordinates and station locations
  // <group>
  Bool           areStationCoordsSet_p;
  String         telescope_p;
  Vector<Double> x_p;
  Vector<Double> y_p;
  Vector<Double> z_p;
  Vector<Double>  diam_p;
  Vector<Double>  offset_p;
  Vector<String> mount_p;
  Vector<String> antName_p;
  String         coordsystem_p;
  MPosition      mRefLocation_p;
  // </group>

  // Info for observed field parameters
  // <group>

  String 	sourceName_p, calCode_p;
  MDirection	sourceDirection_p;
  Quantity      distance_p;

  // </group>


  // Info for spectral window parameters
  // <group>

  // Spectral windows data
  // <group>
  String 	spWindowName_p; 
  Int		nChan_p;
  Quantity     	startFreq_p;
  Quantity     	freqInc_p;
  Quantity     	freqRes_p;
  String     	stokesString_p;   
  // </group>
  // </group>


  // Feed information (there will be much more coming,
  // but we are brain dead at this moment).
  String feedMode_p;
  Int nFeeds_p;
  Bool feedsHaveBeenSet;
  Bool feedsInitialized;

  // Some times which are required for settimes
  // <group>
  Quantity integrationTime_p;
  Bool     useHourAngle_p;
  MEpoch   refTime_p;
  Bool timesHaveBeenSet_p;
  // </group>

  // Some parameters for voltage pattern (vp):
  // <group>
  Bool doVP_p;			// Do we apply VP or not?
  Bool doDefaultVP_p;		// Do we use the default VP for this telescope?
  String vpTableStr_p;		// Otherwise, use the VP specified in this Table
  Quantity  parAngleInc_p;	// Parallactic Angle increment
  Quantity  skyPosThreshold_p;  // a tolerance in the pointing center position
  Float pbLimit_p;              // The PB level (in percentage) after which the PB is assumed to be zero
  BeamSquint::SquintType  squintType_p;	// Control of squint to use
  VPSkyJones* vp_p;		// pointer to VPSkyJones for the sky equation
  VPSkyJones* gvp_p;		// pointer to VPSkyJones for the sky equation
  // </group>

  // Saving some information about the various corrupting terms
  // <group>
  String noisemode_p;
  // </group>

  // Cache the SimulatorProxy
  NewMSSimulator* sim_p;

  // The Jones matrix to hold the antenna pointing offsets and the
  // associated table name.  if applyPointingOffsets is False, only
  // VLA polarization squint will be included in EPJones.  If
  // doPBCorrection is True, the model image will be divided by the
  // primary beam before being used to predict the visibilities.
  // <group>
  EPJones *epJ_p;
  String epJTableName_p;
  Bool applyPointingOffsets_p;
  Bool doPBCorrection_p;
  // </group>
};

} //# NAMESPACE CASA - END

#endif
