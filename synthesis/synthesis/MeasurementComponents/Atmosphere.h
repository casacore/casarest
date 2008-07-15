//# Atmosphere.h: Atmosphere model object for mm radio astron; 
//# Copyright (C) 1996,1997,1998,1999,2000,2001
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
#ifndef SYNTHESIS_ATMOSPHERE_H
#define SYNTHESIS_ATMOSPHERE_H

//#include <casa/Quanta/Quantum.h>

#include <vector>

namespace casa { //# NAMESPACE CASA - BEGIN

//<summary>
// Basic class for atm calc wrapping the fortran code of Pardo's ATM model 
//</summary>


#include <synthesis/MeasurementComponents/Atmlib.h>

/// Atmospheric model type
/*! \enum AtmType
 */
enum AtmType         { TROPICAL=1, 
		       MIDLAT_SUMMER=2,    
		       MIDLAT_WINTER=3, 
		       SUBARCTIC_SUMMER=4, 
		       SUBARCTIC_WINTER=5};

/// Temperature type
/*! \enum TemperatureType
 */
enum TemperatureType { BLACKBODY=1, 
		       RAYLEIGH_JEANS=2};


/// Atmospheric profile of layers [1:natmlayers]
/*! \struct Profile
 - thickness[natmlayers]   - in [m]
 - temperature[natmlayers] - in [K]
 - water[natmlayers]       - in [K/m^3]
 - pressure[natmlayers]    - in [Pa]
 - O3[natmlayers]          - in [m^-3]
 - CO[natmlayers]          - in [m^-3]
 - N2O[natmlayers]         - in [m^-3]
*/
struct Profile {
  std::vector<double> thickness_m;    ///< thickness[natmlayers]
  std::vector<double> temperature_m;  ///< temperature[natmlayers]
  std::vector<double> water_m;        ///< water[natmlayers]
  std::vector<double> pressure_m;     ///< pressure[natmlayers]
  std::vector<double> O3_m;           ///< O3[natmlayers]
  std::vector<double> CO_m;           ///< CO[natmlayers]
  std::vector<double> N2O_m;          ///< N2O[natmlayers]
};

/// Absorption coefficients.
/*! \struct AbsCoeff
 * Absoprtion coefficients for each band, for each channel of band,
 * for each atmospheric layer.  Units are [m^-1].
 */
struct AbsCoeff {
  std::vector<double>  kH2OLines_m;  ///< kH2OLines[natmlayers][ndata][nbands]
  std::vector<double>  kH2OCont_m;   ///< kH2OCont[natmlayers][ndata][nbands]
  std::vector<double>  kO2_m;        ///< kO2[natmlayers][ndata][nbands]
  std::vector<double>  kDryCont_m;   ///< kDryCont[natmlayers][ndata][nbands]
  std::vector<double>  kO3_m;        ///< kO3[natmlayers][ndata][nbands]
  std::vector<double>  kCO_m;        ///< kCO[natmlayers][ndata][nbands]
  std::vector<double>  kN2O_m;       ///< kN2O[natmlayers][ndata][nbands]
};

/// Derivative of absorption coefficients.
/*! \struct AbsCoeffDer
 * Derivative of absoprtion coefficients for each band,
 * for each channel of band, for each atmospheric layer.
 * Units are [m^-1].
 */
struct AbsCoeffDer {
  std::vector<double>  kH2OLinesDer_m; ///< kH2OLinesDer[natmlayers][ndata][nbands]
  std::vector<double>  kH2OContDer_m;  ///< kH2OContDer[natmlayers][ndata][nbands]
  std::vector<double>  kO2Der_m;       ///< kO2Der[natmlayers][ndata][nbands]
  std::vector<double>  kDryContDer_m;  ///< kDryContDer[natmlayers][ndata][nbands]
  std::vector<double>  kO3Der_m;       ///< kO3Der[natmlayers][ndata][nbands]
  std::vector<double>  kCODer_m;       ///< kCODer[natmlayers][ndata][nbands]
  std::vector<double>  kN2ODer_m;      ///< kN2ODer[natmlayers][ndata][nbands]
};

/// Opacity structure
/*! \struct Opacity
 - dryOpacity[nbands]
 - wetOpacity[nbands]   (per millimeter of precipitable water vapor)
*/
struct Opacity {
  std::vector<double>   dryOpacity_m;  ///< dry opacity
  std::vector<double>   wetOpacity_m;  ///< wet opacity in mm of PVW
};

/// Opacity structure for each frequency channel of each band
/*! \struct OpacitySpec
 - dryOpacitySpec[nbands][ndata]
 - wetOpacitySpec[nbands][ndata] (per mm of PWV)
*/
struct OpacitySpec {
  std::vector< std::vector<double> >  dryOpacitySpec_m; ///< dry opacity spec
  std::vector< std::vector<double> >  wetOpacitySpec_m; ///< wet opacity spec [mm^-1 H2O]
};

/// Phase-delay factor structure
/*! \struct PhaseFactor
-  dispPhase[nbands]  in [deg/m]
-  nonDispPhase[nbands] in [deg/m]
*/
struct PhaseFactor {
  std::vector<double>   dispPhase_m;      ///< dispersive phase delay in [deg/m]
  std::vector<double>   nonDispPhase_m;   ///< dispersive phase delay in [deg/m]
};

/// Water vapor results retrieved from radiometric measurements
/*! \struct WaterVaporFit
  - precWater;
  - sigmaFit;
*/
struct WaterVaporFit {
  double precWater_m;   ///< precipitable water
  double sigmaFit_m;    ///< quality of fit
};


/// C++ interface to FORTRAN ATM library of Juan R. Padro
/*! \class Atmosphere
 * The class Atmosphere is the C++ interface to the ATM library
 * used by ALMA TELCAL.
 */
class Atmosphere {

 public:



  //-------------------------------------------------------------------
  /**
   * Constructor:
   *  setup things in the atmosphere that are not supposed to change fast
   *  with time and water content.

   *  Builds the atmospheric profile with a guessed startup water content.
   *  \par
   *  Call fortran subroutine ATM_telluric()
   *
   *  @param  altitude    - at site in [m]
   *  @param  temperature - at site in [K]
   *  @param  pressure    - at site in [Pa]
   *  @param  maxAltitude - to top of modelled atmosphere in [m]
   *  @param  humidity    - percentage humidity used to guess water
   *  @param  dTem_dh     - change of T with height in [K/m]
   *  @param  dP          - initial pressure step (P[1]-P[0]) in [Pa]
   *  @param  dPm         - pressure multiplicative factor for steps : P[i+1]-P[i] = dPm * ( P[i]-P[i-1]) 
   *  @param  h0          - scale height for water (exp distribution) in [m]
   *  @param  atmtype     - enumerated types are:
   *   - 1: tropical
   *   - 2: mid latitude summer
   *   - 3: mid latitude winter
   *   - 4: subarctic summer
   *   - 5: subarctic winter
   *
   * @exception -              invalid parameters. 
   * @exception -              too many atmosphere layers  
   */
  //-------------------------------------------------------------------

  Atmosphere(const double altitude, 
	     const double temperature,
	     const double pressure,
	     const double maxAltitude,
	     const double humidity,
	     const double dTem_dh,
	     const double dP,
	     const double dPm,
	     const double h0,
	     const AtmType atmtype) ;


  //---------------------------------------------------------------------
  /// Destructor
  //---------------------------------------------------------------------

  ~Atmosphere();


  //---------------------------------------------------------------------
  // Methods
  //---------------------------------------------------------------------


  //-------------------------------------------------------------------
  // getStartupWaterContent()
  // ------------------------
  /**
   *  Get the guessed startup water content computed by initAtmosphere()
   * @return 
   *  precWater  - guessed precipitable water content in [m]
   */
  //-------------------------------------------------------------------
  double getStartupWaterContent() const;


  //-------------------------------------------------------------------
  // getProfiles()
  //-------------------------------------------------------------------
  /**
   * Get atmospheric profile
   * @return Profile - atmospheric profile structure
   */
  Profile getProfile() const;

  //-------------------------------------------------------------------
  /**
   * Get atmospheric profile
   * @param[out] p - Profile structure
   */
  void  getProfile(Profile & p) const;
  //-------------------------------------------------------------------
  

  //-------------------------------------------------------------------
  //  initWindow()
  //  -----------
  //
  /**
   *    Define a spectral window, compute absorption and emission
   *    coefficients for this window, using the above atmospheric
   *    parameters.
   *
   *    Call fortran subroutine INI_telluric()
   * @param[in]    nbands           - number of bands
   * @param[in]    fCenter[nbands]  - (sky) frequencies in [Hz]
   * @param[in]    fWidth[nbands]   - frequency widths in [Hz]
   * @param[in]    fRes[nbands]     - resolution inside band in [Hz]
   * @exception - invalid parameters
   */
  void initWindow(const int nbands,
		  const double fCenter[],
		  const double fWidth[],
		  const double fRes[]);


  //-------------------------------------------------------------------
  // getNdata()
  // -----------
  /**
   *  Return the number of channels of ith band
   *  @param[in]     iband - identifier of band 
   *  @return        ndata - number of channels
   *  @exception - invalid parameters
   */
  int getNdata(const int iband) const;


  //-------------------------------------------------------------------
  //  getOpacity()
  //  ------------
  /**  Get the integrated optical depth on each band
   *  @return Opacity structure
   */
  Opacity getOpacity() const;

  /**  Get the integrated optical depth on each band
   *  @param[out] opacity - opacity structure
   */
  void getOpacity(Opacity &opacity) const ;

  //-------------------------------------------------------------------
  //  getOpacitySpec()
  //  ----------------
  /**  Get the integrated optical depth for each channel of each band
   *  @return OpacitySpec structure
   */
  OpacitySpec getOpacitySpec() const;

  /**  Get the integrated optical depth for each channel of each band
   *  @param[out] o - OpacitySpec structure
   */
  void getOpacitySpec(OpacitySpec &o) const;


  //-------------------------------------------------------------------
  //  getAbsCoeff()
  //  -------------
  /**
   * Get absorption coefficients
   * @return AbsCoeff structure
   */
  AbsCoeff getAbsCoeff() const;

  //-------------------------------------------------------------------
  /**
   * Get absorption coefficients
   * @param[out] a - AbsCoeff structure
   */
  void  getAbsCoeff(AbsCoeff &a) const;

  //-------------------------------------------------------------------
  // getAbsCoeffDer
  // --------------
  /**
   * Get absorption coefficients
   * @return AbsCoeffDer - AbsCoeffDer structure
   */
  AbsCoeffDer getAbsCoeffDer() const;

  //-------------------------------------------------------------------
  /**
   * Get absorption coefficients
   * @param[out] a - AbsCoeff structure
   */
  void  getAbsCoeffDer(AbsCoeffDer & a ) const;

         
  //-------------------------------------------------------------------
  // getPhaseFactor()
  // ----------------
  //
  /**  
   *   Get dispersive and non-dispersive phase delay factor 
   *   @return PhaseFactor structure
   */
  PhaseFactor getPhaseFactor() const;

  /**   Get dispersive and non-dispersive phase delay factor 
   *   @param[out] phaseFactor - PhaseFactor structure
   */
  void getPhaseFactor(PhaseFactor &phaseFactor) const ;


  //-------------------------------------------------------------------
  // computeSkyBrightness()
  // ----------------------
  /**
   * Compute atmospheric brightness by integrating the transfer equation.
   * Call  fortran subroutine SPE_telluric()
   *
   * 
   * @param[in] airMass   - Air mass
   * @param[in] tbgr      - Temperature of cosmic background in [K]
   * @param[in] precWater - Precipitable water content in [m]
   *
   *  @exception -  no spectral band initialized
   *  @exception -  invalid parameters
   */
  void computeSkyBrightness(const double airMass,
			    const double tbgr,
			    const double precWater);



  //-------------------------------------------------------------------
  // getSkyBrightness()
  // ----------------------
  //
  /**
   *  Get sky brightness computed by method computeSkyBrightness()
   *
   * @param[in] iopt - TemperatureType type.  Enumerated values are:
   *  - 1:  return blackbody temperature
   *  - 2:  return Rayleigh Jeans temperature
   * @return
   *    Tspec[nbands]
   */
  std::vector<double> getSkyBrightness(const TemperatureType iopt);



  //-------------------------------------------------------------------
  // getSkyBrightnessSpec()
  // ----------------------
  /**
   *  Get sky brightness by integrating the transfer equation for each channel
   * @param[in] iopt -  TemperatureType type.  Enumerated values are one of:
   *    -     1:  return blackbody temperature
   *    -     2:  return Rayleigh Jeans temperature
   * @return
   *    double  Tspec[nbands][ndata]
   */
  std::vector< std::vector<double> >  getSkyBrightnessSpec(const TemperatureType iopt);


 
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// setSkyCoupling() - set the sky coupling
  void   setSkyCoupling(const float c) {skyCoupling_m = c ;}

  /// getSkyCoupling() - get the sky coupling
  float  getSkyCoupling()              {return skyCoupling_m;};


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 private:

  /// Inaccessible copy constructor and assignment operator
  Atmosphere(const Atmosphere& atm);
  Atmosphere& operator=(const Atmosphere& atm);
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  /// data

  int npp_m ,nbands_m;
  AtmType atmtype_m;

  double t0_m, dt_m, p0_m, rh_m, h0_m, dp_m, dp1_m, alti_m, top_m, water0_m;
  double atmProfile1_m[3][NMAX_ATMLAYERS];
  double atmProfile2_m[2][NMAX_ATMLAYERS]; 
  double minorComp_m[3][NMAX_ATMLAYERS];

  int    nwvrbands_m, nastrobands_m, nextrabands_m, iniband_m[NMAX_BANDS];
  int    nbd_m[NMAX_BANDS], npointsSpec_m[NMAX_BANDS];
  int    nData_m[NMAX_BANDS];
  double freq0_m[NMAX_BANDS], vif_m[NMAX_BANDS], res_m[NMAX_BANDS];
  double width_m[NMAX_BANDS], sdb_m[NMAX_BANDS];
  double sidebandGain_m[NMAX_BANDS];

  double xspec_m[NMAX_BD][NMAX_DATOS][NMAX_BANDS];

  double *kmain_p;
  double *kmainDer_p;
  double *kminor_p;
  double *kminorDer_p;

  //  double kmain[4][NMAX_BD][NMAX_ATMLAYERS][NMAX_DATOS][NMAX_BANDS];
  //  double kmain_der[4][NMAX_BD][NMAX_ATMLAYERS][NMAX_DATOS][NMAX_BANDS];
  //  double kminor[NMAX_MINOR][NMAX_BD][NMAX_ATMLAYERS][NMAX_DATOS][NMAX_BANDS];
  //  double kminor_der[NMAX_MINOR][NMAX_BD][NMAX_ATMLAYERS][NMAX_DATOS][NMAX_BANDS];

  double   dlnFactorSpec_m[NMAX_BD][NMAX_DATOS][NMAX_ASTROBAND];
  double   dldFactorSpec_m[NMAX_BD][NMAX_DATOS][NMAX_ASTROBAND];
  double   dlnFactor_m[NMAX_ASTROBAND];
  double   dldFactor_m[NMAX_ASTROBAND];
  double   dryOpSpec_m[NMAX_BD][NMAX_DATOS][NMAX_ASTROBAND];
  double   wetOpSpec_m[NMAX_BD][NMAX_ASTROBAND][NMAX_DATOS];  
  double   dryOp_m[NMAX_ASTROBAND];
  double   wetOp_m[NMAX_ASTROBAND];

  int    notherchan_m, opt_m;
  double rapport_m, deltaT_m, airmass_m, tbgr_m;
  double tSpec_m[3][NMAX_BD][NMAX_DATOS][NMAX_BANDS];


  double h2olines_m, h2ocont_m, o2lines_m, drycont_m, minorlines_m;

  double skyCoupling_m, spillOverTemp_m;


/*
!  Private data used by ATM_telluric()
! 
!               t0_m       Ambient temperature at the site (K)
!               dt_m       tropospheric lapse rate in K/km
!               p0_m       Ground pressure at the site (mb)
!               rh_m       relative humidity at the site (%)    
!                             (It is used only to make an estimate 
!                             of the water vapor column, first guess)
!               h0_m       scale height of water vapor distribution (km)
!               dp_m       Pressure basic step (mb)
!               dp1_m      multiplicative factor for presure steps.
!                                 Example of pressure parameters: 
!                                 P_m=550; DP_m: 10; DP1_m: 
!                                 1.2 ==> The pressure levels will 
!                                 then be 550, 560, 572, 586.4, ....
!               alti_m       Altitude of the site (km)
!               top_m        Top of atmospheric profile (km)
!               atmtype       Atmosphere type (1: Tropical, 2: Midlat. Summer,
!                             3: Midlat. Winter, 4: Subarctic Summer, 
!                             5: Subarctic Winter
!               water0_m      1st guess water vapor column
!               npp_m           Total number of layers in the output 
!                             atmospherice profiles
!               atmProfile1_m  (1:npp,1) -> thickness of layer (m)
!                             (1:npp,2) -> Temp. of layer (K)  
!                             (1:npp,3) -> delta_T layer (K)   
!               atmProfile2_m  (1:npp,1) -> water vapor kg/m**3  
!                             (1:npp,2) -> pressure in mb   
!               minorComp_m    (1,1:npp) --> O3  (2,1:npp) --> CO  
!                             (3,1:npp) --> N2O   [all in m**-3]
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Private data used by INI_telluric()
!
!               nwvrbands_m      Number of radiometric channels
!               nastrobands_m    Number of correlator (inteferometer) bands
!               nextrabands_m    Number of other BANDS defined by the user
!               iniband_m        (1:wvrbands+nastrobands+nextrabands)
!                                if 0 ---> the band corresponding to 
!                                         the index has to be initialized
!                                if =/=0 ---> the band corresponding to the
!                                   index does not need to be initialized
!               freq0_m          Local Oscillator (after multiplication) 
!                              frequency of BAND 
!                              (1:nwvrbands+nastrobands+nextrabands) GHz
!               vif_m            Intermediate Frequency of BAND 
!                              (1:nwvrbands+nastrobands+nextrabands) GHz
!               res_m            Freq. resolution inside BAND 
!                              (1:nwvrbands+nastrobands+nextrabands) GHz
!               width_m          Width of BAND 
!                              (1:nwvrbands+nastrobands+nextrabands) GHz
!               sdb_m            for each BAND, if sdb=1 ==> 1 is lower side 
!                                            band and 2 is upper side band
!                                            if sdb=1 ==> 1 is upper side 
!                                            band and 2 is lower side band
!               sidebandGain  For each BAND, verifies that the sum of the 
!                              gains of the two sidebands is 1.0.
!               wh2o           First gues water vapor column returned by 
!                              ATM_telluric (mm)
!               npp_m            Total number of layers in the output 
!                              atmospherice profiles
!               atmProfile1_m   (1:npp,1) -> thinkness of layer (m)  
!                              (1:npp,2) -> Temp. of layer (K)  
!                              (1:npp,3) -> delta_T layer (K)   
!               atmProfile2_m   (1:npp,1) -> water vapor kg/m**3  
!                              (1:npp,2) -> pression in mb   
!               minorComp_m     (1:npp,1) --> O3   [m**-3]  
!                              (1:npp,2) --> CO   [m**-3]
!                              (1:npp,3) --> N2O  [m**-3]
!               nbd            Number of sidebands in each band  (1:nmax_bands)   
!               xspec_M          individual frequencies considered in each 
!                              band and sideband (GHz)
!               npointsSpec   Number of individual frequencies considered 
!                              inside each band of each channel (1:nmax_bands)   
!               kmain_p          absorption coefficient of 1: H2O lines , 
!                              2: H2O cont , 3: O2 lines , 
!                              4: Dry continuum (last index)
!                              per channel (1st index),  individual frequency 
!                              inside channel (2nd index), atmospheric layer 
!                              (3rd index), sideband (4th index).
!               kmainDer_p      derivative of absorption coefficient respect to 
!                              ground temperature for 1: H2O lines , 
!                              2: H2O cont , 3: O2 lines ,4: Dry continuum 
!                              (last index) per channel (1st index), 
!                              individual frequency inside channel (2nd 
!                              index), atmospheric layer (3rd index), 
!                              sideband (4th index).
!               kminor_p         absorption coefficient of minor gases (5th 
!                              index): 1: O3, 2: CO, 3: N2O,... 
!                              per channel (1st index),  individual frequency 
!                              inside channel (2nd index), atmospheric layer 
!                              (3rd index), sideband (4th index).
!               kminorDer_p     derivative of absorption coefficient of minor 
!                              gases respect to ground temperature (5th 
!                              index): 1: O3, 2: CO, 3: N2O,... per channel 
!                              (1st index), individual frequency inside channel 
!                              (2nd index), atmospheric layer (3rd index), 
!                              sideband (4th index).
!             dlnFactorSpec_m    Non-dispersive phase delay factor (deg/mm)
!                              for each individual frequency and sideband 
!                              in correlator (inteferometer) band.
!             dldFactorSpec_m    Dispersive phase delay factor (deg/mm)
!                              for each individual frequency and sideband 
!                              in correlator (inteferometer) band.
!               dlnFactor_m     average non-dispersive phase in each inter-
!                              ferometer band (deg/mm) [defined: (1:nastrobands)]
!               dldFactor_m     average dispersive phase in each interferometer 
!                              band (deg/mm) [defined: (1:nastrobands)]
!               dryOpSpec_m    Dry opacity (np) for each individual frequency 
!                              and sideband in correlator (inteferometer) band.
!               wetOpSpec_m    "wet" opacity (np) for each individual frequency 
!                              and sideband in correlator (inteferometer) band.
!               dryOp_m         average dry opacity in each interferometer band 
!                              (np) [defined: (1:nastrobands)]
!               wetOp_m         average "wet" opacity in each interferometer 
!                              band (np/mm) [defined: (1:nastrobands)]
!
!
! Private Data associated to SPE_telluric()
!               rapport_m          Ratio between precipitable water vapor column 
!                                and water vapor column used when initializing
!                                absorption coefficiens.
!               deltaT_m        difference respect to average ground 
!                                temperaute.
!               airm_m             Air mass
!               tbgr_m             Temperature of cosmic background  
!               xSpec_m           individual frequencies considered in each 
!                               channel
!               skyCoupling_m     Sky coupling
!               tSpillover_m      Spill over Temperature
!               tSpec_m          TEBB (K) (last index = 1) of the sky for extrabands 
!                               T_RJ (K) (last index = 2) of the sky for extrabands 
!                               (1.0/(exp(h_div_k*v/T_EBB)-1.0)) (last index = 3) of 
!                                    the sky for extrabands 
!                               extrabands: notherchan+1 to notherchan+nextrabands 
!                               (notherchan+1:notherchan+nextrabands,
!                                1:npoints_spec(channel),1:nbd(channel))
*/

};

} //# NAMESPACE CASA - END

#endif
