//# Atmosphere.cc: Implementation of Atmosphere.h
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


#include <iostream>
#include <vector>

#include <synthesis/MeasurementComponents/Atmosphere.h>
#include <casa/Logging/LogIO.h>
#include <casa/OS/File.h>
#include <casa/BasicSL/String.h>

//------------------------------------------------------------------------
#define atm_telluric atm_telluric__
#define ini_telluric ini_telluric__
#define spe_telluric spe_telluric__

#if defined(AIPS_USEATM)
extern "C" {
  void atm_telluric(double *, double *, double *,
		    double *, double *, double *, 
		    double *, double *, double *, 
		    int *,          // atmtype
		    double *, int *, 
		    double *, double *, double *);
  void ini_telluric(int *,    int *,    int *,    int*,
		    double *, double *, double *, double *, 
		    double *, double *, double *, int*,
		    double *, double *, double *,
		    int *,    double *, int *,
		    double *, double *, double *, double *, 
		    double *, double *, double *, double *, 
		    double *, double *, double *, double *,
		    int *);
  void spe_telluric(int *,    int  *,   int *,
		    double *, double *,
		    double *, double *, double *,
		    double *, double *, double *, double *,
		    int *,    int  *,    //npoints_spec, nbd
		    double *,            // xspec
		    double *, double *,  // skycoupling,tspillover
		    double *            // tspec
		    //int *                // opt
		    );

};
#endif
//------------------------------------------------------------------------
namespace casa {
Atmosphere::Atmosphere(const double altitude, 
		       const double temperature,
		       const double pressure,
		       const double maxAltitude,
		       const double humidity,
		       const double dTempdH,
		       const double dP,
		       const double dPm,
		       const double h0,
		       const AtmType atmtype) {

  /// Copy parameters in private data area
  alti_m = altitude/1000.;   //    unit ATM_telluric: km
  t0_m = temperature;        //    unit ATM_telluric: K
  p0_m = pressure/100.;      //    unit ATM_telluric: mb (100 Pascal equals 1 millibars)
  top_m= maxAltitude/1000.;  //    unit ATM_telluric: km
  dt_m = dTempdH * 1000.;    //    unit ATM_telluric: K/km
  rh_m = humidity; 
  dp_m = dP/100.;
  dp1_m= dPm;
  h0_m = h0/1000.;           //    unit ATM_telluric: km
  atmtype_m = atmtype;
  int type = atmtype;
  skyCoupling_m = 1.;

  /// Allocate large arrays
  kmain_p      = new double [4*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS];
  kmainDer_p   = new double [4*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS];
  kminor_p     = new double [NMAX_MINOR*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS];
  kminorDer_p  = new double [NMAX_MINOR*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS];


#if defined(AIPS_USEATM)
  /// call ATM_telluric()
  atm_telluric(&t0_m,     //    temperature (k)
	       &dt_m,     //    delta t (k/km)
	       &p0_m,     //    pressure (mb)
	       &rh_m,     //    relative humidity (%)
	       &h0_m,     //    scale height (km)
	       &dp_m,     //    pressure step (mb)
	       &dp1_m,    //    multi. factor for pressure step 
	       &alti_m,   //    altitude (km)
	       &top_m,    //    top of atmophere (km)
	       &type,     //    atmosphere type
	       &water0_m, //    first guess water
	       &npp_m,    //    nb of layers in atm model
	       &atmProfile1_m[0][0],
	       &atmProfile2_m[0][0],
	       &minorComp_m[0][0]);


// Atmosphere::Atmosphere(Atmosphere::MODEL atmosmodel){

#endif


}
//------------------------------------------------------------------------
Atmosphere::~Atmosphere(){
  delete [] kmain_p;
  delete [] kmainDer_p;
  delete [] kminor_p;
  delete [] kminorDer_p;


  kmain_p     =0;
  kmainDer_p  =0;
  kminor_p    =0;
  kminorDer_p =0;

}
//-------------------------------------------------------------------
double Atmosphere::getStartupWaterContent() const{
  return water0_m/1000;
}
//-------------------------------------------------------------------
Profile Atmosphere::getProfile() const {
  Profile p;
  getProfile(p);
  return p;
}
//-------------------------------------------------------------------
void   Atmosphere::getProfile(Profile &p) const {

  //  p.thickness.resize(npp_m);
  //  p.thickness = valarray<double> ((double * const)atmProfile1_m, npp_m);
  p.thickness_m.clear();
  p.thickness_m.insert(p.thickness_m.begin(),
  		     (double *)atmProfile1_m,(double *)atmProfile1_m+npp_m);

  p.temperature_m.clear();
  p.temperature_m.insert(p.temperature_m.begin(),
  		       (double *)(atmProfile1_m+1),(double *)(atmProfile1_m+1)+npp_m);

  p.water_m.clear();
  p.water_m.insert(p.water_m.begin(),
  		 (double *)atmProfile2_m,(double *)atmProfile2_m+npp_m);

  p.pressure_m.clear();
  p.pressure_m.insert(p.pressure_m.begin(),
		(double *)(atmProfile2_m+1),(double *)(atmProfile2_m+1)+npp_m);

  p.O3_m.clear();
  p.O3_m.insert(p.O3_m.begin(),
  	      (double *)minorComp_m,(double *)minorComp_m+npp_m);

  p.CO_m.clear();
  p.CO_m.insert(p.CO_m.begin(),
  	      (double *)(minorComp_m+1),(double *)(minorComp_m+1)+npp_m); 

  p.N2O_m.resize(npp_m);
  p.N2O_m.insert(p.N2O_m.begin(),
 	       (double *)(minorComp_m+2),(double *)(minorComp_m+2)+npp_m);
  return;
  
}
//-------------------------------------------------------------------

void Atmosphere::initWindow(const int nbands,
			    const double fCenter[],
			    const double fWidth[],
			    const double fRes[]) {

  /// Prepare parameters for INI_telluric()
  int ier       = 0;
  nbands_m      = nbands;
  nwvrbands_m   = 0;
  nastrobands_m = nbands;
  nextrabands_m = 0;
  for (int i=0;i<nbands;i++) {
    iniband_m[i] = 0;
    //    freq0_m[i]   = (fCenter[i] - fWidth[i]/2.)/1.e9;
    freq0_m[i]   = fCenter[i]/1.e9; //    in atm freq0 is at the center of band
    sdb_m[i]     = 2;               //     upper side band
    vif_m[i]     = 0;
    res_m[i]     = fRes[i]/1.e9;
    width_m[i]   = fWidth[i]/1.e9;
    sidebandGain_m[i] = 1.0;
    nData_m[i]     = int(fWidth[i]/fRes[i]+0.5);
    //    printf("%d nData_m=%d res_m=%f width=%f fWidth=%g fCenter=%g freq0=%g\n",i, nData_m[i],res_m[i],width_m[i],fWidth[i],fCenter[i],freq0_m[i]);
  }
 
  //  cout<<"nwvrbands="<<nwvrbands_m<<" water0="<<water0_m<<endl;
  //  cout<<"Atmtransm water airmass freq " <<*water<<" "<<*airmass<<" "<<*freq<<endl;

  /// Call INI_telluric()
#if defined(AIPS_USEATM)
  ini_telluric( &nwvrbands_m,                  // IN Nb of radiometric channels
		&nastrobands_m,                // IN Nb of interferometer bands
		&nextrabands_m,                // IN Nb of other bands   
		&iniband_m[0],                 // IN                     
		&freq0_m[0],                   // IN   
		&vif_m[0],                     // IN   
		&res_m[0],                     // IN   
		&width_m[0],                   // IN   
		&sdb_m[0],                     // IN   
		&sidebandGain_m[0],            // IN   
		&water0_m,                     // IN   
		&npp_m,                        // IN   
		&atmProfile1_m[0][0],          // IN   
		&atmProfile2_m[0][0],          // IN   
		&minorComp_m[0][0],            // IN   
		&nbd_m[0],                     // OUT
		&xspec_m[0][0][0],             // OUT
		&npointsSpec_m[0],             // OUT
		//&kmain[0][0][0][0][0],       // OUT
		//&kmain_der[0][0][0][0][0],   // OUT
		//&kminor[0][0][0][0][0],      // OUT
		//&kminor_der[0][0][0][0][0],  // OUT
		kmain_p,                       // OUT
		kmainDer_p,                    // OUT
		kminor_p,                      // OUT
		kminorDer_p,                   // OUT
		&dlnFactorSpec_m[0][0][0],     // OUT
		&dldFactorSpec_m[0][0][0],     // OUT
		&dlnFactor_m[0],               // OUT  
		&dldFactor_m[0],               // OUT  
		&dryOpSpec_m[0][0][0],         // OUT  
		&wetOpSpec_m[0][0][0],         // OUT
		&dryOp_m[0],                   // OUT
		&wetOp_m[0],                   // OUT
		&ier);                         // OUT
#endif
}

//-------------------------------------------------------------------

int Atmosphere::getNdata(const int iband) const {
  return nData_m[iband];
}


//-------------------------------------------------------------------
Opacity Atmosphere::getOpacity() const {
  Opacity opa;
  getOpacity(opa);
  return opa;
}  
//-------------------------------------------------------------------
void Atmosphere::getOpacity(Opacity &o) const {
  o.dryOpacity_m.clear();
  o.dryOpacity_m.insert(o.dryOpacity_m.begin(),
		      (double *)dryOp_m,(double *)dryOp_m+nbands_m);
  o.wetOpacity_m.clear();
  o.wetOpacity_m.insert(o.wetOpacity_m.begin(),
		      (double *)wetOp_m,(double *)wetOp_m+nbands_m); 
  return;
}
//-------------------------------------------------------------------


OpacitySpec Atmosphere::getOpacitySpec() const {
  std::vector< std::vector<double> > dryOpa(nbands_m);
  std::vector< std::vector<double> > wetOpa(nbands_m);
  
  for (int iband=0;iband<nbands_m;iband++) {
    dryOpa[iband].resize(nData_m[iband]);
    wetOpa[iband].resize(nData_m[iband]);
    for (int idat=0; idat<nData_m[iband];idat++) {
      dryOpa[iband][idat] = dryOpSpec_m[0][idat][iband];
      wetOpa[iband][idat] = wetOpSpec_m[0][idat][iband];
      // PATCH DB: 29-09-2004 
      // take value from 1st channel  because wetOpaSpec is nul 
      wetOpa[iband][idat] = wetOpSpec_m[0][0][iband];
      //=================================================================
    }
  }
  OpacitySpec opa;
  opa.dryOpacitySpec_m = dryOpa;
  opa.wetOpacitySpec_m = wetOpa;
  return opa;
}  

// Call by reference method added by R. Rusk
void Atmosphere::getOpacitySpec(OpacitySpec &opa) const {
  std::vector< std::vector<double> > dryOpa(nbands_m);
  std::vector< std::vector<double> > wetOpa(nbands_m);
  
  for (int iband=0;iband<nbands_m;iband++) {
    dryOpa[iband].resize(nData_m[iband]);
    wetOpa[iband].resize(nData_m[iband]);
    for (int idat=0; idat<nData_m[iband];idat++) {
      dryOpa[iband][idat] = dryOpSpec_m[0][idat][iband];
      wetOpa[iband][idat] = wetOpSpec_m[0][idat][iband];
      // PATCH DB: 29-09-2004 
      // take value from 1st channel  because wetOpaSpec is nul 
      wetOpa[iband][idat] = wetOpSpec_m[0][0][iband];
      //=================================================================
    }
  }
  opa.dryOpacitySpec_m = dryOpa;
  opa.wetOpacitySpec_m = wetOpa;

}  
//-------------------------------------------------------------------

AbsCoeff Atmosphere::getAbsCoeff() const {
  AbsCoeff absCoeff;
  getAbsCoeff(absCoeff);
  return   absCoeff;
}
//-------------------------------------------------------------------
void  Atmosphere::getAbsCoeff(AbsCoeff &a) const {

  a.kH2OLines_m.resize(npp_m*NMAX_DATOS*nbands_m);
  a.kH2OCont_m.resize (npp_m*NMAX_DATOS*nbands_m);
  a.kO2_m.resize      (npp_m*NMAX_DATOS*nbands_m);
  a.kDryCont_m.resize (npp_m*NMAX_DATOS*nbands_m);
  a.kO3_m.resize      (npp_m*NMAX_DATOS*nbands_m);
  a.kCO_m.resize      (npp_m*NMAX_DATOS*nbands_m);
  a.kN2O_m.resize     (npp_m*NMAX_DATOS*nbands_m);
 
  int h2oLinesStart  = 0;
  int h2oContStart   = NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;
  int o2Start        = 2*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;
  int dryContStart   = 3*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;
  int o3Start        = 0;
  int coStart        = NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;;
  int n2oStart       = 2*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;

  int k =0;
  for (int ilayer=0; ilayer<npp_m; ilayer++) {
    int i1=NMAX_DATOS*NMAX_BANDS; 
    for (int idat=0;   idat<NMAX_DATOS; idat++) {
      int i2 = i1+idat*NMAX_BANDS;
      for (int iband=0; iband<nbands_m; iband++) {
	a.kH2OLines_m[k] = kmain_p[h2oLinesStart+ iband+i2];
	a.kH2OCont_m[k]  = kmain_p[h2oContStart + iband+i2];
	a.kO2_m[k]       = kmain_p[o2Start      + iband+i2];
	a.kDryCont_m[k]  = kmain_p[dryContStart + iband+i2];
	a.kO3_m[k]       = kminor_p[o3Start     + iband+i2];
	a.kCO_m[k]       = kminor_p[coStart     + iband+i2];
	a.kN2O_m[k]      = kminor_p[n2oStart    + iband+i2];
	k+=1;
      }
    }
  }
}
//-------------------------------------------------------------------
AbsCoeffDer Atmosphere::getAbsCoeffDer() const {
  AbsCoeffDer absCoeffDer;
  getAbsCoeffDer(absCoeffDer);
  return   absCoeffDer;
}
//-------------------------------------------------------------------
void  Atmosphere::getAbsCoeffDer(AbsCoeffDer &a) const {

  a.kH2OLinesDer_m.resize(npp_m*NMAX_DATOS*nbands_m);
  a.kH2OContDer_m.resize (npp_m*NMAX_DATOS*nbands_m);
  a.kO2Der_m.resize      (npp_m*NMAX_DATOS*nbands_m);
  a.kDryContDer_m.resize (npp_m*NMAX_DATOS*nbands_m);
  a.kO3Der_m.resize      (npp_m*NMAX_DATOS*nbands_m);
  a.kCODer_m.resize      (npp_m*NMAX_DATOS*nbands_m);
  a.kN2ODer_m.resize     (npp_m*NMAX_DATOS*nbands_m);
 
  int h2oLinesStart  = 0;
  int h2oContStart   = NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;
  int o2Start        = 2*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;
  int dryContStart   = 3*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;
  int o3Start        = 0;
  int coStart        = NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;;
  int n2oStart       = 2*NMAX_BD*NMAX_ATMLAYERS*NMAX_DATOS*NMAX_BANDS;

  int k =0;
  for (int ilayer=0; ilayer<npp_m; ilayer++) {
    int i1=NMAX_DATOS*NMAX_BANDS; 
    for (int idat=0;   idat<NMAX_DATOS; idat++) {
      int i2 = i1+idat*NMAX_BANDS;
      for (int iband=0; iband<nbands_m; iband++) {
	a.kH2OLinesDer_m[k] = kmainDer_p[h2oLinesStart+ iband+i2];
	a.kH2OContDer_m[k]  = kmainDer_p[h2oContStart + iband+i2];
	a.kO2Der_m[k]       = kmainDer_p[o2Start      + iband+i2];
	a.kDryContDer_m[k]  = kmainDer_p[dryContStart + iband+i2];
	a.kO3Der_m[k]       = kminorDer_p[o3Start     + iband+i2];
	a.kCODer_m[k]       = kminorDer_p[coStart     + iband+i2];
	a.kN2ODer_m[k]      = kminorDer_p[n2oStart    + iband+i2];
	k+=1;
      }
    }
  }
}


//-------------------------------------------------------------------
PhaseFactor Atmosphere::getPhaseFactor() const {
  PhaseFactor phaseFactor;
  getPhaseFactor(phaseFactor);
  return phaseFactor;
}  
//-------------------------------------------------------------------
void Atmosphere::getPhaseFactor(PhaseFactor &f) const {
  f.dispPhase_m.clear();
  f.dispPhase_m.insert(f.dispPhase_m.begin(),
		      (double *)dlnFactor_m,(double *)dlnFactor_m+nbands_m);
  f.nonDispPhase_m.clear();
  f.nonDispPhase_m.insert(f.nonDispPhase_m.begin(),
		      (double *)dldFactor_m,(double *)dldFactor_m+nbands_m); 
  return;
}
//-------------------------------------------------------------------



void Atmosphere::computeSkyBrightness(const double airMass,
			  const double tbgr, const double precWater) {

  /// Initialize some parameters for SPE_telluric()
  notherchan_m  = 0;
  nextrabands_m = nastrobands_m;
  rapport_m     = (precWater*1000)/water0_m;
  deltaT_m      = 0.; 
  airmass_m     = airMass;
  tbgr_m        = tbgr;
  //    opt_m         = 2;                   // Rayleigh Jeans temperature
  spillOverTemp_m   = 275.;
  
  /// Call SPE_telluric()
#if defined(AIPS_USEATM)
  spe_telluric( &notherchan_m,             // IN   
		&nextrabands_m,            // IN   
		&npp_m,                    // IN   
                &atmProfile1_m[0][0],      // IN   
                &rapport_m,                //  IN   
                &deltaT_m,                 // IN   
                &airmass_m,                // IN   
                &tbgr_m,                   // IN   
		//   &kmain[0][0][0][0][0],          // IN   
		//   &kmain_der[0][0][0][0][0],      // IN   
		//   &kminor[0][0][0][0][0],         // IN   
		//   &kminor_der[0][0][0][0][0],     // IN   
		kmain_p,                    // IN   
                kmainDer_p,                 // IN   
                kminor_p,                   // IN   
                kminorDer_p,                // IN   
		&npointsSpec_m[0],          //  IN   
                &nbd_m[0],                  // IN   
                &xspec_m[0][0][0],          // IN   
		&skyCoupling_m,               // IN
		&spillOverTemp_m,                // IN
		////                &tSpec_m[0][0][0],       // OUT 
		&tSpec_m[0][0][0][0]       // OUT  
		//&opt_m                    // IN  
		);
#endif

}


//-------------------------------------------------------------------
std::vector<double> Atmosphere::getSkyBrightness(const TemperatureType iopt) {

  /// Compute mean on bands
  std::vector<double> tempBand(nbands_m);

  for (int iband=0;iband<nbands_m;iband++) {
    tempBand[iband] = 0;
    double weight = 0;
    for (int idat=0; idat<nData_m[iband];idat++) {
      tempBand[iband] += tSpec_m[iopt-1][0][idat][iband];
      weight += 1;
    }
    if (weight>0) { tempBand[iband] /= weight; }
  }
  return tempBand;
}

//-------------------------------------------------------------------
std::vector< std::vector<double> >  Atmosphere::getSkyBrightnessSpec(const TemperatureType iopt) {

  std::vector< std::vector<double> > tspec(nbands_m);
  
  for (int iband=0;iband<nbands_m;iband++) {
    tspec[iband].resize(nData_m[iband]);
    for (int idat=0; idat<nData_m[iband];idat++) {
      tspec[iband][idat] = tSpec_m[iopt-1][0][idat][iband];
    }
  }
  return tspec;
}

} //# NAMESPACE CASA - END

