#include <stdio.h>
    
#include <ionosphere/Ionosphere/RINEX.h>
#include <ionosphere/Ionosphere/GPSEphemeris.h>
#include <ionosphere/Ionosphere/GPSDCB.h>

#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/casa/Exceptions/Error.h>

#include <casacore/casa/namespace.h>

int main (void)
{
  RINEX::debug_level=2;
  RINEXSat::debug_level=2;
  GPSEphemeris::debug_level=2;
  
  try {
// read group delay data
    GPSDCB tgd("jpl_tgd.dta");
// read ephemeris
//    GPSEphemeris eph("test.orb");                                               
// read RINEX
    RINEX rinex("test.rnx");

// compute TECs
    Vector<Double> mjd,tec,stec,stec30;
    Vector<Int> sat,domain;
    rinex.getTEC(mjd,sat,tec,stec,stec30,domain,tgd);
    
// print 'em
    cerr<<"======================================== Computed TECs\n";
    for( uInt i=0; i<mjd.nelements(); i++ ) 
       cerr<<MVTime(mjd(i))<<" "<<sat(i)<<" "<<tec(i)<<" "<<stec(i)
          <<" "<<stec30(i)<<endl;
  }
  
  catch( AipsError x )
  {
    cerr << "AipsError: " << x.getMesg() << endl;
    return 1;
  }
  return 0;
}
