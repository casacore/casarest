#ifndef MSVIS_ROWNR_H
#define MSVIS_ROWNR_H

#if CASACORE_MAJOR_VERSION<3 || (CASACORE_MAJOR_VERSION==3 && CASACORE_MINOR_VERSION<4)
  typedef unsigned int rownr_t;
#else
  typedef casacore::rownr_t rownr_t;
#endif

#endif
