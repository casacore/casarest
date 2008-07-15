/*****************************************************************************
  WARNING: In case this include file has to be modified, the atmlib.inc file 
           has to be modified accordingly!
*****************************************************************************/
//#include "arrayconfig.h"
#if !defined(ATMLIB_H) || !defined(NMAX_DATPOS)
#define ATMLIB_H

#define NMAX_DATOS        128        //
#define NMAX_BD           2          // max nb of sidebands (2: signal + image)
#define NMAX_ATMLAYERS    40         // max nb of layers in atmosph. model
#define NMAX_MINOR        3          // Max number of minor components
#define NMAX_ASTROBAND    3          // max nb of astro. bands
#define NMAX_NWVRBANDS    5          // max nb of radiometric bands
#define NMAX_NEXTRABANDS  (NMAX_ASTROBAND+NMAX_NWVRBANDS)
#define NMAX_BANDS        (NMAX_NEXTRABANDS+NMAX_ASTROBAND+NMAX_NWVRBANDS)   
#define NUM_MAX_LINES_CATA 63        // max nb of lines de raies for minor  composants mineurs
#define NMAX_SUBCATALOG    100       // Max. number of subcatalogs for minor 
                                     // components lines
#define NMAX_PAR           1   
#endif
