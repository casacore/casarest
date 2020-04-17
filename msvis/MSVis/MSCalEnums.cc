//# MSCalEnums.cc: Implementation of MSCalEnums.h
//# Copyright (C) 1996,1997,1998,2000,2001,2002
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
//# $Id$
//----------------------------------------------------------------------------

#include <msvis/MSVis/MSCalEnums.h>

namespace casacore { //# NAMESPACE CASACORE - BEGIN

//----------------------------------------------------------------------------

// Static data initialization
std::map <Int, String> MSCalEnums::theirFieldMap;
std::map <Int, DataType> MSCalEnums::theirTypeMap;

//----------------------------------------------------------------------------

void MSCalEnums::initMaps ()
{
// Initialize the static map containing the field names.
// Skip this step if already initialized.
//
  if (theirFieldMap.empty()) {
    theirFieldMap.insert (std::make_pair(ANTENNA1, "ANTENNA1"));
    theirFieldMap.insert (std::make_pair(ANTENNA2, "ANTENNA2"));
    theirFieldMap.insert (std::make_pair(FEED1, "FEED1"));
    theirFieldMap.insert (std::make_pair(FEED2, "FEED2"));
    theirFieldMap.insert (std::make_pair(PULSAR_BIN, "PULSAR_BIN"));
    theirFieldMap.insert (std::make_pair(SCAN_NUMBER, "SCAN_NUMBER"));
    theirFieldMap.insert (std::make_pair(TIME, "TIME"));
    theirFieldMap.insert (std::make_pair(TIME_EXTRA_PREC, "TIME_EXTRA_PREC"));
    theirFieldMap.insert (std::make_pair(INTERVAL, "INTERVAL"));
    theirFieldMap.insert (std::make_pair(ARRAY_ID, "ARRAY_ID"));
    theirFieldMap.insert (std::make_pair(PROCESSOR_ID, "PROCESSOR_ID"));
    theirFieldMap.insert (std::make_pair(FIELD_ID, "FIELD_ID"));
    theirFieldMap.insert (std::make_pair(OBSERVATION_ID, "OBSERVATION_ID"));
    theirFieldMap.insert (std::make_pair(PULSAR_GATE_ID, "PULSAR_GATE_ID"));
    theirFieldMap.insert (std::make_pair(SPECTRAL_WINDOW_ID, "SPECTRAL_WINDOW_ID"));
    theirFieldMap.insert (std::make_pair(PHASE_ID, "PHASE_ID"));
    theirFieldMap.insert (std::make_pair(STATE_ID, "STATE_ID"));

    theirFieldMap.insert (std::make_pair(FREQ_GROUP, "FREQ_GROUP"));
    theirFieldMap.insert (std::make_pair(FREQ_GROUP_NAME, "FREQ_GROUP_NAME"));
    theirFieldMap.insert (std::make_pair(FIELD_NAME, "FIELD_NAME"));
    theirFieldMap.insert (std::make_pair(FIELD_CODE, "FIELD_CODE"));
    theirFieldMap.insert (std::make_pair(SOURCE_NAME, "SOURCE_NAME"));
    theirFieldMap.insert (std::make_pair(SOURCE_CODE, "SOURCE_CODE"));
    theirFieldMap.insert (std::make_pair(CALIBRATION_GROUP, "CALIBRATION_GROUP"));

    theirFieldMap.insert (std::make_pair(GAIN, "GAIN"));
    theirFieldMap.insert (std::make_pair(REF_ANT, "REF_ANT"));
    theirFieldMap.insert (std::make_pair(REF_FEED, "REF_FEED")); 
    theirFieldMap.insert (std::make_pair(REF_RECEPTOR, "REF_RECEPTOR"));
    theirFieldMap.insert (std::make_pair(REF_FREQUENCY, "REF_FREQUENCY"));
    theirFieldMap.insert (std::make_pair(MEAS_FREQ_REF, "MEAS_FREQ_REF"));
    theirFieldMap.insert (std::make_pair(REF_DIRECTION, "REF_DIRECTION"));
    theirFieldMap.insert (std::make_pair(MEAS_DIR_REF, "MEAS_DIR_REF"));
    theirFieldMap.insert (std::make_pair(POINTING_OFFSET, "POINTING_OFFSET"));
    theirFieldMap.insert (std::make_pair(MEAS_POINTING, "MEAS_POINTING"));
    theirFieldMap.insert (std::make_pair(CAL_DESC_ID, "CAL_DESC_ID"));
    theirFieldMap.insert (std::make_pair(CAL_HISTORY_ID, "CAL_HISTORY_ID"));
    
    theirFieldMap.insert (std::make_pair(TOTAL_SOLUTION_OK, "TOTAL_SOLUTION_OK"));
    theirFieldMap.insert (std::make_pair(TOTAL_FIT, "TOTAL_FIT"));
    theirFieldMap.insert (std::make_pair(TOTAL_FIT_WEIGHT, "TOTAL_FIT_WEIGHT"));
    theirFieldMap.insert (std::make_pair(SOLUTION_OK, "SOLUTION_OK"));
    theirFieldMap.insert (std::make_pair(FIT, "FIT"));
    theirFieldMap.insert (std::make_pair(FIT_WEIGHT, "FIT_WEIGHT"));
    theirFieldMap.insert (std::make_pair(FLAG, "FLAG"));
    theirFieldMap.insert (std::make_pair(SNR, "SNR"));
    
    theirFieldMap.insert (std::make_pair(NUM_SPW, "NUM_SPW"));
    theirFieldMap.insert (std::make_pair(NUM_CHAN, "NUM_CHAN"));
    theirFieldMap.insert (std::make_pair(NUM_RECEPTORS, "NUM_RECEPTORS"));
    theirFieldMap.insert (std::make_pair(N_JONES, "N_JONES"));
    theirFieldMap.insert (std::make_pair(CHAN_FREQ, "CHAN_FREQ"));
    theirFieldMap.insert (std::make_pair(CHAN_WIDTH, "CHAN_WIDTH")); 
    theirFieldMap.insert (std::make_pair(CHAN_RANGE, "CHAN_RANGE"));
    theirFieldMap.insert (std::make_pair(JONES_TYPE, "JONES_TYPE"));
    theirFieldMap.insert (std::make_pair(POLARIZATION_TYPE, "POLARIZATION_TYPE"));
    theirFieldMap.insert (std::make_pair(MS_NAME, "MS_NAME"));
    
    theirFieldMap.insert (std::make_pair(CAL_PARMS, "CAL_PARMS"));
    theirFieldMap.insert (std::make_pair(CAL_TABLES, "CAL_TABLES"));
    theirFieldMap.insert (std::make_pair(CAL_SELECT, "CAL_SELECT"));
    theirFieldMap.insert (std::make_pair(CAL_NOTES, "CAL_NOTES"));
    
    theirFieldMap.insert (std::make_pair(CAL_DESC, "CAL_DESC"));
    theirFieldMap.insert (std::make_pair(CAL_HISTORY, "CAL_HISTORY"));
    
    theirFieldMap.insert (std::make_pair(ROT_MEASURE, "ROT_MEASURE"));
    theirFieldMap.insert (std::make_pair(ROT_MEASURE_ERROR, "ROT_MEASURE_ERROR"));
    theirFieldMap.insert (std::make_pair(IONOSPH_TEC, "IONOSPH_TEC"));
    theirFieldMap.insert (std::make_pair(IONOSPH_TEC_ERROR, "IONOSPH_TEC_ERROR"));

    theirFieldMap.insert (std::make_pair(PHASE_OFFSET, "PHASE_OFFSET"));
    theirFieldMap.insert (std::make_pair(SB_DELAY, "SB_DELAY"));
    theirFieldMap.insert (std::make_pair(DELAY_RATE, "DELAY_RATE"));

    theirFieldMap.insert (std::make_pair(POLY_TYPE, "POLY_TYPE"));
    theirFieldMap.insert (std::make_pair(POLY_MODE, "POLY_MODE"));
    theirFieldMap.insert (std::make_pair(SCALE_FACTOR, "SCALE_FACTOR"));
    theirFieldMap.insert (std::make_pair(VALID_DOMAIN, "VALID_DOMAIN"));
    theirFieldMap.insert (std::make_pair(N_POLY_AMP, "N_POLY_AMP"));
    theirFieldMap.insert (std::make_pair(N_POLY_PHASE, "N_POLY_PHASE"));
    theirFieldMap.insert (std::make_pair(POLY_COEFF_AMP, "POLY_COEFF_AMP"));
    theirFieldMap.insert (std::make_pair(POLY_COEFF_PHASE, "POLY_COEFF_PHASE"));
    theirFieldMap.insert (std::make_pair(PHASE_UNITS, "PHASE_UNITS"));

    theirFieldMap.insert (std::make_pair(SIDEBAND_REF, "SIDEBAND_REF"));

    theirFieldMap.insert (std::make_pair(N_KNOTS_AMP, "N_KNOTS_AMP"));
    theirFieldMap.insert (std::make_pair(N_KNOTS_PHASE, "N_KNOTS_PHASE"));
    theirFieldMap.insert (std::make_pair(SPLINE_KNOTS_AMP, "SPLINE_KNOTS_AMP"));
    theirFieldMap.insert (std::make_pair(SPLINE_KNOTS_PHASE, "SPLINE_KNOTS_PHASE"));
  };

// Initialize the static map containing the basic field data types
// Skip this step if already initialized.
//
  if (theirTypeMap.empty()) {
    theirTypeMap.insert (std::make_pair(ANTENNA1, TpInt));
    theirTypeMap.insert (std::make_pair(ANTENNA2, TpInt));
    theirTypeMap.insert (std::make_pair(FEED1, TpInt));
    theirTypeMap.insert (std::make_pair(FEED2, TpInt));
    theirTypeMap.insert (std::make_pair(PULSAR_BIN, TpInt));
    theirTypeMap.insert (std::make_pair(SCAN_NUMBER, TpInt));
    theirTypeMap.insert (std::make_pair(TIME, TpDouble));
    theirTypeMap.insert (std::make_pair(TIME_EXTRA_PREC, TpDouble));
    theirTypeMap.insert (std::make_pair(INTERVAL, TpDouble));
    theirTypeMap.insert (std::make_pair(ARRAY_ID, TpInt));
    theirTypeMap.insert (std::make_pair(PROCESSOR_ID, TpInt));
    theirTypeMap.insert (std::make_pair(FIELD_ID, TpInt));
    theirTypeMap.insert (std::make_pair(OBSERVATION_ID, TpInt));
    theirTypeMap.insert (std::make_pair(PULSAR_GATE_ID, TpInt));
    theirTypeMap.insert (std::make_pair(SPECTRAL_WINDOW_ID, TpInt));
    theirTypeMap.insert (std::make_pair(PHASE_ID, TpInt));
    theirTypeMap.insert (std::make_pair(STATE_ID, TpInt));

    theirTypeMap.insert (std::make_pair(FREQ_GROUP, TpInt));
    theirTypeMap.insert (std::make_pair(FREQ_GROUP_NAME, TpString));
    theirTypeMap.insert (std::make_pair(FIELD_NAME, TpString));
    theirTypeMap.insert (std::make_pair(FIELD_CODE, TpString));
    theirTypeMap.insert (std::make_pair(SOURCE_NAME, TpString));
    theirTypeMap.insert (std::make_pair(SOURCE_CODE, TpString));
    theirTypeMap.insert (std::make_pair(CALIBRATION_GROUP, TpInt));

    theirTypeMap.insert (std::make_pair(GAIN, TpComplex));
    theirTypeMap.insert (std::make_pair(REF_ANT, TpInt));
    theirTypeMap.insert (std::make_pair(REF_FEED, TpInt)); 
    theirTypeMap.insert (std::make_pair(REF_RECEPTOR, TpInt));
    theirTypeMap.insert (std::make_pair(REF_FREQUENCY, TpDouble));
    theirTypeMap.insert (std::make_pair(MEAS_FREQ_REF, TpInt));
    theirTypeMap.insert (std::make_pair(REF_DIRECTION, TpDouble));
    theirTypeMap.insert (std::make_pair(MEAS_DIR_REF, TpInt));
    theirTypeMap.insert (std::make_pair(CAL_DESC_ID, TpInt));
    theirTypeMap.insert (std::make_pair(CAL_HISTORY_ID, TpInt));
    
    theirTypeMap.insert (std::make_pair(TOTAL_SOLUTION_OK, TpBool));
    theirTypeMap.insert (std::make_pair(TOTAL_FIT, TpFloat));
    theirTypeMap.insert (std::make_pair(TOTAL_FIT_WEIGHT, TpFloat));
    theirTypeMap.insert (std::make_pair(SOLUTION_OK, TpBool));
    theirTypeMap.insert (std::make_pair(FIT, TpFloat));
    theirTypeMap.insert (std::make_pair(FIT_WEIGHT, TpFloat));
    theirTypeMap.insert (std::make_pair(FLAG, TpBool));
    theirTypeMap.insert (std::make_pair(SNR, TpFloat));
    
    theirTypeMap.insert (std::make_pair(NUM_SPW, TpInt));
    theirTypeMap.insert (std::make_pair(NUM_CHAN, TpInt));
    theirTypeMap.insert (std::make_pair(NUM_RECEPTORS, TpInt));
    theirTypeMap.insert (std::make_pair(N_JONES, TpInt));
    theirTypeMap.insert (std::make_pair(CHAN_FREQ, TpDouble));
    theirTypeMap.insert (std::make_pair(CHAN_WIDTH, TpDouble)); 
    theirTypeMap.insert (std::make_pair(CHAN_RANGE, TpInt));
    theirTypeMap.insert (std::make_pair(JONES_TYPE, TpString));
    theirTypeMap.insert (std::make_pair(POLARIZATION_TYPE, TpString));
    theirTypeMap.insert (std::make_pair(MS_NAME, TpString));
    
    theirTypeMap.insert (std::make_pair(CAL_PARMS, TpString));
    theirTypeMap.insert (std::make_pair(CAL_TABLES, TpString));
    theirTypeMap.insert (std::make_pair(CAL_SELECT, TpString));
    theirTypeMap.insert (std::make_pair(CAL_NOTES, TpString));
    
    theirTypeMap.insert (std::make_pair(CAL_DESC, TpTable));
    theirTypeMap.insert (std::make_pair(CAL_HISTORY, TpTable));
    
    theirTypeMap.insert (std::make_pair(ROT_MEASURE, TpFloat));
    theirTypeMap.insert (std::make_pair(ROT_MEASURE_ERROR, TpFloat));
    theirTypeMap.insert (std::make_pair(IONOSPH_TEC, TpFloat));
    theirTypeMap.insert (std::make_pair(IONOSPH_TEC_ERROR, TpFloat));

    theirTypeMap.insert (std::make_pair(PHASE_OFFSET, TpFloat));
    theirTypeMap.insert (std::make_pair(SB_DELAY, TpFloat));
    theirTypeMap.insert (std::make_pair(DELAY_RATE, TpFloat));

    theirTypeMap.insert (std::make_pair(POLY_TYPE, TpString));
    theirTypeMap.insert (std::make_pair(POLY_MODE, TpString));
    theirTypeMap.insert (std::make_pair(SCALE_FACTOR, TpComplex));
    theirTypeMap.insert (std::make_pair(VALID_DOMAIN, TpDouble));
    theirTypeMap.insert (std::make_pair(N_POLY_AMP, TpInt));
    theirTypeMap.insert (std::make_pair(N_POLY_PHASE, TpInt));
    theirTypeMap.insert (std::make_pair(POLY_COEFF_AMP, TpDouble));
    theirTypeMap.insert (std::make_pair(POLY_COEFF_PHASE, TpDouble));
    theirTypeMap.insert (std::make_pair(PHASE_UNITS, TpString));

    theirTypeMap.insert (std::make_pair(SIDEBAND_REF, TpComplex));

    theirTypeMap.insert (std::make_pair(N_KNOTS_AMP, TpInt));
    theirTypeMap.insert (std::make_pair(N_KNOTS_PHASE, TpInt));
    theirTypeMap.insert (std::make_pair(SPLINE_KNOTS_AMP, TpDouble));
    theirTypeMap.insert (std::make_pair(SPLINE_KNOTS_PHASE, TpDouble));
  };

};

//----------------------------------------------------------------------------

String MSCalEnums::fieldName (Int enumField)
{
// Static function to look up the field name:
// Inputs:
//    enumField   Int     Field enumeration.
// Outputs:
//    fieldName   String  Field name.
// Exceptions:
//    Exception if invalid field enumeration.
//
  // Initialize map if empty
  if (theirFieldMap.empty()) initMaps();
  
  // Return the column name
  return theirFieldMap[enumField];
};

//----------------------------------------------------------------------------

Block<String> MSCalEnums::fieldNames (const Vector<Int>& enumFields)
{
// Static function to look up a set of field names:
// Inputs:
//    enumFields  const Vector<Int>&     Field enumerations.
// Outputs:
//    fieldNames  Block<String>          Field names.
// Exceptions:
//    Exception if invalid field enumeration.
//
  // Return the column names
  uInt nFields = enumFields.nelements();
  Block<String> names(nFields);
  for (uInt i=0; i < nFields; i++) {
    names[i] = fieldName (enumFields(i));
  };
  return names;
};

//----------------------------------------------------------------------------

DataType MSCalEnums::basicType (Int enumField)
{
// Static function to look up the basic field data type:
// Inputs:
//    enumField   Int        Field enumeration.
// Outputs:
//    basicType   DataType   Basic data type
// Exceptions:
//    Exception if invalid field enumeration.
//
  // Initialize map if empty
  if (theirTypeMap.empty()) initMaps();
  
  // Return the column name
  return theirTypeMap.at(enumField);
};

//----------------------------------------------------------------------------







} //# NAMESPACE CASACORE - END

