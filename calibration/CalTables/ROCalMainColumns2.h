//# CalMainColumns.h: Calibration table cal_main column access
//# Copyright (C) 1996,1997,1998,2001,2002,2003
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
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#ifndef CALIBRATION_ROCALMAINCOLUMNS2_H
#define CALIBRATION_ROCALMAINCOLUMNS2_H

#include <casacore/casa/aips.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/measures/TableMeasures/TableMeasColumn.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/TableMeasures/ScalarQuantColumn.h>
#include <calibration/CalTables/CalTable2.h>
#include <msvis/MSVis/MSCalEnums.h>

namespace casacore { //# NAMESPACE CASACORE - BEGIN

  template<class T>
  class ROCalMainColumns2
  {
  public:
    // Construct from a calibration table
    ROCalMainColumns2 (const CalTable2& calTable);
    
    // Default destructor
    virtual ~ROCalMainColumns2() {};
    
    // Read-only column accessors
    const ScalarColumn<Double>& time() const {return time_p;};
    const ScalarMeasColumn<MEpoch>& timeMeas() const {return timeMeas_p;};
    const ScalarColumn<Double>& timeEP() const {return timeEP_p;};
    const ScalarQuantColumn<Double>& timeEPQuant() const 
    {return timeEPQuant_p;};
    const ScalarColumn<Double>& interval() const {return interval_p;};
    const ScalarQuantColumn<Double>& intervalQuant() const
    {return intervalQuant_p;};
    const ScalarColumn<Int>& antenna1() const {return antenna1_p;};
    const ScalarColumn<Int>& feed1() const {return feed1_p;};
    const ScalarColumn<Int>& fieldId() const {return fieldId_p;};
    const ScalarColumn<Int>& arrayId() const {return arrayId_p;};
    const ScalarColumn<Int>& obsId() const {return obsId_p;};
    const ScalarColumn<Int>& scanNo() const {return scanNo_p;};
    const ScalarColumn<Int>& processorId() const {return processorId_p;};
    const ScalarColumn<Int>& stateId() const {return stateId_p;};
    const ScalarColumn<Int>& phaseId() const {return phaseId_p;};
    const ScalarColumn<Int>& pulsarBin() const {return pulsarBin_p;};
    const ScalarColumn<Int>& pulsarGateId() const {return pulsarGateId_p;};
    const ScalarColumn<Int>& freqGrp() const {return freqGrp_p;};
    const ScalarColumn<String>& freqGrpName() const {return freqGrpName_p;};
    const ScalarColumn<String>& fieldName() const {return fieldName_p;};
    const ScalarColumn<String>& fieldCode() const {return fieldCode_p;};
    const ScalarColumn<String>& sourceName() const {return sourceName_p;};
    const ScalarColumn<String>& sourceCode() const {return sourceCode_p;};
    const ScalarColumn<Int>& calGrp() const {return calGrp_p;};
    //const ArrayColumn<Complex>& gain() const {return gain_p;};
    const ArrayColumn<T>& gain() const {return gain_p;};
    const ArrayColumn<Float>& solvePar() const {return solvePar_p;};
    const ArrayColumn<Int>& refAnt() const {return refAnt_p;};
    const ArrayColumn<Int>& refFeed() const {return refFeed_p;};
    const ArrayColumn<Int>& refReceptor() const {return refReceptor_p;};
    const ArrayColumn<Double>& refFreq() const {return refFreq_p;};
    const ArrayMeasColumn<MFrequency>& refFreqMeas() const 
    {return refFreqMeas_p;};
    const ScalarColumn<Int>& measFreqRef() const {return measFreqRef_p;};
    const ArrayColumn<Double>& refDir() const {return refDir_p;};
    const ArrayMeasColumn<MDirection>& refDirMeas() const 
    {return refDirMeas_p;};
    const ScalarColumn<Int>& measDirRef() const {return measDirRef_p;};
    const ScalarColumn<Int>& calDescId() const {return calDescId_p;};
    const ScalarColumn<Int>& calHistoryId() const {return calHistoryId_p;};
    
  protected:
    // Prohibit public use of the null constructor, which
    // does not produce a usable object.
    ROCalMainColumns2() {};
    
    // Return a CalTable as a Table reference. Utilizes friendship
    // relationship with class CalTable.
    const Table& asTable(const CalTable2& calTable) 
    {return calTable.calMainAsTable();}
    
    // Attach a table column accessor
    void attach (const CalTable2& calTable, TableColumn& tabCol, 
		 MSCalEnums::colDef colEnum, const Bool& optional = False);
    void attach (const CalTable2& calTable, 
		 ArrayMeasColumn<MEpoch>& tabCol, 
		 MSCalEnums::colDef colEnum, const Bool& optional = False);
    void attach (const CalTable2& calTable, 
		 ArrayMeasColumn<MFrequency>& tabCol, 
		 MSCalEnums::colDef colEnum, const Bool& optional = False);
    void attach (const CalTable2& calTable, 
		 ArrayMeasColumn<MDirection>& tabCol, 
		 MSCalEnums::colDef colEnum, const Bool& optional = False);
    void attach (const CalTable2& calTable, ScalarMeasColumn<MEpoch>& tabCol, 
		 MSCalEnums::colDef colEnum, const Bool& optional = False);
    void attach (const CalTable2& calTable, ScalarQuantColumn<Double>& tabCol, 
		 MSCalEnums::colDef colEnum, const Bool& optional = False);
    
  private:
    // Prohibit copy constructor and assignment operator 
    ROCalMainColumns2 (const ROCalMainColumns2&);
    ROCalMainColumns2& operator= (const ROCalMainColumns2&);
    
    // Private column accessors
    ScalarColumn<Double> time_p;
    ScalarMeasColumn<MEpoch> timeMeas_p;
    ScalarColumn<Double> timeEP_p;
    ScalarQuantColumn<Double> timeEPQuant_p;
    ScalarColumn<Double> interval_p;
    ScalarQuantColumn<Double> intervalQuant_p;
    ScalarColumn<Int> antenna1_p;
    ScalarColumn<Int> feed1_p;
    ScalarColumn<Int> fieldId_p;
    ScalarColumn<Int> arrayId_p;
    ScalarColumn<Int> obsId_p;
    ScalarColumn<Int> scanNo_p;
    ScalarColumn<Int> processorId_p;
    ScalarColumn<Int> stateId_p;
    ScalarColumn<Int> phaseId_p;
    ScalarColumn<Int> pulsarBin_p;
    ScalarColumn<Int> pulsarGateId_p;
    ScalarColumn<Int> freqGrp_p;
    ScalarColumn<String> freqGrpName_p;
    ScalarColumn<String> fieldName_p;
    ScalarColumn<String> fieldCode_p;
    ScalarColumn<String> sourceName_p;
    ScalarColumn<String> sourceCode_p;
    ScalarColumn<Int> calGrp_p;
    //ArrayColumn<Complex> gain_p;
    ArrayColumn<T> gain_p;
    ArrayColumn<Float> solvePar_p;
    ArrayColumn<Int> refAnt_p;
    ArrayColumn<Int> refFeed_p;
    ArrayColumn<Int> refReceptor_p;
    ArrayColumn<Double> refFreq_p;
    ArrayMeasColumn<MFrequency> refFreqMeas_p;
    ScalarColumn<Int> measFreqRef_p;
    ArrayColumn<Double> refDir_p;
    ArrayMeasColumn<MDirection> refDirMeas_p;
    ScalarColumn<Int> measDirRef_p;
    ScalarColumn<Int> calDescId_p;
    ScalarColumn<Int> calHistoryId_p;
  };
  
}


#ifndef AIPS_NO_TEMPLATE_SRC
#include <calibration/CalTables/ROCalMainColumns2.tcc>
#endif

#endif
