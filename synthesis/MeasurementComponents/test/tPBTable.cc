//# tPBTable.cc: This program tests the PB Table concept
//# Copyright (C) 1998,1999,2000,2001
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
//# $Id$

//# Includes
#include <casacore/casa/aips.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures.h>
#include <casacore/coordinates/Coordinates.h>
#include <casacore/casa/Arrays/Matrix.h>

#include <synthesis/MeasurementComponents/VPSkyJones.h>
#include <synthesis/MeasurementComponents/PBMath.h>
#include <casacore/tables/Tables.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableLock.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ScaRecordColDesc.h>
#include <casacore/tables/Tables/ScaRecordColData.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/DataMan/StManAipsIO.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/ExprNodeSet.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/ArrayUtil.h>
#include <casacore/casa/IO/ArrayIO.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/casa/Arrays/Slice.h>
#include <casacore/casa/Utilities/Sort.h>
#include <casacore/casa/OS/RegularFile.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/casa/Quanta/QuantumHolder.h>
#include <casacore/measures/Measures/MeasureHolder.h>

#include <casacore/casa/iostream.h>
#include <casacore/casa/stdio.h>

#include <casacore/casa/BasicSL/String.h>

#include <casacore/casa/namespace.h>
void buildTable(Table&);

main()
{
  try {
    cout << "Tests PB Table" << endl;
    cout << "--------------------------------------" << endl;


    // Build the table description.
    TableDesc td("", "1", TableDesc::Scratch);
    td.comment() = "A test of PB Table";
    td.addColumn (ScalarColumnDesc<String>("telescope","Which telescope is the PB valide for?"));
    td.addColumn (ScalarColumnDesc<Int>("antenna","Which antenna is this PB valid for? -1=all"));
    td.addColumn (ScalarRecordColumnDesc("pbdescription", "record containing PB description"));
 
    // Now create a new table from the description.
    // Use copy constructor to test if it works fine.
    // (newtab and newtabcp have the same underlying object).
    SetupNewTable newtab("PBTABLE.data", td, Table::New);
    Table tab(newtab);


    buildTable(tab);

    MeasurementSet ms("3C273XC1.ms", Table::Update);
    VPSkyJones myVPSet(ms, tab, Quantity(0.1, "deg"), BeamSquint::GOFIGURE);

  } catch (AipsError x) {
    cout << x.getMesg() << endl;
  } 
  
  exit(0);
};




void buildTable(Table& tab) 
{
  ScalarColumn<String> telCol(tab, "telescope");
  ScalarColumn<Int> antCol(tab, "antenna");
  ScalarColumn<TableRecord> recCol(tab, "pbdescription");

  uInt ii=0;
  
  {
    String telescope("VLA");
    Int antenna = -1;
    TableRecord record;
    
    String error;                   
    record.define (RecordFieldId("name"), "GAUSS");
    {
      Record subrec;                
      Quantity halfwidth(0.5, "deg");       
      if (!QuantumHolder(halfwidth).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("halfwidth"), subrec);
    }
    {
      Record subrec;                
      Quantity reffreq(1.0, "GHz");       
      if (!QuantumHolder(reffreq).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("reffreq"), subrec);
    }
    {
      Record subrec;                
      Quantity maxrad(1.0, "deg");       
      if (!QuantumHolder(maxrad).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("maxrad"), subrec);
    }
    {
      Record subrec;                
      MDirection  squintdir(Quantity(0.0, "'"), Quantity(0.0, "'"),
			    MDirection::Ref(MDirection::AZEL));
      
      if (!MeasureHolder(squintdir).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("squintdir"), subrec);
    }
    {
      Record subrec;                
      Quantity squintreffreq(1.0, "GHz");       
      if (!QuantumHolder(squintreffreq).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("squintreffreq"), subrec);
    }
    {
      Bool use=False;
      record.define(RecordFieldId("usesymmetricbeam"), use);
    }
    {
      Bool is=False;
      record.define(RecordFieldId("isthisvp"), is);
    }
    
    tab.addRow();
    telCol.put(ii, telescope);
    antCol.put(ii, antenna);
    recCol.put(ii, record);
    ii++;
  }
  
  
  {
    String telescope("GMRT");
    Int antenna = -1;
    TableRecord record;
    
    String error;                   
    record.define (RecordFieldId("name"), "GAUSS");
    {
      Record subrec;                
      Quantity halfwidth(0.2, "deg");       
      if (!QuantumHolder(halfwidth).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("halfwidth"), subrec);
    }
    {
      Record subrec;                
      Quantity reffreq(1.0, "GHz");       
      if (!QuantumHolder(reffreq).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("reffreq"), subrec);
    }
    {
      Record subrec;                
      Quantity maxrad(1.0, "deg");       
      if (!QuantumHolder(maxrad).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("maxrad"), subrec);
    }
    {
      Record subrec;                
      MDirection  squintdir(Quantity(1.0, "'"), Quantity(0.0, "'"),
			    MDirection::Ref(MDirection::AZEL));
      
      if (!MeasureHolder(squintdir).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("squintdir"), subrec);
    }
    {
      Record subrec;                
      Quantity squintreffreq(1.0, "GHz");       
      if (!QuantumHolder(squintreffreq).toRecord(error, subrec)) { 
	cout << error << endl;
      };
      record.defineRecord(RecordFieldId("squintreffreq"), subrec);
    }
    {
      Bool use=False;
      record.define(RecordFieldId("usesymmetricbeam"), use);
    }
    {
      Bool is=False;
      record.define(RecordFieldId("isthisvp"), is);
    }
    
    
    tab.addRow();
    telCol.put(ii, telescope);
    antCol.put(ii, antenna);
    recCol.put(ii, record);
    ii++;
  }



  {
    String telescope("WSRT");
    Int antenna = -1;
    TableRecord record;
    
    String error;                   
    record.define (RecordFieldId("name"), "COMMONPB");
    record.define (RecordFieldId("commonpb"), "WSRT");
    {
      Bool use=False;
      record.define(RecordFieldId("usesymmetricbeam"), use);
    }
    {
      Bool is=False;
      record.define(RecordFieldId("isthisvp"), is);
    }
    
    tab.addRow();
    telCol.put(ii, telescope);
    antCol.put(ii, antenna);
    recCol.put(ii, record);
    ii++;
  }
  

  {
    String telescope("FUTURE TELESCOPE");
    Int antenna = -1;
    TableRecord record;
    
    String error;                   
    record.define (RecordFieldId("name"), "COMMONPB");
    record.define (RecordFieldId("commonpb"), "VLA_INVERSE");
    {
      Bool use=False;
      record.define(RecordFieldId("usesymmetricbeam"), use);
    }
    {
      Bool is=False;
      record.define(RecordFieldId("isthisvp"), is);
    }
    
    tab.addRow();
    telCol.put(ii, telescope);
    antCol.put(ii, antenna);
    recCol.put(ii, record);
    ii++;
  }
  


};


