#ifndef IMAGES_ADIOSIMAGE_TCC
#define IMAGES_ADIOSIMAGE_TCC

#include <casacore/images/Regions/ImageRegion.h>
#include <casacore/images/Regions/RegionHandlerTable.h>
#include <casacore/images/Images/ImageInfo.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/Lattices/LatticeNavigator.h>
#include <casacore/lattices/Lattices/LatticeStepper.h>
#include <casacore/lattices/Lattices/LatticeIterator.h>
#include <casacore/lattices/Lattices/PagedArrIter.h>
#include <casacore/lattices/LEL/LatticeExprNode.h>
#include <casacore/lattices/LEL/LatticeExpr.h>
#include <casacore/lattices/LRegions/LatticeRegion.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogMessage.h>

#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Arrays/LogiArray.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/OS/Path.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/TableLock.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableInfo.h>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Quanta/UnitMap.h>

#include <casacore/casa/iostream.h>
#include <casacore/casa/sstream.h>


template <class T> 
ADIOSImage<T>::ADIOSImage (
  const casacore::TiledShape& shape, 
  const casacore::CoordinateSystem& coordinateInfo, 
  const casacore::String& filename, 
  const float tolerance,
  const casacore::String& configname, 
  const casacore::rownr_t rowNumber)
: /* casacore::ImageInterface<T>(casacore::RegionHandlerTable(getTable, this)), */
  itsTolerance(tolerance), itsConfig(configname), itsRow(rowNumber),
  itsIOComms(MPI_COMM_SELF)
{
  this->setup(shape,rowNumber,filename, coordinateInfo);
}

template <class T> 
ADIOSImage<T>::ADIOSImage (
  MPI_Comm comms,
  const casacore::TiledShape& shape, 
  const casacore::CoordinateSystem& coordinateInfo, 
  const casacore::String& filename, 
  const float tolerance,
  const casacore::String& configname, 
  const casacore::rownr_t rowNumber)
: /* casacore::ImageInterface<T>(casacore::RegionHandlerTable(getTable, this)), */
  itsTolerance(tolerance),itsConfig(configname),itsRow(rowNumber),
  itsIOComms(comms)
{
  this->setup(shape,rowNumber,filename,coordinateInfo);
}

template <class T>
ADIOSImage<T>::ADIOSImage (
  const casacore::String& filename,
  const float tolerance,
  const casacore::String& configname, 
  const casacore::rownr_t rowNumber)
: /* casacore::ImageInterface<T>(casacore::RegionHandlerTable(getTable, this)), */
  itsTolerance(tolerance), itsConfig(configname),itsRow(rowNumber),
  itsIOComms(MPI_COMM_SELF) 
{
  this->setup(filename);
}

template <class T>
ADIOSImage<T>::ADIOSImage (
  MPI_Comm comms,
  const casacore::String& filename,
  const float tolerance,
  const casacore::String& configname, 
  const casacore::rownr_t rowNumber)
: /* casacore::ImageInterface<T>(casacore::RegionHandlerTable(getTable, this)), */
  itsTolerance(tolerance), itsConfig(configname),itsRow(rowNumber),
  itsIOComms(comms)
{
  this->setup(filename);
}

template <class T> 
ADIOSImage<T>::ADIOSImage (const ADIOSImage<T>& other)
: casacore::ImageInterface<T>(other),
  itsDataColumn(other.itsDataColumn),
  itsMaskColumn(other.itsMaskColumn),
  itsRow(other.itsRow)
{
}

template <class T>
void ADIOSImage<T>::makeNewTable(const casacore::TiledShape& shape, casacore::rownr_t rowNumber, const casacore::String& filename)
{
  const casacore::IPosition latShape = shape.shape();

  const casacore::uInt ndim = latShape.nelements();

  casacore::TableDesc description;
  description.addColumn(casacore::ArrayColumnDesc<T>("map",
                                                      casacore::String("version 4.0"),
                                                      latShape, casacore::ColumnDesc::FixedShape));
  description.addColumn(casacore::ArrayColumnDesc<bool>("mask",
                                                      casacore::String("version 4.0"),
                                                      latShape, casacore::ColumnDesc::FixedShape));

  casacore::SetupNewTable newtab(filename, description, casacore::Table::New);


  std::shared_ptr<casacore::Adios2StMan> adios2StMan;
  if (itsConfig == "") {
        std::string engineType = "";
        std::map<std::string, std::string> engineParams;
        std::vector<std::map<std::string, std::string>> transportParams;
        std::map<std::string, std::string> m;
        std::map<std::string, std::string> m2;
        std::vector<std::map<std::string, std::string>> v;
        if ( itsTolerance > 0.000001 ) {
            std::string accuracy = std::to_string(itsTolerance);
            m["Variable"] = "map";
            m["Operator"] = "mgard";
            m["accuracy"] = accuracy;
            m["s"] = "0.0";
            m["mode"] = "ABS";

            m2["Variable"] = "mask";
            m2["Operator"] = "mgard";
            m["accuracy"] = accuracy;
            m["s"] = "0.0";
            m["mode"] = "ABS";
            
            v.push_back(m);
            v.push_back(m2);
            if ( itsIOComms != MPI_COMM_SELF ) {
                adios2StMan.reset(new casacore::Adios2StMan(itsIOComms,engineType,
                                                        engineParams,transportParams,v));
            } else {
                adios2StMan.reset(new casacore::Adios2StMan(engineType, engineParams,transportParams,v));
            }
            newtab.bindColumn("map", *adios2StMan);
            newtab.bindColumn("mask", *adios2StMan);
        } else {
            m["Variable"] = "map";
            m2["Variable"] = "mask";
            v.push_back(m);
            v.push_back(m2);
            if ( itsIOComms != MPI_COMM_SELF ) {
                adios2StMan.reset(new casacore::Adios2StMan(itsIOComms,engineType,
                                                        engineParams,transportParams,v));
            } else {
                adios2StMan.reset(new casacore::Adios2StMan(engineType, engineParams,
                                                            transportParams,v));
            }

            newtab.bindColumn("map", *adios2StMan);
            newtab.bindColumn("mask", *adios2StMan);
        }
  } else {
    // invoke itsConfiguration based call 
    casacore::Adios2StMan::from_config_t from_config {};
    if ( itsIOComms != MPI_COMM_SELF ) {
        adios2StMan.reset(new casacore::Adios2StMan(itsIOComms,itsConfig, from_config));
    } else {
        adios2StMan.reset(new casacore::Adios2StMan(itsConfig, from_config));
    }
    newtab.bindColumn("map", *adios2StMan);
    newtab.bindColumn("mask", *adios2StMan);
  }
  
  itsTable = casacore::Table(itsIOComms, newtab);
  itsDataColumn = casacore::ArrayColumn<T>(itsTable,"map");
  itsMaskColumn = casacore::ArrayColumn<bool>(itsTable,"mask");
  const casacore::rownr_t rows = itsTable.nrow();
  if ((rowNumber + 1) > rows) {
    itsTable.addRow(rowNumber - rows + 1);
    for (casacore::rownr_t row = rows; row <= rowNumber - rows; row++){
      itsDataColumn.setShape(row, latShape);
      itsMaskColumn.setShape(row, latShape);
    }
  } 
}

template <class T>
casacore::Bool ADIOSImage<T>::setUnits(const casacore::Unit& newUnits)
{
  this->setUnitMember (newUnits);
  casacore::Table& tab = table();
  if (tab.keywordSet().isDefined("units")) {
    tab.rwKeywordSet().removeField("units");
  }
  if (!tab.isWritable()){
    tab.reopenRW();
  }
  tab.rwKeywordSet().define("units", newUnits.getName());
  return casacore::True;
}

template <class T>
casacore::Bool ADIOSImage<T>::setImageInfo(const casacore::ImageInfo& info)
{
  // Set imageinfo in base class.
  casacore::Bool ok = casacore::ImageInterface<T>::setImageInfo(info);
  if (ok) {
    // Make persistent in table keywords.
    casacore::Table& tab = table();
    // Delete existing one if there.
    if (tab.keywordSet().isDefined("imageinfo")) {
      tab.rwKeywordSet().removeField("imageinfo");
    }
    // Convert info to a record and save as keyword.
    casacore::TableRecord rec;
    casacore::String error;
    if (this->imageInfo().toRecord(error, rec)) {
        tab.rwKeywordSet().defineRecord("imageinfo", rec);
    } else {
      // Could not convert to record.
      casacore::LogIO os;
      os << casacore::LogIO::SEVERE << "Error saving ImageInfo in image " << name()
          << "; " << error << casacore::LogIO::POST;
      ok = casacore::False;
    }
  }
  return ok;
}

template<class T> 
casacore::Bool ADIOSImage<T>::setMiscInfo (const casacore::RecordInterface& newInfo)
{
  this->setMiscInfoMember(newInfo);
  casacore::Table& tab = table();
  if (tab.keywordSet().isDefined("miscinfo")) {
    tab.rwKeywordSet().removeField("miscinfo");
  }
  tab.rwKeywordSet().defineRecord("miscinfo", newInfo);
  return casacore::True;
}

template<class T> 
void ADIOSImage<T>::attach_logtable()
{
  // Open logtable as readonly if main table is not writable.
  casacore::Table& tab = table();
  this->setLogMember(casacore::LoggerHolder (name() + "/logtable", tab.isWritable()));
  // Insert the keyword if possible and if it does not exist yet.
  if (tab.isWritable()  &&  ! tab.keywordSet().isDefined ("logtable")) {
    tab.rwKeywordSet().defineTable("logtable", casacore::Table(name() + "/logtable"));
  }
}

template<class T> 
void ADIOSImage<T>::setTableType()
{
  casacore::TableInfo& info(table().tableInfo());
  const casacore::String reqdType = info.type(casacore::TableInfo::PAGEDIMAGE);
  if (info.type() != reqdType) {
    info.setType(reqdType);
  }
  const casacore::String reqdSubType = info.subType(casacore::TableInfo::PAGEDIMAGE);
  if (info.subType() != reqdSubType) {
    info.setSubType(reqdSubType);
  }
}

template <class T>
void ADIOSImage<T>::restoreAll (const casacore::TableRecord& rec)
{
  // Restore the coordinate system from the record
  // MV: I am not sure whether we're supposed to take the ownership of the returned raw pointer, but the original code 
  // written with raw pointers did the equivalent thing
  std::unique_ptr<casacore::CoordinateSystem> restoredCoords(casacore::CoordinateSystem::restore(rec, "coords"));
  AlwaysAssert(restoredCoords, casacore::AipsError);
  this->setCoordsMember(*restoredCoords);
  // Restore the image info.
  restoreImageInfo (rec);
  // Restore the units.
  restoreUnits (rec);
  // Restore the miscinfo.
  restoreMiscInfo (rec);
}

template<class T>
void ADIOSImage<T>::restoreImageInfo (const casacore::TableRecord& rec)
{
  if (rec.isDefined("imageinfo")) {
    casacore::String error;
    casacore::ImageInfo info;
    casacore::Bool ok = info.fromRecord (error, rec.asRecord("imageinfo"));
    if (ok) {
      this->setImageInfoMember(info);
    } else {
      casacore::LogIO os;
      os << casacore::LogIO::WARN << "Failed to restore the ImageInfo in image " << name()
         << "; " << error << casacore::LogIO::POST;
    }
  }
}

template<class T> 
void ADIOSImage<T>::restoreUnits (const casacore::TableRecord& rec)
{
  casacore::Unit retval;
  casacore::String unitName;
  if (rec.isDefined("units")) {
    if (rec.dataType("units") != casacore::TpString) {
      casacore::LogIO os;
      os << casacore::LogOrigin("ADIOSImage<T>", "units()", WHERE)
	 << "'units' keyword in image table is not a string! Units not restored." 
         << casacore::LogIO::SEVERE << casacore::LogIO::POST;
    } else {
      rec.get("units", unitName);
    }
  }
  if (! unitName.empty()) {
    // OK, non-empty unit, see if it's valid, if not try some known things to
    // make a valid unit out of it.
    if (! casacore::UnitVal::check(unitName)) {
      // Beam and Pixel are the most common undefined units
      casacore::UnitMap::putUser("Pixel",casacore::UnitVal(1.0),"Pixel unit");
      casacore::UnitMap::putUser("Beam",casacore::UnitVal(1.0),"Beam area");
    }
    if (! casacore::UnitVal::check(unitName)) {
      // OK, maybe we need FITS
      casacore::UnitMap::addFITS();
    }
    if (!casacore::UnitVal::check(unitName)) {
      // I give up!
      casacore::LogIO os;
      casacore::UnitMap::putUser(unitName, casacore::UnitVal(1.0, casacore::UnitDim::Dnon), unitName);
      os << casacore::LogIO::WARN << "FITS unit \"" << unitName
         << "\" unknown to CASA - will treat it as non-dimensional."
	 << casacore::LogIO::POST;
      retval.setName(unitName);
      retval.setValue(casacore::UnitVal(1.0, casacore::UnitDim::Dnon));
    } else {
      retval = casacore::Unit(unitName);
    }
  }
  this->setUnitMember(retval);
}

template<class T> 
void ADIOSImage<T>::restoreMiscInfo (const casacore::TableRecord& rec)
{
  if (rec.isDefined("miscinfo")  &&
      rec.dataType("miscinfo") == casacore::TpRecord) {
    this->setMiscInfoMember (rec.asRecord ("miscinfo"));
  }
}

template<class T>
casacore::Table& ADIOSImage<T>::getTable (void* imagePtr, casacore::Bool writable)
{
  ADIOSImage<T>* im = static_cast<ADIOSImage<T>*>(imagePtr);
  AlwaysAssert(im != NULL, casacore::AipsError);
  return im->itsTable;
}

template <class T> 
casacore::Bool ADIOSImage<T>::setCoordinateInfo (const casacore::CoordinateSystem& coords)
{
  casacore::Bool ok = casacore::ImageInterface<T>::setCoordinateInfo(coords);
  if (ok) {
    casacore::Table& tab = table();
    // Update the coordinates
    if (tab.keywordSet().isDefined("coords")) {
    tab.rwKeywordSet().removeField("coords");
    }
    if (!(this->coordinates().save(tab.rwKeywordSet(), "coords"))) {
      casacore::LogIO os;
os << casacore::LogIO::SEVERE << "Error saving coordinates in image " << name()
          << casacore::LogIO::POST;
ok = casacore::False;
    }
  }
  return ok;
}

template <class T> 
casacore::IPosition ADIOSImage<T>::shape() const
{
  return itsDataColumn.shape(itsRow);
}

template <class T> 
casacore::String ADIOSImage<T>::name(casacore::Bool stripPath) const 
{
  return itsTable.tableName();
}

template <class T> 
casacore::Bool ADIOSImage<T>::ok() const
{
  return (itsDataColumn.ndim(itsRow) == this->coordinates().nPixelAxes());
}

template <class T> 
casacore::Bool ADIOSImage<T>::doGetSlice(casacore::Array<T>& buffer, const casacore::Slicer& theSlice)
{
  itsDataColumn.getSlice(itsRow, theSlice, buffer, casacore::True);
  return casacore::False;
}

template <class T>
casacore::Bool ADIOSImage<T>::doGetMaskSlice(casacore::Array<bool>& mask, const casacore::Slicer& theSlice)
{
  itsMaskColumn.getSlice(itsRow, theSlice, mask, casacore::True);
  return casacore::False;
}

template <class T> 
void ADIOSImage<T>::doPutSlice(const casacore::Array<T>& sourceBuffer, const casacore::IPosition& where, const casacore::IPosition& stride)
{
  casacore::Array<bool> maskArr = ! casacore::isNaN(sourceBuffer);
  
  const casacore::uInt arrDim = sourceBuffer.ndim();
  const casacore::uInt latDim = this->ndim();
  AlwaysAssert(arrDim <= latDim, casacore::AipsError);
  if (arrDim == latDim) {
    casacore::Slicer section(where, sourceBuffer.shape(), stride, casacore::Slicer::endIsLength); 
    itsDataColumn.putSlice (itsRow, section, sourceBuffer);
    itsMaskColumn.putSlice (itsRow, section, maskArr);
  } else {
    casacore::Array<T> degenerateArr(sourceBuffer.addDegenerate(latDim-arrDim));
    casacore::Slicer section(where, degenerateArr.shape(), stride, casacore::Slicer::endIsLength); 
    itsDataColumn.putSlice (itsRow, section, degenerateArr);
    casacore::Array<bool> degenerateMaskArr(maskArr.addDegenerate(latDim-arrDim));
    itsMaskColumn.putSlice (itsRow, section, degenerateMaskArr);
  } 
}

template <class T>
void ADIOSImage<T>::doPutMaskSlice(const casacore::Array<bool> &mask, const casacore::IPosition& where, const casacore::IPosition& stride)
{
  const casacore::uInt arrDim = mask.ndim();
  const casacore::uInt latDim = this->ndim();
  AlwaysAssert(arrDim <= latDim, casacore::AipsError);
  if (arrDim == latDim) {
    casacore::Slicer section(where, mask.shape(), stride, casacore::Slicer::endIsLength);
    itsMaskColumn.putSlice (itsRow, section, mask);
  } else {
    casacore::Array<bool> degenerateArr(mask.addDegenerate(latDim-arrDim));
    casacore::Slicer section(where, degenerateArr.shape(), stride, casacore::Slicer::endIsLength);
    itsMaskColumn.putSlice (itsRow, section, degenerateArr);
  }
}

template<class T>
const casacore::LatticeRegion* ADIOSImage<T>::getRegionPtr() const
{
  throw (AipsError ("ADIOSImage::getRegionPtr - not implemented"));
  return nullptr;
}

template <class T> 
casacore::ImageInterface<T>* ADIOSImage<T>::cloneII() const
{
  return new ADIOSImage<T> (*this);
}

template<class T> 
void ADIOSImage<T>::resize (const casacore::TiledShape& newShape)
{
  if (newShape.shape().nelements() != this->coordinates().nPixelAxes()) {
    throw(casacore::AipsError("ADIOSImage<T>::resize: coordinate info is "
		    "the incorrect shape."));
  }
  casacore::IPosition tileShape = newShape.tileShape();
  itsDataColumn.setShape (itsRow, newShape.shape(), tileShape);
  itsMaskColumn.setShape (itsRow, newShape.shape(), tileShape);
}

template<class T>
casacore::String ADIOSImage<T>::imageType() const
{
  return className();
}

template<class T>
void ADIOSImage<T>::reopenRW()
{
  if (!itsTable.isWritable()) {
    itsTable.reopenRW();
    itsDataColumn = casacore::ArrayColumn<T>(itsTable, "map");
    itsMaskColumn = casacore::ArrayColumn<T>(itsTable, "mask");
  }
}

template<class T>
void ADIOSImage<T>::reopenColumn()
{
  itsDataColumn = casacore::ArrayColumn<T>(itsTable, "map");
  itsMaskColumn = casacore::ArrayColumn<T>(itsTable, "mask");
}

template<class T>
casacore::Bool ADIOSImage<T>::isPaged() const
{
  return casacore::True;
}
template<class T>
casacore::Lattice<casacore::Bool>& ADIOSImage<T>::pixelMask()
{
  throw (AipsError ("ADIOSImage::pixelMask - not implemented"));
}

template<class T>
const casacore::Lattice<casacore::Bool>& ADIOSImage<T>::pixelMask() const
{
  throw (AipsError ("ADIOSImage::pixelMask - not implemented"));
}

template <class T>
int ADIOSImage<T>::rank() const
{
#ifdef HAVE_MPI
  int size;
  int result = MPI_Comm_size(itsIOComms, &size);
  DebugAssert(result == MPI_SUCCESS, casacore::AipsError);
  int rank;
  result = MPI_Comm_rank(itsIOComms, &rank);
  DebugAssert(result == MPI_SUCCESS, casacore::AipsError);

  return rank;
#else
  return 0;
#endif
}

template <class T>
void ADIOSImage<T>::setup(const casacore::TiledShape& shape, 
                          const casacore::rownr_t rowNumber, 
                          const casacore::String& filename,
                          const casacore::CoordinateSystem& coordinateInfo)
{
  makeNewTable(shape, rowNumber, filename);
  const int r = rank();
  if ( itsIOComms == MPI_COMM_SELF || r == 0 ) {
    attach_logtable();
    setTableType();
  }
  AlwaysAssert(setCoordinateInfo(coordinateInfo), casacore::AipsError);
}

template <class T>
void ADIOSImage<T>::setup(const casacore::String& filename)
{
  itsTable = casacore::Table(filename,casacore::Table::TableOption::Old);
  itsDataColumn = casacore::ArrayColumn<T>(itsTable, "map");
  itsMaskColumn = casacore::ArrayColumn<bool>(itsTable, "mask");
  if ( itsIOComms == MPI_COMM_SELF ) {
    attach_logtable();
  }
  restoreAll (itsTable.keywordSet());
}
#endif
