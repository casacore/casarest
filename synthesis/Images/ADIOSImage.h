#ifndef IMAGES_ADIOSIMAGE_H
#define IMAGES_ADIOSIMAGE_H

#ifdef CASACORE_HAS_ADIOS2

//# Includes
#include <casacore/casa/aips.h>
#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/ImageAttrHandlerCasa.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/Utilities/DataType.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/DataMan/Adios2StMan.h>

//# Forward Declarations
#include <casacore/casa/iosfwd.h>

// we need to be able to compile without MPI (although this class is not expected to be
// useful in such circumstances)
#ifdef HAVE_MPI

#include <mpi.h>

#else

#ifdef MPI_COMM_SELF
  #error "MPI_COMM_SELF appears to be defined without MPI available - this shouldn't have happened!"
#endif
typedef int MPI_Comm;
// the following can be any number
#define MPI_COMM_SELF ((MPI_Comm)0x44000001)
#endif

namespace casacore {

/// @details this class is a subclass of the casacore ImageInterface. It uses
///          casacore adios2 to read and write image cube from/to file. The 
///          class utilises the Adios2StMan as a Data Manager and it has two
///          columns named (1) "map" which stores the pixel value of the image cube
///          and (2) "mask" that keeps the mask of the image. These column names are
///          used as ADIOS variables. Also, this class does not support the retrival
///          of the mask of the image cube because this operation can cause excessive
///          memory usage on a compute node. However, the mask get/put operation per
///          image plane is supported.
template <class T> class ADIOSImage : public casacore::ImageInterface<T>
{
public:

    /// @brief create adios image accessor
    /// @param[in] mapShape of the image
    /// @param[in] coordinateInfo coordinate of the image
    /// @param[in] nameOfNewFile filenamee
    /// @param[in] tolerance used for compression
    /// @param[in] configname name of configuration file
    /// @param[in] rowNumber number of rows
    ADIOSImage(const casacore::TiledShape& mapShape,
               const casacore::CoordinateSystem& coordinateInfo,
               const casacore::String& nameOfNewFile,
               float tolerance = 0.0,
               casacore::String configname = "", 
               casacore::uInt rowNumber = 0);

    /// @brief create adios image accessor
    /// @param[in] comms MPI communicator
    /// @param[in] mapShape of the image
    /// @param[in] coordinateInfo coordinate of the image
    /// @param[in] nameOfNewFile filenamee
    /// @param[in] tolerance used for compression
    /// @param[in] configname name of configuration file
    /// @param[in] rowNumber number of rows
    ADIOSImage(MPI_Comm comms,
               const casacore::TiledShape& mapShape,
               const casacore::CoordinateSystem& coordinateInfo,
               const casacore::String& nameOfNewFile,
               float tolerance = 0.0,
               casacore::String configname = "",
               casacore::uInt rowNumber = 0);

    /// @brief create adios image accessor
    /// @param[in] filename name of the file to be created
    /// @param[in] tolerance used for compression
    /// @param[in] configname name of configuration file
    /// @param[in] spec mask value
    /// @param[in] rowNumber number of rows
    explicit ADIOSImage(const casacore::String &filename,
                        float tolerance = 0.0,
                        casacore::String configname = "",
                        casacore::uInt rowNumber = 0);

    /// @brief create adios image accessor
    /// @param[in] comms MPI communicator
    /// @param[in] filename name of the file to be created
    /// @param[in] tolerance used for compression
    /// @param[in] configname name of configuration file
    /// @param[in] spec mask value
    /// @param[in] rowNumber number of rows
    ADIOSImage(MPI_Comm comms,
               const casacore::String &filename,
               float tolerance = 0.0,
               casacore::String configname = "",
               casacore::uInt rowNumber = 0);

    /// @brief copy constructor
    /// @param[in] other an instance to copy from
    ADIOSImage (const ADIOSImage<T>& other);

    casacore::Table& table() { return itsTable; }

    static casacore::String className()
    { return "ADIOSImage"; }

    /// @brief Create a table that contains an array of columns to be managed by Adios2StMan
    /// @param[in] shape shape of the array/data
    /// @param[in] rowNumber number of rows
    /// @param[in] filename name of the file
    void makeNewTable (const casacore::TiledShape& shape, casacore::uInt rowNumber, casacore::String filename);

    ///@brief override base abstract methods
    casacore::ImageInterface<T>* cloneII() const override;
    casacore::String imageType() const override;
    void resize (const casacore::TiledShape& newShape) override;
    const casacore::LatticeRegion* getRegionPtr() const override;
    casacore::Bool setImageInfo(const casacore::ImageInfo& info) override;
    casacore::Bool setCoordinateInfo (const casacore::CoordinateSystem& coords);
    casacore::IPosition shape() const;
    ///@brief get a slice of the image cube
    casacore::Bool doGetSlice (casacore::Array<T>& buffer, const casacore::Slicer& theSlice);
    ///@brief get a slice of the mask.
    casacore::Bool doGetMaskSlice (casacore::Array<bool>& buffer, const casacore::Slicer& theSlice);
    ///@brief put a slice of data into the image cube
    void doPutSlice (const casacore::Array<T>& sourceBuffer, const casacore::IPosition& where, 
                     const casacore::IPosition& stride);
    void doPutMaskSlice(const casacore::Array<bool> &mask, const casacore::IPosition& where,
                        const casacore::IPosition& stride);
    casacore::Bool setUnits (const casacore::Unit& newUnits) override;
    casacore::String name(casacore::Bool stripPath=casacore::False) const override;

    ///@brief override base virtual methods
    casacore::Bool ok() const;
    casacore::Bool setMiscInfo (const casacore::RecordInterface& newInfo) override;
    casacore::Bool isPaged() const;
 
    ///@brief reopen the table column
    void reopenColumn();

    // masking operation on the whole image cube is not supported
    bool isMasked() const override { return false; }
    bool hasPixelMask() const override { return false; }
    virtual const Lattice<casacore::Bool>& pixelMask() const override;
    virtual Lattice<casacore::Bool>& pixelMask() override;

    static casacore::Table& getTable(void* imagePtr, casacore::Bool writable);
// MV: the code seems to build file with the following methods declared private (as they were in the
// PagedImage interface where it presumably was copied from). Unless there is a good reason not to, it
// would be better to keep them private. I won't insist on further design improvements for code which 
// goes into casarest/casacore.
private:
    //rewritten PagedImage private functions
    // @see casacore documentation for more details
    void attach_logtable();
    void setTableType();
    void restoreAll (const casacore::TableRecord& rec);
    void restoreImageInfo (const casacore::TableRecord& rec);
    void restoreUnits (const casacore::TableRecord& rec);
    void restoreMiscInfo (const casacore::TableRecord& rec);

    void reopenRW();

    ///@brief tolerance used for compression
    float itsTolerance;
    casacore::ArrayColumn<T> itsDataColumn;
    casacore::ArrayColumn<bool> itsMaskColumn;
    casacore::Table itsTable;
    
    const casacore::uInt itsRow;
    casacore::String itsConfig;
    // @brief MPI communicator the default MPI_COMM_SELF is assigned in the constructor (unless proper communicator is provided)
    MPI_Comm itsIOComms;
    int rank() const; 
    void setup(const casacore::TiledShape& shape, casacore::uInt rowNumber, 
               const casacore::String& filename, const casacore::CoordinateSystem& coordinateInfo);
    void setup(const casacore::String& filename);
};
// The tcc include must be here and not after the namespace
// otherwise the template definition in the tcc file must have 
// askap::imageaccess::ADIOSImage<T>::blah(...) instead of
// ADIOSImage<T>::blah(...)
#include <synthesis/Images//ADIOSImage.tcc>
}
#else
#warning "Using ADIOSImage requires building casacore with ADIOS2 support. It also requires building casarest with -DHAVE_ADIOS2=YES"
#endif // CASACORE_HAS_ADIOS2

#endif // IMAGES_ADIOSIMAGE_H
