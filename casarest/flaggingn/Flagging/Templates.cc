#include <flagging/Flagging/RFCubeLattice.tcc>
#include <flagging/Flagging/RFChunkStats.tcc>
#include <flagging/Flagging/RFFlagCube.tcc>

namespace casa {

  template class RFCubeLattice<uInt>;
  template class RFCubeLattice<Float>;
  template class RFCubeLatticeIterator<uInt>;
  template class RFCubeLatticeIterator<Float>;

  template RFlagWord RFChunkStats::getCorrMask (const Vector<int>&);

  template LogicalArray maskBits (const Array<uInt>&, const uInt&);
}
