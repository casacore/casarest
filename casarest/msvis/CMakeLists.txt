set(msvis_LIB_SRCS
  MSVis/AsynchronousTools.cc
  MSVis/CalVisBuffer.cc
  MSVis/MSAnalysis.cc
  MSVis/MsAverager.cc
  MSVis/MSCalEnums.cc
  MSVis/MSContinuumSubtractor.cc
  MSVis/MSUVWGenerator.cc
  MSVis/SelectAverageSpw.cc
  MSVis/SimpleSubMS.cc
  MSVis/StokesVector.cc
  MSVis/SubMS.cc
  MSVis/UtilJ.cc
  MSVis/VBContinuumSubtractor.cc
  MSVis/VisBuffAccumulator.cc
  MSVis/VisBuffer.cc
  MSVis/VisBufferAsync.cc
  MSVis/VisBufferAsyncWrapper.cc
  MSVis/VisBufferUtil.cc
  MSVis/VisBuffGroupAcc.cc
  MSVis/VisibilityIterator.cc
  MSVis/VisibilityIteratorAsync.cc
  MSVis/VisImagingWeight.cc
  MSVis/VisIterator.cc
  MSVis/VisSet.cc
  MSVis/VisSetUtil.cc
  MSVis/VisChunkAverager.cc
  MSVis/VLAT.cc
)


add_library(casa_msvis ${msvis_LIB_SRCS})
install(TARGETS casa_msvis DESTINATION ${LIBRARY_INSTALL_DIR})
target_link_libraries(casa_msvis ${OTHER_LIBRARIES})
set(MSVIS_LIBRARIES casa_msvis ${OTHER_LIBRARIES})


install (FILES
MSVis.h
DESTINATION include/casarest/msvis
)

install (FILES
MSVis/AsynchronousTools.h
MSVis/CalVisBuffer.h
MSVis/MSAnalysis.h
MSVis/MSCalEnums.h
MSVis/MSContinuumSubtractor.h
MSVis/MSMoments.h
MSVis/MSUVWGenerator.h
MSVis/MsAverager.h
MSVis/SelectAverageSpw.h
MSVis/SimpleSubMS.h
MSVis/StokesVector.h
MSVis/SubMS.h
MSVis/UtilJ.h
MSVis/VBContinuumSubtractor.h
MSVis/VLAT.h
MSVis/VisBuffAccumulator.h
MSVis/VisBuffGroupAcc.h
MSVis/VisBuffer.h
MSVis/VisBufferAsync.h
MSVis/VisBufferAsyncWrapper.h
MSVis/VisBufferUtil.h
MSVis/VisChunkAverager.h
MSVis/VisImagingWeight.h
MSVis/VisIterator.h
MSVis/VisSet.h
MSVis/VisSetUtil.h
MSVis/VisTimeAverager.h
MSVis/VisibilityIterator.h
MSVis/VisibilityIteratorAsync.h
MSVis/VisibilityProcessing.h
MSVis/MSMoments.tcc
MSVis/SubMS.tcc
MSVis/VisBuffer.tcc
MSVis/VisibilityIterator.tcc
DESTINATION include/casarest/msvis/MSVis
)
