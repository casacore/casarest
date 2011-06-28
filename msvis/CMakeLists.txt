set(msvis_LIB_SRCS
    MSVis/CalVisBuffer.cc
    MSVis/MsAverager.cc
    MSVis/MSCalEnums.cc
    MSVis/MSContinuumSubtractor.cc
    MSVis/MSUVWGenerator.cc
    MSVis/SelectAverageSpw.cc
    MSVis/StokesVector.cc
    MSVis/SubMS.cc
    MSVis/VisBuffAccumulator.cc
    MSVis/VisBuffer.cc
    MSVis/VisBufferUtil.cc
    MSVis/VisBuffGroupAcc.cc
    MSVis/VisibilityIterator.cc
    MSVis/VisImagingWeight.cc
    MSVis/VisSet.cc
    MSVis/VisSetUtil.cc
    MSVis/VisTimeAverager.cc
)


add_library(casa_msvis ${msvis_LIB_SRCS})
install(TARGETS casa_msvis DESTINATION ${LIBRARY_INSTALL_DIR})
target_link_libraries(casa_msvis ${OTHER_LIBRARIES})
set(MSVIS_LIBRARIES casa_msvis ${OTHER_LIBRARIES})


install (FILES
MSVis.h
MSVis/CalVisBuffer.h
MSVis/MSCalEnums.h
MSVis/MSContinuumSubtractor.h
MSVis/MSUVWGenerator.h
MSVis/MsAverager.h
MSVis/SelectAverageSpw.h
MSVis/StokesVector.h
MSVis/SubMS.h
MSVis/VisBuffAccumulator.h
MSVis/VisBuffGroupAcc.h
MSVis/VisBuffer.h
MSVis/VisBufferUtil.h
MSVis/VisImagingWeight.h
MSVis/VisSet.h
MSVis/VisSetUtil.h
MSVis/VisTimeAverager.h
MSVis/VisibilityIterator.h
MSVis/SubMS.tcc
DESTINATION include/casarest/msvis/MSVis
)
