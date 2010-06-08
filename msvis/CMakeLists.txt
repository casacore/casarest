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
