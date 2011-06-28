set(flagging_LIB_SRCS
    Flagging/DDMapper.cc
    Flagging/ExtendFlagger.cc
    Flagging/Flagger.cc
    Flagging/FlagIndex.cc
    Flagging/RedFlagger.cc
    Flagging/RFAApplyFlags.cc
    Flagging/RFADiffBase.cc
    Flagging/RFAFlagCubeBase.cc
    Flagging/RFAFlagExaminer.cc
    Flagging/RFAMedianClip.cc
    Flagging/RFANewMedianClip.cc
    Flagging/RFASelector2.cc
    Flagging/RFASpectralRej.cc
    Flagging/RFATimeFreqCrop.cc
    Flagging/RFAUVBinner.cc
    Flagging/RFChunkStats.cc
    Flagging/RFDataMapper.cc
    Flagging/RFDebugPlot.cc
    Flagging/RFFlagCube.cc
    Flagging/RFRowClipper.cc
    Flagging/Templates.cc
)


add_library(casa_flagging ${flagging_LIB_SRCS})
install(TARGETS casa_flagging DESTINATION ${LIBRARY_INSTALL_DIR})
target_link_libraries(casa_flagging casa_msvis ${OTHER_LIBRARIES})
set(FLAGGING_LIBRARIES casa_flagging casa_msvis ${OTHER_LIBRARIES})


install (FILES
Flagging.h
Flagging/DDMapper.h
Flagging/ExtendFlagger.h
Flagging/FlagIndex.h
Flagging/Flagger.h
Flagging/RFAApplyFlags.h
Flagging/RFABase.h
Flagging/RFADiffBase.h
Flagging/RFAFlagCubeBase.h
Flagging/RFAFlagExaminer.h
Flagging/RFAMedianClip.h
Flagging/RFANewMedianClip.h
Flagging/RFASelector.h
Flagging/RFASpectralRej.h
Flagging/RFATimeFreqCrop.h
Flagging/RFAUVBinner.h
Flagging/RFChunkStats.h
Flagging/RFCommon.h
Flagging/RFCubeLattice.h
Flagging/RFDataMapper.h
Flagging/RFDebugPlot.h
Flagging/RFFlagCube.h
Flagging/RFRowClipper.h
Flagging/RedFlagger.h
Flagging/RFASelector.tcc
Flagging/RFChunkStats.tcc
Flagging/RFCubeLattice.tcc
Flagging/RFFlagCube.tcc
DESTINATION include/casarest/msvis/MSVis
)
