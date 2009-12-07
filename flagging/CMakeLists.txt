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
