set(ionosphere_LIB_SRCS
    Ionosphere/GPSEphemeris.cc
    Ionosphere/GPSGroupDelay.cc
    Ionosphere/IonosphDataGPS.cc
    Ionosphere/Ionosphere.cc
    Ionosphere/IonosphModelPIM.cc
    Ionosphere/RINEX.cc
    Ionosphere/RINEXSat.cc
)


add_library(casa_ionosphere ${ionosphere_LIB_SRCS})
install(TARGETS casa_ionosphere DESTINATION ${LIBRARY_INSTALL_DIR})
target_link_libraries(casa_ionosphere ${OTHER_LIBRARIES})
set(IONOSPHERE_LIBRARIES casa_ionosphere ${OTHER_LIBRARIES})
