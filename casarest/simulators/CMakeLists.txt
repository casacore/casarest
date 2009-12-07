set(simulators_LIB_SRCS
    Simulators/SimArray.cc
    Simulators/SimDataDesc.cc
    Simulators/SimDoppler.cc
    Simulators/SimFeed.cc
    Simulators/SimField.cc
    Simulators/SimObservations.cc
    Simulators/SimpleSimulator.cc
    Simulators/SimPol.cc
    Simulators/SimScan.cc
    Simulators/SimSource.cc
    Simulators/SimSpWin.cc
    Simulators/SimTelescope.cc
)


add_library(casa_simulators ${simulators_LIB_SRCS})
install(TARGETS casa_simulators DESTINATION ${LIBRARY_INSTALL_DIR})
target_link_libraries(casa_simulators ${OTHER_LIBRARIES})
set(SIMULATORS_LIBRARIES casa_simulators ${OTHER_LIBRARIES})
