project(casarest)
enable_language(CXX Fortran)
cmake_minimum_required(VERSION 3.3.0)
include(CheckCXXCompilerFlag)
#list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# rpath setup for libraries

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# HDF5 uses an ENV variable.
if(HDF5_ROOT_DIR)
  set(ENV{HDF5_ROOT} ${HDF5_ROOT_DIR})
endif(HDF5_ROOT_DIR)

# project version
set( ${PROJECT_NAME}_MAJOR_VERSION 1 )
set( ${PROJECT_NAME}_MINOR_VERSION 9 )
set( ${PROJECT_NAME}_PATCH_LEVEL 0 )

if (UseCasaNamespace)
    add_definitions (-DUseCasaNamespace)
    message (STATUS "Using namespace casa. This will be deprecated at some point.")
else ()
    message (STATUS "Using namespace casacore.")
endif ()

# Make use of test programs
include (CTest)

# Determine which external packages to use.
# dependencies
find_package(CasaCore 3.4 REQUIRED)
find_package(CFITSIO 3.030 REQUIRED) # Should pad to three decimal digits
find_package(WCSLIB 4.7 REQUIRED)    # needed for CASA
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread system)
find_package(HDF5)
if(HDF5_FOUND)
    add_definitions(-DHAVE_HDF5)
else(HDF5_FOUND)
    message(STATUS "  HDF5 not used")
    set(HDF5_LIBRARIES )
    set(HDF5_INCLUDE_DIRS )
endif(HDF5_FOUND)

# options and defaults
set( BUILD_SHARED_LIBS TRUE )

if(NOT LIBRARY_INSTALL_DIR )
    set( LIBRARY_INSTALL_DIR lib )
endif(NOT LIBRARY_INSTALL_DIR )

if(NOT INCLUDE_INSTALL_DIR )
    set( INCLUDE_INSTALL_DIR "include/casarest" )
endif(NOT INCLUDE_INSTALL_DIR )

#add_definitions(-DAIPS_NO_TEMPLATE_SRC -DCASA_STANDALONE)
add_definitions(-DCASA_STANDALONE)


# -------- library building ----------------------------
INCLUDE_DIRECTORIES(${CASACORE_INCLUDE_DIR} 
    ${CASACORE_INCLUDE_DIR}
    ${CFITSIO_INCLUDE_DIR}
    ${WCSLIB_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(OTHER_LIBRARIES ${CASACORE_LIBRARIES}
    ${WCSLIB_LIBRARIES}
    ${CFITSIO_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${Boost_LIBRARIES}
)

add_subdirectory(components)
add_subdirectory(msvis)
add_subdirectory(calibration)
add_subdirectory(synthesis)
if(BUILD_ALL)
##    message(STATUS "  Also building flagging and simulators")
    message(STATUS "  flagging and simulators are currently not buildable")
# flagging and simulators do not build well at the moment
##    add_subdirectory(flagging)
##    add_subdirectory(simulators)
# Ionosphere requires PIM software which is not available
##  add_subdirectory(ionosphere)
endif(BUILD_ALL)
