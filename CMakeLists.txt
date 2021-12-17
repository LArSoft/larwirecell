# larwirecell/CMakeLists.txt
#
# David Adams
# June 2016
#
# CMake instructions for the larwirecell package.
#
# ======================================================================
# standalone (without mrb) build instructions:
#  cd .../path/to/build/directory
#  source .../path/to/larwirecell/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(larwirecell LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
set(${PROJECT}_CMAKE_PROJECT_VERSION_STRING 09.05.02-rc0)
cet_cmake_env()

set_install_root()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package(canvas_root_io REQUIRED)
find_package(art REQUIRED PUBLIC)
find_package(art_root_io REQUIRED PUBLIC)
find_package (Eigen3 3.3...3.3.9 REQUIRED NO_MODULE)
find_package(ROOT COMPONENTS Core Geom XMLIO Gdml REQUIRED PUBLIC)

find_ups_product( larcore )
find_ups_product( lardata )
find_ups_product( larevt )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( larcoreobj )
find_ups_product( larcorealg )

find_ups_product( jsoncpp )
find_ups_product( gojsonnet )
find_ups_product( spdlog )
find_ups_product( wirecell )


cet_find_library( WIRECELL_SIGPROC_LIB NAMES WireCellSigProc PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_GEN_LIB NAMES WireCellGen PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_IFACE_LIB NAMES WireCellIface PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_UTIL_LIB NAMES WireCellUtil PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
cet_find_library( WIRECELL_APPS_LIB NAMES WireCellApps PATHS ENV WIRECELL_LIB NO_DEFAULT_PATH )
set(WIRECELL_LIBS WIRECELL_APPS_LIB WIRECELL_SIGPROC_LIB WIRECELL_IFACE_LIB WIRECELL_UTIL_LIB WIRECELL_GEN_LIB)

cet_find_library( JSONCPP NAMES jsoncpp PATHS ENV JSONCPP_LIB NO_DEFAULT_PATH )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# sources
add_subdirectory(larwirecell)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

