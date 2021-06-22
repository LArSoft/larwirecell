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

cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(larwirecell VERSION 09.02.06 LANGUAGES CXX)

find_package(cetmodules REQUIRED)
find_package(canvas_root_io)
find_package(art)

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcore )
find_ups_product( lardata )
find_ups_product( lardataalg )
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( art )
find_package(art_root_io)

find_ups_product( jsoncpp )
find_ups_product( jsonnet )
find_ups_product( eigen )
find_ups_product( spdlog )
find_ups_product( wirecell )


cet_find_library( WIRECELL_SIGPROC_LIB NAMES WireCellSigProc PATHS ENV WIRECELL_LIB  )
cet_find_library( WIRECELL_GEN_LIB NAMES WireCellGen PATHS ENV WIRECELL_LIB  )
cet_find_library( WIRECELL_IFACE_LIB NAMES WireCellIface PATHS ENV WIRECELL_LIB  )
cet_find_library( WIRECELL_UTIL_LIB NAMES WireCellUtil PATHS ENV WIRECELL_LIB  )
cet_find_library( WIRECELL_APPS_LIB NAMES WireCellApps PATHS ENV WIRECELL_LIB  )
set(WIRECELL_LIBS ${WIRECELL_APPS_LIB} ${WIRECELL_SIGPROC_LIB} ${WIRECELL_IFACE_LIB} ${WIRECELL_UTIL_LIB} ${WIRECELL_GEN_LIB})

cet_find_library( JSONCPP NAMES jsoncpp PATHS ENV JSONCPP_LIB  )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

# sources
add_subdirectory(larwirecell)
add_subdirectory(cmakeConfig)
# ups - table and config files
if(WANT_UPS)
add_subdirectory(ups)
endif()

# packaging utility
#include(UseCPack)

