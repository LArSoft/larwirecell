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

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(larwirecell VERSION 09.14.01 LANGUAGES CXX)

# cetbuildtools contains our cmake modules

include(CetCMakeEnv)
cet_cmake_env()

cet_cmake_module_directories(Modules)

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(art REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core REQUIRED EXPORT)

find_package(larcore REQUIRED EXPORT)
find_package(lardata REQUIRED EXPORT)
find_package(larevt REQUIRED EXPORT)
find_package(lardataobj REQUIRED EXPORT)
find_package(lardataalg REQUIRED EXPORT)
find_package(larcoreobj REQUIRED EXPORT)
find_package(larcorealg REQUIRED EXPORT)

find_package(WireCell REQUIRED EXPORT)

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(CetMake)
include(BuildPlugins)

add_subdirectory(Modules)

# sources
add_subdirectory(larwirecell)

# ups - table and config files

# packaging utility

cet_cmake_config()
