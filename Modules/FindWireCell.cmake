#[============================================================[.rst:
FindWireCell
============
#]============================================================]
find_program(WIRE-CELL NAMES wire-cell HINTS ENV WIRECELL_FQ_DIR)
mark_as_advanced(WIRE-CELL)

set(_fwc_libs Ress Util Iface Apps Hio Pgraph Root Sig Tbb Aux Gen Img SigProc Sio)
set(_fwc_transitive_deps_Apps WireCell::Iface WireCell::Util)
set(_fwc_transitive_deps_Aux WireCell::Iface WireCell::Util Boost::headers Eigen3::Eigen)
set(_fwc_transitive_deps_Gen WireCell::Aux WireCell::Iface WireCell::Util Boost::headers Eigen3::Eigen)
set(_fwc_transitive_deps_Hio WireCell::Iface WireCell::Util h5cpp)
set(_fwc_transitive_deps_Iface WireCell::Util Boost::graph Boost::headers)
set(_fwc_transitive_deps_Img WireCell::Aux WireCell::Iface WireCell::Util)
set(_fwc_transitive_deps_Pgraph WireCell::Iface WireCell::Util Boost::headers)
set(_fwc_transitive_deps_Ress jsoncpp_lib jsonnet_lib Eigen3::Eigen)
set(_fwc_transitive_deps_Root WireCell::Iface WireCell::Util)
set(_fwc_transitive_deps_Sig WireCell::Iface WireCell::Util)
set(_fwc_transitive_deps_SigProc WireCell::Aux WireCell::Iface WireCell::Util)
set(_fwc_transitive_deps_Sio WireCell::Aux WireCell::Iface WireCell::Util Boost::iostreams)
set(_fwc_transitive_deps_Tbb WireCell::Aux WireCell::Iface WireCell::Util Boost::headers TBB::TBB)
set(_fwc_transitive_deps_Util Boost::date_time Boost::exception Boost::filesystem Boost::iostreams Boost::stacktrace_basic Eigen3::Eigen jsoncpp_lib jsonnet_lib spdlog::spdlog ${CMAKE_DL_LIBS})
set(_fwc_deps Boost Eigen3 jsoncpp jsonnet ROOT spdlog TBB)
set(_fwc_fp_Boost_args COMPONENTS graph headers date_time exception filesystem iostreams stacktrace_basic)
set(_fwc_fp_ROOT_args COMPONENTS Core Hist RIO Tree)

unset(_fwc_fphsa_extra_args)
if (WIRE-CELL)
  if (NOT WireCell_FOUND)
    execute_process(COMMAND ${WIRE-CELL} -v
      OUTPUT_VARIABLE WireCell_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if (NOT WireCell_VERSION)
      string(REGEX REPLACE "^v?([0-9._-]+).*$" "\\1" WireCell_VERSION "$ENV{WIRECELL_VERSION}")
      string(REGEX REPLACE "[_-]" "." WireCell_VERSION "${WireCell_VERSION}")
    endif()
    unset(WireCell_LIBRARIES)
    find_path(WireCell_INCLUDE_DIR WireCellApps/Main.h PATH_SUFFIXES include)
    mark_as_advanced(WireCell_INCLUDE_DIR)
    if (WireCell_INCLUDE_DIR)
      foreach (_fwc_lib IN LISTS _fwc_libs)
        find_library(WireCell_${_fwc_lib}_LIBRARY NAMES WireCell${_fwc_lib})
        mark_as_advanced(WireCell_${_fwc_lib}_LIBRARY)
        if (WireCell_${_fwc_lib}_LIBRARY)
          list(APPEND WireCell_LIBRARIES "WireCell::${_fwc_lib}")
        endif()
      endforeach()
    endif()
  endif()

  list(TRANSFORM _fwc_deps APPEND _FOUND
    OUTPUT_VARIABLE _fwc_fphsa_extra_required_vars)

  if (WireCell_FOUND OR (WireCell_LIBRARIES AND WireCell_INCLUDE_DIR))
    unset(_fwc_missing_deps)
    foreach (_fwc_dep IN LISTS _fwc_deps)
      get_property(_fwc_${_fwc_dep}_alreadyTransitive GLOBAL PROPERTY
        _CMAKE_${_fwc_dep}_TRANSITIVE_DEPENDENCY)
      find_package(${_fwc_dep} ${_fwc_fp_${_fwc_dep}_args} QUIET)
      if (NOT DEFINED cet_${_fwc_dep}_alreadyTransitive OR cet_${_fwc_dep}_alreadyTransitive)
        set_property(GLOBAL PROPERTY _CMAKE_${_fwc_dep}_TRANSITIVE_DEPENDENCY TRUE)
      endif()
      unset(_fwc_${_fwc_dep}_alreadyTransitive)
      if (NOT ${_fwc_dep}_FOUND)
        list(APPEND _fwc_missing_deps ${_fwc_dep})
      endif()
    endforeach()
    if (NOT "${_fwc_missing_deps}" STREQUAL "")
      set(_fwc_fphsa_extra_args
        REASON_FAILURE_MESSAGE "missing dependencies: ${_fwc_missing_deps}"
      )
      unset(_fwc_missing_deps)
    endif()
  endif()
else()
  set(_fwc_fphsa_extra_args "could not find executable \"wire-cell\"")
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(WireCell
  REQUIRED_VARS WIRE-CELL WireCell_INCLUDE_DIR WireCell_LIBRARIES
  ${_fwc_fphsa_extra_required_vars}
  VERSION_VAR WireCell_VERSION
  ${_fwc_fphsa_extra_args}
)
unset(_fwc_fphsa_extra_required_vars)
unset(_fwc_fphsa_extra_args)

if (WireCell_FOUND)
  foreach (_fwc_lib IN LISTS _fwc_libs)
    if (WireCell_${_fwc_lib}_LIBRARY AND NOT TARGET WireCell::${_fwc_lib})
      add_library(WireCell::${_fwc_lib} SHARED IMPORTED)
      set_target_properties(WireCell::${_fwc_lib} PROPERTIES
        IMPORTED_NO_SONAME TRUE
        IMPORTED_LOCATION "${WireCell_${_fwc_lib}_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${WireCell_INCLUDE_DIR}")
      foreach (_fwc_tdep IN LISTS _fwc_transitive_deps_${_fwc_lib})
        set_property(TARGET WireCell::${_fwc_lib}
          APPEND PROPERTY INTERFACE_LINK_LIBRARIES "${_fwc_tdep}")
      endforeach()
      unset(_fwc_tdep)
      unset(_fwc_transitive_deps_${_fwc_lib})
    endif()
  endforeach()
  if (NOT TARGET WireCell::IComponent)
    add_library(WireCell::IComponent INTERFACE IMPORTED
      INTERFACE_INCLUDE_DIRECTORIES "${WireCell_INCLUDE_DIR};${EIGEN3_INCLUDE_DIRS}")
  endif()
endif()

unset(_fwc_lib)
unset(_fwc_libs)
