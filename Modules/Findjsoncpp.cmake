if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FOUND)
  set(_cet_fjcpp_required ${${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED})
  set(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED)
  set(_cet_fjcpp_quietly ${${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY})
  find_package(${CMAKE_FIND_PACKAGE_NAME} CONFIG QUIET)
  if (${CMAKE_FIND_PACKAGE_NAME}_FOUND)
    set(_cet_fjcpp_config_mode CONFIG_MODE)
  else()
    if (DEFINED ENV{JSONCPP_FQ_DIR} AND EXISTS "$ENV{JSONCPP_FQ_DIR}/lib/pkgconfig/jsoncpp.pc")
      set(ENV{PKG_CONFIG_PATH} "$ENV{JSONCPP_FQ_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    endif()
    include(CetFindPkgConfigPackage)
    cet_find_pkg_config_package(jsoncpp QUIET)
  endif()
  set(${CMAKE_FIND_PACKAGE_NAME}_FIND_REQUIRED ${_cet_fjcpp_required})
  set(${CMAKE_FIND_PACKAGE_NAME}_FIND_QUIETLY ${_cet_fjcpp_quietly})
  unset(_cet_fjcpp_required)
  unset(_cet_fjcpp_quietly)
endif()
if (NOT ${CMAKE_FIND_PACKAGE_NAME}_FOUND)
  set(${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIR $ENV{JSONCPP_INC})
  find_library(${CMAKE_FIND_PACKAGE_NAME}_LIBRARY NAME ${CMAKE_FIND_PACKAGE_NAME}
    HINTS ENV JSONCPP_LIB NO_DEFAULT_PATH)
  mark_as_advanced(${CMAKE_FIND_PACKAGE_NAME}_LIBRARY)
  if (${CMAKE_FIND_PACKAGE_NAME}_LIBRARY)
    string(REGEX REPLACE "^v?([0-9._-]+).*$" "\\1" ${CMAKE_FIND_PACKAGE_NAME}_VERSION "$ENV{WIRECELL_VERSION}")
    string(REGEX REPLACE "[_-]" "." ${CMAKE_FIND_PACKAGE_NAME}_VERSION "${${CMAKE_FIND_PACKAGE_NAME}_VERSION}")
  endif()
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(${CMAKE_FIND_PACKAGE_NAME}
  ${_cet_fjcpp_config_mode}
  REQUIRED_VARS ${CMAKE_FIND_PACKAGE_NAME}_LIBRARIES
  ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIRS
  VERSION_VAR ${CMAKE_FIND_PACKAGE_NAME}_VERSION
  )

unset(_cet_fjcpp_config_mode)

if (${CMAKE_FIND_PACKAGE_NAME}_FOUND AND NOT TARGET jsoncpp_lib)
  if (TARGET PkgConfig::${CMAKE_FIND_PACKAGE_NAME})
    add_library(jsoncpp_lib ALIAS PkgConfig::${CMAKE_FIND_PACKAGE_NAME})
  elseif (${CMAKE_FIND_PACKAGE_NAME}_LIBRARY AND
      ${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIR)
    add_library(jsoncpp_lib UNKNOWN IMPORTED)
    set_target_properties(jsoncpp_lib PROPERTIES
      IMPORTED_LOCATION "${${CMAKE_FIND_PACKAGE_NAME}_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${${CMAKE_FIND_PACKAGE_NAME}_INCLUDE_DIR}")
  endif()
endif()
