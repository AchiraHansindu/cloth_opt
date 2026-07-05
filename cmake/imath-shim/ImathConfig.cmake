# Shim ImathConfig.cmake
#
# Ubuntu 22.04's libopenvdb-dev (8.1) is built against the old IlmBase 2.5
# "Half" library, but the bundled deps/openvdb/cmake/FindOpenVDB.cmake calls
# find_package(Imath REQUIRED CONFIG), expecting the standalone Imath 3.x
# package (ImathConfig.cmake), which Ubuntu 22.04 does not ship.
#
# This shim satisfies that find_package by exposing the installed IlmBase 2.5
# Half/Imath libraries (libilmbase-dev) as the Imath::Imath target that
# FindOpenVDB.cmake links OpenVDB against.

if(NOT TARGET Imath::Imath)
  find_library(_IMATH_SHIM_HALF NAMES Half Half-2_5)
  find_library(_IMATH_SHIM_IMATH NAMES Imath Imath-2_5)
  find_path(_IMATH_SHIM_INCLUDE NAMES OpenEXR/half.h)

  add_library(Imath::Imath INTERFACE IMPORTED)
  set_target_properties(Imath::Imath PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${_IMATH_SHIM_INCLUDE};${_IMATH_SHIM_INCLUDE}/OpenEXR"
    INTERFACE_LINK_LIBRARIES "${_IMATH_SHIM_HALF};${_IMATH_SHIM_IMATH}"
  )
endif()

set(Imath_FOUND TRUE)
