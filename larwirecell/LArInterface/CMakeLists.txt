# larwirecell/LArInterface/CMakeLists.txt
#
# David Adams
# June 2016
#
# CMake instructions for larsoft service interface classes.

art_make(
  MODULE_LIBRARIES
    lardataobj_RawData
    larevt_CalibrationDBI_IOVData
    larcorealg_Geometry
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    ${WIRECELL_LIBS}
    ROOT::Core
)

install_headers()
install_fhicl()
