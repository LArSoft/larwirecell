# Build the libWireCellLarsoft WCT plugin library.
# See README for some details.

file(GLOB wcls_components_sources "*.cxx")

# This is what's needed for building he WireCellLarsoft WCT component
# plugin library.
art_make_library(
  LIBRARY_NAME WireCellLarsoft
  SOURCE ${wcls_components_sources}
  LIBRARIES
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Persistency_Provenance
    art::Utilities
    ${JSONCPP}
    ROOT::Core
    ${WIRECELL_LIBS}
    canvas
    cetlib_except
    larcorealg_Geometry
    lardataalg_DetectorInfo
    lardataobj_RawData
    lardataobj_RecoBase
    lardataobj_Simulation
    larevt_CalibrationDBI_IOVData
)
