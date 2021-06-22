# This is what's needed for building he WireCellLarsoft WCT component
# plugin library.
art_make(
  MODULE_LIBRARIES
    lardataobj_RawData
    lardataobj_RecoBase
    ${WIRECELL_LIBS}
)

add_subdirectory(Geometry)

install_headers()
