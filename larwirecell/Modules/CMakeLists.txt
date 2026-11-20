cet_build_plugin(BlipMaker art::EDProducer
  LIBRARIES PRIVATE
  lardataobj::Simulation)

cet_build_plugin(EventButcher art::EDProducer
  LIBRARIES PRIVATE
  lardataobj::RecoBase
  lardataobj::RawData
  canvas::canvas
  art::Persistency_Common
  )

cet_build_plugin(WireCellToolkit art::SharedProducer
  LIBRARIES PRIVATE
  art::Utilities
  )

cet_build_plugin(WireSchemaFileWriter art::EDAnalyzer
  LIBRARIES PRIVATE
  larwirecell::WireCellLarsoft
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  WireCell::Util
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  )

add_subdirectory(Geometry)

install_headers()
install_source()
