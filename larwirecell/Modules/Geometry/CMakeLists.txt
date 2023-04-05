cet_build_plugin(CTreeGeometry art::EDAnalyzer
  LIBRARIES PRIVATE
  lardata::Utilities
  lardataobj::Simulation
  lardataobj::RecoBase
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardata::DetectorPropertiesService
  lardata::DetectorClocksService
  art::Framework_Principal
  canvas::canvas
  )

install_source()
