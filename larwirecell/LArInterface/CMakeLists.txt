# larwirecell/LArInterface/CMakeLists.txt
#
# David Adams
# June 2016
#
# CMake instructions for larsoft service interface classes.

cet_build_plugin(WireCellNoiseFilter art::EDProducer
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  lardata::DetectorClocksService
  lardata::DetectorPropertiesService
  larevt::ChannelStatusService
  larevt::ChannelStatusProvider
  larevt::DetPedestalService
  larevt::DetPedestalProvider
  larevt::ElectronicsCalibService
  larevt::ElectronicsCalibProvider
  lardataobj::RawData
  art::Framework_Principal
  WireCell::Iface
  WireCell::SigProc
  )

install_headers()
install_fhicl()
install_source()
