# Build the libWireCellLarsoft WCT plugin library.
# See README for some details.

# This is what's needed for building the WireCellLarsoft WCT component
# plugin library.
cet_make_library(LIBRARY_NAME WireCellLarsoft
  SOURCE
  ChannelNoiseDB.cxx
  ChannelSelectorDB.cxx
  CookedFrameSink.cxx
  CookedFrameSource.cxx
  DepoFluxWriter.cxx
  FrameSaver.cxx
  LazyFrameSource.cxx
  MultiChannelNoiseDB.cxx
  OpFlashSource.cxx
  RawFrameSource.cxx
  SimDepoSetSource.cxx
  SimDepoSource.cxx
  WireSchemaBuilder.cxx
  # vvv needs WireSchema::correct() from WCT > 0.37.1 vvv
  # WireSchemaService.cxx
  # WireSchemaSource.cxx
  # ^^^
  # vvv obsolete vvv
  DepoSetSimChannelSink.cxx
  SimChannelSink.cxx
  # ^^^ obsolete ^^^
  LIBRARIES
  PUBLIC
  larwirecell::IArtEventVisitor
  larcoreobj::headers
  canvas::canvas
  WireCell::SigProc
  WireCell::Iface
  WireCell::Util
  lardataobj::Simulation
  PRIVATE
  larevt::CalibrationDBI_IOVData
  larevt::ChannelStatusService
  larevt::ChannelStatusProvider
  larevt::DetPedestalService
  larevt::DetPedestalProvider
  larevt::ElectronicsCalibService
  larevt::ElectronicsCalibProvider
  lardata::DetectorPropertiesService
  larcore::Geometry_Geometry_service
  larcore::ServiceUtil
  lardataobj::RecoBase
  lardataobj::RawData
  art::Framework_Core
  art::Framework_Principal
  larcorealg::Geometry
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  WireCell::Gen
  ROOT::Core
)
