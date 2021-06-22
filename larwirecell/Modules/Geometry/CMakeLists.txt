#add_subdirectory(tools)

art_make( MODULE_LIBRARIES larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardataobj_Simulation
                           lardata_Utilities
                           larevt_Filters
                           lardataobj_RawData
                           larevt_CalibrationDBI_Providers
                           art::Framework_Core
                           art::Framework_Principal
                           art::Framework_Services_Registry
                           art::Framework_Services_System_FileCatalogMetadata_service art::Framework_Services_System_FloatingPointControl_service art::Framework_Services_System_TriggerNamesService_service
                           art::Framework_Services_System_FloatingPointControl_service
                           art_root_io::TFileService_service_plugin
                           art_root_io::tfile_support 
                           art_Persistency_Common
                           art_Persistency_Provenance
                           art_Utilities
                           canvas
                           messagefacility::MF_MessageLogger
                           messagefacility::headers
                           fhiclcpp::fhiclcpp fhiclcpp::types
                           ${CETLIB}
                           cetlib_except
                           ROOT::Core
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ROOT::MathCore
                           ROOT::Core
                           
        )

#install_headers()
# install_fhicl()
install_source()
