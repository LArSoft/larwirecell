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
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
                           canvas::canvas
                           messagefacility::MF_MessageLogger
                           ${MF_UTILITIES}
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
                           cetlib_except::cetlib_except
                           ROOT::Geom
                           ROOT::XMLIO
                           ROOT::Gdml
                           ${ROOT_BASIC_LIB_LIST}
        )

#install_headers()
# install_fhicl()
install_source()
