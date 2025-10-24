
cet_make_library(LIBRARY_NAME WireCellAIML
    SOURCE
        Labelling2D.cxx
    LIBRARIES
    PUBLIC
        WireCell::Util
        WireCell::Aux
        WireCell::Iface
        lardataobj::Simulation
        art::Framework_Core
        art::Framework_Principal
    PRIVATE
        larsim::MCCheater_ParticleInventoryService_service
)
