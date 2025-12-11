

cet_make_library(LIBRARY_NAME WireCellAIML
    SOURCE
        Labelling2D.cxx
        Truth2h5.cxx
    LIBRARIES
    PUBLIC
        ${HDF5_LIBRARIES}
        WireCell::Util
        WireCell::Aux
        WireCell::Iface
        lardataobj::Simulation
        nusimdata::SimulationBase
        art::Framework_Core
        art::Framework_Principal
        larsim::MCCheater_ParticleInventoryService_service
)
