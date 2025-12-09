
cet_make_library(LIBRARY_NAME WireCellQLMatch
    SOURCE
        Opflash.cxx
        QLMatching.cxx
        TimingTPCBundle.cxx
        util.cxx
    LIBRARIES
    PUBLIC 
        larsim::PhotonPropagation
        art_plugin_support::toolMaker
        WireCell::Util
        WireCell::Aux
        WireCell::Clus
        WireCell::Iface
)
install_fhicl()