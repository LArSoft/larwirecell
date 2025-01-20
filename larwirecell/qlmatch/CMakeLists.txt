
cet_make_library(LIBRARY_NAME WireCellQLMatch
    SOURCE
        Opflash.cxx
        QLMatching.cxx
        TimingTPCBundle.cxx
        util.cxx
    LIBRARIES
    PUBLIC
        WireCell::Util
        WireCell::Aux
        WireCell::Clus
        WireCell::Iface
)
