
cet_make_library(LIBRARY_NAME WireCellQLMatch
    SOURCE
        QLMatching.cxx
        util.cxx
    LIBRARIES
    PUBLIC
        WireCell::Util
        WireCell::Aux
        WireCell::Iface
)
