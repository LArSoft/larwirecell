add_subdirectory(Interfaces)
add_subdirectory(Components)
# Services needs WireSchema::correct() from WCT > 0.37.1
# add_subdirectory(Services)
add_subdirectory(Tools)
add_subdirectory(Modules)
add_subdirectory(artcfg)
add_subdirectory(qlmatch)
add_subdirectory(aiml)

# deprecated:
add_subdirectory(LArInterface)
