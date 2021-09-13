
include_directories($ENV{EIGEN_INC})

add_subdirectory(Interfaces)
add_subdirectory(Components)
add_subdirectory(Tools)
add_subdirectory(Modules)
add_subdirectory(artcfg)

# deprecated:
add_subdirectory(LArInterface)
