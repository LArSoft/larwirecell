# Interface target for includers of IArtEventVisitor.h
cet_make_library(LIBRARY_NAME IArtEventVisitor
  INTERFACE SOURCE IArtEventVisitor.h
  LIBRARIES INTERFACE WireCell::IComponent)

install_headers()
install_source()
