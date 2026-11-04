cet_build_plugin(WireSchemaService art::service
  LIBRARIES PRIVATE
  larwirecell::WireCellLarsoft
  WireCell::Util
  )

install_headers()
install_source()
