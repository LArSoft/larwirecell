// art-service registration shim for wcls::WireSchemaService.
//
// The implementation (constructor, build(), wire_schema_store()) lives in
// larwirecell/Components/WireSchemaService.cxx so that it is part of
// libWireCellLarsoft and can be called directly by the WCT forwarder component
// wcls::WireSchemaSource.  This file only registers the service with art.
#include "larwirecell/Services/WireSchemaService.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

DEFINE_ART_SERVICE(wcls::WireSchemaService)
