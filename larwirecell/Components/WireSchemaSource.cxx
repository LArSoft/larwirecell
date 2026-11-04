#include "larwirecell/Components/WireSchemaSource.h"

#include "larwirecell/Services/WireSchemaService.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(wclsWireSchemaSource,
                 wcls::WireSchemaSource,
                 WireCell::IWireSchema,
                 WireCell::IConfigurable)

using namespace WireCell;

wcls::WireSchemaSource::WireSchemaSource() {}
wcls::WireSchemaSource::~WireSchemaSource() {}

WireCell::Configuration wcls::WireSchemaSource::default_configuration() const
{
    // Nothing to configure: the geometry is owned by the art service.  The
    // grouping/correction/validate knobs live on the WireSchemaService FHiCL.
    return Configuration{};
}

void wcls::WireSchemaSource::configure(const WireCell::Configuration& /*cfg*/)
{
    // Pull the store from the art service.  This runs during WCT configuration
    // (eg from AnodePlane::configure via the "uses" ordering), by which point
    // the geometry/wire-readout art services are available.
    m_store = art::ServiceHandle<wcls::WireSchemaService const>()->wire_schema_store();
}

const WireSchema::Store& wcls::WireSchemaSource::wire_schema_store() const
{
    return m_store;
}
