#include "larwirecell/Services/WireSchemaService.h"

#include "larwirecell/Components/WireSchemaBuilder.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "WireCellUtil/WireSchema.h"

wcls::WireSchemaService::WireSchemaService(fhicl::ParameterSet const& pset,
                                           art::ActivityRegistry&)
    : fGrouping(pset.get<std::string>("grouping", "tpcset"))
    , fCorrection(pset.get<int>("correction",
                      static_cast<int>(WireCell::WireSchema::Correction::pitch)))
    , fValidate(pset.get<bool>("validate", true))
{
}

void wcls::WireSchemaService::build() const
{
    if (m_built) return;

    auto const* geom = art::ServiceHandle<geo::Geometry const>().get();
    auto const& wires = art::ServiceHandle<geo::WireReadout const>()->Get();

    GroupingPolicy policy;
    policy.grouping = grouping_from_string(fGrouping);
    policy.correction = static_cast<WireCell::WireSchema::Correction>(fCorrection);
    policy.validate = fValidate;

    mf::LogInfo("WireSchemaService")
        << "building WCT wire schema from detector '" << geom->DetectorName()
        << "' (" << wires.Nchannels() << " channels), grouping=" << fGrouping;

    m_store = build_wire_schema(*geom, wires, policy);
    m_built = true;

    mf::LogInfo("WireSchemaService")
        << "built " << m_store.anodes().size() << " anodes, "
        << m_store.faces().size() << " faces, " << m_store.planes().size()
        << " planes, " << m_store.wires().size() << " wires";
}

const WireCell::WireSchema::Store& wcls::WireSchemaService::wire_schema_store() const
{
    build();
    return m_store;
}
