#include "larwirecell/Components/DriftPlanesSource.h"

#include "larwirecell/Services/WireSchemaService.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <cmath>

WIRECELL_FACTORY(wclsDriftPlanesSource,
                 wcls::DriftPlanesSource,
                 WireCell::IDriftPlanes,
                 WireCell::IConfigurable)

using namespace WireCell;

wcls::DriftPlanesSource::DriftPlanesSource() {}
wcls::DriftPlanesSource::~DriftPlanesSource() {}

WireCell::Configuration wcls::DriftPlanesSource::default_configuration() const
{
    Configuration cfg;
    // Distance from the collection plane to the response plane, against the
    // drift.  Transcribe it from the "origin" of the field response.
    cfg["response_distance"] = 10 * units::cm;
    // Faces whose planes all agree within this are one xregion.
    cfg["tolerance"] = 0.1 * units::mm;
    return cfg;
}

void wcls::DriftPlanesSource::configure(const WireCell::Configuration& cfg)
{
    const double response_distance = get<double>(cfg, "response_distance", 10 * units::cm);
    const double tolerance = get<double>(cfg, "tolerance", 0.1 * units::mm);

    const auto& dfs = art::ServiceHandle<wcls::WireSchemaService const>()->drift_faces();

    m_faces.clear();
    m_xregions = Json::arrayValue;
    std::vector<FaceDrift> distinct;
    for (const auto& [anode, faces] : dfs) {
        Configuration jfaces = Json::arrayValue;
        for (const auto& fd : faces) {
            // "Against the drift" is from the planes toward the cathode.
            const double response = fd.collection - fd.drift_sign * response_distance;
            if (fd.drift_sign * (fd.anode - response) <= 0) {
                raise<ValueError>("anode %d: anode bound at x=%f mm is not on the wire "
                                  "side of the response plane at x=%f mm",
                                  anode, fd.anode, response);
            }
            if (fd.drift_sign * (response - fd.cathode) <= 0) {
                raise<ValueError>("anode %d: response plane at x=%f mm is not between "
                                  "the cathode at x=%f mm and the wires",
                                  anode, response, fd.cathode);
            }
            Configuration jface;
            jface["anode"] = fd.anode;
            jface["response"] = response;
            jface["cathode"] = fd.cathode;
            jfaces.append(jface);

            bool seen = false;
            for (const auto& d : distinct) {
                if (std::abs(d.anode - fd.anode) < tolerance &&
                    std::abs(d.collection - fd.collection) < tolerance &&
                    std::abs(d.cathode - fd.cathode) < tolerance) {
                    seen = true;
                    break;
                }
            }
            if (!seen) {
                distinct.push_back(fd);
                m_xregions.append(jface);
            }
        }
        m_faces[anode] = jfaces;
    }

    auto log = Log::logger("wcls");
    log->info("DriftPlanesSource: {} anodes, {} xregions, response distance {} mm",
              m_faces.size(), m_xregions.size(), response_distance);
    for (const auto& jxr : m_xregions) {
        log->info("DriftPlanesSource: xregion anode={} response={} cathode={} mm",
                  jxr["anode"].asDouble(), jxr["response"].asDouble(),
                  jxr["cathode"].asDouble());
    }
}

WireCell::Configuration wcls::DriftPlanesSource::faces(int anode_ident) const
{
    auto it = m_faces.find(anode_ident);
    if (it == m_faces.end()) {
        raise<KeyError>("DriftPlanesSource: no faces for anode %d", anode_ident);
    }
    return it->second;
}

WireCell::Configuration wcls::DriftPlanesSource::xregions() const { return m_xregions; }
