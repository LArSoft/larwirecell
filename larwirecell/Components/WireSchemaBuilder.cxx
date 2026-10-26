#include "larwirecell/Components/WireSchemaBuilder.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"

#include <map>

using namespace WireCell;
using namespace WireCell::WireSchema;

namespace {

    // Convert a LArSoft point (centimeters, world frame) into a WCT Point
    // (WCT base length unit is mm; units::cm == 10).
    inline WireCell::Point to_wct(const geo::Point_t& p)
    {
        return WireCell::Point(p.X() * units::cm, p.Y() * units::cm, p.Z() * units::cm);
    }

    // Precomputed mapping from LArSoft structural IDs to WCT (anode,face) idents
    // according to the grouping policy.  Channel/plane/wire idents never depend
    // on the policy (they come straight from the channel map), only the anode
    // and face organisation does.
    struct Grouper {
        wcls::Grouping grouping;

        // kTPCset: TPCsetID -> anode ident ; TPCID -> face index
        std::map<readout::TPCsetID, int> tpcset_anode;
        std::map<geo::TPCID, int> tpc_face;

        // kPerTPC / kPairTPC: TPCID -> running global TPC index
        std::map<geo::TPCID, int> tpc_index;

        Grouper(const geo::GeometryCore& geom,
                const geo::WireReadoutGeom& wires,
                wcls::Grouping g)
            : grouping(g)
        {
            // Visit (cryostat, TPCset, TPC) in LArSoft order to assign dense,
            // reproducible running idents.
            int anode_running = 0;
            int tpc_running = 0;
            for (unsigned int c = 0; c < geom.Ncryostats(); ++c) {
                readout::CryostatID cid(c);
                const unsigned int nsets = wires.NTPCsets(cid);
                for (unsigned int s = 0; s < nsets; ++s) {
                    readout::TPCsetID sid(cid, s);
                    if (!wires.HasTPCset(sid)) continue;
                    tpcset_anode[sid] = anode_running++;
                    const auto tpcs = wires.TPCsetToTPCs(sid);
                    int face = 0;
                    for (const auto& tid : tpcs) {
                        tpc_face[tid] = face++;
                        tpc_index[tid] = tpc_running++;
                    }
                }
            }
        }

        // Return (anode_ident, face_ident) for a wire's TPC.
        std::pair<int, int> anode_face(const geo::WireReadoutGeom& wires,
                                       const geo::WireID& wid) const
        {
            const geo::TPCID tid = wid.asTPCID();
            if (grouping == wcls::Grouping::kTPCset) {
                const readout::TPCsetID sid = wires.TPCtoTPCset(tid);
                return {tpcset_anode.at(sid), tpc_face.at(tid)};
            }
            const int gtpc = tpc_index.at(tid);
            if (grouping == wcls::Grouping::kPerTPC) {
                return {gtpc, 0};
            }
            // kPairTPC: legacy wrapped two-face convention.
            return {gtpc / 2, (1 + gtpc) % 2};
        }
    };

}  // namespace

wcls::Grouping wcls::grouping_from_string(const std::string& name)
{
    if (name == "tpcset")   return wcls::Grouping::kTPCset;
    if (name == "per-tpc")  return wcls::Grouping::kPerTPC;
    if (name == "pair-tpc") return wcls::Grouping::kPairTPC;
    THROW(WireCell::ValueError()
          << WireCell::errmsg{"unknown wire-schema grouping policy: '" + name
                              + "' (want tpcset|per-tpc|pair-tpc)"});
}

WireCell::WireSchema::Store
wcls::build_wire_schema(const geo::GeometryCore& geom,
                        const geo::WireReadoutGeom& wires,
                        const wcls::GroupingPolicy& policy)
{
    auto log = Log::logger("wire-schema");

    Grouper grouper(geom, wires, policy.grouping);

    auto db = std::make_shared<StoreDB>();
    const int det_ident = 0;

    // Iterate every readout channel.  ChannelToWire() yields one WireID per
    // physical wire segment that feeds that channel: a single entry for an
    // un-wrapped (collection or independent-CRU) channel, and multiple entries
    // -- one per TPC/face segment -- for a wrapped induction channel.  The
    // position in that list is the segment number (0 == the segment nearest the
    // electronics, per the channel map's own ordering).
    const raw::ChannelID_t nchan = wires.Nchannels();
    size_t nwires = 0;
    for (raw::ChannelID_t ch = 0; ch < nchan; ++ch) {
        if (!wires.HasChannel(ch)) continue;
        const auto segwires = wires.ChannelToWire(ch);
        int segment = 0;
        for (const geo::WireID& wid : segwires) {
            const auto [aid, fid] = grouper.anode_face(wires, wid);
            const int pid = wid.Plane;          // 0=U,1=V,2=collection
            const int wire_ident = wid.Wire;    // LArSoft wire-in-plane index

            const geo::WireGeo& wg = wires.Wire(wid);

            Wire w;
            w.ident = wire_ident;
            w.channel = static_cast<int>(ch);
            w.segment = segment;
            w.tail = to_wct(wg.GetStart());
            w.head = to_wct(wg.GetEnd());

            Plane& plane = get_append(*db, pid, fid, aid, det_ident);
            db->wires.push_back(w);
            plane.wires.push_back(db->wires.size() - 1);

            ++segment;
            ++nwires;
        }
    }

    log->debug("wire-schema: built {} anodes, {} faces, {} planes, {} wires from {} channels",
               db->anodes.size(), db->faces.size(), db->planes.size(), nwires, nchan);

    // Make the raw store compliant with the wire-schema conventions "by
    // construction" using the exact same fixers wcwires applies on load.
    WireCell::WireSchema::correct(*db, policy.correction);

    Store store(db);

    if (policy.validate) {
        // Throws WireCell::ValueError on a convention violation.
        WireCell::WireSchema::validate(store, policy.repsilon, false);
        log->debug("wire-schema: validated OK");
    }

    return store;
}
