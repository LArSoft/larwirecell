/** Build a Wire-Cell Toolkit wire-schema Store directly from the LArSoft
 *  geometry + wire-readout (channel map) services.
 *
 *  This is the single, shared "dump transform" used by both:
 *
 *   - the WireSchemaFileWriter_module (Option 1): an art job that writes a
 *     WCT "wires file" (.json.bz2) as a precursor step; and
 *
 *   - the WireSchemaService / wclsWireSchemaSource component (Option 3): an
 *     in-process IWireSchema provider that hands WCT its geometry on the fly.
 *
 *  Because the anode/face/plane grouping, the channel numbers and the wire
 *  segmentation are all read from the SAME `geo::WireReadoutGeom` the rest of
 *  the art job uses, the resulting Store is consistent with the LArSoft channel
 *  map by construction -- regardless of which concrete mapper (WireReadoutCRUGeom,
 *  CRPWireReadoutGeom, ColdBoxWireReadoutGeom, ...) is installed.
 *
 *  The logic here replaces, and unifies, the several per-detector transforms in
 *  wire-cell-python (wirecell.util.wires.{multitpc,dunevd,gdml}).  The one piece
 *  that genuinely varies between detectors -- how LArSoft TPCs are grouped into
 *  WCT anodes/faces -- is captured by GroupingPolicy below rather than by forked
 *  code paths.
 *
 *  See util/docs/wire-schema.org for the "legal" conventions the output obeys
 *  (ordering, wire/pitch direction, segments).  Compliance is achieved by
 *  running WireCell::WireSchema::correct() and asserted by validate().
 */
#ifndef LARWIRECELL_COMPONENTS_WIRESCHEMABUILDER_H
#define LARWIRECELL_COMPONENTS_WIRESCHEMABUILDER_H

#include "WireCellUtil/WireSchema.h"

#include <map>
#include <string>
#include <vector>

namespace geo {
    class GeometryCore;
    class WireReadoutGeom;
}

namespace wcls {

    /// How LArSoft (cryostat, TPC) elements map onto WCT (anode, face).
    ///
    /// In all cases the WCT plane ident is the LArSoft plane index within its
    /// TPC (0,1,2 == U,V,collection) and the WCT channel is exactly the LArSoft
    /// offline channel, so the channel numbering is identical to the installed
    /// channel map for every policy.  The policy only decides the anode/face
    /// *organisation* of those wires.
    enum class Grouping {
        /// Anode == LArSoft TPCset (ie one CRP); its faces are the TPCs of that
        /// set in LArSoft order.  This is the recommended, mapper-faithful
        /// choice: for WireReadoutCRUGeom (TPCset==TPC) it yields single-face
        /// anodes with no wrapping; for CRP/ColdBox mappers it yields
        /// multi-face anodes whose shared induction channels wrap across faces.
        kTPCset,

        /// Anode == LArSoft TPC, one face per anode (face ident 0).  Matches the
        /// wire-cell-python "dunevd 3view" (single-face) convention.
        kPerTPC,

        /// Anode == TPC/2, face == (1+TPC)%2.  Reproduces the legacy
        /// wire-cell-python "multitpc"/"dunevd" wrapped two-face convention
        /// (eg the shipped dunevd10kt_3view_30deg_v6 file).  Provided for
        /// backward compatibility; prefer kTPCset for new files.
        kPairTPC,
    };

    struct GroupingPolicy {
        Grouping grouping{Grouping::kTPCset};

        /// Correction level applied to the built store (see WireSchema::Correction).
        /// Default "pitch" == the full set wcwires applies.
        WireCell::WireSchema::Correction correction{WireCell::WireSchema::Correction::pitch};

        /// If true, throw if the corrected store fails validate().
        bool validate{true};

        /// Relative epsilon passed to validate().
        double repsilon{1e-6};
    };

    /// Parse a grouping policy name, one of "tpcset", "per-tpc", "pair-tpc".
    /// Throws WireCell::ValueError on an unknown name.
    Grouping grouping_from_string(const std::string& name);

    /// Build and return a corrected, validated WCT wire-schema Store from the
    /// given LArSoft geometry and wire-readout.  Throws WireCell::ValueError on
    /// a validation failure (when policy.validate is true).
    WireCell::WireSchema::Store
    build_wire_schema(const geo::GeometryCore& geom,
                      const geo::WireReadoutGeom& wires,
                      const GroupingPolicy& policy = {});

    /// The drift-axis geometry of one WCT anode face, from its LArSoft TPC.
    /// Locations are X coordinates in WCT units.
    struct FaceDrift {
        double collection{0};  ///< center of the collection plane
        double anode{0};       ///< anode bound, built coplanar with collection
        double cathode{0};     ///< center of the active volume face opposite the planes
        int drift_sign{0};     ///< +1 (-1) if electrons drift toward +X (-X)
    };

    /// Anode ident -> one FaceDrift per face, in the face order of the Store.
    using DriftFaces = std::map<int, std::vector<FaceDrift>>;

    /// Return the drift geometry of the faces of a store made by
    /// build_wire_schema() from the same geometry and grouping.  Throws
    /// WireCell::ValueError if a TPC does not drift along X.
    DriftFaces build_drift_faces(const geo::GeometryCore& geom,
                                 const geo::WireReadoutGeom& wires,
                                 const WireCell::WireSchema::Store& store,
                                 Grouping grouping);

}  // namespace wcls

#endif
