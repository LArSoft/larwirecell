/** WireSchemaService (Option 3, data owner)
 *
 *  An art service that builds -- on first use -- a Wire-Cell Toolkit wire-schema
 *  Store from the LArSoft geometry + wire-readout (channel map) services of the
 *  running job, and hands it out via wire_schema_store().
 *
 *  It is paired with the trivial WCT component wcls::WireSchemaSource
 *  (larwirecell/Components/WireSchemaSource) which implements WCT's IWireSchema
 *  by forwarding to this service.  The split exists because art and WCT each
 *  insist on constructing and owning their own objects: art owns this (heavy)
 *  service; WCT owns the (trivial) forwarder it instantiates from Jsonnet.
 *
 *  The resulting geometry is consistent-by-construction with whatever channel
 *  map the job loaded -- no wires file, no file/map mismatch possible.
 *
 *  It also provides the drift-axis geometry of each anode face (collection
 *  plane, cathode, drift direction), which the WCT component
 *  wcls::DriftPlanesSource turns into the "anode", "response" and "cathode"
 *  planes of AnodePlane and Drifter.
 *
 *  FHiCL (services block):
 *    WireSchemaService: {
 *      grouping:    "tpcset"   # tpcset (default) | per-tpc | pair-tpc
 *      correction:  4          # 1=load 2=order 3=direction 4=pitch
 *      validate:    true
 *      anode_offset: 0.0       # cm from the collection plane, positive
 *                              # against the drift (into the drift volume)
 *    }
 */
#ifndef LARWIRECELL_SERVICES_WIRESCHEMASERVICE_H
#define LARWIRECELL_SERVICES_WIRESCHEMASERVICE_H

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "larwirecell/Components/WireSchemaBuilder.h"

#include "WireCellUtil/WireSchema.h"

#include <mutex>
#include <string>

namespace art {
    class ActivityRegistry;
}

namespace wcls {

    class WireSchemaService {
      public:
        WireSchemaService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

        /// Return the (lazily built) wire-schema Store.  Safe to call once the
        /// geometry/wire-readout art services are available (ie from begin-job
        /// onward), which is the case when WCT configures its AnodePlane.
        const WireCell::WireSchema::Store& wire_schema_store() const;

        /// Return the (lazily built) drift geometry of the faces of the store,
        /// with the anode bound displaced from the collection plane by the
        /// configured anode_offset.
        const DriftFaces& drift_faces() const;

      private:
        void build() const;

        std::string fGrouping;
        int fCorrection;
        bool fValidate;
        double fAnodeOffset;  // WCT units

        mutable WireCell::WireSchema::Store m_store;
        mutable DriftFaces m_drift_faces;
        mutable std::once_flag m_built;
    };

}  // namespace wcls

// SHARED so a WCT-only job can run art (and so TbbFlow) with many threads; the
// one lazy build is guarded by std::call_once.
DECLARE_ART_SERVICE(wcls::WireSchemaService, SHARED)

#endif
