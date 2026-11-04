/** DriftPlanesSource (WCT forwarder)
 *
 *  A WCT component implementing IDriftPlanes from the face drift geometry of
 *  the art-owned wcls::WireSchemaService, in the same way WireSchemaSource
 *  provides IWireSchema.  Name it in the "drift_planes" parameter of each
 *  AnodePlane and of the Drifter (and list it in their "uses"):
 *
 *     { type: "wclsDriftPlanesSource", name: "",
 *       data: { response_distance: 10*wc.cm } }
 *
 *  Per face, from the live LArSoft geometry:
 *
 *  - anode :: the collection plane, displaced by the service's anode_offset.
 *  - response :: response_distance from the collection plane, against the
 *    drift.  It must match the "origin" of the field response in use, which
 *    is NOT read here and must be transcribed by hand.
 *  - cathode :: the face of the TPC active volume opposite the planes.
 *
 *  The xregions are the faces with distinct planes, within "tolerance".
 */
#ifndef LARWIRECELL_COMPONENTS_DRIFTPLANESSOURCE_H
#define LARWIRECELL_COMPONENTS_DRIFTPLANESSOURCE_H

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDriftPlanes.h"

#include <map>

namespace wcls {

    class DriftPlanesSource : public WireCell::IDriftPlanes,
                              public WireCell::IConfigurable {
      public:
        DriftPlanesSource();
        virtual ~DriftPlanesSource();

        // IDriftPlanes
        virtual WireCell::Configuration faces(int anode_ident) const;
        virtual WireCell::Configuration xregions() const;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

      private:
        std::map<int, WireCell::Configuration> m_faces;
        WireCell::Configuration m_xregions;
    };

}  // namespace wcls

#endif
