/** WireSchemaSource (Option 3, WCT forwarder)
 *
 *  A trivial, default-constructible WCT component that implements WireCell's
 *  IWireSchema by forwarding to the art-owned wcls::WireSchemaService.  It is a
 *  drop-in replacement for Gen::WireSchemaFile: in Jsonnet declare it like
 *
 *     { type: "wclsWireSchemaSource", name: "", data: {} }
 *
 *  and point AnodePlane.data.wire_schema at "wclsWireSchemaSource:" (with that
 *  object in AnodePlane's "uses" list).  No wires file is involved; the geometry
 *  comes from the live LArSoft channel map via the service.
 */
#ifndef LARWIRECELL_COMPONENTS_WIRESCHEMASOURCE_H
#define LARWIRECELL_COMPONENTS_WIRESCHEMASOURCE_H

#include "WireCellIface/IWireSchema.h"
#include "WireCellIface/IConfigurable.h"

namespace wcls {

    class WireSchemaSource : public WireCell::IWireSchema,
                             public WireCell::IConfigurable {
      public:
        WireSchemaSource();
        virtual ~WireSchemaSource();

        // IWireSchema
        virtual const WireCell::WireSchema::Store& wire_schema_store() const;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

      private:
        WireCell::WireSchema::Store m_store;
    };

}  // namespace wcls

#endif
