/** A WCT component which is a sink of "cooked" frames which saves
 * frames into an art::Event that it visits.
 *
 * Cooked means that some processing of the frame has occurred and the
 * frame is saved into the art::Event as recob::Wires.
 */

#ifndef LARWIRECELL_COMPONENTS_COOKEDFRAMESINK
#define LARWIRECELL_COMPONENTS_COOKEDFRAMESINK

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameSink.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include <string>
#include <vector>

namespace wcls {
  class CookedFrameSink : public IArtEventVisitor,
                          public WireCell::IFrameSink,
                          public WireCell::IConfigurable {
  public:
    CookedFrameSink();
    virtual ~CookedFrameSink();

    /// IArtEventVisitor
    virtual void produces(art::ProducesCollector& collector);
    virtual void visit(art::Event& event);

    /// IFrameSink
    virtual bool operator()(const WireCell::IFrame::pointer& frame);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

  private:
    WireCell::IFrame::pointer m_frame;
    WireCell::IAnodePlane::pointer m_anode;
    std::vector<std::string> m_frame_tags;
    int m_nticks;
  };
}

#endif
