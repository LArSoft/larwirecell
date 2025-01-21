/** A WCT component which is a source of "cooked" frames which it
 * produces by also being an art::Event visitor.
 *
 * Cooked means that the waveforms are taken from the art::Event as a
 * labeled std::vector<recob::Wire> collection.
 */

#ifndef LARWIRECELL_COMPONENTS_COOKEDFRAMESOURCE
#define LARWIRECELL_COMPONENTS_COOKEDFRAMESOURCE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameSource.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include "canvas/Utilities/InputTag.h"

#include <deque>
#include <string>
#include <vector>

namespace wcls {
  class CookedFrameSource : public IArtEventVisitor,
                            public WireCell::IFrameSource,
                            public WireCell::IConfigurable {
  public:
    CookedFrameSource();
    virtual ~CookedFrameSource();

    /// IArtEventVisitor
    virtual void visit(art::Event& event);

    /// IFrameSource
    virtual bool operator()(WireCell::IFrame::pointer& frame);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

  private:
    std::deque<WireCell::IFrame::pointer> m_frames;
    art::InputTag m_inputTag;
    double m_tick;
    int m_nticks;
    std::vector<std::string> m_frame_tags;
  };

}

#endif
