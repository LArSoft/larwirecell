/** A WCT component which is a source of raw frames which it produces
 * by also being an art::Event visitor.
 *
 * Raw means that the waveforms are taken from the art::Event as a
 * labeled std::vector<raw::RawDigit> collection.
 */

#ifndef LARWIRECELL_COMPONENTS_RAWFRAMESOURCE
#define LARWIRECELL_COMPONENTS_RAWFRAMESOURCE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameSource.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include "canvas/Utilities/InputTag.h"

#include <deque>
#include <string>
#include <vector>

namespace wcls {
  class RawFrameSource : public IArtEventVisitor,
                         public WireCell::IFrameSource,
                         public WireCell::IConfigurable {
  public:
    RawFrameSource();
    virtual ~RawFrameSource();

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
