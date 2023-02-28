/** A WCT component which is a source of frames which it produces
 * by also being an art::Event visitor.
 *
 * Lazy means that there is a delay between conversion of the short
 * int samples of the raw::RawDigit and the float samples of
 * IFrame/ITrace.  This can help memory usage if a subset of the
 * frame is processed serially.
 */

#ifndef LARWIRECELL_COMPONENTS_LAZYFRAMESOURCE
#define LARWIRECELL_COMPONENTS_LAZYFRAMESOURCE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameSource.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include "canvas/Utilities/InputTag.h"

#include <deque>
#include <string>
#include <vector>

namespace wcls {
  class LazyFrameSource : public IArtEventVisitor,
                          public WireCell::IFrameSource,
                          public WireCell::IConfigurable {
  public:
    LazyFrameSource();
    virtual ~LazyFrameSource();

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
