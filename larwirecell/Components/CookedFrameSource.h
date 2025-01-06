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
#include "WireCellUtil/Logging.h"
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
    double m_tick;
    int m_nticks;
    double m_scale{50}; // scale up input recob::Wire by this factor
    std::vector<std::string> m_frame_tags;
    std::vector<std::string> m_recobwire_tags;
    std::vector<std::string> m_trace_tags;
    std::vector<std::string> m_summary_tags;
    std::vector<std::string> m_input_mask_tags;
    std::vector<std::string> m_output_mask_tags;
    WireCell::Log::logptr_t l;
  };

}

#endif
