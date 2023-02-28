#ifndef LARWIRECELL_COMPONENTS_CHANNELNOISEDB
#define LARWIRECELL_COMPONENTS_CHANNELNOISEDB

#include "WireCellSigProc/OmniChannelNoiseDB.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

namespace wcls {

  class ChannelNoiseDB : public WireCell::SigProc::OmniChannelNoiseDB, public IArtEventVisitor {
  public:
    ChannelNoiseDB();
    virtual ~ChannelNoiseDB();

    /// IArtEventVisitor.
    //
    // Note: we don't actually poke at the event but use this
    // entry to refresh info from services in case they change
    // event-to-event.
    virtual void visit(art::Event& event);

    /// IConfigurable.
    //
    // Defer default to parent.  By default this class does not
    // overriding.

    // Defer to parent but override if asked to.
    virtual void configure(const WireCell::Configuration& config);

    /// All IChannelNoiseDatabase interface is defered to parent.

  private:
    // The policy for overriding parent class information.
    enum OverridePolicy_t {
      kNothing = 0, // this class leaves info untouched.
      kReplace,     // this class replaces all parent info with what is given by service.
      kUnion        // this class adds to parent only for channels given by service.
    };

    OverridePolicy_t parse_policy(const WireCell::Configuration& jpol);

    OverridePolicy_t m_bad_channel_policy;
    OverridePolicy_t m_misconfig_channel_policy;
    double m_fgstgs[4];
  };

} // wcls

#endif /* LARWIRECELL_COMPONENTS_CHANNELNOISEDB */

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
