#ifndef LARWIRECELL_COMPONENTS_CHANNELSELECTORDB
#define LARWIRECELL_COMPONENTS_CHANNELSELECTORDB

#include "WireCellIface/IChannelNoiseDatabase.h"
#include "larwirecell/Components/ChannelNoiseDB.h"

namespace wcls {

  class ChannelSelectorDB : public ChannelNoiseDB {
  public:
    ChannelSelectorDB();
    virtual ~ChannelSelectorDB();

    /// IArtEventVisitor.
    // Note: we don't actually poke at the event but use this
    // entry to refresh info from services in case they change
    // event-to-event.
    virtual void visit(art::Event& event);

    /// IConfigurable.
    virtual void configure(const WireCell::Configuration& config);

    virtual channel_group_t bad_channels() const { return m_bad_channels; }
    virtual channel_group_t miscfg_channels() const { return m_miscfg_channels; }

  private:
    std::string m_type;
    channel_group_t m_bad_channels;
    channel_group_t m_miscfg_channels;
  };

} // wcls

#endif /* LARWIRECELL_COMPONENTS_CHANNELSELECTORDB */

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
