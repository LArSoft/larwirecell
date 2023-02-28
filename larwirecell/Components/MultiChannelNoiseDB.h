#ifndef LARWIRECELL_COMPONENTS_MULTICHANNELNOISEDB
#define LARWIRECELL_COMPONENTS_MULTICHANNELNOISEDB

#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IConfigurable.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

namespace wcls {

  class MultiChannelNoiseDB : public WireCell::IChannelNoiseDatabase,
                              public IArtEventVisitor,
                              public WireCell::IConfigurable {
  public:
    MultiChannelNoiseDB();
    virtual ~MultiChannelNoiseDB();

    /// IArtEventVisitor.
    //
    // Note: we don't actually poke at the event but use this
    // entry to refresh info from services in case they change
    // event-to-event.
    virtual void visit(art::Event& event);

    /// IConfigurable.
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& cfg);

    // IChannelNoiseDatabase

    // Big Fat Pimpl!
    virtual double sample_time() const { return m_pimpl->sample_time(); }
    virtual double nominal_baseline(int channel) const
    {
      return m_pimpl->nominal_baseline(channel);
    }
    virtual double gain_correction(int channel) const { return m_pimpl->gain_correction(channel); }
    virtual double response_offset(int channel) const { return m_pimpl->response_offset(channel); }
    virtual double min_rms_cut(int channel) const { return m_pimpl->min_rms_cut(channel); }
    virtual double max_rms_cut(int channel) const { return m_pimpl->max_rms_cut(channel); }
    virtual int pad_window_front(int channel) const { return m_pimpl->pad_window_front(channel); }
    virtual int pad_window_back(int channel) const { return m_pimpl->pad_window_back(channel); }
    virtual float coherent_nf_decon_limit(int channel) const
    {
      return m_pimpl->coherent_nf_decon_limit(channel);
    }
    virtual float coherent_nf_decon_lf_cutoff(int channel) const
    {
      return m_pimpl->coherent_nf_decon_lf_cutoff(channel);
    }
    virtual float coherent_nf_adc_limit(int channel) const
    {
      return m_pimpl->coherent_nf_adc_limit(channel);
    }
    virtual float coherent_nf_decon_limit1(int channel) const
    {
      return m_pimpl->coherent_nf_decon_limit1(channel);
    }
    virtual float coherent_nf_protection_factor(int channel) const
    {
      return m_pimpl->coherent_nf_protection_factor(channel);
    }
    virtual float coherent_nf_min_adc_limit(int channel) const
    {
      return m_pimpl->coherent_nf_min_adc_limit(channel);
    }
    virtual float coherent_nf_roi_min_max_ratio(int channel) const
    {
      return m_pimpl->coherent_nf_roi_min_max_ratio(channel);
    }
    virtual const filter_t& rcrc(int channel) const { return m_pimpl->rcrc(channel); }
    virtual const filter_t& config(int channel) const { return m_pimpl->config(channel); }
    virtual const filter_t& noise(int channel) const { return m_pimpl->noise(channel); }
    virtual const filter_t& response(int channel) const { return m_pimpl->response(channel); }
    virtual std::vector<channel_group_t> coherent_channels() const
    {
      return m_pimpl->coherent_channels();
    }
    virtual channel_group_t bad_channels() const { return m_pimpl->bad_channels(); }
    virtual channel_group_t miscfg_channels() const { return m_pimpl->miscfg_channels(); }

  private:
    // little helper struct
    struct SubDB {
      std::function<bool(art::Event& event)> check;
      WireCell::IChannelNoiseDatabase::pointer chndb;
      IArtEventVisitor::pointer visitor;
      SubDB(std::function<bool(art::Event& event)> f,
            WireCell::IChannelNoiseDatabase::pointer d,
            IArtEventVisitor::pointer v)
        : check(f), chndb(d), visitor(v)
      {}
    };
    typedef std::vector<SubDB> rulelist_t;
    rulelist_t m_rules;
    WireCell::IChannelNoiseDatabase::pointer m_pimpl;
    IArtEventVisitor::pointer m_pimpl_visitor;
  };

} // wcls

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
