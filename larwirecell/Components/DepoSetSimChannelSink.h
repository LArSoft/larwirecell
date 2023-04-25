/**
 * Same as SimChannelSink but takes in DepoSet other than Depo
 */

#ifndef LARWIRECELL_COMPONENTS_DEPOSETSIMCHANNELSINK
#define LARWIRECELL_COMPONENTS_DEPOSETSIMCHANNELSINK

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDepoSetFilter.h"
#include "WireCellIface/IRandom.h"
#include "WireCellUtil/Pimpos.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

namespace wcls {

  class DepoSetSimChannelSink : public IArtEventVisitor,
                                public WireCell::IDepoSetFilter,
                                public WireCell::IConfigurable {

  public:
    /// IArtEventVisitor
    virtual void produces(art::ProducesCollector& collector);
    virtual void visit(art::Event& event);

    /// IDepoSetFilter
    virtual bool operator()(const WireCell::IDepoSet::pointer& indepos,
                            WireCell::IDepoSet::pointer& outdepos);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

  private:
    // WireCell::IDepoSet::pointer m_depo;
    // WireCell::IAnodePlane::pointer m_anode;
    std::vector<WireCell::IAnodePlane::pointer> m_anodes; // multiple volumes
    WireCell::IRandom::pointer m_rng;

    std::map<unsigned int, sim::SimChannel> m_mapSC;

    void save_as_simchannel(const WireCell::IDepo::pointer& depos);

    std::string m_artlabel;
    double m_readout_time;
    double m_tick;
    double m_start_time;
    double m_nsigma;
    double m_drift_speed;
    double m_u_to_rp;
    double m_v_to_rp;
    double m_y_to_rp;
    double m_u_time_offset;
    double m_v_time_offset;
    double m_y_time_offset;
    double m_g4_ref_time;
    bool m_use_energy;
    bool m_use_extra_sigma; // extra smearing from signal processing
  };
}

#endif
