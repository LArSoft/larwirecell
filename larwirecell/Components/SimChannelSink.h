/** A WCT component that is a sink of sim::SimChannel
 * (converted from WireCell::IDepo) into an art::Event.
 *
 * Contact brussell@yale.edu for comments/questions.
 *
 * Modified by Wenqiang Gu (wgu@bnl.gov), 9/16/2019
 * A generic SimChannel saver for multiple volumes
 */

#ifndef LARWIRECELL_COMPONENTS_SIMCHANNELSINK
#define LARWIRECELL_COMPONENTS_SIMCHANNELSINK

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDepoFilter.h"
#include "WireCellIface/IRandom.h"
#include "WireCellUtil/Pimpos.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

namespace wcls {

  class SimChannelSink : public IArtEventVisitor,
                         public WireCell::IDepoFilter,
                         public WireCell::IConfigurable {

  public:
    SimChannelSink();
    // virtual ~SimChannelSink();

    /// IArtEventVisitor
    virtual void produces(art::ProducesCollector& collector);
    virtual void visit(art::Event& event);

    /// IDepoFilter
    virtual bool operator()(const WireCell::IDepo::pointer& indepo,
                            WireCell::IDepo::pointer& outdepo);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

  private:
    WireCell::IDepo::pointer m_depo;
    // WireCell::IAnodePlane::pointer m_anode;
    std::vector<WireCell::IAnodePlane::pointer> m_anodes; // multiple volumes
    WireCell::IRandom::pointer m_rng;

    std::map<unsigned int, sim::SimChannel> m_mapSC;

    void save_as_simchannel(const WireCell::IDepo::pointer& depo);

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

    // double Pi = 3.141592653589;
    // WireCell::Pimpos *uboone_u;
    // WireCell::Pimpos *uboone_v;
    // WireCell::Pimpos *uboone_y;
    // //WireCell::Pimpos *pimpos; // unused
  };
}

#endif
