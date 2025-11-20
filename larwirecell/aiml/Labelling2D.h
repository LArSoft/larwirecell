#ifndef WIRECELL_AIML_LABELLING2D
#define WIRECELL_AIML_LABELLING2D

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"
#include "lardataobj/Simulation/SimChannel.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace WireCell::AIML {
  class Labelling2D : public Aux::Logger,
                      public wcls::IArtEventVisitor,
                      public IFrameFilter,
                      public IConfigurable {
  public:
    Labelling2D();
    ~Labelling2D() override;

    // IFrameFilter
    bool operator()(const input_pointer& in, output_pointer& out) override;

    // IConfigurable
    void configure(const WireCell::Configuration& config) override;
    WireCell::Configuration default_configuration() const override;

    // IArtEventVisitor
    void visit(art::Event& event) override;

  private:
    using channel_index_t = std::unordered_map<unsigned int, std::size_t>;
    using TrackChargeInfo = std::pair<int, double>; // {trackID, charge}
    using PIDChargeInfo = std::pair<int, double>;   // {PID, merged_charge}

    const sim::SimChannel* find_simchannel(unsigned int channel) const;
    int select_track_id(const sim::SimChannel& sc, int tdc_begin, int tdc_end) const;
    void cache_simchannels(const std::vector<sim::SimChannel>& simchs);
    void clear_cache();
    void populate_trackid_pid_map();
    int pid_from_track(int track_id) const;

    // New functions for charge-based labeling with rebin
    std::vector<TrackChargeInfo> extract_track_charges(const sim::SimChannel& sc,
                                                       int tdc_begin, int tdc_end) const;
    std::pair<int, int> select_top2_track_ids(const sim::SimChannel& sc,
                                              int tdc_begin, int tdc_end) const;
    std::pair<int, int> select_top2_pids(const sim::SimChannel& sc,
                                         int tdc_begin, int tdc_end) const;

    WireCell::IAnodePlane::pointer m_anode;
    std::string m_anode_tn;
    std::string m_reco_tag;
    std::string m_output_trace_tag_trackid;
    std::string m_output_trace_tag_pid;
    std::string m_output_trace_tag_trackid_1st;
    std::string m_output_trace_tag_trackid_2nd;
    std::string m_output_trace_tag_pid_1st;
    std::string m_output_trace_tag_pid_2nd;
    std::vector<std::string> m_frame_tags;
    std::string m_simchannel_label;
    int m_default_label;
    int m_tdc_offset;
    double m_min_charge;
    bool m_copy_input_traces;
    int m_rebin_time_tick;

    std::vector<sim::SimChannel> m_simchannels;
    channel_index_t m_channel_index;
    std::unordered_map<int, int> m_trackid_to_pid;
  };
}

#endif
