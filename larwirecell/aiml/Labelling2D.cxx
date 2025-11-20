#include "Labelling2D.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <algorithm>
#include <cmath>

WIRECELL_FACTORY(Labelling2D,
                 WireCell::AIML::Labelling2D,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

using namespace WireCell;

using WireCell::Aux::SimpleFrame;
using WireCell::Aux::SimpleTrace;

namespace {
  template <typename T>
  inline T abs_charge(T value)
  {
    return std::fabs(value);
  }
} // namespace

AIML::Labelling2D::Labelling2D()
  : Aux::Logger("Labelling2D", "aiml")
  , m_anode(nullptr)
  , m_reco_tag("reco")
  , m_output_trace_tag_trackid("trackid")
  , m_output_trace_tag_pid("pid")
  , m_output_trace_tag_trackid_1st("trackid_1st")
  , m_output_trace_tag_trackid_2nd("trackid_2nd")
  , m_output_trace_tag_pid_1st("pid_1st")
  , m_output_trace_tag_pid_2nd("pid_2nd")
  , m_default_label(0)
  , m_tdc_offset(0)
  , m_min_charge(0.0)
  , m_copy_input_traces(false)
  , m_rebin_time_tick(1)
{
  m_frame_tags.push_back("truth");
}

AIML::Labelling2D::~Labelling2D() = default;

Configuration AIML::Labelling2D::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = m_anode_tn;
  cfg["reco_tag"] = m_reco_tag;
  cfg["output_trace_tag_trackid"] = m_output_trace_tag_trackid;
  cfg["output_trace_tag_pid"] = m_output_trace_tag_pid;
  cfg["output_trace_tag_trackid_1st"] = m_output_trace_tag_trackid_1st;
  cfg["output_trace_tag_trackid_2nd"] = m_output_trace_tag_trackid_2nd;
  cfg["output_trace_tag_pid_1st"] = m_output_trace_tag_pid_1st;
  cfg["output_trace_tag_pid_2nd"] = m_output_trace_tag_pid_2nd;
  cfg["default_label"] = m_default_label;
  cfg["tdc_offset"] = m_tdc_offset;
  cfg["min_charge"] = m_min_charge;
  cfg["copy_input_traces"] = m_copy_input_traces;
  cfg["rebin_time_tick"] = m_rebin_time_tick;
  cfg["simchannel"]["label"] = m_simchannel_label;
  cfg["frame_tags"] = Json::arrayValue;
  for (auto const& tag : m_frame_tags) {
    cfg["frame_tags"].append(tag);
  }
  return cfg;
}

void AIML::Labelling2D::configure(const Configuration& cfg)
{
  m_anode_tn = get(cfg, "anode", m_anode_tn);
  if (!m_anode_tn.empty()) {
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
  }

  m_reco_tag = get(cfg, "reco_tag", m_reco_tag);
  m_output_trace_tag_trackid = get(cfg, "output_trace_tag_trackid", m_output_trace_tag_trackid);
  m_output_trace_tag_pid = get(cfg, "output_trace_tag_pid", m_output_trace_tag_pid);
  m_output_trace_tag_trackid_1st = get(cfg, "output_trace_tag_trackid_1st", m_output_trace_tag_trackid_1st);
  m_output_trace_tag_trackid_2nd = get(cfg, "output_trace_tag_trackid_2nd", m_output_trace_tag_trackid_2nd);
  m_output_trace_tag_pid_1st = get(cfg, "output_trace_tag_pid_1st", m_output_trace_tag_pid_1st);
  m_output_trace_tag_pid_2nd = get(cfg, "output_trace_tag_pid_2nd", m_output_trace_tag_pid_2nd);
  const int configured_default = get(cfg, "default_label", m_default_label);
  if (configured_default != 0) {
    log->warn("Labelling2D overrides configured default_label {} with 0", configured_default);
  }
  m_default_label = 0;
  m_tdc_offset = get(cfg, "tdc_offset", m_tdc_offset);
  m_min_charge = get(cfg, "min_charge", m_min_charge);
  m_copy_input_traces = get(cfg, "copy_input_traces", m_copy_input_traces);
  m_rebin_time_tick = get(cfg, "rebin_time_tick", m_rebin_time_tick);
  if (m_rebin_time_tick <= 0) {
    log->warn("Labelling2D rebin_time_tick must be > 0, setting to 1");
    m_rebin_time_tick = 1;
  }

  m_simchannel_label = get(cfg, "simchannel_label", m_simchannel_label);

  if (cfg.isMember("frame_tags")) {
    m_frame_tags.clear();
    for (auto const& tag : cfg["frame_tags"]) {
      m_frame_tags.push_back(tag.asString());
    }
  }

  clear_cache();
}

void AIML::Labelling2D::visit(art::Event& event)
{
  log->debug("Labelling2D::visit: event: {}", event.event());
  clear_cache();

  if (m_simchannel_label.empty()) {
    log->debug("SimChannel label not configured; skipping visit");
    return;
  }

  art::Handle<std::vector<sim::SimChannel>> handle;
  if (!event.getByLabel(art::InputTag{m_simchannel_label}, handle)) {
    log->warn("Labelling2D failed to fetch SimChannel with label '{}'", m_simchannel_label);
    return;
  }

  cache_simchannels(*handle);
  populate_trackid_pid_map();
  log->debug("Labelling2D cached {} SimChannels and {} track->pid entries",
             m_channel_index.size(),
             m_trackid_to_pid.size());
}

bool AIML::Labelling2D::operator()(const input_pointer& in, output_pointer& out)
{
  out = nullptr;

  if (!in) {
    return true;
  }

  auto reco_traces = Aux::tagged_traces(in, m_reco_tag);
  if (reco_traces.empty()) {
    log->warn("Labelling2D: no traces tagged '{}' in frame {}", m_reco_tag, in->ident());
    out = in;
    return true;
  }

  auto traces_buffer = std::make_shared<ITrace::vector>();
  if (m_copy_input_traces) {
    auto in_traces = in->traces();
    if (in_traces) {
      traces_buffer->insert(traces_buffer->end(), in_traces->begin(), in_traces->end());
    }
  }

  IFrame::trace_list_t trackid_indices;
  IFrame::trace_list_t pid_indices;
  IFrame::trace_list_t trackid_1st_indices;
  IFrame::trace_list_t trackid_2nd_indices;
  IFrame::trace_list_t pid_1st_indices;
  IFrame::trace_list_t pid_2nd_indices;

  trackid_indices.reserve(reco_traces.size());
  pid_indices.reserve(reco_traces.size());
  trackid_1st_indices.reserve(reco_traces.size());
  trackid_2nd_indices.reserve(reco_traces.size());
  pid_1st_indices.reserve(reco_traces.size());
  pid_2nd_indices.reserve(reco_traces.size());

  for (auto const& trace : reco_traces) {
    if (!trace) {
      continue;
    }

    const int chid = trace->channel();
    const sim::SimChannel* sc = find_simchannel(chid);
    // if (!sc) {
    //   log->debug("Labelling2D: no SimChannel found for channel {}", chid);
    // } else {
    //   const auto & idemap = sc->TDCIDEMap();
    //   log->debug("Labelling2D: found SimChannel for channel {} with {} TDC entries",
    //              chid,
    //              idemap.size());
    // }

    const auto& reco_charge = trace->charge();

    // Original traces at sample level
    SimpleTrace* label_trace = new SimpleTrace(chid, trace->tbin(), reco_charge.size());
    SimpleTrace* pid_trace = new SimpleTrace(chid, trace->tbin(), reco_charge.size());
    auto& label_values = label_trace->charge();
    auto& pid_values = pid_trace->charge();
    std::fill(label_values.begin(), label_values.end(), static_cast<float>(m_default_label));
    std::fill(pid_values.begin(), pid_values.end(), 0.0f);

    // Rebinned traces
    std::size_t rebinned_size = (reco_charge.size() + m_rebin_time_tick - 1) / m_rebin_time_tick;
    SimpleTrace* trackid_1st_trace = new SimpleTrace(chid, trace->tbin(), rebinned_size);
    SimpleTrace* trackid_2nd_trace = new SimpleTrace(chid, trace->tbin(), rebinned_size);
    SimpleTrace* pid_1st_trace = new SimpleTrace(chid, trace->tbin(), rebinned_size);
    SimpleTrace* pid_2nd_trace = new SimpleTrace(chid, trace->tbin(), rebinned_size);

    auto& trackid_1st_values = trackid_1st_trace->charge();
    auto& trackid_2nd_values = trackid_2nd_trace->charge();
    auto& pid_1st_values = pid_1st_trace->charge();
    auto& pid_2nd_values = pid_2nd_trace->charge();

    std::fill(trackid_1st_values.begin(), trackid_1st_values.end(), static_cast<float>(m_default_label));
    std::fill(trackid_2nd_values.begin(), trackid_2nd_values.end(), static_cast<float>(m_default_label));
    std::fill(pid_1st_values.begin(), pid_1st_values.end(), 0.0f);
    std::fill(pid_2nd_values.begin(), pid_2nd_values.end(), 0.0f);

    if (sc) {
      const int base_tbin = trace->tbin();

      // Process original samples
      for (std::size_t isample = 0; isample < reco_charge.size(); ++isample) {
        if (!(abs_charge(reco_charge[isample]) > m_min_charge)) {
          continue;
        }
        const int tdc_begin = m_tdc_offset + base_tbin + static_cast<int>(isample);
        const int track_id = select_track_id(*sc, tdc_begin, tdc_begin + 1);
        if (track_id == m_default_label) {
          continue;
        }
        label_values[isample] = static_cast<float>(track_id);
        pid_values[isample] = static_cast<float>(pid_from_track(track_id));
      }

      // Process rebinned samples
      for (std::size_t ibin = 0; ibin < rebinned_size; ++ibin) {
        std::size_t start_sample = ibin * m_rebin_time_tick;
        std::size_t end_sample = std::min((ibin + 1) * m_rebin_time_tick, reco_charge.size());

        // Check if this rebinned bin has any significant charge
        bool has_charge = false;
        for (std::size_t isample = start_sample; isample < end_sample; ++isample) {
          if (abs_charge(reco_charge[isample]) > m_min_charge) {
            has_charge = true;
            break;
          }
        }

        if (!has_charge) {
          continue;
        }

        // For the rebinned bin, use the span of TDCs
        const int tdc_begin = m_tdc_offset + base_tbin + static_cast<int>(start_sample);
        const int tdc_end = m_tdc_offset + base_tbin + static_cast<int>(end_sample);

        auto top2_tracks = select_top2_track_ids(*sc, tdc_begin, tdc_end);
        int track_id_1st = top2_tracks.first;
        int track_id_2nd = top2_tracks.second;

        auto top2_pids = select_top2_pids(*sc, tdc_begin, tdc_end);
        int pid_1st = top2_pids.first;
        int pid_2nd = top2_pids.second;

        trackid_1st_values[ibin] = static_cast<float>(track_id_1st);
        trackid_2nd_values[ibin] = static_cast<float>(track_id_2nd);
        pid_1st_values[ibin] = static_cast<float>(pid_1st);
        pid_2nd_values[ibin] = static_cast<float>(pid_2nd);
      }
    }

    trackid_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(label_trace));
    pid_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(pid_trace));

    trackid_1st_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(trackid_1st_trace));
    trackid_2nd_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(trackid_2nd_trace));
    pid_1st_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(pid_1st_trace));
    pid_2nd_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(pid_2nd_trace));
  }

  ITrace::shared_vector traces_out = traces_buffer;

  auto sframe =
    new SimpleFrame(in->ident(), in->time(), traces_out, in->tick(), in->masks());

  for (auto const& tag : m_frame_tags) {
    sframe->tag_frame(tag);
  }

  if (!m_output_trace_tag_trackid.empty()) {
    sframe->tag_traces(m_output_trace_tag_trackid, trackid_indices);
  }
  if (!m_output_trace_tag_pid.empty()) {
    sframe->tag_traces(m_output_trace_tag_pid, pid_indices);
  }
  if (!m_output_trace_tag_trackid_1st.empty()) {
    sframe->tag_traces(m_output_trace_tag_trackid_1st, trackid_1st_indices);
  }
  if (!m_output_trace_tag_trackid_2nd.empty()) {
    sframe->tag_traces(m_output_trace_tag_trackid_2nd, trackid_2nd_indices);
  }
  if (!m_output_trace_tag_pid_1st.empty()) {
    sframe->tag_traces(m_output_trace_tag_pid_1st, pid_1st_indices);
  }
  if (!m_output_trace_tag_pid_2nd.empty()) {
    sframe->tag_traces(m_output_trace_tag_pid_2nd, pid_2nd_indices);
  }

  out = IFrame::pointer(sframe);
  log->debug("output frame: {}", Aux::taginfo(out));
  return true;
}

const sim::SimChannel* AIML::Labelling2D::find_simchannel(unsigned int channel) const
{
  auto it = m_channel_index.find(channel);
  if (it == m_channel_index.end()) {
    return nullptr;
  }
  return &m_simchannels[it->second];
}

int AIML::Labelling2D::select_track_id(const sim::SimChannel& sc, int tdc_begin, int tdc_end) const
{
  if (tdc_end <= tdc_begin) {
    tdc_end = tdc_begin + 1;
  }

  const double start = static_cast<double>(tdc_begin);
  const double stop = static_cast<double>(tdc_end);
  auto matches = sc.TrackIDEs(start, stop);
  if (matches.empty()) {
    return m_default_label;
  }

  const sim::TrackIDE* best = nullptr;
  double best_charge = 0.0;
  for (auto const& match : matches) {
    const double charge = abs_charge(match.numElectrons);
    if (!best || charge > best_charge) {
      best = &match;
      best_charge = charge;
    }
  }
  if (!best) {
    return m_default_label;
  }
  return best->trackID;
}

std::vector<AIML::Labelling2D::TrackChargeInfo> AIML::Labelling2D::extract_track_charges(
  const sim::SimChannel& sc, int tdc_begin, int tdc_end) const
{
  if (tdc_end <= tdc_begin) {
    tdc_end = tdc_begin + 1;
  }

  const double start = static_cast<double>(tdc_begin);
  const double stop = static_cast<double>(tdc_end);
  auto matches = sc.TrackIDEs(start, stop);

  std::vector<TrackChargeInfo> track_charges;
  for (auto const& match : matches) {
    const double charge = abs_charge(match.numElectrons);
    if (charge > 0.0) {
      track_charges.emplace_back(match.trackID, charge);
    }
  }

  // Sort by charge in descending order
  std::sort(track_charges.begin(), track_charges.end(),
            [](const TrackChargeInfo& a, const TrackChargeInfo& b) {
              return a.second > b.second;
            });

  return track_charges;
}

std::pair<int, int> AIML::Labelling2D::select_top2_track_ids(const sim::SimChannel& sc,
                                                              int tdc_begin, int tdc_end) const
{
  auto track_charges = extract_track_charges(sc, tdc_begin, tdc_end);

  int top1_track_id = m_default_label;
  int top2_track_id = m_default_label;

  if (!track_charges.empty()) {
    top1_track_id = track_charges[0].first;
  }
  if (track_charges.size() > 1) {
    top2_track_id = track_charges[1].first;
  }

  return std::make_pair(top1_track_id, top2_track_id);
}

std::pair<int, int> AIML::Labelling2D::select_top2_pids(const sim::SimChannel& sc,
                                                        int tdc_begin, int tdc_end) const
{
  auto track_charges = extract_track_charges(sc, tdc_begin, tdc_end);

  // Merge charge contributions by PID
  std::unordered_map<int, double> pid_merged_charges;
  for (auto const& track_charge : track_charges) {
    int track_id = track_charge.first;
    double charge = track_charge.second;
    int pid = pid_from_track(track_id);

    // Accumulate charge for this PID
    if (pid_merged_charges.find(pid) == pid_merged_charges.end()) {
      pid_merged_charges[pid] = charge;
    } else {
      pid_merged_charges[pid] += charge;
    }
  }

  // Convert map to vector of PIDChargeInfo and sort by charge
  std::vector<PIDChargeInfo> pid_charges;
  for (auto const& pid_charge : pid_merged_charges) {
    pid_charges.emplace_back(pid_charge.first, pid_charge.second);
  }

  std::sort(pid_charges.begin(), pid_charges.end(),
            [](const PIDChargeInfo& a, const PIDChargeInfo& b) {
              return a.second > b.second;
            });

  int top1_pid = 0;
  int top2_pid = 0;

  if (!pid_charges.empty()) {
    top1_pid = pid_charges[0].first;
  }
  if (pid_charges.size() > 1) {
    top2_pid = pid_charges[1].first;
  }

  return std::make_pair(top1_pid, top2_pid);
}

void AIML::Labelling2D::cache_simchannels(const std::vector<sim::SimChannel>& simchs)
{
  m_simchannels = simchs;
  m_channel_index.clear();
  m_channel_index.reserve(m_simchannels.size());
  m_trackid_to_pid.clear();
  for (std::size_t idx = 0; idx < m_simchannels.size(); ++idx) {
    m_channel_index.emplace(m_simchannels[idx].Channel(), idx);
  }
}

void AIML::Labelling2D::populate_trackid_pid_map()
{
  m_trackid_to_pid.clear();
  if (m_simchannels.empty()) {
    return;
  }

  try {
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    for (auto const& sc : m_simchannels) {
      for (auto const& tdc_entry : sc.TDCIDEMap()) {
        for (auto const& ide : tdc_entry.second) {
          const int track_id = ide.trackID;
          if (track_id == 0) {
            continue;
          }
          if (m_trackid_to_pid.find(track_id) != m_trackid_to_pid.end()) {
            continue;
          }
          int pid = 0;
          try {
            auto const particle = pi_serv->TrackIdToParticle_P(track_id);
            if (particle) {
              pid = particle->PdgCode();
            }
          }
          catch (const cet::exception& ex) {
            log->debug("Labelling2D: failed to fetch MCParticle for track {}: {}", track_id, ex.what());
          }
          m_trackid_to_pid.emplace(track_id, pid);
        }
      }
    }
  }
  catch (const cet::exception& ex) {
    log->warn("Labelling2D: ParticleInventoryService unavailable: {}", ex.what());
  }
}

int AIML::Labelling2D::pid_from_track(int track_id) const
{
  if (track_id == 0) {
    return 0;
  }
  auto it = m_trackid_to_pid.find(track_id);
  if (it == m_trackid_to_pid.end()) {
    return 0;
  }
  return it->second;
}

void AIML::Labelling2D::clear_cache()
{
  m_simchannels.clear();
  m_channel_index.clear();
  m_trackid_to_pid.clear();
}
