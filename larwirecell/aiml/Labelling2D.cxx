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
  , m_default_label(0)
  , m_tdc_offset(0)
  , m_min_charge(0.0)
  , m_copy_input_traces(false)
{
  m_frame_tags.push_back(m_output_trace_tag_trackid);
}

AIML::Labelling2D::~Labelling2D() = default;

Configuration AIML::Labelling2D::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = m_anode_tn;
  cfg["reco_tag"] = m_reco_tag;
  cfg["output_trace_tag_trackid"] = m_output_trace_tag_trackid;
  cfg["output_trace_tag_pid"] = m_output_trace_tag_pid;
  cfg["default_label"] = m_default_label;
  cfg["tdc_offset"] = m_tdc_offset;
  cfg["min_charge"] = m_min_charge;
  cfg["copy_input_traces"] = m_copy_input_traces;
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
  const int configured_default = get(cfg, "default_label", m_default_label);
  if (configured_default != 0) {
    log->warn("Labelling2D overrides configured default_label {} with 0", configured_default);
  }
  m_default_label = 0;
  m_tdc_offset = get(cfg, "tdc_offset", m_tdc_offset);
  m_min_charge = get(cfg, "min_charge", m_min_charge);
  m_copy_input_traces = get(cfg, "copy_input_traces", m_copy_input_traces);

  m_simchannel_label = get(cfg, "simchannel_label", m_simchannel_label);

  m_frame_tags.clear();
  if (cfg.isMember("frame_tags")) {
    for (auto const& tag : cfg["frame_tags"]) {
      m_frame_tags.push_back(tag.asString());
    }
  }

  if (m_frame_tags.empty() && !m_output_trace_tag_trackid.empty()) {
    m_frame_tags.push_back(m_output_trace_tag_trackid);
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
  trackid_indices.reserve(reco_traces.size());
  pid_indices.reserve(reco_traces.size());

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

    SimpleTrace* label_trace = new SimpleTrace(chid, trace->tbin(), reco_charge.size());
    SimpleTrace* pid_trace = new SimpleTrace(chid, trace->tbin(), reco_charge.size());
    auto& label_values = label_trace->charge();
    auto& pid_values = pid_trace->charge();
    std::fill(label_values.begin(), label_values.end(), static_cast<float>(m_default_label));
    std::fill(pid_values.begin(), pid_values.end(), 0.0f);

    if (sc) {
      const int base_tbin = trace->tbin();
      for (std::size_t isample = 0; isample < reco_charge.size(); ++isample) {
        if (m_min_charge > 0.0 && abs_charge(reco_charge[isample]) < m_min_charge) {
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
    }

    trackid_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(label_trace));
    pid_indices.push_back(
      static_cast<IFrame::trace_list_t::value_type>(traces_buffer->size()));
    traces_buffer->push_back(ITrace::pointer(pid_trace));
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
