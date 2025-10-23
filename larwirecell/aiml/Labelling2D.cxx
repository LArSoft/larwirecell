#include "Labelling2D.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"

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
  , m_output_tag("truth")
  , m_default_label(-1)
  , m_tdc_offset(0)
  , m_min_charge(0.0)
  , m_copy_input_traces(false)
{
  m_frame_tags.push_back(m_output_tag);
}

AIML::Labelling2D::~Labelling2D() = default;

Configuration AIML::Labelling2D::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = m_anode_tn;
  cfg["reco_tag"] = m_reco_tag;
  cfg["output_tag"] = m_output_tag;
  cfg["default_label"] = m_default_label;
  cfg["tdc_offset"] = m_tdc_offset;
  cfg["min_charge"] = m_min_charge;
  cfg["copy_input_traces"] = m_copy_input_traces;
  cfg["simchannel"]["label"] = m_simchannel_label;
  cfg["frame_tags"] = Json::arrayValue;
  for (std::size_t ind = 0; ind < m_frame_tags.size(); ++ind) {
    cfg["frame_tags"][ind] = m_frame_tags[ind];
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
  m_output_tag = get(cfg, "output_tag", m_output_tag);
  m_default_label = get(cfg, "default_label", m_default_label);
  m_tdc_offset = get(cfg, "tdc_offset", m_tdc_offset);
  m_min_charge = get(cfg, "min_charge", m_min_charge);
  m_copy_input_traces = get(cfg, "copy_input_traces", m_copy_input_traces);

  if (cfg.isMember("simchannel")) {
    const auto& scfg = cfg["simchannel"];
    if (scfg.isMember("module_label")) {
      const std::string module = scfg["module_label"].asString();
      const std::string instance = get(scfg, "product_instance", std::string());
      const std::string process = get(scfg, "process", std::string());
      art::InputTag tag(module, instance, process);
      m_simchannel_label = tag.encode();
    }
    else {
      m_simchannel_label = get(scfg, "label", m_simchannel_label);
    }
  }
  else {
    m_simchannel_label = get(cfg, "simchannel_label", m_simchannel_label);
  }

  m_frame_tags.clear();
  if (cfg.isMember("frame_tags")) {
    for (auto const& tag : cfg["frame_tags"]) {
      m_frame_tags.push_back(tag.asString());
    }
  }

  if (m_frame_tags.empty() && !m_output_tag.empty()) {
    m_frame_tags.push_back(m_output_tag);
  }

  clear_cache();
}

void AIML::Labelling2D::visit(art::Event& event)
{
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
  log->debug("Labelling2D cached {} SimChannels", m_channel_index.size());
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

  ITrace::shared_vector traces_out;
  if (m_copy_input_traces) {
    auto in_traces = in->traces();
    if (in_traces) {
      traces_out = std::make_shared<ITrace::vector>(*in_traces);
    }
    else {
      traces_out = std::make_shared<ITrace::vector>();
    }
  }
  else {
    traces_out = std::make_shared<ITrace::vector>();
  }

  IFrame::trace_list_t label_indices;
  label_indices.reserve(reco_traces.size());

  for (auto const& trace : reco_traces) {
    if (!trace) {
      continue;
    }

    const int chid = trace->channel();
    const sim::SimChannel* sc = find_simchannel(chid);
    const auto& reco_charge = trace->charge();

    SimpleTrace* label_trace = new SimpleTrace(chid, trace->tbin(), reco_charge.size());
    auto& label_values = label_trace->charge();
    std::fill(label_values.begin(), label_values.end(), static_cast<float>(m_default_label));

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
      }
    }

    label_indices.push_back(traces_out->size());
    traces_out->push_back(ITrace::pointer(label_trace));
  }

  auto sframe =
    new SimpleFrame(in->ident(), in->time(), traces_out, in->tick(), in->masks());

  for (auto const& tag : m_frame_tags) {
    sframe->tag_frame(tag);
  }

  if (!m_output_tag.empty()) {
    sframe->tag_traces(m_output_tag, label_indices);
  }

  out = IFrame::pointer(sframe);
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

  auto best = std::max_element(matches.begin(),
                               matches.end(),
                               [](const sim::IDE& lhs, const sim::IDE& rhs) {
                                 return abs_charge(lhs.numElectrons) < abs_charge(rhs.numElectrons);
                               });

  if (best == matches.end()) {
    return m_default_label;
  }
  return best->trackID;
}

void AIML::Labelling2D::cache_simchannels(const std::vector<sim::SimChannel>& simchs)
{
  m_simchannels = simchs;
  m_channel_index.clear();
  m_channel_index.reserve(m_simchannels.size());
  for (std::size_t idx = 0; idx < m_simchannels.size(); ++idx) {
    m_channel_index.emplace(m_simchannels[idx].Channel(), idx);
  }
}

void AIML::Labelling2D::clear_cache()
{
  m_simchannels.clear();
  m_channel_index.clear();
}
