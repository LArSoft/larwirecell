#include "CookedFrameSource.h"
#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "lardataobj/RecoBase/Wire.h"
 
#include "TTimeStamp.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Waveform.h"

WIRECELL_FACTORY(wclsCookedFrameSource,
                 wcls::CookedFrameSource,
                 wcls::IArtEventVisitor,
                 WireCell::IFrameSource)

using namespace wcls;
using namespace WireCell;
using WireCell::Aux::SimpleFrame;
using WireCell::Aux::SimpleTrace;

CookedFrameSource::CookedFrameSource() : m_nticks(0), l(Log::logger("io")) {}

CookedFrameSource::~CookedFrameSource() {}

WireCell::Configuration CookedFrameSource::default_configuration() const
{
  Configuration cfg;
  cfg["art_tag"] = ""; // how to look up the cooked digits
  cfg["tick"] = 0.5 * WireCell::units::us;
  cfg["frame_tags"][0] = "orig"; // the tags to apply to this frame
  cfg["nticks"] = m_nticks;      // if nonzero, truncate or zero-pad frame to this number of ticks.
  return cfg;
}

void CookedFrameSource::configure(const WireCell::Configuration& cfg)
{
  m_scale = cfg["scale"].asDouble();
  m_tick = cfg["tick"].asDouble();
  m_nticks = get(cfg, "nticks", m_nticks);


  for (auto recobwire_tag : cfg["recobwire_tags"]) {
    m_recobwire_tags.push_back(recobwire_tag.asString());
  }
  for (auto jtag : cfg["frame_tags"]) {
    m_frame_tags.push_back(jtag.asString());
  }
  for (auto jtag : cfg["trace_tags"]) {
    m_trace_tags.push_back(jtag.asString());
  }
  for (auto summary_tag : cfg["summary_tags"]) {
    m_summary_tags.push_back(summary_tag.asString());
  }
  for (auto mask_tag : cfg["input_mask_tags"]) {
    m_input_mask_tags.push_back(mask_tag.asString());
  }
  for (auto mask_tag : cfg["output_mask_tags"]) {
    m_output_mask_tags.push_back(mask_tag.asString());
  }
  if (m_recobwire_tags.size() != m_summary_tags.size()) {
    raise<ValueError>("m_recobwire_tags.size %d != m_summary_tags.size %d",
                      m_recobwire_tags.size(),
                      m_summary_tags.size());
  }
  if (m_recobwire_tags.size() != m_trace_tags.size()) {
    raise<ValueError>("m_recobwire_tags.size %d != m_trace_tags.size %d",
                      m_recobwire_tags.size(),
                      m_trace_tags.size());
  }
  if (m_input_mask_tags.size() != m_output_mask_tags.size()) {
    raise<ValueError>("m_input_mask_tags.size %d != m_output_mask_tags.size %d",
                      m_input_mask_tags.size(),
                      m_output_mask_tags.size());
  }
}

// this code assumes that the high part of timestamp represents number of seconds from Jan 1st, 1970 and the low part
// represents the number of nanoseconds.
static double tdiff(const art::Timestamp& ts1, const art::Timestamp& ts2)
{
  TTimeStamp tts1(ts1.timeHigh(), ts1.timeLow());
  TTimeStamp tts2(ts2.timeHigh(), ts2.timeLow());
  return tts2.AsDouble() - tts1.AsDouble();
}

static SimpleTrace* make_trace(const recob::Wire& rw, unsigned int nticks_want, const double scale = 50)
{
  // uint
  const raw::ChannelID_t chid = rw.Channel();
  const int tbin = 0;
  const std::vector<float> sig = rw.Signal();

  const float baseline = 0.0;
  unsigned int nsamp = sig.size();
  if (nticks_want > 0) { nsamp = std::min(nsamp, nticks_want); }
  else {
    nticks_want = nsamp;
  }

  auto strace = new SimpleTrace(chid, tbin, nticks_want);
  auto& q = strace->charge();
  for (unsigned int itick = 0; itick < nsamp; ++itick) {
    q[itick] = scale*sig[itick]; // changed Ewerton 2023-10-06 recob::Wire it scaled up by a factor (make it onfigurable!!!)
  }
  for (unsigned int itick = nsamp; itick < nticks_want; ++itick) {
    q[itick] = baseline;
  }
  return strace;
}

void CookedFrameSource::visit(art::Event& event)
{
  l->debug("CookedFrameSource::visit: event: {}", event.event());
  // fixme: want to avoid depending on DetectorPropertiesService for now.
  const double tick = m_tick;
  const double time = tdiff(event.getRun().beginTime(), event.time());

  ITrace::vector* itraces = new ITrace::vector; // will become shared_ptr.
  IFrame::trace_list_t indices;
  std::map<std::string, IFrame::trace_list_t> tag2indices;

  for (auto const& recobwire_tag : m_recobwire_tags) {
    art::Handle<std::vector<recob::Wire>> rwvh;
    art::InputTag recobwire_tag_art(recobwire_tag);
    bool okay = event.getByLabel(recobwire_tag_art, rwvh);
    if (!okay) {
      raise<RuntimeError>("CookedFrameSource failed to get vector<recob::Wire>: %s", recobwire_tag.c_str());
    }
    const std::vector<recob::Wire>& rwv(*rwvh);
    const size_t nchannels = rwv.size();
    l->debug("wcls::CookedFrameSource: got {} {} recob::Wire objects", nchannels, recobwire_tag);
    for (size_t ind = 0; ind < nchannels; ++ind) {
      auto const& rw = rwv.at(ind);
      SimpleTrace* trace = make_trace(rw, m_nticks, m_scale);
      const size_t trace_index = itraces->size();

      indices.push_back(trace_index);
      tag2indices[recobwire_tag].push_back(trace_index);
      itraces->push_back(ITrace::pointer(trace));
    }
  }

  std::map<std::string, art::Handle<std::vector<double>>> tag2summaryh;
  for (auto const& summary_tag : m_summary_tags) {
    if (summary_tag.empty()) {
      continue;
    }
    bool okay = event.getByLabel(summary_tag, tag2summaryh[summary_tag]);
    if (!okay) {
      raise<RuntimeError>("CookedFrameSource failed to get vector<double>: %s", summary_tag.c_str());
    }
  }

  std::map<std::string, art::Handle<std::vector<int>>> tag2maskh;
  for (auto mask_tag : m_input_mask_tags) {
    bool okay = event.getByLabel(mask_tag, tag2maskh[mask_tag]);
    if (!okay) {
      raise<RuntimeError>("CookedFrameSource failed to get vector<int>: %s", mask_tag.c_str());
    }
  }
  std::map<std::string, std::string> tag2outtag;
  for (size_t ind = 0; ind < m_input_mask_tags.size(); ++ind) {
    tag2outtag[m_input_mask_tags[ind]] = m_output_mask_tags[ind];
  }

  Waveform::ChannelMaskMap cmm;

  for (auto [tag, mask] : tag2maskh) {
    Waveform::ChannelMasks& cm = cmm[tag2outtag[tag]];
    size_t nchannels = mask->size() / 3;
    for (size_t i = 0; i < nchannels; i++) {
      size_t ch_idx = 3 * i;
      size_t low_idx = 3 * i + 1;
      size_t up_idx = 3 * i + 2;
      auto cmch = mask->at(ch_idx);
      auto cm_first = mask->at(low_idx);
      auto cm_second = mask->at(up_idx);
      Waveform::BinRange bins(cm_first, cm_second);
      cm[cmch].push_back(bins);
    }
  }

  auto sframe =
    new Aux::SimpleFrame(event.event(), time, ITrace::shared_vector(itraces), tick, cmm);

  for (auto tag : m_frame_tags) {
    sframe->tag_frame(tag);
  }
  for (size_t ind = 0; ind < m_trace_tags.size(); ++ind) {
    auto recobwire_tag = m_recobwire_tags[ind];
    auto trace_tag = m_trace_tags[ind];
    auto summary_tag = m_summary_tags[ind];
    if (tag2summaryh.find(summary_tag) != tag2summaryh.end()) {
      sframe->tag_traces(trace_tag, tag2indices[recobwire_tag], *tag2summaryh[summary_tag]);
    }
    else {
      sframe->tag_traces(trace_tag, tag2indices[recobwire_tag]);
    }
  }

  m_frames.push_back(WireCell::IFrame::pointer(sframe));
  m_frames.push_back(nullptr); // <= put back 2023-10-05
}

bool CookedFrameSource::operator()(WireCell::IFrame::pointer& frame)
{
  frame = nullptr;
  if (m_frames.empty()) { return false; }
  frame = m_frames.front();
  m_frames.pop_front();
  return true;
}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
