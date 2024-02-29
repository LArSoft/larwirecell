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

  // ----- Ewerton: Modify original CookedFrameSource
  // ----- to read in products for wire-cell imaging

  cfg["wiener_inputTag"] = "";                 // art tag for wiener recob::Wire
  cfg["gauss_inputTag"] = "";                  // art tag for gauss recob::Wire
  cfg["badmasks_inputTag"] = Json::arrayValue; // list of art tags for bad channels
  cfg["threshold_inputTag"] = "";              // art tag for thresholds
  cfg["cmm_tag"] = "";                         // for now use a single tag for output cmm
  return cfg;
}

void CookedFrameSource::configure(const WireCell::Configuration& cfg)
{
  m_tick = cfg["tick"].asDouble();
  for (auto jtag : cfg["frame_tags"]) {
    m_frame_tags.push_back(jtag.asString());
  }
  m_nticks = get(cfg, "nticks", m_nticks);

  // check setup to be used
  m_wiener_inputTag = cfg["wiener_inputTag"].asString();
  m_gauss_inputTag = cfg["gauss_inputTag"].asString();
  for (auto badmask : cfg["badmasks_inputTag"]) {
    m_badmasks_inputTag.push_back(badmask.asString());
  }
  m_threshold_inputTag = cfg["threshold_inputTag"].asString();
  m_cmm_tag = cfg["cmm_tag"].asString();

  if (m_wiener_inputTag.empty() && m_gauss_inputTag.empty() && !m_badmasks_inputTag.size() &&
      m_threshold_inputTag.empty() && m_cmm_tag.empty()) { // use standard setup
    l->debug("wcls::CookedFrameSource: using standard setup");
    m_cmm_setup = false;
  }
  else { // use cmm setup
    l->debug("wcls::CookedFrameSource: using cmm setup");
    m_cmm_setup = true;
  }

  if (!m_cmm_setup) {
    const std::string art_tag = cfg["art_tag"].asString();
    if (art_tag.empty()) {
      THROW(ValueError() << errmsg{"wcls::CookedFrameSource requires a source_label"});
    }
    m_inputTag = cfg["art_tag"].asString();
  }
  else {
    if (m_wiener_inputTag.empty())
      THROW(ValueError() << errmsg{
              "wcls::CookedFrameSource requires a source_label for recob::Wire (wiener)"});
    if (m_gauss_inputTag.empty())
      THROW(ValueError() << errmsg{
              "wcls::CookedFrameSource requires a source_label for recob::Wire (gauss)"});
    if (m_badmasks_inputTag.size()) {
      for (auto badmask : m_badmasks_inputTag) {
        if (badmask.empty())
          THROW(ValueError() << errmsg{
                  "wcls::CookedFrameSource cannot use an empty source_label for bad masks"});
      }
    }
    else {
      THROW(
        ValueError() << errmsg{"wcls::CookedFrameSource requires a source_label for bad masks"});
    }
    if (m_threshold_inputTag.empty())
      THROW(
        ValueError() << errmsg{"wcls::CookedFrameSource requires a source_label for threshold"});
    if (m_cmm_tag.empty())
      THROW(ValueError() << errmsg{"wcls::CookedFrameSource requires a label for cmm_tag"});
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

static SimpleTrace* make_trace(const recob::Wire& rw, unsigned int nticks_want)
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
    q[itick] = sig[itick];
  }
  for (unsigned int itick = nsamp; itick < nticks_want; ++itick) {
    q[itick] = baseline;
  }
  return strace;
}

void CookedFrameSource::visit(art::Event& e)
{
  auto const& event = e;
  // fixme: want to avoid depending on DetectorPropertiesService for now.
  const double tick = m_tick;
  const double time = tdiff(event.getRun().beginTime(), event.time());
  if (!m_cmm_setup) { // do standard setup
    art::Handle<std::vector<recob::Wire>> rwvh;
    bool okay = event.getByLabel(m_inputTag, rwvh);
    if (!okay) {
      std::string msg =
        "WireCell::CookedFrameSource failed to get vector<recob::Wire>: " + m_inputTag.encode();
      std::cerr << msg << std::endl;
      THROW(RuntimeError() << errmsg{msg});
    }
    else if (rwvh->size() == 0)
      return;

    const std::vector<recob::Wire>& rwv(*rwvh);
    const size_t nchannels = rwv.size();
    std::cerr << "CookedFrameSource: got " << nchannels << " recob::Wire objects\n";

    WireCell::ITrace::vector traces(nchannels);
    for (size_t ind = 0; ind < nchannels; ++ind) {
      auto const& rw = rwv.at(ind);
      traces[ind] = ITrace::pointer(make_trace(rw, m_nticks));
      if (!ind) { // first time through
        if (m_nticks) {
          std::cerr << "\tinput nticks=" << rw.NSignal() << " setting to " << m_nticks << std::endl;
        }
        else {
          std::cerr << "\tinput nticks=" << rw.NSignal() << " keeping as is" << std::endl;
        }
      }
    }

    auto sframe = new SimpleFrame(event.event(), time, traces, tick);
    for (auto tag : m_frame_tags) {
      //std::cerr << "\ttagged: " << tag << std::endl;
      sframe->tag_frame(tag);
    }
    m_frames.push_back(WireCell::IFrame::pointer(sframe));
    //m_frames.push_back(nullptr); //<- passes empty frame to next module?
  }
  else { // do cmm setup

    art::Handle<std::vector<recob::Wire>> rwvh_wiener;
    art::Handle<std::vector<recob::Wire>> rwvh_gauss;
    std::map<std::string, art::Handle<std::vector<int>>> bad_masks;
    art::Handle<std::vector<double>> threshold;
    bool okay_wiener = event.getByLabel(m_wiener_inputTag, rwvh_wiener);
    bool okay_gauss = event.getByLabel(m_gauss_inputTag, rwvh_gauss);
    for (auto badmask : m_badmasks_inputTag) {
      bool okay = e.getByLabel(badmask, bad_masks[badmask.label()]);
      if (!okay) {
        std::string msg = "wcls::CookedFrameSource failed to get vector<int>: " + badmask.encode() +
                          " (badmask not found)";
        std::cerr << msg << std::endl;
        THROW(RuntimeError() << errmsg{msg});
      }
      else if (bad_masks[badmask.label()]->size() == 0) {
        std::string msg = "wcls::CookedFrameSource bad mask " + badmask.encode() + " is empty";
        std::cerr << msg << std::endl;
        return;
      }
    }
    bool okay_th = e.getByLabel(m_threshold_inputTag, threshold);

    if (!okay_wiener) {
      std::string msg =
        "wcls::CookedFrameSource failed to get vector<recob::Wire>: " + m_wiener_inputTag.encode() +
        " (wire-cell imaging needs both wiener and gauss recob::Wire)";
      std::cerr << msg << std::endl;
      THROW(RuntimeError() << errmsg{msg});
    }
    else if (rwvh_wiener->size() == 0)
      return;

    if (!okay_gauss) {
      std::string msg =
        "wcls::CookedFrameSource failed to get vector<recob::Wire>: " + m_gauss_inputTag.encode() +
        " (wire-cell imaging needs both wiener and gauss recob::Wire)";
      std::cerr << msg << std::endl;
      THROW(RuntimeError() << errmsg{msg});
    }
    else if (rwvh_gauss->size() == 0)
      return;

    if (!okay_th) {
      std::string msg =
        "wcls::CookedFrameSource failed to get vector<double>: " + m_threshold_inputTag.encode() +
        " (wire-cell imaging needs thresholds)";
      std::cerr << msg << std::endl;
      THROW(RuntimeError() << errmsg{msg});
    }
    else if (threshold->size() == 0)
      return;

    const std::vector<recob::Wire>& wiener_rwv(*rwvh_wiener);
    const std::vector<recob::Wire>& gauss_rwv(*rwvh_gauss);
    const size_t nchannels_wiener = wiener_rwv.size();
    const size_t nchannels_gauss = gauss_rwv.size();
    l->debug("wcls::CookedFrameSource: got {} (wiener) recob::Wire objects", nchannels_wiener);
    l->debug("wcls::CookedFrameSource: got {} (gauss) recob::Wire objects", nchannels_gauss);

    ITrace::vector* itraces = new ITrace::vector; // will become shared_ptr.
    IFrame::trace_list_t indices;
    IFrame::trace_list_t wiener_traces, gauss_traces;

    for (size_t ind = 0; ind < nchannels_wiener; ++ind) {
      auto const& rw = wiener_rwv.at(ind);
      SimpleTrace* trace = make_trace(rw, m_nticks);
      const size_t trace_index = itraces->size();

      indices.push_back(trace_index);
      wiener_traces.push_back(trace_index);
      itraces->push_back(ITrace::pointer(trace));

      if (!ind) { // first time through
        if (m_nticks) {
          l->debug("wcls::CookedFrameSource: input nticks= {} setting to {} (wiener)",
                   rw.NSignal(),
                   m_nticks);
        }
        else {
          l->debug("wcls::CookedFrameSource: input nticks= {} keeping as is (wiener)",
                   rw.NSignal());
        }
      }
    }

    for (size_t ind = 0; ind < nchannels_gauss; ++ind) {
      auto const& rw = gauss_rwv.at(ind);
      SimpleTrace* trace = make_trace(rw, m_nticks);
      const size_t trace_index = itraces->size();

      indices.push_back(trace_index);
      gauss_traces.push_back(trace_index);
      itraces->push_back(ITrace::pointer(trace));

      if (!ind) { // first time through
        if (m_nticks) {
          l->debug("wcls::CookedFrameSource: input nticks= {} setting to {} (gauss)\n",
                   rw.NSignal(),
                   m_nticks);
        }
        else {
          l->debug("wcls::CookedFrameSource: input nticks= {} keeping as is (gauss)\n",
                   rw.NSignal());
        }
      }
    }

    Waveform::ChannelMaskMap cmm;
    Waveform::ChannelMasks cm_out;

    for (auto bad_mask : bad_masks) {
      Waveform::ChannelMasks cm;
      size_t nchannels = bad_mask.second->size() / 3;
      for (size_t i = 0; i < nchannels; i++) {
        size_t ch_idx = 3 * i;
        size_t low_idx = 3 * i + 1;
        size_t up_idx = 3 * i + 2;
        auto cmch = bad_mask.second->at(ch_idx);
        auto cm_first = bad_mask.second->at(low_idx);
        auto cm_second = bad_mask.second->at(up_idx);
        Waveform::BinRange bins(cm_first, cm_second);
        cm[cmch].push_back(bins);
      }
      cm_out = Waveform::merge(cm_out, cm);
    }

    cmm[m_cmm_tag] = cm_out; // for now use a single tag for output cmm

    auto sframe =
      new Aux::SimpleFrame(event.event(), time, ITrace::shared_vector(itraces), tick, cmm);
    for (auto tag : m_frame_tags) {
      sframe->tag_frame(tag);
    }

    sframe->tag_traces(m_wiener_tag, wiener_traces, *threshold);
    sframe->tag_traces(m_gauss_tag, gauss_traces);

    m_frames.push_back(WireCell::IFrame::pointer(sframe));
    //m_frames.push_back(nullptr); //<- passes empty frame to next module?
  }
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
