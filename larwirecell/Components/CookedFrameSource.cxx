#include "CookedFrameSource.h"
#include "art/Framework/Principal/Handle.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "lardataobj/RecoBase/Wire.h"

#include "TTimeStamp.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(wclsCookedFrameSource,
                 wcls::CookedFrameSource,
                 wcls::IArtEventVisitor,
                 WireCell::IFrameSource)

using namespace wcls;
using namespace WireCell;
using WireCell::Aux::SimpleFrame;
using WireCell::Aux::SimpleTrace;

CookedFrameSource::CookedFrameSource() : m_nticks(0) {}

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
  const std::string art_tag = cfg["art_tag"].asString();
  if (art_tag.empty()) {
    THROW(ValueError() << errmsg{"WireCell::CookedFrameSource requires a source_label"});
  }
  m_inputTag = cfg["art_tag"].asString();

  m_tick = cfg["tick"].asDouble();
  for (auto jtag : cfg["frame_tags"]) {
    m_frame_tags.push_back(jtag.asString());
  }
  m_nticks = get(cfg, "nticks", m_nticks);
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

  const double time = tdiff(event.getRun().beginTime(), event.time());
  auto sframe = new SimpleFrame(event.event(), time, traces, tick);
  for (auto tag : m_frame_tags) {
    //std::cerr << "\ttagged: " << tag << std::endl;
    sframe->tag_frame(tag);
  }
  m_frames.push_back(WireCell::IFrame::pointer(sframe));
  m_frames.push_back(nullptr);
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
