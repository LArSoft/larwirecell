#include "CookedFrameSink.h"
//#include "art/Framework/Principal/Handle.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "WireCellIface/IFrame.h"
#include "WireCellIface/ITrace.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(wclsCookedFrameSink,
                 wcls::CookedFrameSink,
                 wcls::IArtEventVisitor,
                 WireCell::IFrameSink)

using namespace wcls;
using namespace WireCell;

CookedFrameSink::CookedFrameSink() : m_frame(nullptr), m_nticks(0) {}

CookedFrameSink::~CookedFrameSink() {}

WireCell::Configuration CookedFrameSink::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = "AnodePlane";
  // frames to output
  cfg["frame_tags"][0] = "gauss";
  cfg["frame_tags"][1] = "wiener";
  cfg["nticks"] = m_nticks; // if nonzero, force number of ticks in output waveforms.
  return cfg;
}

void CookedFrameSink::configure(const WireCell::Configuration& cfg)
{
  const std::string anode_tn = cfg["anode"].asString();
  if (anode_tn.empty()) {
    THROW(ValueError() << errmsg{"CookedFrameSink requires an anode plane"});
  }

  // this throws if not found
  m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

  auto jtags = cfg["frame_tags"];
  std::cerr << "CookedFrameSink: saving " << jtags.size() << " tags\n";
  for (auto jtag : jtags) {
    std::string tag = jtag.asString();
    std::cerr << "\t" << tag << "\n";
    m_frame_tags.push_back(tag);
  }
  m_nticks = get(cfg, "nticks", m_nticks);
}

void CookedFrameSink::produces(art::ProducesCollector& collector)
{
  for (auto tag : m_frame_tags) {
    std::cerr << "CookedFrameSink: promising to produce recob::Wires named \"" << tag << "\"\n";
    collector.produces<std::vector<recob::Wire>>(tag);
  }
}

// fixme: this is copied from sigproc/FrameUtil but it's currently private code.
static ITrace::vector tagged_traces(IFrame::pointer frame, IFrame::tag_t tag)
{
  ITrace::vector ret;
  auto const& all_traces = frame->traces();
  for (size_t index : frame->tagged_traces(tag)) {
    ret.push_back(all_traces->at(index));
  }
  if (!ret.empty()) { return ret; }
  auto ftags = frame->frame_tags();
  if (std::find(ftags.begin(), ftags.end(), tag) == ftags.end()) { return ret; }
  return *all_traces; // must make copy of shared pointers
}

void CookedFrameSink::visit(art::Event& event)
{
  if (!m_frame) {
    std::cerr << "CookedFrameSink: I have no frame to save to art::Event\n";
    return;
  }

  std::cerr << "CookedFrameSink: got " << m_frame->traces()->size() << " total traces\n";

  for (auto tag : m_frame_tags) {

    auto traces = tagged_traces(m_frame, tag);
    if (traces.empty()) {
      std::cerr << "CookedFrameSink: no traces for tag \"" << tag << "\"\n";
      continue;
    }

    std::unique_ptr<std::vector<recob::Wire>> outwires(new std::vector<recob::Wire>);

    // what about the frame's time() and ident()?

    for (const auto& trace : traces) {

      const int tbin = trace->tbin();
      const int chid = trace->channel();
      const auto& charge = trace->charge();

      //std::cerr << tag << ": chid=" << chid << " tbin=" << tbin << " nq=" << charge.size() << std::endl;

      // enforce number of ticks if we are so configured.
      size_t ncharge = charge.size();
      int nticks = tbin + ncharge;
      if (m_nticks) { // force output waveform size
        if (m_nticks < nticks) { ncharge = m_nticks - tbin; }
        nticks = m_nticks;
      }
      recob::Wire::RegionsOfInterest_t roi(nticks);
      roi.add_range(tbin, charge.begin(), charge.begin() + ncharge);

      // FIXME: the current assumption in this code is that LS channel
      // numbers are identified with WCT channel IDs.
      // Fact: the plane view for the ICARUS induction-1 is "geo::kY",
      // instead of "geo::kU"
      auto const& gc = *lar::providerFrom<geo::Geometry>();
      auto view = gc.View(chid);

      // what about those pesky channel map masks?
      // they are dropped for now.

      outwires->emplace_back(recob::Wire(roi, chid, view));
    }
    std::cerr << "CookedFrameSink saving " << outwires->size() << " recob::Wires named \"" << tag
              << "\"\n";
    event.put(std::move(outwires), tag);
  }

  m_frame = nullptr;
}

bool CookedFrameSink::operator()(const WireCell::IFrame::pointer& frame)
{
  // set an IFrame based on last visited event.
  m_frame = frame;
  return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
