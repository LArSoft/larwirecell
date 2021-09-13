#include "RawFrameSource.h"
#include "art/Framework/Principal/Handle.h"

// for tick
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include "TTimeStamp.h"


#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(wclsRawFrameSource, wcls::RawFrameSource,
		 wcls::IArtEventVisitor, WireCell::IFrameSource)


using namespace wcls;
using namespace WireCell;

RawFrameSource::RawFrameSource()
    : m_nticks(0)
{
}

RawFrameSource::~RawFrameSource()
{
}


WireCell::Configuration RawFrameSource::default_configuration() const
{
    Configuration cfg;
    cfg["art_tag"] = "";        // how to look up the raw digits
    cfg["tick"] = 0.5*WireCell::units::us;
    cfg["frame_tags"][0] = "orig"; // the tags to apply to this frame
    cfg["nticks"] = m_nticks; // if nonzero, truncate or baseline-pad frame to this number of ticks.
    return cfg;
}

void RawFrameSource::configure(const WireCell::Configuration& cfg)
{
    const std::string art_tag = cfg["art_tag"].asString();
    if (art_tag.empty()) {
        THROW(ValueError() << errmsg{"WireCell::RawFrameSource requires a source_label"});
    }
    m_inputTag = cfg["art_tag"].asString();

    m_tick = cfg["tick"].asDouble();
    for (auto jtag : cfg["frame_tags"]) {
        m_frame_tags.push_back(jtag.asString());
    }
    m_nticks = get(cfg, "nticks", m_nticks);
}



// is this the right way to diff an art::Timestamp?
static
double tdiff(const art::Timestamp& ts1, const art::Timestamp& ts2)
{
    TTimeStamp tts1(ts1.timeHigh(), ts1.timeLow());
    TTimeStamp tts2(ts2.timeHigh(), ts2.timeLow());
    return tts2.AsDouble() - tts1.AsDouble();
}


static
SimpleTrace* make_trace(const raw::RawDigit& rd, unsigned int nticks_want)
{
    const int chid = rd.Channel();
    const int tbin = 0;
    const raw::RawDigit::ADCvector_t& adcv = rd.ADCs();

    short baseline = 0;
    unsigned int nadcs = adcv.size();
    if (nticks_want > 0) {      // don't want natural input size
        if (nticks_want > nadcs) {
            baseline = Waveform::most_frequent(adcv);
        }
        nadcs = std::min(nadcs, nticks_want);
    }
    else {
        nticks_want = nadcs;
    }

    auto strace = new SimpleTrace(chid, tbin, nticks_want);
    for (unsigned int itick=0; itick < nadcs; ++ itick) {
	strace->charge()[itick] = adcv[itick];
    }
    for (unsigned int itick = nadcs; itick < nticks_want; ++itick) {
	strace->charge()[itick] = baseline;
    }
    return strace;
}

void RawFrameSource::visit(art::Event & event)
{
    // fixme: want to avoid depending on DetectorPropertiesService for now.
    const double tick = m_tick;
    art::Handle< std::vector<raw::RawDigit> > rdvh;
    bool okay = event.getByLabel(m_inputTag, rdvh);
    if (!okay) {
        std::string msg = "WireCell::RawFrameSource failed to get vector<raw::RawDigit>: " + m_inputTag.encode();
        std::cerr << msg << std::endl;
        THROW(RuntimeError() << errmsg{msg});
    }
    else if (rdvh->size() == 0) return;

    const std::vector<raw::RawDigit>& rdv(*rdvh);
    const size_t nchannels = rdv.size();
    std::cerr << "RawFrameSource: got " << nchannels << " raw::RawDigit objects\n";

    WireCell::ITrace::vector traces(nchannels);
    for (size_t ind=0; ind<nchannels; ++ind) {
        auto const& rd = rdv.at(ind);
        traces[ind] = ITrace::pointer(make_trace(rd, m_nticks));
	if (!ind) {
            if (m_nticks) {
                std::cerr
                    << "\tinput nticks=" << rd.ADCs().size() << " setting to " << m_nticks
                    << std::endl;
            }
            else {
                std::cerr
                    << "\tinput nticks=" << rd.ADCs().size() << " keeping as is"
                    << std::endl;
            }
	}
    }

    const double time = tdiff(event.getRun().beginTime(), event.time());
    auto sframe = new WireCell::SimpleFrame(event.event(), time, traces, tick);
    for (auto tag : m_frame_tags) {
        //std::cerr << "\ttagged: " << tag << std::endl;
        sframe->tag_frame(tag);
    }
    m_frames.push_back(WireCell::IFrame::pointer(sframe));
    m_frames.push_back(nullptr);
}

bool RawFrameSource::operator()(WireCell::IFrame::pointer& frame)
{
    frame = nullptr;
    if (m_frames.empty()) {
        return false;
    }
    frame = m_frames.front();
    m_frames.pop_front();
    return true;
}
// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
