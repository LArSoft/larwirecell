#include "LazyFrameSource.h"
#include "art/Framework/Principal/Handle.h"

// for tick
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/RawDigit.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"

#include "TTimeStamp.h"


#include "WireCellIface/IFrame.h"

namespace wcls {

    class LazyTrace : public WireCell::ITrace {
        mutable art::Handle< std::vector<raw::RawDigit> > m_rdvh;
        size_t m_index;
        int m_channel;
        
        mutable WireCell::ITrace::ChargeSequence m_charge;


    public:
        LazyTrace(art::Handle< std::vector<raw::RawDigit> > rdvh, size_t index)
            : m_rdvh(rdvh), m_index(index), m_channel(rdvh->at(index).Channel()) {}


	virtual int channel() const { return m_channel; }
	virtual int tbin() const { return 0; }

	virtual const ChargeSequence& charge() const {
            if (m_charge.empty()) {
                auto const& rd = m_rdvh->at(m_index);
                const raw::RawDigit::ADCvector_t& adcv = rd.ADCs();
                //std::cerr << "trace " << m_index << " chan " << m_channel << " with " << adcv.size() << " samples\n";
                m_charge.insert(m_charge.end(), adcv.begin(), adcv.end());
                m_rdvh.clear(); // bye bye
            }
            return m_charge;
        }

    };

    class LazyFrame : public WireCell::IFrame {
        int m_ident;
        double m_time, m_tick;
        tag_list_t m_tags;
        WireCell::ITrace::shared_vector m_traces;
    public:
        LazyFrame(art::Handle< std::vector<raw::RawDigit> > rdvh,
                  int ident, double time, double tick, const tag_list_t& tags)
            : m_ident(ident), m_time(time), m_tick(tick), m_tags(tags.begin(), tags.end()) {
            const size_t nrds = rdvh->size();
            auto* traces = new std::vector<LazyTrace::pointer>(nrds);
            for (size_t ind = 0; ind < nrds; ++ind) {
                traces->at(ind) = std::make_shared<LazyTrace>(rdvh, ind);
            }
            m_traces = WireCell::ITrace::shared_vector(traces);
        }

        virtual ~LazyFrame() { }

        virtual const tag_list_t& frame_tags() const {
            return m_tags;
        }

        virtual const tag_list_t& trace_tags() const {
            static tag_list_t dummy; // this frame doesn't support trace tags
            return dummy;
        }

        virtual const trace_list_t& tagged_traces(const tag_t& tag) const {
            static trace_list_t dummy; // this frame doesn't support trace tags
            return dummy;
        }

        virtual const trace_summary_t& trace_summary(const tag_t& tag) const {
            static trace_summary_t dummy; // this frame doesn't support trace tags
            return dummy;
        }


	virtual WireCell::ITrace::shared_vector traces() const {
            return m_traces;
        }

	virtual int ident() const {
            return m_ident;
        }

	virtual double time() const {
            return m_time;
        }

	virtual double tick() const {
            return m_tick;
        }

        
    };

}



#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(wclsLazyFrameSource, wcls::LazyFrameSource,
		 wcls::IArtEventVisitor, WireCell::IFrameSource)


using namespace wcls;
using namespace WireCell;

LazyFrameSource::LazyFrameSource()
    : m_nticks(0)
{
}

LazyFrameSource::~LazyFrameSource()
{
}


WireCell::Configuration LazyFrameSource::default_configuration() const
{
    Configuration cfg;
    cfg["art_tag"] = "";        // how to look up the raw digits
    cfg["tick"] = 0.5*WireCell::units::us;
    cfg["frame_tags"][0] = "orig"; // the tags to apply to this frame
    cfg["nticks"] = m_nticks; // if nonzero, truncate or baseline-pad frame to this number of ticks.
    return cfg;
}

void LazyFrameSource::configure(const WireCell::Configuration& cfg)
{
    const std::string art_tag = cfg["art_tag"].asString();
    if (art_tag.empty()) {
        THROW(ValueError() << errmsg{"LazyFrameSource requires a source_label"});
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

void LazyFrameSource::visit(art::Event & event)
{
    // fixme: want to avoid depending on DetectorPropertiesService for now.
    const double tick = m_tick;

    art::Handle< std::vector<raw::RawDigit> > rdvh;
    bool okay = event.getByLabel(m_inputTag, rdvh);
    if (!okay) {
        std::string msg = "LazyFrameSource failed to get vector<raw::RawDigit>: " + m_inputTag.encode();
        std::cerr << msg << std::endl;
        THROW(RuntimeError() << errmsg{msg});
    }
    else if (rdvh->size() == 0) return;
    const double time = tdiff(event.getRun().beginTime(), event.time());

    std::cerr << "LazyFrameSource: got " << rdvh->size() << " raw::RawDigit objects\n";

    m_frames.push_back(std::make_shared<LazyFrame>(rdvh, event.event(), time, tick, m_frame_tags));
    m_frames.push_back(nullptr);
}

bool LazyFrameSource::operator()(WireCell::IFrame::pointer& frame)
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
