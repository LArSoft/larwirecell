#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"

namespace butcher {

    struct EventButcherConfig {

	// fixme: provide reasonable defaults

	fhicl::Atom<std::string> inRawTag  {
	    fhicl::Name("inRawTag"),
		fhicl::Comment("Input tag for the raw::RawDigit collection.") };

	fhicl::Atom<std::string> inSigTag  { fhicl::Name("inSigTag"),
		fhicl::Comment("Input tag for the recob::Wire collection.") };

	fhicl::Atom<std::string> outRawTag {
	    fhicl::Name("outRawTag"),
		fhicl::Comment("Label the output raw::RawDigit collection."),
		"truncraw"};

	fhicl::Atom<std::string> outSigTag  {
	    fhicl::Name("outSigTag"),
		fhicl::Comment("Label the output recob::Wire collection."),
		"truncsig"};

	fhicl::Atom<std::string> outAssnTag  {
	    fhicl::Name("outAssnTag"),
		fhicl::Comment("Label the output associations."),
		"rawsigassns" };

	fhicl::Atom<int> ndrop { fhicl::Name("ndrop"),
		fhicl::Comment("Number of ticks from start of waveform to drop"),
		0 };
	fhicl::Atom<int> nkeep { fhicl::Name("nkeep"),
		fhicl::Comment("Number of remaining ticks to keep after initial drop"),
		-1 };

	fhicl::Atom<double> sigscale { fhicl::Name("sigscale"),
		fhicl::Comment("A multiplicative scale factor applied to the output recob::Wires"),
		 1.0 };
    };

    class EventButcher : public art::EDProducer {
    public:

	using Parameters = art::EDProducer::Table<EventButcherConfig>;

        explicit EventButcher(Parameters const& params);
        virtual ~EventButcher();

        void produce(art::Event & evt);

    private:

	const EventButcherConfig m_cfg;
	// inputs
	art::InputTag m_rawtag, m_sigtag;
	// this needs art 2.08
	//art::ProductToken< std::vector<raw::RawDigit> > m_rawtok;
	//art::ProductToken< std::vector<recob::Wire>   > m_sigtok;
    };


}

using namespace std;

butcher::EventButcher::EventButcher(Parameters const& params)
    : EDProducer{params}
    , m_cfg(params())
    , m_rawtag{m_cfg.inRawTag()}
    , m_sigtag{m_cfg.inSigTag()}
//    , m_rawtok{consumes< std::vector<raw::RawDigit> >(m_rawtag)}
//    , m_sigtok{consumes< std::vector<recob::Wire> >(m_sigtag)}
{
    //cerr << "Producing: outraw:"<<m_cfg.outRawTag()<<" outsig:"<<m_cfg.outSigTag()<<" outras:" << m_cfg.outAssnTag() << endl;
    produces< std::vector<raw::RawDigit> >(m_cfg.outRawTag());
    produces< std::vector<recob::Wire> >(m_cfg.outSigTag());
    produces< art::Assns<raw::RawDigit,recob::Wire> >(m_cfg.outAssnTag());
}

butcher::EventButcher::~EventButcher()
{
}


void butcher::EventButcher::produce(art::Event & event)
{

    //cerr <<"In raw=" << m_rawtag << " in sig=" << m_sigtag << "\n";

    art::Handle< std::vector<raw::RawDigit> > raw;
    //event.getByToken(m_rawtok, raw);
    event.getByLabel(m_rawtag, raw);
    const size_t nraw = raw->size();

    art::Handle< std::vector<recob::Wire> > sig;
    //event.getByToken(m_sigtok, sig);
    event.getByLabel(m_sigtag, sig);
    const size_t nsig = sig->size();

    // https://cdcvs.fnal.gov/redmine/projects/art/wiki/The_PtrMaker_utility
    art::PtrMaker<raw::RawDigit> RawPtr(event, m_cfg.outRawTag());
    art::PtrMaker<recob::Wire> SigPtr(event, m_cfg.outSigTag());

    // keep track of raw digit Ptr by its channel id, for later look
    // up in making associations.
    std::unordered_map<raw::ChannelID_t, art::Ptr<raw::RawDigit> > chid2rawptr;

    // raw-signal association
    auto outrsa = std::make_unique< art::Assns<raw::RawDigit,recob::Wire> >();
    auto outraw = std::make_unique< std::vector<raw::RawDigit> >();
    auto outsig = std::make_unique< std::vector<recob::Wire> >();

    const int ndrop = m_cfg.ndrop();
    const int nkeep = m_cfg.nkeep();
    const double sigscale = m_cfg.sigscale();

    // Truncate the raw digits
    for (size_t iraw=0; iraw != nraw; ++iraw) {
	const auto& inrd = raw->at(iraw);
	const auto& inadcs = inrd.ADCs();
	const size_t inlen = inadcs.size();

	const int outlen = std::min(inlen-ndrop, nkeep < 0 ? inlen : nkeep);

	if (outlen <= 0) {
	    continue;		// add a "keep_empty" option and save an empty place holder RawDigit here?
	}

	size_t outind = outraw->size();
	const raw::ChannelID_t chid = inrd.Channel();
	raw::RawDigit::ADCvector_t outadc(inadcs.begin()+ndrop, inadcs.begin()+ndrop+outlen);
	outraw->emplace_back(raw::RawDigit(chid, outlen, outadc, raw::kNone));
	// given the truncationn, this is technically a BOGUS thing to do
	auto& outrd = outraw->back();
	outrd.SetPedestal(inrd.GetPedestal(), inrd.GetSigma());

	chid2rawptr[chid] = RawPtr(outind);
    }


    // Truncate the signal and make assns
    for (size_t isig=0; isig != nsig; ++isig) {
	const auto& inw = sig->at(isig);
	std::vector<float> wave = inw.Signal();
	const size_t inlen = wave.size();

	const int outlen = std::min(inlen-ndrop, nkeep < 0 ? inlen : nkeep);

	if (outlen <= 0) {
	    continue;
	}

	const auto chid = inw.Channel();
	const auto view = inw.View();

	// resparsify
	recob::Wire::RegionsOfInterest_t roi(outlen);
	auto first = wave.begin()+ndrop;
	auto done = wave.begin()+ndrop+outlen;
	auto beg = first;
	while (true) {
	    beg = std::find_if(beg, done, [](float v){return v != 0.0;});
	    if (beg == done) {
		break;
	    }
	    auto end = std::find_if(beg, done, [](float v){return v == 0.0;});

	    std::vector<float> scaled(beg, end);
	    for (int ind=0; ind<end-beg; ++ind) {
		scaled[ind] *= sigscale;
	    }
	    roi.add_range(beg-first, scaled.begin(), scaled.end());
	    beg = end;
	}

	const size_t outind = outsig->size();
	outsig->emplace_back(recob::Wire(roi, chid, view));

	// associate
	auto rawit = chid2rawptr.find(chid);
	if (rawit == chid2rawptr.end()) {
	    continue;		// emit warning about no accompaning raw digit?
	}
	auto const& rawptr = rawit->second;
	auto const sigptr = SigPtr(outind);
	outrsa->addSingle(rawptr, sigptr);
    }

    event.put(std::move(outraw), m_cfg.outRawTag());
    event.put(std::move(outsig), m_cfg.outSigTag());
    event.put(std::move(outrsa), m_cfg.outAssnTag());

}

namespace butcher {
    DEFINE_ART_MODULE(EventButcher)
}
