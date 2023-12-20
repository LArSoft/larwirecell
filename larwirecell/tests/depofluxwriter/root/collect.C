#include <algorithm>
#include <functional>
#include <regex>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include "TFile.h"
#include "TH1D.h"

using Hists = std::vector<TH2D*>;

std::vector<std::string> split(std::string s, const std::string& delim)
{
  std::vector<std::string> chunks;
  size_t pos = 0;
  while ((pos = s.find(delim)) != std::string::npos) {
    std::string token = s.substr(0, pos);
    chunks.push_back(token);
    s.erase(0, pos + delim.length());
  }
  chunks.push_back(s);
  return chunks;
}

template <typename T>
std::pair<int, int> chan_mm(const std::vector<T>& things)
{
  int min_chan = -1, max_chan = -1;
  for (const auto& thing : things) {
    const int chid = thing.Channel();
    if (min_chan < 0) min_chan = max_chan = chid;
    if (min_chan > chid) min_chan = chid;
    if (max_chan < chid) max_chan = chid;
  }
  return std::make_pair(min_chan, max_chan);
}

template <typename T>
const auto& getv(gallery::Event& ev, const art::InputTag& tag)
{
  const auto& v = *ev.getValidHandle<std::vector<T>>(tag);
  if (v.empty()) { throw std::runtime_error("failed to get full vector at " + tag.instance()); }
  return v;
}

template <typename T>
const auto& getv(gallery::Event& ev, const std::string& label)
{
  art::InputTag tag{label};
  return getv<T>(ev, tag);
}

struct Collect {

  gallery::Event ev;

  // User options.  Require plane channels to be contiguous,
  // monotonic in each plane and in apa-major
  // order. P:[wN,wN+1,...], apaN:[U,V,W], [apa1,apa2,...]
  const std::vector<int> chaninplanes{800, 800, 960};
  const int nticks{6000};
  const int napas{12};

  Collect(const std::string artfile, const std::string exclude = "") : ev({artfile})
  {
    for (auto e : split(exclude, ",")) {
      excludes_re.push_back(std::regex(e));
    }
  }

private:
  // Hardcoded.
  const int nplanes{3};
  int _nchanperapa = 0;
  std::vector<std::regex> excludes_re;

public:
  // Methods:
  //
  // one_hist_*(tag) assume product is per-apa and return nplanes per-plane histograms.
  //
  // hist_*(tag) assume product is all apa and return napas*nplanes histograms.

  int nchanperapa()
  {
    if (!_nchanperapa) {
      for (auto n : chaninplanes) {
        _nchanperapa += n;
      }
    }
    return _nchanperapa;
  }

  // Return relative number of channel in its apa
  int apa_channel(int chid) { return chid % nchanperapa(); }

  // Return relative number of plane holding of channel in {0,1,2}.
  int channel_plane(int chid)
  {
    int rchid = apa_channel(chid);
    for (int ind = 0; ind < nplanes; ++ind) {
      if (rchid < chaninplanes[ind]) { return ind; }
      rchid -= chaninplanes[ind];
    }
    return -1;
  }

  // Return index of hist holding channel in apa-major ordering.
  int index_channel(int chid)
  {
    int apa = chid / nchanperapa();
    int lch = chid % nchanperapa();

    for (int ind = 0; ind < nplanes; ++ind) {
      if (lch < chaninplanes[ind]) return apa * nplanes + ind;
      lch -= chaninplanes[ind];
    }
    return -1;
  };

  // Reimplement to provide sorting and excluding.
  template <typename T>
  std::vector<art::InputTag> vtags()
  {
    std::vector<art::InputTag> stags;
    for (auto t : ev.getInputTags<std::vector<T>>()) {
      bool skip = false;
      for (const auto& re : excludes_re) {
        std::smatch sre;
        std::string inst = t.instance();
        if (std::regex_match(inst, sre, re)) {
          skip = true;
          break;
        }
      }
      if (skip) {
        std::cerr << "skip excluded tag: " << t << "\n";
        continue;
      }
      stags.push_back(t);
    }

    std::sort(
      stags.begin(), stags.end(), [](const art::InputTag& a, const art::InputTag& b) -> bool {
        return a.instance() < b.instance();
      });
    return stags;
  }

  // Return napa * nplanes histograms.
  Hists make_hists(const std::string& name, const std::string& title)
  {
    Hists hists;
    const std::string uvw = "UVW";
    int min_chan = 0;
    for (int apa = 0; apa < napas; ++apa) {
      for (size_t ind = 0; ind < nplanes; ++ind) {
        int nchans = chaninplanes[ind];
        std::string hname = Form("%s_%c_%d", name.c_str(), uvw[ind], apa);
        std::string htitle = Form("%s - %c-plane, APA %d", title.c_str(), uvw[ind], apa);

        auto hist = new TH2D(
          hname.c_str(), htitle.c_str(), nticks, 0, nticks, nchans, min_chan, min_chan + nchans);
        std::cerr << hname << " ticks=" << nticks << " chans=" << nchans << "\n";
        hists.push_back(hist);
        min_chan += nchans;
      }
    }
    return hists;
  }
  Hists make_hists(const art::InputTag& tag)
  {
    return make_hists("h_" + tag.instance(), tag.instance());
  }

  // Despite the name, return 3 hists, one per plane.  It assumes to
  // be called no more than once for a given tag.
  template <typename T>
  Hists make_one_hist(const art::InputTag& tag)
  {
    const auto& things = getv<T>(ev, tag);

    std::string name = "h_" + tag.instance();
    std::string title = tag.instance();

    Hists hists;
    const std::string uvw = "UVW";

    auto cmm = chan_mm(things);

    const int apa = index_channel(cmm.first) / nplanes;
    int min_chan = apa * nchanperapa();

    std::cerr << "cmm:[" << cmm.first << "," << cmm.second << "] "
              << " apa_chan=[" << apa_channel(cmm.first) << "," << apa_channel(cmm.second) << "]"
              << " chan_pln=[" << channel_plane(cmm.first) << "," << channel_plane(cmm.second)
              << "]"
              << " ind_chan=[" << index_channel(cmm.first) << "," << index_channel(cmm.second)
              << "]"
              << " apa=" << apa << " min_chan=" << min_chan << "\n";

    for (int ind = 0; ind < 3; ++ind) {
      int nchans = chaninplanes[ind];
      std::string hname = Form("%s_%c", name.c_str(), uvw[ind]);
      std::string htitle = Form("%s - %c-plane, APA %d", title.c_str(), uvw[ind], apa);

      auto hist = new TH2D(
        hname.c_str(), htitle.c_str(), nticks, 0, nticks, nchans, min_chan, min_chan + nchans);
      std::cerr << hname << " ticks=" << nticks << " chans=" << nchans << "\n";
      hists.push_back(hist);
      min_chan += nchans;
    }
    return hists;
  }

  Hists one_hist_rawdigit(const art::InputTag& tag)
  {
    std::cerr << "rawdigit: label=" << tag.label() << " instance=" << tag.instance()
              << " process=" << tag.process() << "\n";

    const auto& things = getv<raw::RawDigit>(ev, tag);
    auto hists = make_one_hist<raw::RawDigit>(tag);

    for (const auto& thing : things) {
      const int chid = thing.Channel();
      const int pind = channel_plane(chid);
      if (pind < 0) { continue; }
      auto hist = hists[pind];
      const auto& adcs = thing.ADCs();
      const size_t nadcs = adcs.size();
      for (size_t ind = 0; ind < nadcs; ++ind) {
        hist->Fill(ind + 0.5, chid + 0.5, adcs[ind]);
      }
    }
    return hists;
  }

  // Return vector of ADC histograms, one per apa per plane.
  // label is Event key, name/title are prefixes.
  Hists hist_rawdigit(const std::string& label,
                      const std::string& name = "h_adc",
                      const std::string& title = "ADC")
  {
    Hists hists;

    const auto& rds = getv<raw::RawDigit>(ev, label);
    if (rds.empty()) { return hists; }

    hists = make_hists(name, title);

    for (const auto& rd : rds) {
      const int chid = rd.Channel();
      const int pind = index_channel(chid);
      if (pind < 0) { continue; }

      const auto& adcs = rd.ADCs();
      const size_t nadcs = adcs.size();

      auto hist = hists[pind];
      for (size_t ind = 0; ind < nadcs; ++ind) {
        hist->Fill(ind + 0.5, chid + 0.5, adcs[ind]);
      }
    }
    return hists;
  }

  Hists one_hist_recobwire(const art::InputTag& tag)
  {
    std::cerr << "recobwire: label=" << tag.label() << " instance=" << tag.instance()
              << " process=" << tag.process() << "\n";

    const auto& things = getv<recob::Wire>(ev, tag);
    auto hists = make_one_hist<recob::Wire>(tag);

    for (const auto& thing : things) {
      const int chid = thing.Channel();
      const int pind = channel_plane(chid);
      if (pind < 0) { continue; }
      auto hist = hists[pind];
      auto sig = thing.Signal();
      const size_t ntot = std::min(sig.size(), (size_t)nticks);
      for (size_t ind = 0; ind < ntot; ++ind) {
        hist->Fill(ind + 0.5, chid + 0.5, sig[ind]);
      }
    }
    return hists;
  }

  // Return vector of ADC histograms, one per apa per plane.
  // label is Event key, name/title are prefixes.
  Hists hist_recobwire(const std::string& label,
                       const std::string& name = "h_sig",
                       const std::string& title = "SIG")
  {
    Hists hists;

    const auto& things = getv<recob::Wire>(ev, label);
    if (things.empty()) { return hists; }

    hists = make_hists(name, title);

    for (const auto& thing : things) {
      const int chid = thing.Channel();
      const int pind = index_channel(chid);
      if (pind < 0) { continue; }
      auto hist = hists[pind];

      auto sig = thing.Signal();
      const size_t ntot = std::min(sig.size(), (size_t)nticks);
      for (size_t ind = 0; ind < ntot; ++ind) {
        hist->Fill(ind + 0.5, chid + 0.5, sig[ind]);
      }
    }
    return hists;
  }

  Hists one_hist_simchannel(const art::InputTag& tag)
  {
    const auto& things = getv<sim::SimChannel>(ev, tag);
    auto hists = make_one_hist<sim::SimChannel>(tag);

    for (const auto& thing : things) {
      int chid = thing.Channel();
      const int pind = channel_plane(chid);
      if (pind < 0) { continue; }
      auto hist = hists[pind];

      for (const auto& ts : thing.TDCIDEMap()) {
        const auto tdc = ts.first;
        for (const auto& ide : ts.second) {
          hist->Fill(tdc + 0.5, chid + 0.5, ide.numElectrons);
        }
      }
    }
    return hists;
  }

  Hists hist_simchannel(const art::InputTag& tag)
  {
    std::cerr << "hist_simchannel(" << tag << ")\n";

    Hists hists;

    const auto& scs = getv<sim::SimChannel>(ev, tag);

    hists = make_hists(tag);

    for (const auto& sc : scs) {
      int chid = sc.Channel();
      const int pind = index_channel(chid);
      if (pind < 0) { continue; }
      auto hist = hists[pind];

      for (const auto& ts : sc.TDCIDEMap()) {
        const auto tdc = ts.first;
        for (const auto& ide : ts.second) {
          hist->Fill(tdc + 0.5, chid + 0.5, ide.numElectrons);
        }
      }
    }
    return hists;
  }

}; // Collect

const std::string default_exclude = "dnnsp.*";

void collect(const std::string& artfile = "single_gen_dunefd.root",
             const std::string& histfile = "single_gen_dunefd_hist.root",
             const std::string& exclude = default_exclude);
void collect(const std::string& artfile, const std::string& histfile, const std::string& exclude)
{
  Collect col(artfile, exclude);

  auto ofile = TFile::Open(histfile.c_str(), "recreate");

  for (auto& tag : col.vtags<raw::RawDigit>()) {
    auto hists = col.one_hist_rawdigit(tag);
    ofile->Write();
    for (auto h : hists)
      delete h;
  }

  for (auto& tag : col.vtags<recob::Wire>()) {
    auto hists = col.one_hist_recobwire(tag);
    ofile->Write();
    for (auto h : hists)
      delete h;
  }

  for (auto& tag : col.vtags<sim::SimChannel>()) {
    auto hists = col.hist_simchannel(tag);
    ofile->Write();
    for (auto h : hists)
      delete h;
  }
}

/* Sim job
PROCESS NAME | MODULE LABEL.. | PRODUCT INSTANCE NAME................ | DATA PRODUCT TYPE.................................................... | .SIZE
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCPlaneUInner | std::vector<sim::SimEnergyDeposit>................................... | ....0
SinglesGen.. | generator..... | ..................................... | std::vector<simb::MCTruth>........................................... | ....1
SinglesGen.. | largeant...... | ..................................... | sim::ParticleAncestryMap............................................. | ....-
SinglesGen.. | largeant...... | ..................................... | std::vector<simb::MCParticle>........................................ | ...11
SinglesGen.. | largeant...... | ..................................... | art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo> | ...11
SinglesGen.. | rns........... | ..................................... | std::vector<art::RNGsnapshot>........................................ | ....3
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCInner...... | std::vector<sim::SimEnergyDeposit>................................... | ....6
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCActiveOuter | std::vector<sim::SimEnergyDeposit>................................... | ....0
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCActiveInner | std::vector<sim::SimEnergyDeposit>................................... | 11687
SinglesGen.. | IonAndScint... | ..................................... | std::vector<sim::SimEnergyDeposit>................................... | 11693
SinglesGen.. | tpcrawdecoder. | simpleSC............................. | std::vector<sim::SimChannel>......................................... | ..925
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCPlaneZInner | std::vector<sim::SimEnergyDeposit>................................... | ....0
SinglesGen.. | TriggerResults | ..................................... | art::TriggerResults.................................................. | ....1
SinglesGen.. | tpcrawdecoder. | daq.................................. | std::vector<raw::RawDigit>........................................... | 30720
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCPlaneVInner | std::vector<sim::SimEnergyDeposit>................................... | ....0
 */

/* Sig job
PROCESS NAME | MODULE LABEL.. | PRODUCT INSTANCE NAME................ | DATA PRODUCT TYPE.................................................... | .SIZE
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCPlaneUInner | std::vector<sim::SimEnergyDeposit>................................... | ....0
SinglesGen.. | generator..... | ..................................... | std::vector<simb::MCTruth>........................................... | ....1
SinglesGen.. | largeant...... | ..................................... | sim::ParticleAncestryMap............................................. | ....-
SinglesGen.. | largeant...... | ..................................... | std::vector<simb::MCParticle>........................................ | ...15
SinglesGen.. | tpcrawdecoder. | gauss................................ | std::vector<recob::Wire>............................................. | 30720
SinglesGen.. | largeant...... | ..................................... | art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo> | ...15
SinglesGen.. | rns........... | ..................................... | std::vector<art::RNGsnapshot>........................................ | ....3
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCInner...... | std::vector<sim::SimEnergyDeposit>................................... | ....0
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCActiveOuter | std::vector<sim::SimEnergyDeposit>................................... | ....0
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCActiveInner | std::vector<sim::SimEnergyDeposit>................................... | 11792
SinglesGen.. | IonAndScint... | ..................................... | std::vector<sim::SimEnergyDeposit>................................... | 11792
SinglesGen.. | tpcrawdecoder. | wiener............................... | std::vector<recob::Wire>............................................. | 30720
SinglesGen.. | tpcrawdecoder. | dnnsp................................ | std::vector<recob::Wire>............................................. | 30720
SinglesGen.. | tpcrawdecoder. | simpleSC............................. | std::vector<sim::SimChannel>......................................... | .2249
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCPlaneZInner | std::vector<sim::SimEnergyDeposit>................................... | ....0
SinglesGen.. | TriggerResults | ..................................... | art::TriggerResults.................................................. | ....1
SinglesGen.. | largeant...... | LArG4DetectorServicevolTPCPlaneVInner | std::vector<sim::SimEnergyDeposit>................................... | ....0
 */
