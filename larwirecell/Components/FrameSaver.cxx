/*
 */

#include "FrameSaver.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/ITrace.h"
#include "WireCellUtil/NamedFactory.h"

// it would be nice to remove this dependency since it is needed only
// to add bogus information to raw::RawDigit.
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"

#include <algorithm>
#include <map>

WIRECELL_FACTORY(wclsFrameSaver, wcls::FrameSaver, wcls::IArtEventVisitor, WireCell::IFrameFilter)

using namespace wcls;
using namespace WireCell;

FrameSaver::FrameSaver() : m_frame(nullptr), m_nticks(0) {}

FrameSaver::~FrameSaver() {}

WireCell::Configuration FrameSaver::default_configuration() const
{
  Configuration cfg;
  cfg["anode"] = "AnodePlane";

  // If true, truncate frame waveforms and save as raw::RawDigit,
  // else leave as floating point recob::Wire
  cfg["digitize"] = false;

  // If true, save recob::Wires as sparse (true zero suppressed).
  // Note, this sparsification is performed independently from if
  // the input IFrame itself is sparse or not.
  cfg["sparse"] = true;

  // if FALSE (DEFAULT BEHAVIOUR) continue saving RawDigit frame
  // if true, save an empty frame (RawDigit), used for saving CMM only
  cfg["skip_frame"] = false;

  // Provide a configurable translation layer between WCT Plane View IDs
  // and those from larsoft. Default is to assume they're the same,
  // except for in certain cases (such as the vertical drift geometry)
  // where the map is provided in the jsonnet/json
  cfg["plane_map"][std::to_string((int)WireCell::kUlayer)] = (int)geo::kU;
  cfg["plane_map"][std::to_string((int)WireCell::kVlayer)] = (int)geo::kV;
  cfg["plane_map"][std::to_string((int)WireCell::kWlayer)] = (int)geo::kW;

  // Transform channels to have differe plane IDs.
  // This function is provided in case the plane IDs have been
  // overwritten by WCT.
  cfg["channels_transform"]["to_U"] = Json::arrayValue;
  cfg["channels_transform"]["to_V"] = Json::arrayValue;
  cfg["channels_transform"]["to_W"] = Json::arrayValue;

  // If digitize, raw::RawDigit has slots for pedestal mean and
  // sigma.  Legacy/obsolete code stuff unrelated values into these
  // slots.  If pedestal_mean is a number, it will be stuffed.  If
  // it is the string "fiction" then the DetPedestalService will be
  // used to set some unrelated pedestal value.
  cfg["pedestal_mean"] = 0.0;
  // Stuff some value into pedestal sigma.  This value has nothing
  // to do with the produced NF'ed waveforms.
  cfg["pedestal_sigma"] = 0.0;

  // frames to output, if any
  cfg["frame_tags"] = Json::arrayValue;
  cfg["frame_scale"] = 1.0; // multiply this number to all
    // waveform samples.  If list, then
    // one number per frame tags.
  // If nonzero, force number of ticks in output waveforms.
  // If zero, use whatever input data has.
  // If -1, use value as per LS's detector properties service.
  cfg["nticks"] = m_nticks;

  // Summaries to output, if any
  cfg["summary_tags"] = Json::arrayValue;
  cfg["summary_scale"] = 1.0;
  cfg["summary_suffix"] = "summary";
  // Sumaries are *per trace* quantities coming in but it is likely
  // that some (most?) consumers of the output will expect *per
  // channel* quantities.  Aggregating by channel requires some
  // operator to be applied to the sequence of summary values from a
  // common channel.  The operator is defined as an object keyed by
  // the summary tag.  Values may be "set" which will simply save a
  // summary value to the output element (last one from a channel
  // wins) or "sum" (the default) which will add up the values.
  cfg["summary_operator"] = Json::objectValue;

  // Names of channel mask maps to save, if any.
  cfg["chanmaskmaps"] = Json::arrayValue;
  cfg["cmm_masks_suffix"] = "masks";
  cfg["cmm_channels_suffix"] = "channels";

  return cfg;
}

static float summary_sum(const std::vector<float>& tsvals)
{
  return std::accumulate(tsvals.begin(), tsvals.end(), 0);
}
static float summary_set(const std::vector<float>& tsvals)
{
  if (tsvals.empty()) { return 0.0; }
  return tsvals.back();
}

void FrameSaver::configure(const WireCell::Configuration& cfg)
{
  // Populate the channel2layer transform
  std::map<int, std::string> channel2layer;
  std::vector<std::string> tranforms = {"to_U", "to_V", "to_W"};
  for (const auto& key: tranforms) {
      int layerValue = (key == "to_U") ? WireCell::kUlayer :
                       (key == "to_V") ? WireCell::kVlayer :
                       WireCell::kWlayer;
      for (const auto& ch: cfg["channels_transform"][key]) {
        channel2layer[ch.asInt()] = std::to_string(layerValue);
      }
  }

  const std::string anode_tn = cfg["anode"].asString();
  if (anode_tn.empty()) { THROW(ValueError() << errmsg{"FrameSaver requires an anode plane"}); }

  WireCell::IAnodePlane::pointer anode = Factory::find_tn<IAnodePlane>(anode_tn);
  for (auto chid : anode->channels()) {

    auto wpid = anode->resolve(chid);
    geo::View_t view;


    // Use configurable translation between WCT and larsoft
    // plane view IDs. Relevant especially for VD 3 view
    // since the 2nd induction plane is actually labelled
    // kY in larsoft vs kV in WCT
    // Unless otherwise specified, this map amounts to
    // kU->kU, kV->kV, kW->kW
    std::string wct_layer = std::to_string((int)wpid.layer());
    // Overwrite wct_layer if necessary
    if (!channel2layer.empty() && channel2layer.count(chid)) {
        wct_layer = channel2layer[chid];
    }
    view = (geo::View_t)(cfg["plane_map"][wct_layer].asInt());

    m_chview[chid] = view;
  }

  m_digitize = get(cfg, "digitize", false);
  m_sparse = get(cfg, "sparse", true);
  m_skipframe = get(cfg, "skip_frame", false);

  m_cmms = cfg["chanmaskmaps"];
  m_cmm_masks_suffix = cfg["cmm_masks_suffix"].asString();
  m_cmm_channels_suffix = cfg["cmm_channels_suffix"].asString();

  m_pedestal_mean = cfg["pedestal_mean"];
  m_pedestal_sigma = get(cfg, "pedestal_sigma", 0.0);

  if (!cfg["frame_scale"].isNull()) {
    m_frame_scale.clear();
    auto jscale = cfg["frame_scale"];
    m_frame_tags.clear();
    auto jtags = cfg["frame_tags"];
    for (auto jtag : jtags) {
      std::string tag = jtag.asString();

      // get any waveform scaling for this frame
      const int ind = m_frame_tags.size();
      double scale = 1.0;
      //std::cerr << "JSCALE:" << jscale << std::endl;
      if (!jscale.isNull()) {
        if (jscale.isArray()) { scale = jscale[ind].asDouble(); }
        else {
          scale = jscale.asDouble();
        }
      }

      m_frame_scale.push_back(scale);
      m_frame_tags.push_back(tag);
    }
    m_nticks = get(cfg, "nticks", m_nticks);
  }

  auto jso = cfg["summary_operator"];

  if (!cfg["summary_tags"].isNull()) {
    m_summary_scale.clear();
    auto jscale = cfg["summary_scale"];
    m_summary_tags.clear();
    auto jtags = cfg["summary_tags"];
    m_summary_suffix.clear();
    m_summary_suffix = cfg["summary_suffix"].asString();
    for (auto jtag : jtags) {
      std::string tag = jtag.asString();

      const int ind = m_summary_tags.size();
      double scale = 1.0;
      if (!jscale.isNull()) {
        if (jscale.isArray()) { scale = jscale[ind].asDouble(); }
        else {
          scale = jscale.asDouble();
        }
      }

      m_summary_tags.push_back(tag);
      m_summary_scale.push_back(scale);
      auto so = get<std::string>(jso, tag, "sum");
      if (so == "set") { m_summary_operators[tag] = summary_set; }
      else {
        m_summary_operators[tag] = summary_sum;
      }
    }
  }
}

void FrameSaver::produces(art::ProducesCollector& collector)
{
  for (auto tag : m_frame_tags) {
    if (!m_digitize && !m_skipframe) {
      std::cerr << "wclsFrameSaver: promising to produce recob::Wires named \"" << tag << "\"\n";
      collector.produces<std::vector<recob::Wire>>(tag);
    }
    else if (!m_skipframe) {
      std::cerr << "wclsFrameSaver: promising to produce raw::RawDigits named \"" << tag << "\"\n";
      collector.produces<std::vector<raw::RawDigit>>(tag);
    }
  }
  for (auto tag : m_summary_tags) {
    std::cerr << "wclsFrameSaver: promising to produce channel summary named \"" << tag << "\"\n";
    collector.produces<std::vector<double>>(tag + m_summary_suffix);
  }
  for (auto cmm : m_cmms) {
    const std::string cmm_name = cmm.asString();
    std::cerr << "wclsFrameSaver: promising to produce channel masks named \"" << cmm_name
              << "\"\n";
    collector.produces<channel_list>(cmm_name + m_cmm_channels_suffix);
    collector.produces<channel_masks>(cmm_name + m_cmm_masks_suffix);
  }
}

static void tagged_traces(IFrame::pointer frame, std::string tag, ITrace::vector& ret)
{
  auto const& all_traces = frame->traces();
  const auto& ttinds = frame->tagged_traces(tag);
  if (ttinds.size()) {
    for (size_t index : ttinds) {
      ret.push_back(all_traces->at(index));
    }
    return;
  }
  auto ftags = frame->frame_tags();
  if (std::find(ftags.begin(), ftags.end(), tag) == ftags.end()) { return; }
  ret.insert(ret.begin(), all_traces->begin(), all_traces->end());
}

typedef std::unordered_map<int, ITrace::vector> traces_bychan_t;
static void traces_bychan(ITrace::vector& traces, traces_bychan_t& ret)
{
  for (auto trace : traces) {
    int chid = trace->channel();
    ret[chid].push_back(trace);
  }
}

// Issolate some silly legacy shenanigans to keep the rest of the code
// blissfully ignorant of the evilness this implies.
struct PU {
  Json::Value pu;

  PU(Json::Value pu) : pu(pu) {}

  float operator()(int chid)
  {
    if (pu.isNumeric()) { return pu.asFloat(); }
    if (pu.asString() == "fiction") {
      art::ServiceHandle<lariov::DetPedestalService const> dps;
      const auto& pv = dps->GetPedestalProvider();
      return pv.PedMean(chid);
    }
    return 0.0;
  }
};

void FrameSaver::save_as_raw(art::Event& event)
{
  int nticks_want = m_nticks;
  if (nticks_want < 0) {
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);
    nticks_want = detProp.NumberTimeSamples();
  }

  size_t nftags = m_frame_tags.size();
  for (size_t iftag = 0; iftag < nftags; ++iftag) {
    std::string ftag = m_frame_tags[iftag];
    std::cerr << "wclsFrameSaver: saving raw::RawDigits tagged \"" << ftag << "\"\n";

    double scale = m_frame_scale[iftag];

    ITrace::vector traces;
    tagged_traces(m_frame, ftag, traces);
    traces_bychan_t bychan;
    traces_bychan(traces, bychan);

    std::unique_ptr<std::vector<raw::RawDigit>> out(new std::vector<raw::RawDigit>);

    for (auto chv : m_chview) {
      const int chid = chv.first;
      const auto& traces = bychan[chid];
      const size_t ntraces = traces.size();

      int tbin = 0;
      std::vector<float> charge;
      if (ntraces) {
        auto trace = traces[0];
        tbin = trace->tbin();
        charge = trace->charge();
      }
      // charge may be empty here

      // enforce number of ticks if we are so configured.
      size_t ncharge = charge.size();
      int nticks = tbin + ncharge;
      if (nticks_want) { // force output waveform size
        if (nticks_want < nticks) { ncharge = nticks_want - tbin; }
        nticks = nticks_want;
      }
      raw::RawDigit::ADCvector_t adcv(nticks);
      for (size_t ind = 0; ind < ncharge; ++ind) {
        adcv[tbin + ind] = scale * charge[ind]; // scale + truncate/redigitize
      }
      out->emplace_back(raw::RawDigit(chid, nticks, adcv, raw::kNone));
      if (m_pedestal_mean.asString() == "native") {
        short baseline = Waveform::most_frequent(adcv);
        out->back().SetPedestal(baseline, m_pedestal_sigma);
      }
      else {
        PU pu(m_pedestal_mean);
        out->back().SetPedestal(pu(chid), m_pedestal_sigma);
      }
    }
    event.put(std::move(out), ftag);
  }
}

void FrameSaver::save_as_cooked(art::Event& event)
{
  int nticks_want = m_nticks;
  if (nticks_want < 0) {
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);
    nticks_want = detProp.NumberTimeSamples();
    std::cerr << "wclsFrameSaver saving cooked to " << nticks_want << " ticks\n";
  }

  size_t nftags = m_frame_tags.size();
  for (size_t iftag = 0; iftag < nftags; ++iftag) {
    std::string ftag = m_frame_tags[iftag];

    double scale = m_frame_scale[iftag];

    ITrace::vector traces;
    tagged_traces(m_frame, ftag, traces);
    if (traces.empty()) {
      std::cerr << "wclsFrameSaver: no traces tagged \"" << ftag << "\"\n";
      // fall through loop so we put (empty) outwires
    }
    else {
      std::cerr << "wclsFrameSaver: saving " << traces.size() << " traces tagged \"" << ftag
                << "\"\n";
    }

    traces_bychan_t bychan;
    traces_bychan(traces, bychan);

    std::unique_ptr<std::vector<recob::Wire>> outwires(new std::vector<recob::Wire>);

    double total_charge = 0.0;
    int total_samples = 0;

    for (auto chv : m_chview) {
      const int chid = chv.first;
      const auto& traces = bychan[chid];

      recob::Wire::RegionsOfInterest_t rois(nticks_want);

      for (const auto& trace : traces) {
        const int tbin = trace->tbin();
        const auto& charge = trace->charge();

        auto beg = charge.begin();
        const auto first = beg;
        auto end = charge.end();
        if (nticks_want) { // user set waveform size
          if (tbin >= nticks_want) { beg = end; }
          else {
            int backup = tbin + charge.size() - nticks_want;
            if (backup > 0) { end -= backup; }
          }
        }
        if (beg >= end) {
          std::cerr << "wclsFrameSaver: no samples within desired window for channel " << chid
                    << "\n";
          continue;
        }
        if (!m_sparse) {
          std::vector<float> scaled(beg, end);
          for (size_t i = 0; i < scaled.size(); ++i)
            scaled[i] *= scale;
          // prefer combine_range() but it segfaults.
          rois.add_range(tbin, scaled.begin(), scaled.end());
          continue;
        }
        // sparsify trace whether or not it may itself already
        // represents a sparse ROI
        while (true) {
          beg = std::find_if(beg, end, [](float v) { return v != 0.0; });
          if (beg == end) { break; }
          auto mid = std::find_if(beg, end, [](float v) { return v == 0.0; });
          std::vector<float> scaled(beg, mid);
          for (int ind = 0; ind < mid - beg; ++ind) {
            scaled[ind] *= scale;
            total_charge += scaled[ind];
            ++total_samples;
          }
          rois.add_range(tbin + beg - first, scaled.begin(), scaled.end());
          beg = mid;
        }
      }

      const geo::View_t view = chv.second;
      outwires->emplace_back(recob::Wire(rois, chid, view));
    }
    std::cerr << "FrameSaver: q=" << total_charge << " n=" << total_samples << " tag=" << ftag
              << "\n";
    event.put(std::move(outwires), ftag);
  } // loop over tags
}

void FrameSaver::save_summaries(art::Event& event)
{
  const int ntags = m_summary_tags.size();
  if (0 == ntags) {
    return; // no tags
  }

  const size_t nchans = m_chview.size();

  // for each summary
  for (int tag_ind = 0; tag_ind < ntags; ++tag_ind) {
    // The scale set for the tag.
    const double scale = m_summary_scale[tag_ind];

    std::unique_ptr<std::vector<double>> outsum(new std::vector<double>(nchans, 0.0));

    // The "summary" and "traces" vectors of the same tag are
    // synced, element-by-element.  Each element corresponds to
    // one trace (ROI).  No particular order or correlation by
    // channel exists, and that's what the rest of this code
    // creates.
    auto tag = m_summary_tags[tag_ind];
    const auto& summary = m_frame->trace_summary(tag);
    ITrace::vector traces;
    tagged_traces(m_frame, tag, traces);
    const size_t ntraces = traces.size();

    std::unordered_map<int, std::vector<float>> bychan;
    for (size_t itrace = 0; itrace < ntraces; ++itrace) {
      const int chid = traces[itrace]->channel();
      const double summary_value = summary[itrace];
      bychan[chid].push_back(summary_value);
    }
    auto oper = m_summary_operators[tag];

    size_t chanind = 0;
    for (auto chv : m_chview) {
      const int chid = chv.first;
      const float val = oper(bychan[chid]);
      outsum->at(chanind) = val * scale;
      ++chanind;
    }
    event.put(std::move(outsum), tag + m_summary_suffix);
  }
}

void FrameSaver::save_cmms(art::Event& event)
{
  if (m_cmms.isNull()) { return; }
  if (!m_cmms.isArray()) {
    std::cerr
      << "wclsFrameSaver: wrong type for configuration array of channel mask maps to save\n";
    return;
  }
  for (auto jcmm : m_cmms) {
    std::string name = jcmm.asString();
    std::unique_ptr<channel_list> out_list(new channel_list);
    std::unique_ptr<channel_masks> out_masks(new channel_masks);

    auto cmm = m_frame->masks();
    auto it = cmm.find(name);
    if (it == cmm.end()) {
      std::cerr << "wclsFrameSaver: failed to find requested channel masks \"" << name << "\"\n";
    }
    else {
      for (auto cmit : it->second) { // int->vec<pair<int,int>>
        out_list->push_back(cmit.first);
        for (auto be : cmit.second) {
          out_masks->push_back(cmit.first);
          out_masks->push_back(be.first);
          out_masks->push_back(be.second);
        }
      }
    }

    if (out_list->empty()) {
      std::cerr << "wclsFrameSaver: found empty channel masks for \"" << name << "\"\n";
    }
    event.put(std::move(out_list), name + m_cmm_channels_suffix);
    event.put(std::move(out_masks), name + m_cmm_masks_suffix);
  }
}

void FrameSaver::save_empty(art::Event& event)
{
  // art (apparently?) requires something to be saved if a produces() is promised.
  std::cerr << "wclsFrameSaver: saving empty frame to art::Event\n";

  for (auto ftag : m_frame_tags) {
    if (m_digitize) {
      std::unique_ptr<std::vector<raw::RawDigit>> out(new std::vector<raw::RawDigit>);
      event.put(std::move(out), ftag);
    }
    else {
      std::unique_ptr<std::vector<recob::Wire>> outwires(new std::vector<recob::Wire>);
      event.put(std::move(outwires), ftag);
    }
  }

  for (auto stag : m_summary_tags) {
    std::unique_ptr<std::vector<double>> outsum(new std::vector<double>);
    event.put(std::move(outsum), stag + m_summary_suffix);
  }

  for (auto jcmm : m_cmms) {
    std::string name = jcmm.asString();
    std::unique_ptr<channel_list> out_list(new channel_list);
    std::unique_ptr<channel_masks> out_masks(new channel_masks);
    event.put(std::move(out_list), name + m_cmm_channels_suffix);
    event.put(std::move(out_masks), name + m_cmm_masks_suffix);
  }
}

void FrameSaver::visit(art::Event& event)
{
  if (!m_frame) {
    save_empty(event);
    return;
  }

  if (m_digitize && !m_skipframe) { save_as_raw(event); }
  else if (!m_digitize && !m_skipframe) {
    save_as_cooked(event);
  }

  save_summaries(event);

  save_cmms(event);

  m_frame = nullptr; // done with stashed frame
}

bool FrameSaver::operator()(const WireCell::IFrame::pointer& inframe,
                            WireCell::IFrame::pointer& outframe)
{
  // set an IFrame based on last visited event.
  outframe = inframe;
  if (inframe) {
    if (m_frame) {
      std::cerr
        << "wclsFrameSaver: warning: dropping prior frame.  Fixme to handle queue of frames.\n";
    }
    // else {
    //     std::cerr << "wclsFrameSaver got frame\n";
    // }
    m_frame = inframe;
  }
  // else {
  //     std::cerr << "wclsFrameSaver sees EOS\n";
  // }
  return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
