#include "DepoFluxWriter.h"

#include "WireCellIface/IFieldResponse.h"

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

// fixme: this is an "illegal" header as gen is not a public
// dependency but a provider of a WCT plugin.  To fix this,
// GaussianDiffusion should move into aux....
#include "WireCellGen/GaussianDiffusion.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "DebugDumper.h" // for debug

WIRECELL_FACTORY(wclsDepoFluxWriter,
                 wcls::DepoFluxWriter,
                 wcls::IArtEventVisitor,
                 WireCell::IDepoSetFilter)

using namespace wcls;
using namespace WireCell;

void DepoFluxWriter::produces(art::ProducesCollector& collector)
{
  // fixme: need to extend WireCell_toolkit module to for consumes().
  collector.produces<std::vector<sim::SimChannel>>(m_simchan_label);
}

WireCell::Configuration DepoFluxWriter::default_configuration() const
{
  Configuration cfg;

  // Accept array of strings or single string
  cfg["anodes"] = Json::arrayValue;
  cfg["field_response"] = Json::stringValue;

  // time binning
  const double default_tick = 0.5 * units::us;
  cfg["tick"] = default_tick;
  cfg["window_start"] = 0;
  cfg["window_duration"] = 8096 * default_tick;

  cfg["nsigma"] = m_nsigma;

  cfg["reference_time"] = 0;
  cfg["time_offsets"] = Json::arrayValue;

  cfg["energy"] = m_energy;

  cfg["smear_long"] = 0.0;
  cfg["smear_tran"] = 0.0;

  // input from art::Event
  cfg["sed_label"] = m_sed_label;

  // output to art::Event
  cfg["simchan_label"] = m_simchan_label;
  cfg["sparse"] = true;

  // Provide file name into which validation text is dumped.
  cfg["debug_file"] = m_debug_file;

  return cfg;
}

void DepoFluxWriter::configure(const WireCell::Configuration& cfg)
{
  // Response plane
  auto frname = cfg["field_response"].asString();
  auto ifr = Factory::find_tn<IFieldResponse>(frname);
  const auto& fr = ifr->field_response();
  m_speed = fr.speed;
  m_origin = fr.origin;


  // Anode planes.
  Configuration janode_names = cfg["anodes"];
  if (janode_names.isString()) {
    janode_names = Json::arrayValue;
    janode_names.append(cfg["anodes"]);
  }
  m_anodes.clear();
  for (const auto& jname : janode_names) {
    auto anode = Factory::find_tn<IAnodePlane>(jname.asString());
    m_anodes.push_back(anode);
  }
  if (m_anodes.empty()) {
    THROW(ValueError() << errmsg{"DepoFluxWriter: requires one or more anode planes"});
  }

  // input/output designations
  m_simchan_label = get(cfg, "simchan_label", m_simchan_label);
  m_sed_label = get(cfg, "sed_label", m_sed_label);
  m_debug_file = get(cfg, "debug_file", m_debug_file);

  // sparsity 
  m_sparse = get(cfg, "sparse", true);

  // time-binning
  const double wtick = get(cfg, "tick", 0.5 * units::us);
  const double wstart = cfg["window_start"].asDouble();
  const double wduration = cfg["window_duration"].asDouble();
  const int nwbins = wduration / wtick;
  m_tbins = Binning(nwbins, wstart, wstart + nwbins * wtick);

  // Gaussian cut-off.
  m_nsigma = get(cfg, "nsigma", m_nsigma);

  m_energy = get(cfg, "energy", m_energy);
  m_smear_long = get(cfg, "smear_long", m_smear_long);
  m_smear_tran.clear();
  auto jst = cfg["smear_tran"];
  if (jst.isDouble()) {
    const double s = jst.asDouble();
    m_smear_tran = {s, s, s};
  }
  else if (jst.size() == 3) {
    for (const auto& js : jst) {
      m_smear_tran.push_back(js.asDouble());
    }
  }

  // Last, per-plane, arbitrary time offsets in terms of ticks.
  const double reftime = get(cfg, "reference_time", 0);
  std::vector<double> time_offsets = {-reftime, -reftime, -reftime};
  auto jto = cfg["time_offsets"];
  if (jto.isArray()) {
    if (jto.size() == 3) {
      for (int ind = 0; ind < 3; ++ind) {
        time_offsets[ind] += jto[ind].asDouble();
      }
    }
    else if (!jto.empty()) {
      THROW(ValueError() << errmsg{"DepoFluxWriter: time_offsets must be empty or be a 3-array"});
    }
  }
  m_tick_offsets.clear();
  for (const auto& to : time_offsets) {
    m_tick_offsets.push_back(to / wtick);
  }
}

IAnodeFace::pointer DepoFluxWriter::find_face(const IDepo::pointer& depo)
{
  for (auto anode : m_anodes) {
    for (auto face : anode->faces()) {
      auto bb = face->sensitive();
      if (bb.inside(depo->pos())) { return face; }
    }
  }
  return nullptr;
}

void DepoFluxWriter::visit(art::Event& event)
{
  art::Handle<std::vector<sim::SimEnergyDeposit>> sedvh;
  if (not m_sed_label.empty()) {
    bool okay = event.getByLabel(m_sed_label, sedvh);
    if (!okay) {
      std::string msg =
        "DepoFluxWriter failed to get sim::SimEnergyDeposit from art label: " + m_sed_label;
      std::cerr << msg << std::endl;
      THROW(WireCell::RuntimeError() << WireCell::errmsg{msg});
    }
    sed_dumper(event, m_sed_label, m_debug_file, "DepoFluxWriter ");
    depo_dumper(m_depos, m_simchan_label, m_debug_file, "DepoFluxWriter ");
  }

  std::map<unsigned int, sim::SimChannel> simchans;

  // if dense, initialize SimChannels
  if (!m_sparse){
    for (auto& anode : m_anodes){
      for (auto& channel : anode->channels()) {
        simchans.try_emplace(channel, sim::SimChannel(channel));
      }
    }
  }

  for (const auto& depo : m_depos) {
    if (!depo) continue;

    auto face = find_face(depo);
    if (!face) continue;

    // Depo is at response plane.  Find its time at the collection
    // plane assuming it were to continue along a uniform field.
    // After this, all times are nominal up until we add arbitrary
    // time offsets in delivering electrons to the SimChannel
    const double nominal_depo_time = depo->time() + m_origin / m_speed;

    // Allow for extra smear in time.
    double sigma_L = depo->extent_long(); // [length]
    if (m_smear_long) {
      const double extra = m_smear_long * m_tbins.binsize() * m_speed;
      sigma_L = sqrt(sigma_L * sigma_L + extra * extra);
    }
    Gen::GausDesc time_desc(nominal_depo_time, sigma_L / m_speed);

    // Check if patch is outside time binning
    {
      double nmin_sigma = time_desc.distance(m_tbins.min());
      double nmax_sigma = time_desc.distance(m_tbins.max());

      double eff_nsigma = depo->extent_long() > 0 ? m_nsigma : 0;
      if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) { continue; }
    }

    // Location of original depo.
    WireCell::IDepo::pointer orig = depo_chain(depo).back();
    const double xyz_cm[3] = {
      orig->pos().x() / units::cm, orig->pos().y() / units::cm, orig->pos().z() / units::cm};

    int trackID, origTrackID = -999; // kBogusI
    double energy = m_energy;
    if (depo->prior()) {
      trackID = depo->prior()->id();
      if (!energy) { energy = depo->prior()->energy(); }
    }
    else {
      trackID = depo->id();
      if (!energy) { energy = depo->energy(); }
    }
    if (sedvh.isValid()) { // IDepo::id() is index
      const auto& sed = sedvh->at(trackID);
      trackID = sed.TrackID();
      origTrackID = sed.OrigTrackID();
    }

    // Tabulate depo flux for wire regions from each plane
    for (auto plane : face->planes()) {
      int iplane = plane->planeid().index();
      if (iplane < 0) continue;

      const Pimpos* pimpos = plane->pimpos();
      auto& wires = plane->wires();
      auto wbins = pimpos->region_binning(); // wire binning

      double sigma_T = depo->extent_tran();
      if (m_smear_tran.size()) {
        const double extra = m_smear_tran[iplane] * wbins.binsize();
        sigma_T = sqrt(sigma_T * sigma_T + extra * extra);
      }

      const double center_pitch = pimpos->distance(depo->pos());
      Gen::GausDesc pitch_desc(center_pitch, sigma_T);
      {
        double nmin_sigma = pitch_desc.distance(wbins.min());
        double nmax_sigma = pitch_desc.distance(wbins.max());

        double eff_nsigma = depo->extent_tran() > 0 ? m_nsigma : 0;
        if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) { break; }
      }

      auto gd = std::make_shared<Gen::GaussianDiffusion>(depo, time_desc, pitch_desc);
      gd->set_sampling(m_tbins, wbins, m_nsigma, 0, 1);

      const auto patch = gd->patch();
      const int poffset_bin = gd->poffset_bin();
      const int toffset_bin = gd->toffset_bin();
      const int np = patch.rows();
      const int nt = patch.cols();

      const int min_imp = 0;
      const int max_imp = wbins.nbins();

      for (int pbin = 0; pbin != np; pbin++) {
        const int abs_pbin = pbin + poffset_bin;
        if (abs_pbin < min_imp || abs_pbin >= max_imp) continue;

        auto iwire = wires[abs_pbin];
        const int channel = iwire->channel();

        auto scit = simchans.find(channel);
        sim::SimChannel& sc =
          (scit == simchans.end()) ? (simchans[channel] = sim::SimChannel(channel)) : scit->second;

        for (int tbin = 0; tbin != nt; tbin++) {
          const double charge = std::abs(patch(pbin, tbin));
          if (charge < 1.0) { continue; }

          const unsigned int tdc = tbin + toffset_bin + m_tick_offsets[iplane];

          sc.AddIonizationElectrons(
            trackID, tdc, charge, xyz_cm, energy * abs(charge / depo->charge()), origTrackID);
        } // tbins
      }   // pbins
    }     // plane
  }       // depo
  m_depos.clear();

  auto out = std::make_unique<std::vector<sim::SimChannel>>();
  for (const auto& scit : simchans) {
    out->push_back(scit.second);
  }
  event.put(std::move(out), m_simchan_label);
}

bool DepoFluxWriter::operator()(const WireCell::IDepoSet::pointer& indepos,
                                WireCell::IDepoSet::pointer& outdepos)
{
  outdepos = indepos;

  auto depos = indepos->depos();
  m_depos.insert(m_depos.end(), depos->begin(), depos->end());

  return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
