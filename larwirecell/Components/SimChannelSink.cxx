#include "SimChannelSink.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "WireCellGen/GaussianDiffusion.h"
#include "WireCellIface/IDepo.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

WIRECELL_FACTORY(wclsSimChannelSink,
                 wcls::SimChannelSink,
                 wcls::IArtEventVisitor,
                 WireCell::IDepoFilter)

using namespace wcls;
using namespace WireCell;

SimChannelSink::SimChannelSink() : m_depo(nullptr)
{
  m_mapSC.clear();
}

WireCell::Configuration SimChannelSink::default_configuration() const
{
  Configuration cfg;
  cfg["anodes_tn"] = Json::arrayValue;
  cfg["anode"] = "AnodePlane"; // either `anodes_tn` or `anode`
  cfg["rng"] = "Random";
  cfg["artlabel"] = "simpleSC";
  cfg["start_time"] = -1.6 * units::ms;
  cfg["readout_time"] = 4.8 * units::ms;
  cfg["tick"] = 0.5 * units::us;
  cfg["nsigma"] = 3.0;
  cfg["drift_speed"] = 1.098 * units::mm / units::us;
  cfg["uboone_u_to_rp"] = 94 * units::mm;
  cfg["uboone_v_to_rp"] = 97 * units::mm;
  cfg["uboone_y_to_rp"] = 100 * units::mm;
  cfg["u_time_offset"] = 0.0 * units::us;
  cfg["v_time_offset"] = 0.0 * units::us;
  cfg["y_time_offset"] = 0.0 * units::us;
  cfg["g4_ref_time"] = -4050 * units::us; // uboone: -4050us, pdsp: -250us
  cfg["use_energy"] = false;
  cfg["use_extra_sigma"] = false;
  return cfg;
}

void SimChannelSink::configure(const WireCell::Configuration& cfg)
{
  m_anodes.clear();
  auto anodes_tn = cfg["anodes_tn"];
  for (auto anode_tn : anodes_tn) {
    auto anode = Factory::find_tn<IAnodePlane>(anode_tn.asString());
    m_anodes.push_back(anode);
  }

  //  Initialize a SimChannel map. Two assumptions are followed in
  //  larsim::BackTracker::FindCimChannel()
  //  1) fill all channels even with empty SimChannel
  //  2) SimChannels are sorted in channel number
  //  In SimChannelSink::visit(), the map is reinitialized insted of
  //  calling map::clear()
  for (auto& anode : m_anodes) {
    for (auto& channel : anode->channels()) {
      // std::cout << "-> channel: " << channel << std::endl;
      m_mapSC.emplace(channel, sim::SimChannel(channel));
    }
  }
  // std::cout << "m_mapSC.size(): " << m_mapSC.size() << std::endl;

  if (m_anodes.empty()) {
    const std::string anode_tn = cfg["anode"].asString();
    if (anode_tn.empty()) {
      THROW(ValueError() << errmsg{"SimChannelSink requires an anode plane"});
    }
    auto anode = Factory::find_tn<IAnodePlane>(anode_tn);
    m_anodes.push_back(anode);
  }

  const std::string rng_tn = cfg["rng"].asString();
  if (rng_tn.empty()) { THROW(ValueError() << errmsg{"SimChannelSink requires a noise source"}); }
  m_rng = Factory::find_tn<IRandom>(rng_tn);

  m_artlabel = cfg["artlabel"].asString();
  m_artlabel = get(cfg, "artlabel", m_artlabel);
  m_start_time = get(cfg, "start_time", -1.6 * units::ms);
  m_readout_time = get(cfg, "readout_time", 4.8 * units::ms);
  m_tick = get(cfg, "tick", 0.5 * units::us);
  m_nsigma = get(cfg, "nsigma", 3.0);
  m_drift_speed = get(cfg, "drift_speed", 1.098 * units::mm / units::us);
  m_u_to_rp = get(cfg, "uboone_u_to_rp", 94 * units::mm); // compatible interface
  m_v_to_rp = get(cfg, "uboone_v_to_rp", 97 * units::mm);
  m_y_to_rp = get(cfg, "uboone_y_to_rp", 100 * units::mm);
  // compatible interface for protoDUNE-SP
  if (cfg.isMember("u_to_rp")) m_u_to_rp = get(cfg, "u_to_rp", 90.58 * units::mm);
  if (cfg.isMember("v_to_rp")) m_v_to_rp = get(cfg, "v_to_rp", 95.29 * units::mm);
  if (cfg.isMember("y_to_rp")) m_y_to_rp = get(cfg, "y_to_rp", 100.0 * units::mm);

  m_u_time_offset = get(cfg, "u_time_offset", 0.0 * units::us);
  m_v_time_offset = get(cfg, "v_time_offset", 0.0 * units::us);
  m_y_time_offset = get(cfg, "y_time_offset", 0.0 * units::us);
  m_g4_ref_time = get(cfg, "g4_ref_time", -4050 * units::us);
  m_use_energy = get(cfg, "use_energy", false);
  m_use_extra_sigma = get(cfg, "use_extra_sigma", false);
}

void SimChannelSink::produces(art::ProducesCollector& collector)
{
  collector.produces<std::vector<sim::SimChannel>>(m_artlabel);
}

void SimChannelSink::save_as_simchannel(const WireCell::IDepo::pointer& depo)
{
  // Binning tbins(m_readout_time/m_tick, m_start_time, m_start_time+m_readout_time);

  /*  Start the gate ealier for the depos between the response
   *  plane and the anode plane. Those depos are anti-drifted
   *  to the reponse plane, so the start time is earlier.
   *  c.f. jsonnet config in wirecell toolkit: params.sim.ductor
   */
  double response_plane = 10.0 * units::cm;
  double response_time_offset = response_plane / m_drift_speed;
  int response_nticks = (int)(response_time_offset / m_tick);
  Binning tbins(m_readout_time / m_tick + response_nticks,
                m_start_time - response_time_offset,
                m_start_time + m_readout_time);

  if (!depo) return;

  int ctr = 0;
  while (ctr < 1) {
    ctr++;
    //      if(ctr % 10000==0){ std::cout<<"COUNTER "<<ctr<<std::endl;}
    for (auto anode : m_anodes) {
      for (auto face : anode->faces()) {
        auto boundbox = face->sensitive();
        if (!boundbox.inside(depo->pos())) continue;

        // int plane = -1;
        // for(Pimpos* pimpos : {uboone_u, uboone_v, uboone_y}){
        for (auto plane : face->planes()) {
          // plane++;
          int iplane = plane->planeid().index();
          if (iplane < 0) continue;
          const Pimpos* pimpos = plane->pimpos();
          auto& wires = plane->wires();

          const double center_time = depo->time();
          const double center_pitch = pimpos->distance(depo->pos());

          double sigma_L = depo->extent_long();
          if (m_use_extra_sigma) {
            int nrebin = 1;
            double time_slice_width = nrebin * m_drift_speed * m_tick; // units::mm
            double add_sigma_L =
              1.428249 * time_slice_width / nrebin / (m_tick / units::us); // units::mm
            sigma_L =
              sqrt(pow(depo->extent_long(), 2) + pow(add_sigma_L, 2)); // / time_slice_width;
          }
          Gen::GausDesc time_desc(center_time, sigma_L / m_drift_speed);
          // Gen::GausDesc time_desc(center_time, depo->extent_long() / m_drift_speed);
          {
            double nmin_sigma = time_desc.distance(tbins.min());
            double nmax_sigma = time_desc.distance(tbins.max());

            double eff_nsigma = depo->extent_long() / m_drift_speed > 0 ? m_nsigma : 0;
            if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) { break; }
          }

          // auto ibins = pimpos->impact_binning();
          auto wbins = pimpos->region_binning(); // wire binning

          double sigma_T = depo->extent_tran();
          if (m_use_extra_sigma) {
            double add_sigma_T = wbins.binsize();
            if (iplane == 0)
              add_sigma_T *= (0.402993 * 0.3);
            else if (iplane == 1)
              add_sigma_T *= (0.402993 * 0.5);
            else if (iplane == 2)
              add_sigma_T *= (0.188060 * 0.2);
            sigma_T = sqrt(pow(depo->extent_tran(), 2) + pow(add_sigma_T, 2)); // / wbins.binsize();
          }
          Gen::GausDesc pitch_desc(center_pitch, sigma_T);
          // Gen::GausDesc pitch_desc(center_pitch, depo->extent_tran());
          {
            double nmin_sigma = pitch_desc.distance(wbins.min());
            double nmax_sigma = pitch_desc.distance(wbins.max());

            double eff_nsigma = depo->extent_tran() > 0 ? m_nsigma : 0;
            if (nmin_sigma > eff_nsigma || nmax_sigma < -eff_nsigma) { break; }
          }

          auto gd = std::make_shared<Gen::GaussianDiffusion>(depo, time_desc, pitch_desc);
          gd->set_sampling(tbins, wbins, m_nsigma, 0, 1);

          double xyz[3];
          int id = -10000;
          double energy = 100.0;
          if (depo->prior()) {
            id = depo->prior()->id();
            if (m_use_energy) { energy = depo->prior()->energy(); }
          }
          else {
            id = depo->id();
            if (m_use_energy) { energy = depo->energy(); }
          }

          const auto patch = gd->patch();
          const int poffset_bin = gd->poffset_bin();
          const int toffset_bin = gd->toffset_bin();
          const int np = patch.rows();
          const int nt = patch.cols();

          int min_imp = 0;
          int max_imp = wbins.nbins();

          for (int pbin = 0; pbin != np; pbin++) {
            int abs_pbin = pbin + poffset_bin;
            if (abs_pbin < min_imp || abs_pbin >= max_imp) continue;

            auto iwire = wires[abs_pbin];
            int channel = iwire->channel();
            // int channel = abs_pbin;
            // if(plane == 1){ channel = abs_pbin+2400; }
            // if(plane == 2){ channel = abs_pbin+4800; }

            auto channelData = m_mapSC.find(channel);
            sim::SimChannel& sc = (channelData == m_mapSC.end()) ?
                                    (m_mapSC[channel] = sim::SimChannel(channel)) :
                                    channelData->second;

            for (int tbin = 0; tbin != nt; tbin++) {
              int abs_tbin = tbin + toffset_bin;
              double charge = patch(pbin, tbin);
              double tdc = tbins.center(abs_tbin);

              // double wire_response_offset = iwire->center().x() - pimpos->origin().x();
              if (iplane == 0) {
                tdc = tdc + (m_u_to_rp / m_drift_speed) + m_u_time_offset;
                // xyz[0] = depo->pos().x()/units::cm - 94*units::mm/units::cm; // m_u_to_rp/units::cm;
              }
              if (iplane == 1) {
                tdc = tdc + (m_v_to_rp / m_drift_speed) + m_v_time_offset;
                // xyz[0] = depo->pos().x()/units::cm - 97*units::mm/units::cm; // m_v_to_rp/units::cm;
              }
              if (iplane == 2) {
                tdc = tdc + (m_y_to_rp / m_drift_speed) + m_y_time_offset;
                // xyz[0] = depo->pos().x()/units::cm - 100*units::mm/units::cm; // m_y_to_rp/units::cm;
              }
              // xyz[0] = depo->pos().x()/units::cm + wire_response_offset/units::cm;
              WireCell::IDepo::pointer orig = depo_chain(depo).back(); // first depo in the chain
              xyz[0] = orig->pos().x() / units::cm;
              xyz[1] = orig->pos().y() / units::cm;
              xyz[2] = orig->pos().z() / units::cm;

              unsigned int temp_time = (unsigned int)((tdc - m_g4_ref_time) / m_tick);
              charge = abs(charge);
              if (charge > 1) {
                sc.AddIonizationElectrons(
                  id, temp_time, charge, xyz, energy * abs(charge / depo->charge()));
              }
            }
          }
        } // plane
      }   //face
    }     // anode
  }
}

void SimChannelSink::visit(art::Event& event)
{
  std::unique_ptr<std::vector<sim::SimChannel>> out(new std::vector<sim::SimChannel>);

  for (auto& m : m_mapSC) {
    out->emplace_back(m.second);
  }

  event.put(std::move(out), m_artlabel);
  // m_mapSC.clear();
  for (auto& elem : m_mapSC) {
    elem.second = sim::SimChannel(elem.first);
  }
}

bool SimChannelSink::operator()(const WireCell::IDepo::pointer& indepo,
                                WireCell::IDepo::pointer& outdepo)
{
  outdepo = indepo;
  if (indepo) {
    m_depo = indepo;
    save_as_simchannel(m_depo);
  }

  return true;
}
