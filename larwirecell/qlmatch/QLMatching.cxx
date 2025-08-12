#include "QLMatching.h"
#include "util.h"

#include "Opflash.h"
#include "TimingTPCBundle.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellClus/Facade.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Ress.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Units.h"

#include "art/Utilities/make_tool.h"
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larsim/PhotonPropagation/OpticalPathTools/OpticalPath.h"
#include "larsim/PhotonPropagation/SemiAnalyticalModel.h"

WIRECELL_FACTORY(QLMatching,
                 WireCell::QLMatch::QLMatching,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::PointCloud::Facade;

WireCell::QLMatch::QLMatching::QLMatching() : Aux::Logger("QLMatching", "matching") {}

WireCell::QLMatch::QLMatching::~QLMatching() {}

std::vector<std::string> WireCell::QLMatch::QLMatching::input_types()
{
  const std::string tname = std::string(typeid(input_type).name());
  std::vector<std::string> ret(m_multiplicity, tname);
  return ret;
}

void WireCell::QLMatch::QLMatching::configure(const WireCell::Configuration& cfg)
{
  m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());

  // m_VUVEfficiency = get(cfg, "VUVEfficiency", m_VUVEfficiency);
  // m_VISEfficiency = get(cfg, "VISEfficiency", m_VISEfficiency);

  m_inpath = get(cfg, "inpath", m_inpath);
  m_outpath = get(cfg, "outpath", m_outpath);
  m_bee_dir = get(cfg, "bee_dir", m_bee_dir);

  m_pmts = get(cfg, "pmts", m_pmts);
  m_data = get(cfg, "data", m_data);
  m_beamonly = get(cfg, "beamonly", m_beamonly);

  Configuration jch_mask = cfg["ch_mask"];
  for (const auto& jch : jch_mask) {
    m_ch_mask.push_back(jch.asInt());
  }
  m_flash_minPE = get(cfg, "flash_minPE", m_flash_minPE);
  if (m_beamonly) {
    m_flash_mintime = -1e3;
    m_flash_maxtime = 5e3;
  }

  m_QtoL = get(cfg, "QtoL", m_QtoL);
}

WireCell::Configuration WireCell::QLMatch::QLMatching::default_configuration() const
{
  Configuration cfg;

  cfg["inpath"] = m_inpath;
  cfg["outpath"] = m_outpath;
  cfg["bee_dir"] = m_bee_dir;

  cfg["pmts"] = m_pmts;
  cfg["data"] = m_data;
  cfg["beamonly"] = m_beamonly;
  cfg["ch_mask"] = Json::arrayValue;

  cfg["flash_minPE"] = m_flash_minPE;
  cfg["flash_mintime"] = m_flash_mintime;
  cfg["flash_maxtime"] = m_flash_maxtime;
  cfg["QtoL"] = m_QtoL;

  return cfg;
}

bool WireCell::QLMatch::QLMatching::operator()(const input_vector& invec, output_pointer& out)
{
  out = nullptr;
  using WireCell::PointCloud::Facade::float_t;
  // check input size
  if (invec.size() != m_multiplicity) {
    raise<ValueError>("unexpected multiplicity got %d want %d", invec.size(), m_multiplicity);
    return true;
  }

  // boilerplate for EOS handling
  size_t neos = 0;
  for (const auto& in : invec) {
    if (!in) { ++neos; }
  }
  if (neos == invec.size()) {
    // all inputs are EOS, good.
    log->debug("EOS at call {}", m_count++);
    return true;
  }
  if (neos) { raise<ValueError>("missing %d input tensors ", neos); }

  ExecMon em("starting QLMatching");

  // ! start things to move to config block

  std::vector<uint> opdet_mask;
  if (m_pmts)
    opdet_mask = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};

  if (m_data) {
    for (size_t idet = 0; idet < m_ch_mask.size(); idet++) {
      opdet_mask[m_ch_mask[idet]] = 0;
    }
  }

  std::vector<double> m_VUVEfficiency{
    0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.,    0.,
    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.,    0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
    0.021, 0.021, 0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
    0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.,    0.,
    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.,    0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
    0.021, 0.021, 0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
    0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.,    0.,
    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.,    0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
    0.021, 0.021, 0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
    0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.,    0.,
    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.,    0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
    0.021, 0.021, 0.,    0.,    0.,    0.,    0.,    0.,    0.096, 0.096, 0.096, 0.096, 0.096,
    0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.096, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021};
  std::vector<double> m_VISEfficiency{
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.,    0.014, 0.014,
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.152, 0.152, 0.152,
    0.152, 0.152, 0.152, 0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.,
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.,    0.014, 0.014,
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.152, 0.152, 0.152,
    0.152, 0.152, 0.152, 0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.,
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.,    0.014, 0.014,
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.152, 0.152, 0.152,
    0.152, 0.152, 0.152, 0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.,
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.,    0.014, 0.014,
    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.152, 0.152, 0.152,
    0.152, 0.152, 0.152, 0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
    0.,    0.,    0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.104, 0.104, 0.104, 0.104, 0.104,
    0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.,    0.,    0.,    0.,    0.,    0.};

  // ! end things to move to config block

  std::string filename = "semimodel_sbnd.fcl";
  std::string pathvar("FHICL_FILE_PATH");
  const std::string vuv_key = "VUVHits";
  const std::string vis_key = "VISHits";
  cet::filepath_lookup maker(pathvar);
  fhicl::ParameterSet pset = fhicl::ParameterSet::make(filename, maker);
  auto const vuv_pset = pset.get<fhicl::ParameterSet>(vuv_key);
  auto const vis_pset = pset.get<fhicl::ParameterSet>(vis_key);
  auto opticalPath = std::shared_ptr<phot::OpticalPath>(
    std::move(art::make_tool<phot::OpticalPath>(pset.get<fhicl::ParameterSet>("OpticalPathTool"))));
  auto semi_model = std::make_unique<phot::SemiAnalyticalModel>(
    vuv_pset, vis_pset, opticalPath, true, false, false);

  const auto& charge_ts = invec[0];
  const int charge_ident = charge_ts->ident();
  std::string inpath = m_inpath;
  if (inpath.find("%") != std::string::npos) { inpath = String::format(inpath, charge_ident); }

  const auto& charge_tens = *charge_ts->tensors();
  log->debug("charge_tens.size {}", charge_tens.size());
  auto root_live = Aux::TensorDM::as_pctree(charge_tens, inpath + "/live");
  if (!root_live) {
    log->error("Failed to get point cloud tree from \"{}\"", inpath);
    return false;
  }
  log->debug("Got live pctree with {} children", root_live->nchildren());
  log->debug(em("got live pctree"));

  // check Opflahs object
  std::vector<Opflash::pointer> flashes;
  log->debug("checking Opflahs object");
  const auto& tens = invec[1]->tensors();
  if (tens->size() != 1) { raise<ValueError>("Expected 1 tensor, got %d", tens->size()); }
  const auto& ten = tens->at(0);
  if (ten->shape().size() != 2) {
    raise<ValueError>("input tensor dim %d != 2", ten->shape().size());
  }
  const int nrow = ten->shape()[0];
  const int ncol = ten->shape()[1];
  log->debug("nrow {} ncol {}", nrow, ncol);
  const int nchan = ncol - 1;
  // if (nrow < 1) { raise<ValueError>("input tensor nrow %d < 1", nrow); }
  for (int iflash = 0; iflash < nrow; ++iflash) {
    // Opflash flash(ten, iflash, 0.0, nchan);
    Opflash::pointer flash = std::make_shared<Opflash>(ten, iflash, 0.0, nchan);
    if (flash->get_time() < m_flash_mintime || flash->get_time() > m_flash_maxtime) { continue; }
    if (flash->get_total_PE() < m_flash_minPE) { continue; }
    flashes.push_back(flash);
  }

  // check TimingTPCBundle object
  auto grouping = root_live->value.facade<Grouping>();
  std::vector<Cluster*> clusters = grouping->children();
  std::sort(clusters.begin(), clusters.end(), [](const Cluster* cluster1, const Cluster* cluster2) {
    return cluster1->get_length() > cluster2->get_length();
  });
  // create global maps
  std::map<Opflash*, int> global_flash_idx_map;
  std::map<Cluster*, int> global_cluster_idx_map;

  for (size_t i = 0; i < flashes.size(); ++i) {
    global_flash_idx_map[flashes[i].get()] = i;
  }

  for (size_t i = 0; i < clusters.size(); ++i) {
    global_cluster_idx_map[clusters[i]] = i;
  }

  std::vector<TimingTPCBundle::pointer> all_bundles;
  TimingTPCBundleSet pre_bundles;

  uint tpc = m_anode->ident();
  int sign_offset = (tpc == 0) ? -1 : 1;
  double lo_x_bound = (tpc == 0) ? -2000 : 0;
  double hi_x_bound = (tpc == 0) ? 0 : 2000;

  for (size_t idet = 0; idet < opdet_mask.size(); idet++) {
    if ((tpc == 0) && (idet % 2 == 1)) opdet_mask[idet] = 0;
    if ((tpc == 1) && (idet % 2 == 0)) opdet_mask[idet] = 0;
  }

  for (auto flash : flashes) {
    auto flash_time = flash->get_time();
    auto flash_x_offset =
      sign_offset * flash_time * 1.563e-3; // 1.563e-3 is SBND drift velocity in mm/ns

    // per flash mask
    std::vector<uint> flash_opdet_mask = opdet_mask;
    // ! warning: this is a temporary fix to identify simulated saturated PMTs
    for (size_t idet = 0; idet < size_t(flash->get_num_channels()); idet++) {
      auto pe_det = flash->get_PE(idet);
      if ((flash->get_total_PE() > 5000) & (pe_det == 0) & (m_data == false))
        flash_opdet_mask[idet] = 0;
    }

    log->debug("flash time {} flash PE {} flash_x_offset {}",
               int(flash_time) / 100.,
               int(flash->get_total_PE() * 100) / 100.,
               int(flash_x_offset * 100) / 100.);

    for (size_t icluster = 0; icluster < clusters.size(); ++icluster) {
      Cluster* cluster = clusters[icluster];
      // TimingTPCBundle bundle(flash.get(), cluster, flash->get_flash_id(), icluster);
      TimingTPCBundle::pointer bundle =
        std::make_shared<TimingTPCBundle>(flash.get(), cluster, flash->get_flash_id(), icluster);
      all_bundles.push_back(bundle);

      bundle->set_opdet_mask(flash_opdet_mask);

      size_t nopdets = flash->get_num_channels();
      std::vector<double> pred_flash(nopdets, 0.0);

      size_t npt = cluster->nbpoints();
      int npt_outside_drift = 0;
      int npt_outside_bounds = 0;

      bool drifted_outside = false;

      // log->debug("flash {} and cluster {} with {} children", flash->get_flash_id(), icluster, cluster->nchildren());
      std::vector<Blob*> blobs = cluster->children();
      for (auto blob : blobs) {
        auto q = blob->charge() / blob->npoints();
        std::vector<geo_point_t> points = blob->points();

        for (int i = 0; i != blob->npoints(); i++) {
          auto x = points.at(i).x() + flash_x_offset;
          auto y = points.at(i).y();
          auto z = points.at(i).z();

          if (x < lo_x_bound || x > hi_x_bound) {
            npt_outside_drift++;
            continue;
          }
          if (abs(y) > 2000 || z < 0 || z > 5000) {
            npt_outside_bounds++;
            continue;
          }

          if (abs(x) && bundle->get_flag_at_x_boundary() == false)
            bundle->set_flag_close_to_PMT(true);

          if (abs(x) > 1950 && bundle->get_flag_close_to_PMT() == false)
            bundle->set_flag_close_to_PMT(true);

          if (npt_outside_drift > 0.25 * npt) {
            drifted_outside = true;
            break;
          }

          geo::Point_t xyz_cm = {x / 10, y / 10, z / 10};
          std::vector<double> direct_visibilities;
          semi_model->detectedDirectVisibilities(direct_visibilities, xyz_cm);
          std::vector<double> reflected_visibilities;
          semi_model->detectedReflectedVisibilities(reflected_visibilities, xyz_cm);

          // TODO: add opdet channel mask configurable
          for (size_t idet = 0; idet < nopdets; idet++) {
            if (flash_opdet_mask.at(idet) == 0) continue;
            auto dir_vis = direct_visibilities.at(idet);
            auto ref_vis = reflected_visibilities.at(idet);
            auto dir_eff = m_VUVEfficiency.at(idet);
            auto ref_eff = m_VISEfficiency.at(idet);
            pred_flash.at(idet) += q * m_QtoL * dir_vis * dir_eff + q * m_QtoL * ref_vis * ref_eff;
          }
        }
        if (drifted_outside) break;
      }

      // pre-selection stage
      // * if over 20% of the points ar outside the drift volume, skip
      if (drifted_outside) {
        bundle->set_potential_bad_match_flag(true);
        continue;
      }

      bundle->set_pred_flash(pred_flash);
      if (bundle->get_total_pred_light() < 1) continue;
      bundle->examine_bundle();

      if (bundle->get_ks_dis() == 1) {
        bundle->set_potential_bad_match_flag(true);
        continue;
      }
      if (bundle->get_chi2() / bundle->get_ndf() > 1e4) {
        bundle->set_potential_bad_match_flag(true);
        continue;
      }

      log->debug(
        "flash {} and cluster {} meas PE {} pred PE {} npts {} ks_dis {} chi2/ndf {} ndf {}",
        flash->get_flash_id(),
        global_cluster_idx_map[cluster],
        int(flash->get_total_PE() * 100) / 100.,
        int(bundle->get_total_pred_light() * 100) / 100.,
        npt,
        int(bundle->get_ks_dis() * 1000) / 1000.,
        int(bundle->get_chi2() / bundle->get_ndf() * 100) / 100.,
        bundle->get_ndf());

      pre_bundles.insert(bundle);

    } // end first cluster loop

  } // end flash loop
  log->debug("n preselected bundles: {}", pre_bundles.size());

  // * construct maps
  FlashBundlesMap flash_bundles_map;
  ClusterBundlesMap cluster_bundles_map;
  std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer> flash_cluster_bundles_map;

  std::vector<TimingTPCBundle::pointer> consistent_bundles;

  for (auto it = pre_bundles.begin(); it != pre_bundles.end(); ++it) {
    auto bundle = *it;
    auto flash = bundle->get_flash();
    auto cluster = bundle->get_main_cluster();

    if (bundle->get_consistent_flag()) consistent_bundles.push_back(bundle);

    if (flash_bundles_map.find(flash) == flash_bundles_map.end()) {
      std::vector<TimingTPCBundle::pointer> bundle_v;
      bundle_v.push_back(bundle);
      flash_bundles_map[flash] = bundle_v;
    }
    else
      flash_bundles_map[flash].push_back(bundle);
    if (cluster_bundles_map.find(cluster) == cluster_bundles_map.end()) {
      std::vector<TimingTPCBundle::pointer> bundle_v;
      bundle_v.push_back(bundle);
      cluster_bundles_map[cluster] = bundle_v;
    }
    else
      cluster_bundles_map[cluster].push_back(bundle);

    flash_cluster_bundles_map[std::make_pair(flash, cluster)] = bundle;
  } // end pre-selected bundle loop

  TimingTPCBundleSelection to_be_removed;
  for (auto good_bundle : consistent_bundles) {
    auto flash = good_bundle->get_flash();
    auto cluster = good_bundle->get_main_cluster();
    auto flash_bundles = flash_bundles_map[flash];
    auto cluster_bundles = cluster_bundles_map[cluster];

    for (auto bundle : cluster_bundles) {
      if (bundle == good_bundle) continue;
      if (bundle->get_consistent_flag()) continue;
      to_be_removed.push_back(bundle);
    }
  }
  remove_bundle_selection(
    to_be_removed, flash_bundles_map, cluster_bundles_map, flash_cluster_bundles_map);
  remove_bundle_selection(to_be_removed, pre_bundles);

  to_be_removed.clear();

  // set parameters (used for both matching rounds)
  double lambda = 0.1;
  double delta_charge = 0.01;
  double delta_light = 0.025;

  // set "fudge factors" for the total error
  double factor_pe = 0.2;
  double factor_pe_err = 1.0;

  uint nopdet = 0;
  std::vector<int> opdet_idx_v;
  for (size_t idet = 0; idet < opdet_mask.size(); idet++) {
    if (opdet_mask.at(idet) == 1) {
      opdet_idx_v.push_back(int(idet));
      nopdet++;
    }
  }
  log->debug("nopdet {}", nopdet);
  log->debug("opdet_idx_v size {}", opdet_idx_v.size());

  // * first matching round
  {
    uint nbundle = pre_bundles.size();
    uint nflash = flash_bundles_map.size();
    uint ncluster = cluster_bundles_map.size();

    // create map between flash object and flash vector/matrix index
    std::map<Opflash*, int> flash_idx_map;
    // create map between cluster object and cluster vector/matrix index
    std::map<Cluster*, int> cluster_idx_map;

    int cluster_idx = 0;
    int flash_idx = 0;
    for (auto it = cluster_bundles_map.begin(); it != cluster_bundles_map.end(); ++it) {
      auto cluster = it->first;
      auto index = cluster_idx;
      cluster_idx_map[cluster] = index;
      cluster_idx++;
    }
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      auto index = flash_idx;
      flash_idx_map[flash] = index;
      flash_idx++;
    }

    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      auto bundles = it->second;
      for (size_t i = 0; i < bundles.size(); i++) {
        auto bundle = bundles.at(i);
        if (bundle->get_consistent_flag()) {
          log->debug("flash {}, cluster {}, consistent bundle: with ks_dis {}, chi2/ndf {}, ndf "
                     "{}, pred light {}, meas light {}",
                     flash->get_flash_id(),
                     global_cluster_idx_map[bundle->get_main_cluster()],
                     int(bundle->get_ks_dis() * 1000) / 1000.,
                     int(bundle->get_chi2() / bundle->get_ndf() * 100) / 100.,
                     bundle->get_ndf(),
                     int(bundle->get_total_pred_light() * 100) / 100.,
                     int(flash->get_total_PE() * 100) / 100.);
        }
      }
    }

    Ress::vector_t M = Ress::vector_t::Zero(nopdet * nflash);                   // measurement
    Ress::matrix_t P = Ress::matrix_t::Zero(nopdet * nflash, nbundle + nflash); // prediction
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster + nflash);                // measurement flag
    Ress::matrix_t PF =
      Ress::matrix_t::Zero(ncluster + nflash, nbundle + nflash); // prediction flag
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle + nflash);

    log->debug("M dim {}", M.rows());
    log->debug("P dim {} {}", P.rows(), P.cols());
    log->debug("MF dim {} {}", MF.rows(), MF.cols());
    log->debug("PF dim {} {}", PF.rows(), PF.cols());
    log->debug("weights dim {}", weights.rows());

    std::vector<std::pair<Opflash*, Cluster*>> pairs;

    size_t i = 0;  // flash index counter
    size_t ik = 0; // weights index counter
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      auto bundles = it->second;

      for (uint j = 0; j < nopdet; j++) {
        auto opdet_idx = opdet_idx_v.at(j);
        auto pe = flash->get_PE(opdet_idx);
        auto pe_err = sqrt(pow(flash->get_PE(opdet_idx) * factor_pe, 2) +
                           pow(flash->get_PE_err(opdet_idx) * factor_pe_err, 2));

        M(i * nopdet + j) = pe / pe_err;              // measurement term
        P(i * nopdet + j, nbundle + i) = pe / pe_err; // measurement alone term
      }

      for (size_t k = 0; k < bundles.size(); k++) {
        auto bundle = bundles.at(k);
        auto pred_flash = bundle->get_pred_flash();

        for (uint j = 0; j < nopdet; j++) {
          auto opdet_idx = opdet_idx_v.at(j);
          auto pred_pe = pred_flash.at(opdet_idx);
          auto pe_err = sqrt(pow(flash->get_PE(opdet_idx) * factor_pe, 2) +
                             pow(flash->get_PE_err(opdet_idx) * factor_pe_err, 2));
          P(i * nopdet + j, pairs.size()) = pred_pe / pe_err;
        }

        pairs.push_back(std::make_pair(flash, bundle->get_main_cluster()));

        auto meas_pe_tot = flash->get_total_PE();
        auto pred_pe_tot = bundle->get_total_pred_light();
        if (abs(pred_pe_tot - meas_pe_tot) > 0.25 * meas_pe_tot) {
          weights(ik) = 0.25;
          ik++;
        }
        else {
          weights(ik) = 1.0;
          ik++;
        }
      } // loop over bundles in flash
      // MF(ncluster+i) = 0;
      PF(ncluster + i, nbundle + i) = 1. / delta_light;

      flash_idx_map[flash] = nbundle + i;
      i++;
    } // end loop over flashes

    for (uint i = 0; i < nflash; i++) {
      weights(nbundle + i) = 1.0;
    }

    for (uint k = 0; k < ncluster; k++) {
      MF(k) = 1. / delta_charge;
    }

    for (size_t n = 0; n < pairs.size(); n++) {
      auto cluster = pairs.at(n).second;
      PF(cluster_idx_map[cluster], n) = 1. / delta_charge;
    }

    Ress::matrix_t PT = P.transpose();
    Ress::matrix_t PFT = PF.transpose();

    Ress::vector_t y = PT * M + PFT * MF; // predicted x measured + p1/p2 x measured (bi x Mij)
    Ress::matrix_t X = PT * P + PFT * PF; // predicted^2 + p1/p2 x predicted

    Ress::vector_t initial = Ress::vector_t::Zero(nbundle + nflash);
    for (size_t n = 0; n < pairs.size(); n++) {
      auto bundle = flash_cluster_bundles_map[pairs.at(n)];
      auto meas_pe_tot = bundle->get_flash()->get_total_PE();
      auto pred_pe_tot = bundle->get_total_pred_light();
      if (abs(pred_pe_tot - meas_pe_tot) > 0.25 * meas_pe_tot) { initial(n) = 1.0; }
      else //if (bundle->get_flag_close_to_PMT())
        initial(n) = 0.5;
    }
    log->debug("initial dim {}", initial.rows());

    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;

    log->debug("solving");
    Ress::vector_t solution = Ress::solve(X, y, params, initial, weights);

    int n = 0;
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      auto bundles = it->second;
      for (size_t k = 0; k < bundles.size(); k++) {
        auto bundle = bundles.at(k);

        if (solution(n) > 0.05 || m_beamonly)
          log->debug("flash+bundle: flash {}, cluster {}, consistent flag {} solution={}",
                     flash->get_flash_id(),
                     global_cluster_idx_map[bundle->get_main_cluster()],
                     bundle->get_consistent_flag(),
                     solution(n));
        else {
          to_be_removed.push_back(bundle);
        }
        n++;
      }
    }
    int m = 0;
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      if (solution(nbundle + m) != 0)
        log->debug(
          "flash-only: flash {}, solution={}", flash->get_flash_id(), solution(nbundle + m));
      m++;
    }
    remove_bundle_selection(
      to_be_removed, flash_bundles_map, cluster_bundles_map, flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, pre_bundles);
    to_be_removed.clear();
  } // end matching round
  // second matching round
  {
    uint nbundle = pre_bundles.size();
    uint nflash = flash_bundles_map.size();
    uint ncluster = cluster_bundles_map.size();

    // create map between cluster object and cluster vector/matrix index
    // create map between flash object and flash vector/matrix index
    std::map<Cluster*, int> cluster_idx_map;
    std::map<Opflash*, int> flash_idx_map;

    int cluster_idx = 0;
    int flash_idx = 0;
    for (auto it = cluster_bundles_map.begin(); it != cluster_bundles_map.end(); ++it) {
      auto cluster = it->first;
      auto index = cluster_idx;
      cluster_idx_map[cluster] = index;
      cluster_idx++;
    }
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      auto index = flash_idx;
      flash_idx_map[flash] = index;
      flash_idx++;
    }

    Ress::vector_t M = Ress::vector_t::Zero(nopdet * nflash);          // measurement
    Ress::matrix_t P = Ress::matrix_t::Zero(nopdet * nflash, nbundle); // prediction
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster);                // measurement flag
    Ress::matrix_t PF = Ress::matrix_t::Zero(ncluster, nbundle);       // prediction flag
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle);
    std::vector<std::pair<Opflash*, Cluster*>> pairs;

    log->debug("M dim {}", M.rows());
    log->debug("P dim {} {}", P.rows(), P.cols());
    log->debug("MF dim {} {}", MF.rows(), MF.cols());
    log->debug("PF dim {} {}", PF.rows(), PF.cols());
    log->debug("weights dim {}", weights.rows());

    size_t i = 0;  // flash index counter
    size_t ik = 0; // weights index counter
    log->debug("flash_bundles_map size {}", flash_bundles_map.size());
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      auto bundles = it->second;

      for (uint j = 0; j < nopdet; j++) {
        auto opdet_idx = opdet_idx_v.at(j);
        auto pe = flash->get_PE(opdet_idx);
        auto pe_err = sqrt(pow(flash->get_PE(opdet_idx) * factor_pe, 2) +
                           pow(flash->get_PE_err(opdet_idx) * factor_pe_err, 2));

        M(i * nopdet + j) = pe / pe_err; // measurement term
      }

      for (size_t k = 0; k < bundles.size(); k++) {
        auto bundle = bundles.at(k);
        auto pred_flash = bundle->get_pred_flash();

        for (uint j = 0; j < nopdet; j++) {
          auto opdet_idx = opdet_idx_v.at(j);
          auto pred_pe = pred_flash.at(opdet_idx);
          auto pe_err = sqrt(pow(flash->get_PE(opdet_idx) * factor_pe, 2) +
                             pow(flash->get_PE_err(opdet_idx) * factor_pe_err, 2));
          P(i * nopdet + j, pairs.size()) = pred_pe / pe_err;
        }

        pairs.push_back(std::make_pair(flash, bundle->get_main_cluster()));

        auto meas_pe_tot = flash->get_total_PE();
        auto pred_pe_tot = bundle->get_total_pred_light();

        if (abs(pred_pe_tot - meas_pe_tot) > 0.25 * meas_pe_tot) {
          weights(ik) = 0.25;
          ik++;
        }
        else {
          weights(ik) = 1.0;
          ik++;
        }
        // if (bundle->get_flag_close_to_PMT()){
        //   weights(ik) = 0.5;
        //   ik++;
        // }
        // else{
        //   weights(ik) = 1.0;
        //   ik++;
        // }
      } // loop over bundles in flash
      i++;
    } // end loop over flashes

    for (uint k = 0; k < ncluster; k++) {
      MF(k) = 1. / delta_charge;
    }

    for (size_t n = 0; n < pairs.size(); n++) {
      auto cluster = pairs.at(n).second;
      PF(cluster_idx_map[cluster], n) = 1. / delta_charge;
    }

    Ress::matrix_t PT = P.transpose();
    Ress::matrix_t PFT = PF.transpose();

    Ress::vector_t y = PT * M + PFT * MF;
    Ress::matrix_t X = PT * P + PFT * PF;

    Ress::vector_t initial = Ress::vector_t::Zero(nbundle);
    for (size_t n = 0; n < pairs.size(); n++) {
      auto bundle = flash_cluster_bundles_map[pairs.at(n)];
      auto meas_pe_tot = bundle->get_flash()->get_total_PE();
      auto pred_pe_tot = bundle->get_total_pred_light();
      if (abs(pred_pe_tot - meas_pe_tot) > 0.25 * meas_pe_tot) { initial(n) = 1.0; }
      else //if (bundle->get_flag_close_to_PMT())
        initial(n) = 0.5;
    }
    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;

    log->debug("solving");
    Ress::vector_t solution = Ress::solve(X, y, params, initial, weights);
    log->debug("solution size {}", solution.size());
    int n = 0;
    for (auto it = flash_bundles_map.begin(); it != flash_bundles_map.end(); ++it) {
      auto flash = it->first;
      auto bundles = it->second;
      for (size_t k = 0; k < bundles.size(); k++) {
        auto bundle = bundles.at(k);

        if (solution(n) > 0.05 || m_beamonly) {
          log->debug("flash+bundle: flash {}, cluster {} time {} meas PE {}, pred PE {}, "
                     "solution={}, ks_dis {}, chi2/ndf {}, consistent {}",
                     flash->get_flash_id(),
                     global_cluster_idx_map[bundle->get_main_cluster()],
                     int(flash->get_time()) / 1e3,
                     int(flash->get_total_PE() * 100) / 100.,
                     int(bundle->get_total_pred_light() * 100) / 100.,
                     int(solution(n) * 1e4) / 10000.,
                     int(bundle->get_ks_dis() * 1000) / 1000.,
                     int(bundle->get_chi2() / bundle->get_ndf() * 100) / 100.,
                     bundle->get_consistent_flag());
        }
        else {
          to_be_removed.push_back(bundle);
        }
        n++;
      }
    }
    remove_bundle_selection(
      to_be_removed, flash_bundles_map, cluster_bundles_map, flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, pre_bundles);
    to_be_removed.clear();

  } // end second matching round

  // BEE debug direct imaging output and dead blobs
  log->debug("done with matching");
  if (!m_bee_dir.empty()) {
    std::string sub_dir = String::format("%s/%d", m_bee_dir, charge_ident);
    Persist::assuredir(sub_dir);
    QLMatch::dump_bee_3d(
      *root_live.get(),
      String::format("%s/%d-img-apa%d.json", sub_dir, charge_ident, m_anode->ident()));
    // QLMatch::dump_bee_flash(
    //   invec[1], String::format("%s/%d-op-apa%d.json", sub_dir, charge_ident, m_anode->ident()));
    QLMatch::dump_bee_bundle(
      flash_bundles_map,
      global_cluster_idx_map,
      String::format("%s/%d-op-apa%d.json", sub_dir, charge_ident, m_anode->ident()));
  }
  log->debug(em("dump bee"));

  // TODO: actual impl.
  out = invec[0];
  return true;
}

void WireCell::QLMatch::QLMatching::remove_bundle_selection(TimingTPCBundleSelection to_be_removed,
                                                            TimingTPCBundleSet& bundle_set)
{
  for (auto it = to_be_removed.begin(); it != to_be_removed.end(); ++it) {
    auto rm_bundle = *it;
    bundle_set.erase(rm_bundle);
  }
}

void WireCell::QLMatch::QLMatching::remove_bundle_selection(
  TimingTPCBundleSelection to_be_removed,
  FlashBundlesMap& flash_bundles_map,
  ClusterBundlesMap& cluster_bundles_map,
  std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer>& flash_cluster_bundles_map)
{
  for (auto it = to_be_removed.begin(); it != to_be_removed.end(); ++it) {
    auto rm_bundle = *it;
    auto rm_flash = rm_bundle->get_flash();
    auto rm_cluster = rm_bundle->get_main_cluster();

    flash_cluster_bundles_map.erase(std::make_pair(rm_flash, rm_cluster));
    {
      auto& temp_bundles = flash_bundles_map[rm_flash];
      temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), rm_bundle));
      if (temp_bundles.size() == 0) flash_bundles_map.erase(rm_flash);
    }
    {
      auto& temp_bundles = cluster_bundles_map[rm_cluster];
      temp_bundles.erase(find(temp_bundles.begin(), temp_bundles.end(), rm_bundle));
      if (temp_bundles.size() == 0) cluster_bundles_map.erase(rm_cluster);
    }
  }
}