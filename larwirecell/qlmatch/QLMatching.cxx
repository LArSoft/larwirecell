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
#include "WireCellUtil/Units.h"

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

  m_inpath = get(cfg, "inpath", m_inpath);
  m_outpath = get(cfg, "outpath", m_outpath);
  m_bee_dir = get(cfg, "bee_dir", m_bee_dir);
}

WireCell::Configuration WireCell::QLMatch::QLMatching::default_configuration() const
{
  Configuration cfg;
  cfg["inpath"] = m_inpath;
  cfg["outpath"] = m_outpath;
  cfg["bee_dir"] = m_bee_dir;
  return cfg;
}

bool WireCell::QLMatch::QLMatching::operator()(const input_vector& invec, output_pointer& out)
{
  out = nullptr;

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
    flashes.push_back(flash);
    log->debug("flash {} time {} total_PE {} num_channels {}",
               flash->get_flash_id(),
               flash->get_time(),
               flash->get_total_PE(),
               flash->get_num_channels());
  }
  log->debug("flashes.size {}", flashes.size());

  // check TimingTPCBundle object
  auto grouping = root_live->value.facade<Grouping>();
  std::vector<Cluster*> clusters = grouping->children();
  std::vector<TimingTPCBundle::pointer> bundles;
  log->debug("checking TimingTPCBundle object");
  /// FIXME: rm this fake code for testing
  std::vector<double> cos_pe_low(nchan, 0.0);
  std::vector<double> cos_pe_mid(nchan, 1e9);
  for (auto flash : flashes) {
    for (size_t icluster = 0; icluster < clusters.size(); ++icluster) {
      Cluster* cluster = clusters[icluster];
      // TimingTPCBundle bundle(flash.get(), cluster, flash->get_flash_id(), icluster);
      TimingTPCBundle::pointer bundle =
        std::make_shared<TimingTPCBundle>(flash.get(), cluster, flash->get_flash_id(), icluster);
      bundles.push_back(bundle);
      /// FIXME: rm this fake code for testing
      if (icluster==0)
      {
        bundle->set_pred_pmt_light(flash->get_PEs());
      }
      bundle->examine_bundle(cos_pe_low, cos_pe_mid);
      log->debug(
        "bundle flash {} icluster{} time {} total_PE {} ks_dis {} chi2 {} ndf {} consistent_flag {}",
        bundle->get_flash()->get_flash_id(),
        icluster,
        bundle->get_flash()->get_time(),
        bundle->get_flash()->get_total_PE(),
        bundle->get_ks_dis(),
        bundle->get_chi2(),
        bundle->get_ndf(),
        bundle->get_consistent_flag());
    }
  }
  log->debug("bundles.size {}", bundles.size());

  // BEE debug direct imaging output and dead blobs
  if (!m_bee_dir.empty()) {
    std::string sub_dir = String::format("%s/%d", m_bee_dir, charge_ident);
    Persist::assuredir(sub_dir);
    QLMatch::dump_bee_3d(
      *root_live.get(),
      String::format("%s/%d-img-apa%d.json", sub_dir, charge_ident, m_anode->ident()));
    QLMatch::dump_bee_flash(
      invec[1], String::format("%s/%d-op-apa%d.json", sub_dir, charge_ident, m_anode->ident()));
  }
  log->debug(em("dump bee"));

  // TODO: actual impl.
  out = invec[0];
  return true;
}