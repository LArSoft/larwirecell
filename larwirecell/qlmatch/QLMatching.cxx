#include "QLMatching.h"
#include "util.h"

#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMpointtree.h"
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