#include "Truth2h5.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"

#include "TLorentzVector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <hdf5.h>

WIRECELL_FACTORY(Truth2h5,
                 WireCell::AIML::Truth2h5,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

using namespace WireCell;

namespace {
  hid_t neutrino_type()
  {
    static hid_t t = -1;
    if (t >= 0) { return t; }
    t = H5Tcreate(H5T_COMPOUND, sizeof(WireCell::AIML::Truth2h5::NeutrinoInfo));
    H5Tinsert(t, "nu_pdg", HOFFSET(WireCell::AIML::Truth2h5::NeutrinoInfo, pdg), H5T_NATIVE_INT);
    H5Tinsert(t, "nu_ccnc", HOFFSET(WireCell::AIML::Truth2h5::NeutrinoInfo, ccnc), H5T_NATIVE_INT);
    H5Tinsert(
      t, "nu_intType", HOFFSET(WireCell::AIML::Truth2h5::NeutrinoInfo, intType), H5T_NATIVE_INT);
    H5Tinsert(
      t, "nu_energy", HOFFSET(WireCell::AIML::Truth2h5::NeutrinoInfo, energy), H5T_NATIVE_DOUBLE);
    H5Tinsert(
      t, "nu_vertex_x", HOFFSET(WireCell::AIML::Truth2h5::NeutrinoInfo, vx), H5T_NATIVE_DOUBLE);
    H5Tinsert(
      t, "nu_vertex_y", HOFFSET(WireCell::AIML::Truth2h5::NeutrinoInfo, vy), H5T_NATIVE_DOUBLE);
    H5Tinsert(
      t, "nu_vertex_z", HOFFSET(WireCell::AIML::Truth2h5::NeutrinoInfo, vz), H5T_NATIVE_DOUBLE);
    return t;
  }
} // namespace

AIML::Truth2h5::Truth2h5()
  : Aux::Logger("Truth2h5", "aiml")
  , m_mctruth_label("generator")
  , m_output_file("g4-tru.h5")
  , m_file(-1)
{}

AIML::Truth2h5::~Truth2h5()
{
  if (m_file >= 0) {
    H5Fclose(m_file);
    m_file = -1;
  }
}

Configuration AIML::Truth2h5::default_configuration() const
{
  Configuration cfg;
  cfg["mctruth_label"] = m_mctruth_label;
  cfg["output_file"] = m_output_file;
  return cfg;
}

void AIML::Truth2h5::configure(const Configuration& cfg)
{
  m_mctruth_label = get(cfg, "mctruth_label", m_mctruth_label);
  m_output_file = get(cfg, "output_file", m_output_file);
  reset();
  if (m_file >= 0) {
    H5Fclose(m_file);
    m_file = -1;
  }
}

void AIML::Truth2h5::reset()
{
  m_info = NeutrinoInfo{};
}

void AIML::Truth2h5::visit(art::Event& event)
{
  reset();
  if (m_mctruth_label.empty()) {
    log->debug("Truth2h5: MCTruth label not configured");
    return;
  }

  art::Handle<std::vector<simb::MCTruth>> mctruth_handle;
  if (!event.getByLabel(art::InputTag{m_mctruth_label}, mctruth_handle)) {
    log->warn("Truth2h5 failed to fetch MCTruth with label '{}'", m_mctruth_label);
    return;
  }

  if (!mctruth_handle.isValid() || mctruth_handle->empty()) {
    log->warn("Truth2h5 MCTruth handle for '{}' is empty", m_mctruth_label);
    return;
  }

  const auto& mctruth = mctruth_handle->front();
  if (!mctruth.NeutrinoSet()) {
    log->debug("Truth2h5 MCTruth for '{}' has no neutrino set", m_mctruth_label);
    return;
  }

  const auto& nu = mctruth.GetNeutrino();
  const auto& nu_particle = nu.Nu();
  const TLorentzVector& position = nu_particle.Position(0);
  const TLorentzVector& momentum = nu_particle.Momentum(0);

  m_info.valid = true;
  m_info.pdg = nu_particle.PdgCode();
  m_info.ccnc = nu.CCNC();
  m_info.intType = nu.InteractionType();
  m_info.energy = momentum.E();
  m_info.vx = position.X();
  m_info.vy = position.Y();
  m_info.vz = position.Z();
}

bool AIML::Truth2h5::operator()(const input_pointer& in, output_pointer& out)
{
  out = in;
  if (!in) { return true; }

  if (m_info.valid) {
    ensure_file();
    write_event(*in);
  }
  else {
    log->debug("Truth2h5: no neutrino info to write for frame {}", in->ident());
  }
  return true;
}

void AIML::Truth2h5::ensure_file()
{
  if (m_file >= 0) { return; }
  m_file = H5Fopen(m_output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (m_file < 0) {
    m_file = H5Fcreate(m_output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (m_file < 0) { log->error("Truth2h5 failed to open output file {}", m_output_file); }
  }
}

void AIML::Truth2h5::write_event(const IFrame& frame)
{
  if (m_file < 0) { return; }

  const hsize_t dims[1] = {1};
  hid_t dataspace = H5Screate_simple(1, dims, nullptr);
  const hid_t dtype = neutrino_type();
  const std::string dset_name = "/" + std::to_string(frame.ident()) + "/metadata";

  if (dataspace < 0) {
    log->warn("Truth2h5 failed to create dataspace for {}", dset_name);
    return;
  }

  hid_t lcpl = H5Pcreate(H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl, 1);
  hid_t dset =
    H5Dcreate2(m_file, dset_name.c_str(), dtype, dataspace, lcpl, H5P_DEFAULT, H5P_DEFAULT);
  H5Pclose(lcpl);
  if (dset < 0) {
    log->warn("Truth2h5 failed to create dataset {}", dset_name);
    H5Sclose(dataspace);
    return;
  }

  herr_t status = H5Dwrite(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_info);
  if (status < 0) { log->warn("Truth2h5 failed to write dataset {}", dset_name); }

  H5Dclose(dset);
  H5Sclose(dataspace);
}
