/** The WC/LS sim depo set source adapts larsoft SimEnergyDeposit to WCT's IDepoSets.
    Modified from the SimDepoSource for IDepoSet for better efficiency
*/

#include "SimDepoSetSource.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Units.h"

WIRECELL_FACTORY(wclsSimDepoSetSource,
                 wcls::SimDepoSetSource,
                 wcls::IArtEventVisitor,
                 WireCell::IDepoSetSource,
                 WireCell::IConfigurable)

using WireCell::Aux::SimpleDepo;
using WireCell::Aux::SimpleDepoSet;

namespace units = WireCell::units;

namespace wcls {

  // There is more than one way to make ionization electrons.
  // These adapters erase these differences.
  class SimDepoSetSource::DepoAdapter {
  public:
    virtual ~DepoAdapter() {}
    virtual double operator()(const sim::SimEnergyDeposit& sed) const = 0;
  };

  namespace bits {

    // This takes number of electrons directly, and applies a
    // multiplicative scale.
    class ElectronsAdapter : public SimDepoSetSource::DepoAdapter {
      double m_scale;

    public:
      ElectronsAdapter(double scale = 1.0) : m_scale(scale) {}
      virtual ~ElectronsAdapter(){};
      virtual double operator()(const sim::SimEnergyDeposit& sed) const
      {
        return m_scale * sed.NumElectrons();
      }
    };

    // This one takes a recombination model which only requires dE
    // (ie, assumes MIP).
    class PointAdapter : public SimDepoSetSource::DepoAdapter {
      WireCell::IRecombinationModel::pointer m_model;
      double m_scale;

    public:
      PointAdapter(WireCell::IRecombinationModel::pointer model, double scale = 1.0)
        : m_model(model), m_scale(scale)
      {}
      virtual ~PointAdapter() {}
      virtual double operator()(const sim::SimEnergyDeposit& sed) const
      {
        const double dE = sed.Energy() * units::MeV;
        return m_scale * (*m_model)(dE);
      }
    };

    // This one takes a recombination which is a function of both dE and dX.
    class StepAdapter : public SimDepoSetSource::DepoAdapter {
      WireCell::IRecombinationModel::pointer m_model;
      double m_scale;

    public:
      StepAdapter(WireCell::IRecombinationModel::pointer model, double scale = 1.0)
        : m_model(model), m_scale(scale)
      {}
      virtual ~StepAdapter() {}
      virtual double operator()(const sim::SimEnergyDeposit& sed) const
      {
        const double dE = sed.Energy() * units::MeV;
        const double dX = sed.StepLength() * units::cm;
        return m_scale * (*m_model)(dE, dX);
      }
    };
  }
}

using namespace wcls;

SimDepoSetSource::SimDepoSetSource() : m_count(0), m_adapter(nullptr) {}

SimDepoSetSource::~SimDepoSetSource()
{
  if (m_adapter) {
    delete m_adapter;
    m_adapter = nullptr;
  }
}

WireCell::Configuration SimDepoSetSource::default_configuration() const
{
  WireCell::Configuration cfg;

  // An empty model or the special model "electrons" means to take
  // precalcualted input NumElectrons.  Anything else means some
  // recombination model is used.
  cfg["model"] = "";
  // Multiply this number to the number of electrons before forming
  // a WC depo.
  cfg["scale"] = 1.0;

  // For locating input in the art::Event
  cfg["art_tag"] = "";      // eg, "plopper:bogus"
  cfg["assn_art_tag"] = ""; // eg, "largeant"

  return cfg;
}
void SimDepoSetSource::configure(const WireCell::Configuration& cfg)
{
  if (m_adapter) {
    delete m_adapter;
    m_adapter = nullptr;
  }
  double scale = WireCell::get(cfg, "scale", 1.0);

  std::string model_tn = WireCell::get<std::string>(cfg, "model", "");
  std::string model_type = "";
  if (!model_tn.empty()) {
    // FIXME: Right now, the split(..)[0] is just a no-op.  Is the value supposed to be
    // assigned to something?
    (void)WireCell::String::split(model_tn)[0];
  }

  if (model_type == "" or model_type == "electrons") {
    m_adapter = new wcls::bits::ElectronsAdapter(scale);
  }
  else {
    auto model = WireCell::Factory::lookup_tn<WireCell::IRecombinationModel>(model_tn);
    if (!model) {
      std::cerr << "wcls::SimDepoSetSource: unknown recombination model: \"" << model_tn << "\"\n";
      return; // I should throw something here!
    }
    if (model_type == "MipRecombination") {
      m_adapter = new wcls::bits::PointAdapter(model, scale);
    }
    if (model_type == "BirksRecombination" || model_type == "BoxRecombination") {
      m_adapter = new wcls::bits::StepAdapter(model, scale);
    }
  }

  m_inputTag = cfg["art_tag"].asString();
  m_assnTag = cfg["assn_art_tag"].asString();
}

void SimDepoSetSource::visit(art::Event& event)
{
  art::Handle<std::vector<sim::SimEnergyDeposit>> sedvh;

  bool okay = event.getByLabel(m_inputTag, sedvh);
  if (!okay) {
    std::string msg =
      "SimDepoSetSource failed to get sim::SimEnergyDeposit from art tag: " + m_inputTag.encode();
    std::cerr << msg << std::endl;
    THROW(WireCell::RuntimeError() << WireCell::errmsg{msg});
  }
  //else if (sedvh->empty()) return;

  const size_t ndepos = sedvh->size();

  std::cerr << "SimDepoSetSource got " << ndepos << " depos from art tag \"" << m_inputTag
            << "\" returns: " << (okay ? "okay" : "fail") << std::endl;

  if (!m_depos.empty()) {
    std::cerr << "SimDepoSetSource dropping " << m_depos.size() << " unused, prior depos\n";
    m_depos.clear();
  }

  // associate the input SED with the other set of SED (eg, before SCE)
  std::vector<sim::SimEnergyDeposit> assn_sedv;
  if (m_assnTag != "") {
    art::Handle<std::vector<sim::SimEnergyDeposit>> assn_sedvh;
    okay = event.getByLabel(m_assnTag, assn_sedvh);
    if (!okay) {
      std::string msg =
        "SimDepoSetSource failed to get sim::SimEnergyDeposit from art tag: " + m_assnTag.encode();
      std::cerr << msg << std::endl;
      THROW(WireCell::RuntimeError() << WireCell::errmsg{msg});
    }
    else {
      std::cout << "Larwirecell::SimDepoSetSource got " << assn_sedvh->size()
                << " associated depos from " << m_assnTag << std::endl;
      assn_sedv.insert(assn_sedv.end(), assn_sedvh->begin(), assn_sedvh->end());
    }
    // safty check for the associated SED
    if (ndepos != assn_sedv.size()) {
      std::string msg = "Larwirecell::SimDepoSetSource Inconsistent size of SimDepoSetSources";
      std::cerr << msg << std::endl;
      THROW(WireCell::RuntimeError() << WireCell::errmsg{msg});
    }
  }

  for (size_t ind = 0; ind < ndepos; ++ind) {
    auto const& sed = sedvh->at(ind);
    auto pt = sed.MidPoint();
    const WireCell::Point wpt(pt.x() * units::cm, pt.y() * units::cm, pt.z() * units::cm);
    double wt = sed.Time() * units::ns;
    double wq = (*m_adapter)(sed);
    int wid = sed.TrackID();
    int pdg = sed.PdgCode();
    double we = sed.Energy() * units::MeV;

    if (assn_sedv.size() == 0) {
      WireCell::IDepo::pointer depo =
        std::make_shared<SimpleDepo>(wt, wpt, wq, nullptr, 0.0, 0.0, wid, pdg, we);
      m_depos.push_back(depo);
      // std::cerr << ind << ": t=" << wt/units::us << "us,"
      //           << " r=" << wpt/units::cm << "cm, "
      //           << " q=" << wq
      //           << " e=" << we/units::MeV << "\n";
    }
    else {
      auto const& sed1 = assn_sedv.at(ind);
      auto pt1 = sed1.MidPoint();
      const WireCell::Point wpt1(pt1.x() * units::cm, pt1.y() * units::cm, pt1.z() * units::cm);
      double wt1 = sed1.Time() * units::ns;
      double wq1 = (*m_adapter)(sed1);
      int wid1 = sed1.TrackID();
      int pdg1 = sed1.PdgCode();
      double we1 = sed1.Energy() * units::MeV;

      WireCell::IDepo::pointer assn_depo =
        std::make_shared<SimpleDepo>(wt1, wpt1, wq1, nullptr, 0.0, 0.0, wid1, pdg1, we1);

      WireCell::IDepo::pointer depo =
        std::make_shared<SimpleDepo>(wt, wpt, wq, assn_depo, 0.0, 0.0, wid, pdg, we);
      m_depos.push_back(depo);
      // std::cerr << ind << ": t1=" << wt1/units::us << "us,"
      //           << " r1=" << wpt1/units::cm << "cm, "
      //           << " q1=" << wq1
      //           << " e1=" << we1/units::MeV << "\n";
    }
  }

  // empty "ionization": no TPC activity
  if (ndepos == 0) {
    WireCell::Point wpt(0, 0, 0);
    WireCell::IDepo::pointer depo = std::make_shared<SimpleDepo>(0, wpt, 0, nullptr, 0.0, 0.0);
    m_depos.push_back(depo);
  }

  // don't trust user to honor time ordering.
  std::sort(m_depos.begin(), m_depos.end(), WireCell::ascending_time);
  std::cerr << "SimDepoSetSource: ready with " << m_depos.size() << " depos spanning: ["
            << m_depos.front()->time() / units::us << ", " << m_depos.back()->time() / units::us
            << "]us\n";
}

bool SimDepoSetSource::operator()(WireCell::IDepoSet::pointer& out)
{
  if (m_depos.empty()) { return false; }

  out = std::make_shared<SimpleDepoSet>(m_count, m_depos);
  m_depos.clear();
  ++m_count;

  return true;
}
