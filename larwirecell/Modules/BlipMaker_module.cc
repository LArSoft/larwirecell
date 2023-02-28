#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/Simulation/SimEnergyDeposit.h"

namespace bogoblip {
  // struct BlipMakerConfig {
  // };
  class BlipMaker : public art::EDProducer {
  public:
    // using Parameters = art::EDProducer::Table<BlipMakerConfig>;

    // explicit BlipMaker(Parameters const& params);
    explicit BlipMaker(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt);

  private:
    // const BlipMakerConfig m_cfg;
    int m_count;

    std::vector<double> m_head; // start of the ray [cm]
    std::vector<double> m_tail; // end of the ray [cm]
    double m_time;              // depo's time [ns]
    double m_charge;            // charge per step [electron]
    double m_step;              // step length [cm]
  };
}

#include <string>
const std::string instance = "bogus"; // fixme: make configurable

// bogoblip::BlipMaker::BlipMaker(Parameters const& params)
//   : EDProducer{params}
//   , m_cfg(params())
//     , m_count(0)
// {
//     produces< std::vector<sim::SimEnergyDeposit> >(instance);
// }

bogoblip::BlipMaker::BlipMaker(fhicl::ParameterSet const& pset)
  : EDProducer{pset}
  , m_count(0)
  , m_head{pset.get<std::vector<double>>("Head", {100., 0., 0.})}
  , m_tail{pset.get<std::vector<double>>("Tail", {50., 10., 50.})}
  , m_time{pset.get<double>("Time", 0.0)}
  , m_charge{pset.get<double>("ElectronPerCm", 50000)}
  , m_step{pset.get<double>("StepSize", 0.1)}
{
  produces<std::vector<sim::SimEnergyDeposit>>(instance);
}

void bogoblip::BlipMaker::produce(art::Event& event)
{
  ++m_count;

  std::cout << "head: " << m_head.at(0) << " " << m_head.at(1) << " " << m_head.at(2) << std::endl;
  std::cout << "tail: " << m_tail.at(0) << " " << m_tail.at(1) << " " << m_tail.at(2) << std::endl;
  std::cout << "time: " << m_time << " e/cm: " << m_charge << " step: " << m_step << std::endl;

  auto out = std::make_unique<std::vector<sim::SimEnergyDeposit>>();

  int nphotons = 0;
  const int nelepercm = m_charge; // 50000;
  const double mevpercm = 2.0;

  double t0 = 0.;
  double t1 = 0.;
  int trackid = 0;

  // implicit units are cm, ns and MeV.
  sim::SimEnergyDeposit::Point_t start = {
    m_head.at(0), m_head.at(1), m_head.at(2)}; // {100.,0.,0.};
  sim::SimEnergyDeposit::Point_t end = {
    m_tail.at(0), m_tail.at(1), m_tail.at(2)}; // {150.,10.,50.};
  const auto vdiff = end - start;
  const auto vlen = sqrt(vdiff.Mag2());
  const auto vdir = vdiff.unit();

  const double stepsize = m_step; // 0.1; // cm
  const int nsteps = vlen / stepsize;

  // larsoft works in ns
  const double ns = 1.0;
  // const double us = 1000.0*ns;
  // const double ms = 1000.0*us;

  // MB: WCT sim should cut the first, the second should just be on
  // the edge of time acceptance.  The last two bracket the BNB beam
  // gate.
  for (double jump : {m_time * ns} /*{ -1.6*ms, -1*ms, +3125*ns, (3125+1600)*ns, }*/) {
    sim::SimEnergyDeposit::Point_t last = start;
    for (int istep = 1; istep < nsteps; ++istep) {
      const sim::SimEnergyDeposit::Point_t next(start + stepsize * istep * vdir);
      //std::cerr << last << " -> " << next << "\n";
      out->emplace_back(sim::SimEnergyDeposit(nphotons,
                                              stepsize * nelepercm,
                                              1.0, //scintillatin yield ratio
                                              stepsize * mevpercm,
                                              last,
                                              next,
                                              jump + t0,
                                              jump + t1,
                                              trackid));
      last = next;
    }
  }

  std::cerr << "BlipMaker making " << out->size() << " depos to instance: " << instance
            << std::endl;

  event.put(std::move(out), instance);
}

namespace bogoblip {
  DEFINE_ART_MODULE(BlipMaker)
}
