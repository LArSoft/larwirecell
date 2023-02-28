#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/SharedProducer.h"
#include "art/Utilities/make_tool.h"
#include "larwirecell/Interfaces/MainTool.h"

using namespace std;

namespace wcls {

  class WireCellToolkit : public art::SharedProducer {
  public:
    explicit WireCellToolkit(fhicl::ParameterSet const& pset, art::ProcessingFrame const&);
    virtual ~WireCellToolkit();

    void produce(art::Event& evt, art::ProcessingFrame const&);
    void reconfigure(fhicl::ParameterSet const& pset);

  private:
    std::unique_ptr<wcls::MainTool> m_wcls;
  };
}

wcls::WireCellToolkit::WireCellToolkit(fhicl::ParameterSet const& pset, art::ProcessingFrame const&)
  : SharedProducer(pset)
{
  const std::string s{"WCT"};
  serializeExternal(s);
  this->reconfigure(pset);
}
wcls::WireCellToolkit::~WireCellToolkit() {}

void wcls::WireCellToolkit::produce(art::Event& evt, art::ProcessingFrame const&)
{
  m_wcls->process(evt);
}

void wcls::WireCellToolkit::reconfigure(fhicl::ParameterSet const& pset)
{
  auto const& wclsPS = pset.get<fhicl::ParameterSet>("wcls_main");
  m_wcls = art::make_tool<wcls::MainTool>(wclsPS);
  if (!m_wcls) {
    throw cet::exception("WireCellToolkit_module") << "Failed to get Art Tool \"wcls_main\"\n";
  }
  m_wcls->produces(producesCollector());
}

namespace wcls {
  DEFINE_ART_MODULE(WireCellToolkit)
}
