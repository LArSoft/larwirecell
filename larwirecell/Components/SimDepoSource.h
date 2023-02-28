/** A sim depo source provides IDepo objects to WCT from LS simulation objects.

    Multiple strategies are possible for ingesting depositions into
    WCT simulation.

    1) Treat the vector of SimEnergyDeposit in each art::Event in
    whole and independently from those in other art::Event's.  That's
    what this converter does.

    2) Allow for multiple, independent vectors of SimEnergyDeposit
    otherwise as above to each be ingested into WCT on each art::Event
    and let WCT properly mix them (eg with WCT's Gen::DepoMerger
    component).

    3) As above but relax treating each art::Event atomically and
    allow accumulation of depos spanning multiple visit()'s.  This can
    be used to allow art modules to feed WCT to set up a streaming
    simulation such as one might want to do to simulate for DUNE
    triggering studies.  The back end of the WCT simulation will still
    need to spit out discrete frames but it can do so such that their
    boundaries can be stitched together seemlessly.

*/

#ifndef LARWIRECELL_COMPONENTS_SIMDEPOSOURCE
#define LARWIRECELL_COMPONENTS_SIMDEPOSOURCE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDepo.h"
#include "WireCellIface/IDepoSource.h"
#include "canvas/Utilities/InputTag.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include <deque>

namespace wcls {

  namespace bits {
    class DepoAdapter;
  }

  class SimDepoSource : public IArtEventVisitor,
                        public WireCell::IDepoSource,
                        public WireCell::IConfigurable {
  public:
    SimDepoSource();
    virtual ~SimDepoSource();

    /// IArtEventVisitor
    virtual void visit(art::Event& event);

    /// IDepoSource
    virtual bool operator()(WireCell::IDepo::pointer& out);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

  private:
    std::deque<WireCell::IDepo::pointer> m_depos;
    bits::DepoAdapter* m_adapter;

    art::InputTag m_inputTag;
    art::InputTag m_assnTag; // associated input
  };
}
#endif
