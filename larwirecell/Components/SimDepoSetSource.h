/** A sim depo set source provides IDepoSet objects to WCT from LS simulation objects.
    Modified from the SimDepoSource for IDepoSet for better efficiency
*/

#ifndef LARWIRECELL_COMPONENTS_SIMDEPOSETSOURCE
#define LARWIRECELL_COMPONENTS_SIMDEPOSETSOURCE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDepoSet.h"
#include "WireCellIface/IDepoSetSource.h"
#include "canvas/Utilities/InputTag.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

namespace wcls {

  class SimDepoSetSource : public IArtEventVisitor,
                           public WireCell::IDepoSetSource,
                           public WireCell::IConfigurable {
  public:
    class DepoAdapter;

    SimDepoSetSource();
    virtual ~SimDepoSetSource();

    /// IArtEventVisitor
    virtual void visit(art::Event& event);

    /// IDepoSetSource
    virtual bool operator()(WireCell::IDepoSet::pointer& out);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

  private:
    // Count how many we've produced, use this for the depo set ident.
    int m_count{0};

    // Temporary holding of accepted depos.
    WireCell::IDepo::vector m_depos;

    DepoAdapter* m_adapter;

    art::InputTag m_inputTag;
    art::InputTag m_assnTag; // associated input
  };
}
#endif
