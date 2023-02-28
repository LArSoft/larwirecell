/** An art event visitor component takes a reference to an art::Event
 * object.  Components with this interface will typically implement a
 * sink or a source WCT interface.
 *
 * Note, this is a Wire Cell Toolkit Interface class which depends on
 * external types so is not kept in wire-cell-iface.  See that package
 * for in-toolkit WCT interfaces.
 *
 * https://github.com/WireCell/wire-cell-iface/tree/master/inc/WireCellIface
 */

#ifndef LARWIRECELL_INTERFACES_IARTEVENTVISITOR
#define LARWIRECELL_INTERFACES_IARTEVENTVISITOR

#include "WireCellUtil/IComponent.h"

namespace art {
  class Event;
  class EDProducer;
  class ProducesCollector;
}
namespace wcls {
  class IArtEventVisitor : public WireCell::IComponent<IArtEventVisitor> {
  public:
    virtual ~IArtEventVisitor() {}

    /// If data is produced, must implement in order to call:
    ///   collector.produces<DataType>();
    /// If only reading data, implementation is not required.
    virtual void produces(art::ProducesCollector& collector) {}

    /// Implement to visit an Art event.
    virtual void visit(art::Event& event) = 0;
  };
}
#endif
