#ifndef WIRECELLDEV_QLMATCH_QLMATCHING
#define WIRECELLDEV_QLMATCH_QLMATCHING

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetFanin.h"

#include "TimingTPCBundle.h"

using namespace WireCell;
using namespace WireCell::PointCloud::Facade;

namespace WireCell::QLMatch {
  class QLMatching : public Aux::Logger, public ITensorSetFanin, public IConfigurable {
  public:
    QLMatching();
    virtual ~QLMatching();

    // INode, override because we get multiplicity at run time.
    virtual std::vector<std::string> input_types();

    // ITensorSetFanin
    // input: 0: charge, 1: light
    virtual bool operator()(const input_vector& invec, output_pointer& out);

    // IConfigurable
    virtual void configure(const WireCell::Configuration& config);
    virtual WireCell::Configuration default_configuration() const;

  private:
    // Count how many times we are called
    size_t m_count{0};
    // Currently can only be 2, TODO: remove this?
    size_t m_multiplicity{2};

    IAnodePlane::pointer m_anode{nullptr};

    // refer to MultiAlgBlobClustering for the following
    std::string m_inpath{"pointtrees/%d"};
    std::string m_outpath{"pointtrees/%d"};
    std::string m_bee_dir{"data"};

    void remove_bundle_selection(TimingTPCBundleSelection to_be_removed, TimingTPCBundleSet& bundle_set);
    void remove_bundle_selection(TimingTPCBundleSelection to_be_removed, FlashBundlesMap& flash_bundles_map, ClusterBundlesMap& cluster_bundles_map, std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer>& flash_cluster_bundles_map);
  };
}

#endif