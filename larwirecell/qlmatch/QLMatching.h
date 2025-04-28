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
    // QL Matching configs 
    bool   m_pmts{true};
    bool   m_data{true};
    std::vector<int> m_ch_mask{39, 69, 71, 85, 87, 115, 141, 197, 217, 221, 222, 223, 226, 245, 249, 66, 86, 138, 302};
    bool   m_beamonly{false};
    double  m_flash_minPE{500};       // PE 
    double  m_flash_mintime{-1500e3}; // ns
    double  m_flash_maxtime{ 1500e3};  // ns
    double  m_beam_mintime{-5e3}; // ns
    double  m_beam_maxtime{ 5e3}; // ns
    double  m_QtoL{0.5};

    std::vector<double> m_VUVEfficiency{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.035, 0.03, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021};
    std::vector<double> m_VISEfficiency{0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0375, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.039, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    IAnodePlane::pointer m_anode{nullptr};

    // refer to MultiAlgBlobClustering for the following
    std::string m_inpath{"pointtrees/%d"};
    std::string m_outpath{"pointtrees/%d"};
    std::string m_bee_dir{"data"};

    void remove_bundle_selection(TimingTPCBundleSelection to_be_removed, TimingTPCBundleSet& bundle_set);
    void remove_bundle_selection(TimingTPCBundleSelection to_be_removed, FlashBundlesMap& flash_bundles_map, ClusterBundlesMap& cluster_bundles_map, std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer>& flash_cluster_bundles_map);
    void organize_bundles(TimingTPCBundleSelection& results_bundles, std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer>& flash_cluster_bundles_map);
  };
}

#endif