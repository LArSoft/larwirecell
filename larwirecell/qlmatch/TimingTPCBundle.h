#ifndef WIRECELL_QLMATCH_TIMINGTPCBUNDLE
#define WIRECELL_QLMATCH_TIMINGTPCBUNDLE

#include "Opflash.h"
#include "WireCellClus/Facade.h"
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace WireCell::QLMatch {
  typedef std::vector<WireCell::PointCloud::Facade::Cluster*> ClusterSelection;
  class TimingTPCBundle {
  public:
    typedef std::shared_ptr<TimingTPCBundle> pointer;
    using Cluster = WireCell::PointCloud::Facade::Cluster;

    TimingTPCBundle(Opflash* flash,
                    Cluster* main_cluster,
                    int flash_index_id,
                    int cluster_index_id);
    ~TimingTPCBundle();

    void set_flag_close_to_PMT(bool value) { flag_close_to_PMT = value; };
    void set_flag_at_x_boundary(bool value) { flag_at_x_boundary = value; };

    bool get_flag_close_to_PMT() { return flag_close_to_PMT; };
    bool get_flag_at_x_boundary() { return flag_at_x_boundary; };

    std::vector<double>& get_pred_flash() { return pred_flash; };
    void set_pred_flash(const std::vector<double>& values) { pred_flash = values; };
    Opflash* get_flash() { return flash; };
    void set_flash(Opflash* flash1) { flash = flash1; };
    double get_total_pred_light();

    Cluster* get_main_cluster() { return main_cluster; };
    void set_main_cluster(Cluster* cluster) { main_cluster = cluster; };

    Cluster* get_orig_cluster() { return orig_main_cluster; };
    void set_orig_cluster(Cluster* cluster) { orig_main_cluster = cluster; };

    ClusterSelection& get_other_clusters() { return other_clusters; };
    ClusterSelection& get_more_clusters() { return more_clusters; };
    void clear_other_clusters() { other_clusters.clear(); };
    void clear_more_clusters() { more_clusters.clear(); };
    void add_other_cluster(Cluster* cluster) { other_clusters.push_back(cluster); };

    std::vector<uint>& get_opdet_mask() { return opdet_mask; };
    void set_opdet_mask(const std::vector<uint>& mask) { opdet_mask = mask; };

    bool examine_bundle();
    bool examine_bundle(TimingTPCBundle* new_bundle);
    void add_bundle(TimingTPCBundle* new_bundle);

    double get_chi2() { return chi2; };
    void set_chi2(double value) { chi2 = value; };
    int get_ndf() { return ndf; };
    void set_ndf(int value) { ndf = value; };
    double get_ks_dis() { return ks_dis; };
    void set_ks_dis(double value) { ks_dis = value; };
    void set_consistent_flag(bool value) { flag_high_consistent = value; };
    bool get_consistent_flag() { return flag_high_consistent; };

    void set_spec_end_flag(bool value) { flag_spec_end = value; };
    bool get_spec_end_flag() { return flag_spec_end; };
    bool get_potential_bad_match_flag() { return flag_potential_bad_match; };
    void set_potential_bad_match_flag(bool value) { flag_potential_bad_match = value; };

    double get_strength() { return strength; };
    void set_strength(double value) { strength = value; };

    void examine_merge_clusters(double dis_cut = 3.6 * units::cm);

    int get_flash_index_id() { return flash_index_id; };
    int get_cluster_index_id() { return cluster_index_id; };

  private:
    Opflash* flash;
    Cluster* main_cluster;
    Cluster* orig_main_cluster;
    int flash_index_id;
    int cluster_index_id;
    std::vector<uint> opdet_mask;

    bool flag_close_to_PMT;
    bool flag_at_x_boundary;
    bool flag_spec_end;
    bool flag_potential_bad_match;
    bool flag_high_consistent;

    double ks_dis;
    double chi2;
    int ndf;

    double strength;

    int m_nchan;
    std::vector<double> pred_flash; // prediction
    std::vector<Cluster*> other_clusters; // save every other one
    std::vector<Cluster*> more_clusters;  // save ones satisfying the cut
  };

  /// TODO: implement comparison operators
  typedef std::vector<TimingTPCBundle::pointer> TimingTPCBundleSelection;
  typedef std::set<TimingTPCBundle::pointer> TimingTPCBundleSet;
  typedef std::map<Opflash*, TimingTPCBundleSelection> FlashBundlesMap;
  typedef std::map<WireCell::PointCloud::Facade::Cluster*, TimingTPCBundleSelection> ClusterBundlesMap;

} // namespace WCP

#endif
