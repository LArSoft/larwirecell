#include "TimingTPCBundle.h"

using namespace WireCell::QLMatch;
using namespace WireCell::Clus::Facade;

// #include <boost/math/distributions/kolmogorov_smirnov.hpp>

// Calculate KS test statistic manually
// TODO: check that this is correct
double calc_ks_test(const std::vector<double>& measured, const std::vector<double>& predicted)
{
  // Both distributions should be normalized and same length
  size_t n = measured.size();

  // Calculate cumulative distributions
  std::vector<double> cum_measured(n);
  std::vector<double> cum_predicted(n);

  // First value
  cum_measured[0] = measured[0];
  cum_predicted[0] = predicted[0];

  // Calculate cumulative sums
  for (size_t i = 1; i < n; i++) {
    cum_measured[i] = cum_measured[i - 1] + measured[i];
    cum_predicted[i] = cum_predicted[i - 1] + predicted[i];
  }

  // Find maximum absolute difference
  double max_diff = 0.0;
  for (size_t i = 0; i < n; i++) {
    double diff = std::abs(cum_measured[i] - cum_predicted[i]);
    max_diff = std::max(max_diff, diff);
  }

  return max_diff;
}

TimingTPCBundle::TimingTPCBundle(Opflash* flash,
                                 Cluster* main_cluster,
                                 int flash_index_id,
                                 int cluster_index_id)
  : flash(flash)
  , main_cluster(main_cluster)
  , orig_main_cluster(0)
  , flash_index_id(flash_index_id)
  , cluster_index_id(cluster_index_id)
  , flag_close_to_PMT(false)
  , flag_at_x_boundary(false)
  , flag_spec_end(false)
  , flag_potential_bad_match(false)
  , flag_high_consistent(false)
  , ks_dis(1)
  , chi2(0)
  , ndf(0)
  , strength(0)
{
  m_nchan = flash->get_num_channels();
  pred_flash.resize(m_nchan, 0);
  opdet_mask.resize(m_nchan, 0);
}

TimingTPCBundle::~TimingTPCBundle() {}

double TimingTPCBundle::get_total_pred_light()
{
  double sum = 0;
  for (size_t i = 0; i != pred_flash.size(); i++) {
    sum += pred_flash.at(i);
  }
  return sum;
}
bool TimingTPCBundle::examine_bundle(TimingTPCBundle* candidate_bundle)
{
  // Store data in vectors instead of histograms
  std::vector<double> measured_dist(m_nchan);
  std::vector<double> predicted_dist(m_nchan);

  std::vector<double> pe(m_nchan);
  std::vector<double> pe_err(m_nchan);
  std::vector<double> pred_pe(m_nchan);

  double candidate_ks_dis=0;
  double candidate_chi2=0; 

  // Fill the initial data arrays
  for (int i = 0; i != m_nchan; i++) {
    pe[i] = flash->get_PE(i);
    pe_err[i] = flash->get_PE_err(i);
    pred_pe[i] = pred_flash.at(i) + candidate_bundle->get_pred_flash().at(i);
  }

  // Process data and fill distributions
  double total_predicted = 0;
  double total_measured = 0;
  for (int j = 0; j != m_nchan; j++) {
    measured_dist[j] = pe[j];
    predicted_dist[j] = pred_pe[j];
    total_predicted += pred_pe[j];
    total_measured += pe[j];
  }

  // Normalize distributions for KS test
  if (total_predicted > 0) {
    for (int j = 0; j != m_nchan; j++) {
      predicted_dist[j] /= total_predicted;
    }
  }
  if (total_measured > 0) {
    for (int j = 0; j != m_nchan; j++) {
      measured_dist[j] /= total_measured;
    }
  }

  if (total_predicted > 0) {
    candidate_ks_dis = calc_ks_test(measured_dist, predicted_dist);
  }

  ndf = 0;
  int nvalidopdets = 0;
  // Calculate chi-square statistics
  for (int j = 0; j != m_nchan; j++) {
    if (opdet_mask[j] == 0) continue;
    else nvalidopdets++;
    double cur_chi2 = 0;

    // TODO: add config for noise threshold of PE
    if (pe[j] < 1  && pred_pe[j] < 1) {}
    else ndf++;
    // * can add different chisq calculation (or different denominator) for cluster flags
    cur_chi2 = pow(pred_pe[j] - pe[j], 2) / (pe[j] + pow(pe_err[j], 2));
    candidate_chi2 += cur_chi2;
  }

  // if at least one of the metrics has improved
  // and if the overal metrics pass some threshold
  std::cout << "original ks_dis " << ks_dis << " chi2/ndf " << chi2/ndf << std::endl;
  std::cout << "candidate ks_dis " << candidate_ks_dis << " chi2/ndf " << candidate_chi2/ndf << std::endl;
  if ((candidate_ks_dis < ks_dis || candidate_chi2 < chi2) && ((candidate_ks_dis < 0.2) && (candidate_chi2/ndf < 20))){ 
    return true;
  }
  else { 
    return false;
  }
}

// void TimingTPCBundle::examine_merge_clusters(double dis_cut)
// {
//   // int main_cluster_id = main_cluster->get_cluster_id();

//   ClusterSelection merge_clusters;
//   for (size_t i = 0; i != other_clusters.size(); i++) {
//     Cluster* temp_cluster = other_clusters.at(i);

//     double dis_save = 1e9;

//     {
//       Cluster* cluster1 = temp_cluster;
//       Cluster* cluster2 = main_cluster;
//       const Blob* prev_mcell1 = 0;
//       const Blob* prev_mcell2 = 0;
//       const Blob* mcell1 = 0;
//       Point p1; //
//       const Blob* mcell2 = 0;
//       Point p2;

//       mcell1 = *(cluster1->time_blob_map().begin()->second.begin());
//       p1 = {mcell1->center_x(), mcell1->center_y(), mcell1->center_z()};

//       while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2) {
//         prev_mcell1 = mcell1;
//         prev_mcell2 = mcell2;

//         // find the closest point and merged cell in cluster2
//         std::pair<Point, const Blob*> temp_results = cluster2->get_closest_point_blob(p1);
//         p2 = temp_results.first;
//         mcell2 = temp_results.second;
//         // find the closest point and merged cell in cluster1
//         temp_results = cluster1->get_closest_point_blob(p2);
//         p1 = temp_results.first;
//         mcell1 = temp_results.second;
//       }
//       double dis =
//         sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));

//       if (dis < dis_save) { dis_save = dis; }

//       prev_mcell1 = 0;
//       prev_mcell2 = 0;

//       mcell1 = *(cluster1->time_blob_map().rbegin()->second.begin());
//       // p1 = mcell1->center();
//       p1 = {mcell1->center_x(), mcell1->center_y(), mcell1->center_z()};

//       while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2) {
//         prev_mcell1 = mcell1;
//         prev_mcell2 = mcell2;

//         // find the closest point and merged cell in cluster2
//         std::pair<Point, const Blob*> temp_results = cluster2->get_closest_point_blob(p1);
//         p2 = temp_results.first;
//         mcell2 = temp_results.second;
//         // find the closest point and merged cell in cluster1
//         temp_results = cluster1->get_closest_point_blob(p2);
//         p1 = temp_results.first;
//         mcell1 = temp_results.second;
//       }
//       dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));

//       if (dis < dis_save) { dis_save = dis; }
//     }

//     if (dis_save < dis_cut) { 
//       std::cout << "Merging clusters with distance " << dis_save << std::endl;
//       merge_clusters.push_back(temp_cluster); }
//   }
//   if (merge_clusters.size() > 0) {
//     // merge_clusters.push_back(main_cluster);
//     // Cluster* ncluster = new Cluster(main_cluster_id);
//     // for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
//     //   Cluster* ocluster = *(it1);
//     //   SMGCSelection& mcells = ocluster->get_mcells();
//     //   for (auto it2 = mcells.begin(); it2 != mcells.end(); it2++) {
//     //     Blob* mcell = (*it2);
//     //     // std::cout << ocluster->get_cluster_id() << " " << mcell << std::endl;
//     //     int time_slice = mcell->GetTimeSlice();
//     //     ncluster->AddCell(mcell, time_slice);
//     //   }
//     // }
//     auto grouping = main_cluster->grouping();
//     for (auto tobemerger : merge_clusters) {
//       main_cluster->take_children(*tobemerger, true);
//       grouping->destroy_child(tobemerger);
//     }
//     // delete old clusters
//     for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
//       Cluster* ocluster = *(it1);
//       if (ocluster == main_cluster) { delete ocluster; }
//       else {
//         auto it2 = find(other_clusters.begin(), other_clusters.end(), ocluster);
//         if (it2 != other_clusters.end()) { other_clusters.erase(it2); }
//         auto it3 = find(more_clusters.begin(), more_clusters.end(), ocluster);
//         if (it3 != more_clusters.end()) { more_clusters.erase(it3); }
//         // delete ocluster; //! this line causes seg fault
//       }
//     }
//     // main_cluster = ncluster;
//   }
//   std::cout << "examine_merge_clusters end " << std::endl;
// }

void TimingTPCBundle::add_bundle(TimingTPCBundle* candidate_bundle)
{
  // if the candidate bundle is worse than the current one
  if (ks_dis * pow(chi2 / ndf, 0.8) <
      candidate_bundle->get_ks_dis() * pow(candidate_bundle->get_chi2() / candidate_bundle->get_ndf(), 0.8)) {
    other_clusters.push_back(candidate_bundle->get_main_cluster());
    more_clusters.push_back(candidate_bundle->get_main_cluster());

    std::copy(candidate_bundle->get_other_clusters().begin(),
              candidate_bundle->get_other_clusters().end(),
              std::back_inserter(other_clusters));
    std::copy(candidate_bundle->get_more_clusters().begin(),
              candidate_bundle->get_more_clusters().end(),
              std::back_inserter(more_clusters));
  }
  // if the candidate bundle is actually better
  else {
    other_clusters.push_back(main_cluster);
    more_clusters.push_back(main_cluster);

    main_cluster = candidate_bundle->get_main_cluster();
    std::copy(candidate_bundle->get_other_clusters().begin(),
              candidate_bundle->get_other_clusters().end(),
              std::back_inserter(other_clusters));
    std::copy(candidate_bundle->get_more_clusters().begin(),
              candidate_bundle->get_more_clusters().end(),
              std::back_inserter(more_clusters));

    flag_close_to_PMT = candidate_bundle->get_flag_close_to_PMT();
    flag_at_x_boundary = candidate_bundle->get_flag_at_x_boundary();
  }

  std::vector<double>& pes = candidate_bundle->get_pred_flash();
  for (size_t i = 0; i != pred_flash.size(); i++) {
    pred_flash.at(i) += pes.at(i);
  }
  examine_bundle();
}

bool TimingTPCBundle::examine_bundle()
{
  // Store data in vectors instead of histograms
  std::vector<double> measured_dist(m_nchan);
  std::vector<double> predicted_dist(m_nchan);

  std::vector<double> pe(m_nchan);
  std::vector<double> pe_err(m_nchan);
  std::vector<double> pred_pe(m_nchan);

  // Fill initial arrays
  for (int i = 0; i != m_nchan; i++) {
    pe[i] = flash->get_PE(i);
    pe_err[i] = flash->get_PE_err(i);
    pred_pe[i] = pred_flash.at(i);
  }

  // Process data and fill distributions
  double total_predicted = 0;
  double total_measured = 0;
  for (int j = 0; j != m_nchan; j++) {
    measured_dist[j] = pe[j];
    predicted_dist[j] = pred_pe[j];
    total_predicted += pred_pe[j];
    total_measured += pe[j];
  }

  // Normalize distributions for KS test
  if (total_predicted > 0) {
    for (int j = 0; j != m_nchan; j++) {
      predicted_dist[j] /= total_predicted;
    }
  }
  if (total_measured > 0) {
    for (int j = 0; j != m_nchan; j++) {
      measured_dist[j] /= total_measured;
    }
  }

  if (total_predicted > 0) { ks_dis = calc_ks_test(measured_dist, predicted_dist); }

  chi2 = 0;
  ndf = 0;
  // double max_chi2 = 0;
  // int max_bin = -1;
  int nvalidopdets = 0;
  // Calculate chi-square statistics
  for (int j = 0; j != m_nchan; j++) {
    if (opdet_mask[j] == 0)
      continue;
    else
      nvalidopdets++;
    double cur_chi2 = 0;

    // TODO: add config for noise threshold of PE
    if (pe[j] < 1 && pred_pe[j] < 1) {}
    else
      ndf++;
    // * can add different chisq calculation (or different denominator) for cluster flags
    cur_chi2 = pow(pred_pe[j] - pe[j], 2) / (pe[j] + pow(pe_err[j], 2));
    chi2 += cur_chi2;

    // std::cout << "opdet" << j << " " << pe[j] << " " << pred_pe[j] << " " << pe_err[j] << " " << cur_chi2 << std::endl;

    // * can do some special treatment when the worst chi2 is from a channel with no measured light
    // if (cur_chi2 > max_chi2) {
    //   max_chi2 = cur_chi2;
    // max_bin = j;
    // }
  }
  // if (pe[max_bin] == 0 && pred_pe[max_bin] > 0) { chi2 -= max_chi2 - 1; }

  // Check multiple consistency conditions
  flag_high_consistent = false;
  if (ks_dis < 0.06 && ndf >= 3 && chi2 < ndf * nvalidopdets) { flag_high_consistent = true; }

  if (flag_high_consistent) { return true; }
  else {
    return false;
  }
}
