#include "TimingTPCBundle.h"

using namespace WireCell::QLMatch;
using namespace WireCell::PointCloud::Facade;

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
bool TimingTPCBundle::examine_bundle(TimingTPCBundle* bundle,
                                     const std::vector<double>& cos_pe_low,
                                     const std::vector<double>& cos_pe_mid)
{
  // Store data in vectors instead of histograms
  std::vector<double> measured_dist(m_nchan);
  std::vector<double> predicted_dist(m_nchan);

  std::vector<double> pe(m_nchan);
  std::vector<double> pe_err(m_nchan);
  std::vector<double> pred_pe(m_nchan);

  // Fill the initial data arrays
  for (int i = 0; i != m_nchan; i++) {
    pe[i] = flash->get_PE(i);
    pe_err[i] = flash->get_PE_err(i);
    pred_pe[i] = pred_flash.at(i) + bundle->get_pred_flash().at(i);
  }

  // Process the data and fill distributions
  double total_predicted = 0;
  double total_measured = 0;
  for (int j = 0; j != m_nchan; j++) {
    measured_dist[j] = pe[j];
    if ((pred_pe[j] < cos_pe_low[j] || (pred_pe[j] < cos_pe_mid[j] * 1.1 && pe[j] == 0)) &&
        flash->get_type() == 1) {
      pred_pe[j] = 0;
    }
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

  double temp_chi2 = 0;
  double temp_ndf = 0;
  double temp_ks_dis = 1;

  double max_chi2 = 0;
  int max_bin = -1;

  // Calculate KS test using Boost
  if (total_predicted > 0) {
    // boost::math::kolmogorov_smirnov_test(measured_dist, predicted_dist, &temp_ks_dis);
    temp_ks_dis = calc_ks_test(measured_dist, predicted_dist);
  }

  // Calculate chi-square statistics
  for (int j = 0; j != m_nchan; j++) {
    double cur_chi2 = 0;
    if (flag_close_to_PMT) {
      if (pe[j] - pred_pe[j] > 350 && pe[j] > pred_pe[j] * 1.3) {
        cur_chi2 = pow(pred_pe[j] - pe[j], 2) / (pow(pe_err[j], 2) + pow(pe[j] * 0.5, 2));
      }
      else {
        cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
      }
    }
    else {
      cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
    }
    temp_chi2 += cur_chi2;

    if (cur_chi2 > max_chi2) {
      max_chi2 = cur_chi2;
      max_bin = j;
    }

    if (pe[j] == 0 && pred_pe[j] == 0) {}
    else {
      temp_ndf++;
    }
  }

  if (pe[max_bin] == 0 && pred_pe[max_bin] > 0) { temp_chi2 -= max_chi2 - 1; }

  if ((temp_ks_dis < ks_dis + 0.06 && (temp_ks_dis < ks_dis * 1.2) && temp_chi2 < chi2 + ndf * 5 &&
       temp_chi2 < chi2 * 1.21) ||
      (temp_ks_dis < ks_dis && temp_chi2 < chi2 + ndf * 10 && temp_chi2 < chi2 * 1.45)) {
    return true;
  }
  else {
    return false;
  }
}
bool TimingTPCBundle::examine_bundle_rank(TimingTPCBundle* bundle,
                                          const std::vector<double>& cos_pe_low,
                                          const std::vector<double>& cos_pe_mid)
{
  // Store data in vectors instead of histograms
  std::vector<double> measured_dist(m_nchan);
  std::vector<double> predicted_dist(m_nchan);

  std::vector<double> pe(m_nchan);
  std::vector<double> pe_err(m_nchan);
  std::vector<double> pred_pe(m_nchan);

  // Fill the initial data arrays
  for (int i = 0; i != m_nchan; i++) {
    pe[i] = flash->get_PE(i);
    pe_err[i] = flash->get_PE_err(i);
    pred_pe[i] = pred_flash.at(i) + bundle->get_pred_flash().at(i);
  }

  // Process the data and fill distributions
  double total_predicted = 0;
  double total_measured = 0;
  for (int j = 0; j != m_nchan; j++) {
    measured_dist[j] = pe[j];
    if ((pred_pe[j] < cos_pe_low[j] || (pred_pe[j] < cos_pe_mid[j] * 1.1 && pe[j] == 0)) &&
        flash->get_type() == 1) {
      pred_pe[j] = 0;
    }
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

  double temp_chi2 = 0;
  double temp_ndf = 0;
  double temp_ks_dis = 1;

  double max_chi2 = 0;
  int max_bin = -1;

  // Calculate KS test using Boost
  if (total_predicted > 0) {
    // boost::math::kolmogorov_smirnov_test(measured_dist, predicted_dist, &temp_ks_dis);
    temp_ks_dis = calc_ks_test(measured_dist, predicted_dist);
  }

  // Calculate chi-square statistics
  for (int j = 0; j != m_nchan; j++) {
    double cur_chi2 = 0;
    if (flag_close_to_PMT) {
      if (pe[j] - pred_pe[j] > 350 &&
          pe[j] > pred_pe[j] * 1.3) { // if the measurement is much larger than the prediction
        cur_chi2 = pow(pred_pe[j] - pe[j], 2) / (pow(pe_err[j], 2) + pow(pe[j] * 0.5, 2));
      }
      else {
        cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
      }
    }
    else {
      cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
    }
    temp_chi2 += cur_chi2;

    if (cur_chi2 > max_chi2) {
      max_chi2 = cur_chi2;
      max_bin = j;
    }

    if (pe[j] == 0 && pred_pe[j] == 0) {}
    else {
      temp_ndf++;
    }
  }

  if (pe[max_bin] == 0 &&
      pred_pe[max_bin] > 0) // allow one PMT to be inefficient in measurement ...
    temp_chi2 -= max_chi2 - 1;

  // Final comparison with extended conditions
  if ((temp_ks_dis < ks_dis + 0.06 &&
       (temp_ks_dis < ks_dis * 1.2 || temp_ks_dis < 0.05 || temp_ks_dis < ks_dis + 0.03) &&
       temp_chi2 < chi2 + ndf * 5 && temp_chi2 < chi2 * 1.21) ||
      (temp_ks_dis < ks_dis && temp_chi2 < chi2 + ndf * 10 && temp_chi2 < chi2 * 1.45) ||
      (temp_ks_dis * temp_chi2 < ks_dis * chi2)) {
    return true;
  }
  else {
    return false;
  }
}

void TimingTPCBundle::examine_merge_clusters(double dis_cut)
{
  // int main_cluster_id = main_cluster->get_cluster_id();

  ClusterSelection merge_clusters;
  for (size_t i = 0; i != other_clusters.size(); i++) {
    Cluster* temp_cluster = other_clusters.at(i);

    double dis_save = 1e9;

    {
      Cluster* cluster1 = temp_cluster;
      Cluster* cluster2 = main_cluster;
      const Blob* prev_mcell1 = 0;
      const Blob* prev_mcell2 = 0;
      const Blob* mcell1 = 0;
      Point p1; //
      const Blob* mcell2 = 0;
      Point p2;

      mcell1 = *(cluster1->time_blob_map().begin()->second.begin());
      p1 = {mcell1->center_x(), mcell1->center_y(), mcell1->center_z()};

      while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2) {
        prev_mcell1 = mcell1;
        prev_mcell2 = mcell2;

        // find the closest point and merged cell in cluster2
        std::pair<Point, const Blob*> temp_results = cluster2->get_closest_point_blob(p1);
        p2 = temp_results.first;
        mcell2 = temp_results.second;
        // find the closest point and merged cell in cluster1
        temp_results = cluster1->get_closest_point_blob(p2);
        p1 = temp_results.first;
        mcell1 = temp_results.second;
      }
      double dis =
        sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));

      if (dis < dis_save) { dis_save = dis; }

      prev_mcell1 = 0;
      prev_mcell2 = 0;

      mcell1 = *(cluster1->time_blob_map().rbegin()->second.begin());
      // p1 = mcell1->center();
      p1 = {mcell1->center_x(), mcell1->center_y(), mcell1->center_z()};

      while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2) {
        prev_mcell1 = mcell1;
        prev_mcell2 = mcell2;

        // find the closest point and merged cell in cluster2
        std::pair<Point, const Blob*> temp_results = cluster2->get_closest_point_blob(p1);
        p2 = temp_results.first;
        mcell2 = temp_results.second;
        // find the closest point and merged cell in cluster1
        temp_results = cluster1->get_closest_point_blob(p2);
        p1 = temp_results.first;
        mcell1 = temp_results.second;
      }
      dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));

      if (dis < dis_save) { dis_save = dis; }
    }

    if (dis_save < dis_cut) { merge_clusters.push_back(temp_cluster); }
  }

  if (merge_clusters.size() > 0) {
    // merge_clusters.push_back(main_cluster);
    // Cluster* ncluster = new Cluster(main_cluster_id);
    // for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
    //   Cluster* ocluster = *(it1);
    //   SMGCSelection& mcells = ocluster->get_mcells();
    //   for (auto it2 = mcells.begin(); it2 != mcells.end(); it2++) {
    //     Blob* mcell = (*it2);
    //     // std::cout << ocluster->get_cluster_id() << " " << mcell << std::endl;
    //     int time_slice = mcell->GetTimeSlice();
    //     ncluster->AddCell(mcell, time_slice);
    //   }
    // }
    auto grouping = main_cluster->grouping();
    for (auto tobemerger : merge_clusters) {
      main_cluster->take_children(*tobemerger, true);
      grouping->destroy_child(tobemerger);
    }
    // delete old clusters
    for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
      Cluster* ocluster = *(it1);
      if (ocluster == main_cluster) { delete ocluster; }
      else {
        auto it2 = find(other_clusters.begin(), other_clusters.end(), ocluster);
        if (it2 != other_clusters.end()) { other_clusters.erase(it2); }
        auto it3 = find(more_clusters.begin(), more_clusters.end(), ocluster);
        if (it3 != more_clusters.end()) { more_clusters.erase(it3); }
        delete ocluster;
      }
    }
    // main_cluster = ncluster;
  }
}

void TimingTPCBundle::add_bundle(TimingTPCBundle* bundle,
                                 const std::vector<double>& cos_pe_low,
                                 const std::vector<double>& cos_pe_mid)
{

  if (ks_dis * pow(chi2 / ndf, 0.8) / get_total_pred_light() <
      bundle->get_ks_dis() * pow(bundle->get_chi2() / bundle->get_ndf(), 0.8) /
        bundle->get_total_pred_light()) {
    other_clusters.push_back(bundle->get_main_cluster());
    more_clusters.push_back(bundle->get_main_cluster());

    std::copy(bundle->get_other_clusters().begin(),
              bundle->get_other_clusters().end(),
              std::back_inserter(other_clusters));
    std::copy(bundle->get_more_clusters().begin(),
              bundle->get_more_clusters().end(),
              std::back_inserter(more_clusters));
  }
  else {
    other_clusters.push_back(main_cluster);
    more_clusters.push_back(main_cluster);

    main_cluster = bundle->get_main_cluster();
    std::copy(bundle->get_other_clusters().begin(),
              bundle->get_other_clusters().end(),
              std::back_inserter(other_clusters));
    std::copy(bundle->get_more_clusters().begin(),
              bundle->get_more_clusters().end(),
              std::back_inserter(more_clusters));

    flag_close_to_PMT = bundle->get_flag_close_to_PMT();
    flag_at_x_boundary = bundle->get_flag_at_x_boundary();
  }

  std::vector<double>& pes = bundle->get_pred_flash();
  for (size_t i = 0; i != pred_flash.size(); i++) {
    pred_flash.at(i) += pes.at(i);
  }
  examine_bundle();
}

bool TimingTPCBundle::examine_beam_bundle()
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

  // Fill distributions
  double total_predicted = 0;
  double total_measured = 0;
  for (int j = 0; j != m_nchan; j++) {
    measured_dist[j] = pe[j];
    predicted_dist[j] = pred_pe[j];
    total_predicted += pred_pe[j];
    total_measured += pe[j];
  }

  // Normalize distributions for first KS test
  std::vector<double> norm_measured = measured_dist;
  std::vector<double> norm_predicted = predicted_dist;
  if (total_predicted > 0) {
    for (int j = 0; j != m_nchan; j++) {
      norm_predicted[j] /= total_predicted;
    }
  }
  if (total_measured > 0) {
    for (int j = 0; j != m_nchan; j++) {
      norm_measured[j] /= total_measured;
    }
  }

  // First KS test
  // double temp_ks_dis = 1.0;
  // boost::math::kolmogorov_smirnov_test(norm_measured, norm_predicted, &temp_ks_dis);
  double temp_ks_dis = calc_ks_test(norm_measured, norm_predicted);

  // Calculate chi-square statistics
  double temp_chi2 = 0;
  double temp_ndf = 0;
  double max_chi2 = 0;
  int max_bin = -1;

  for (int j = 0; j != m_nchan; j++) {
    double cur_chi2 = 0;
    if (flag_close_to_PMT) {
      if (pe[j] - pred_pe[j] > 350 && pe[j] > pred_pe[j] * 1.3) {
        cur_chi2 = pow(pred_pe[j] - pe[j], 2) / (pow(pe_err[j], 2) + pow(pe[j] * 0.5, 2));
      }
      else {
        cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
      }
    }
    else {
      cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
    }
    temp_chi2 += cur_chi2;

    if (cur_chi2 > max_chi2) {
      max_chi2 = cur_chi2;
      max_bin = j;
    }

    if (pe[j] == 0 && pred_pe[j] == 0) {}
    else {
      temp_ndf++;
    }
  }

  // Create new distributions with max_bin removed for second KS test
  std::vector<double> measured_dist2;
  std::vector<double> predicted_dist2;
  total_predicted = 0;
  total_measured = 0;

  for (int j = 0; j != m_nchan; j++) {
    if (j != max_bin) {
      measured_dist2.push_back(pe[j]);
      predicted_dist2.push_back(pred_pe[j]);
      total_predicted += pred_pe[j];
      total_measured += pe[j];
    }
  }

  // Normalize distributions for second KS test
  if (total_predicted > 0) {
    for (size_t j = 0; j < measured_dist2.size(); j++) {
      predicted_dist2[j] /= total_predicted;
    }
  }
  if (total_measured > 0) {
    for (size_t j = 0; j < predicted_dist2.size(); j++) {
      measured_dist2[j] /= total_measured;
    }
  }

  // Second KS test without max_bin
  // double temp_ks_dis1 = 1.0;
  // boost::math::kolmogorov_smirnov_test(measured_dist2, predicted_dist2, &temp_ks_dis1);
  double temp_ks_dis1 = calc_ks_test(measured_dist2, predicted_dist2);

  if ((temp_ks_dis < 0.1 || temp_ks_dis1 < 0.05) &&
      (temp_chi2 < temp_ndf * 12 || temp_chi2 - max_chi2 < (temp_ndf - 1) * 6)) {
    return true;
  }

  return false;
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
