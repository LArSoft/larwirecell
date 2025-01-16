#include "TimingTPCBundle.h"

using namespace WireCell::QLMatch;
using namespace WireCell::PointCloud::Facade;

#include <boost/math/statistics/kolmogorov_smirnov.hpp>

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
  , ks_dis(1)
  , chi2(0)
  , ndf(0)
  , flag_high_consistent(false)
  , flag_spec_end(false)
  , flag_potential_bad_match(false)
  , strength(0)
{
  m_nchan = flash->get_num_channels();
  pred_pmt_light.resize(m_nchan, 0);
}

TimingTPCBundle::~TimingTPCBundle() {}

double TimingTPCBundle::get_total_pred_light()
{
  double sum = 0;
  for (size_t i = 0; i != pred_pmt_light.size(); i++) {
    sum += pred_pmt_light.at(i);
  }
  return sum;
}
bool TimingTPCBundle::examine_bundle(TimingTPCBundle* bundle,
                                   Double_t* cos_pe_low,
                                   Double_t* cos_pe_mid) {
    // Store data in vectors instead of histograms
    std::vector<double> measured_dist(m_nchan);
    std::vector<double> predicted_dist(m_nchan);
    
    double pe[m_nchan], pe_err[m_nchan];
    double pred_pe[m_nchan];
    
    // Fill the initial data arrays
    for (int i = 0; i != m_nchan; i++) {
        pe[i] = flash->get_PE(i);
        pe_err[i] = flash->get_PE_err(i);
        pred_pe[i] = pred_pmt_light.at(i) + bundle->get_pred_pmt_light().at(i);
    }
    
    // Process the data and fill distributions
    double total_predicted = 0;
    double total_measured = 0;
    for (int j = 0; j != m_nchan; j++) {
        measured_dist[j] = pe[j];
        if ((pred_pe[j] < cos_pe_low[j] || 
            (pred_pe[j] < cos_pe_mid[j] * 1.1 && pe[j] == 0)) &&
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
        boost::math::kolmogorov_smirnov_test(measured_dist, predicted_dist, &temp_ks_dis);
    }
    
    // Calculate chi-square statistics
    for (int j = 0; j != m_nchan; j++) {
        double cur_chi2 = 0;
        if (flag_close_to_PMT) {
            if (pe[j] - pred_pe[j] > 350 &&
                pe[j] > pred_pe[j] * 1.3) {
                cur_chi2 = pow(pred_pe[j] - pe[j], 2) / (pow(pe_err[j], 2) + pow(pe[j] * 0.5, 2));
            } else {
                cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
            }
        } else {
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
        pred_pe[max_bin] > 0) {
        temp_chi2 -= max_chi2 - 1;
    }
    
    if ((temp_ks_dis < ks_dis + 0.06 && (temp_ks_dis < ks_dis * 1.2) && temp_chi2 < chi2 + ndf * 5 &&
         temp_chi2 < chi2 * 1.21) ||
        (temp_ks_dis < ks_dis && temp_chi2 < chi2 + ndf * 10 && temp_chi2 < chi2 * 1.45)) {
        return true;
    } else {
        return false;
    }
}
bool TimingTPCBundle::examine_bundle_rank(TimingTPCBundle* bundle,
                                        Double_t* cos_pe_low,
                                        Double_t* cos_pe_mid)
{
    // Store data in vectors instead of histograms
    std::vector<double> measured_dist(m_nchan);
    std::vector<double> predicted_dist(m_nchan);
    
    double pe[m_nchan], pe_err[m_nchan];
    double pred_pe[m_nchan];

    // Fill the initial data arrays
    for (int i = 0; i != m_nchan; i++) {
        pe[i] = flash->get_PE(i);
        pe_err[i] = flash->get_PE_err(i);
        pred_pe[i] = pred_pmt_light.at(i) + bundle->get_pred_pmt_light().at(i);
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
        boost::math::kolmogorov_smirnov_test(measured_dist, predicted_dist, &temp_ks_dis);
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
  int main_cluster_id = main_cluster->get_cluster_id();

  ClusterSelection merge_clusters;
  for (size_t i = 0; i != other_clusters.size(); i++) {
    Cluster* temp_cluster = other_clusters.at(i);

    double dis_save = 1e9;

    {
      Cluster* cluster1 = temp_cluster;
      Cluster* cluster2 = main_cluster;
      SlimMergeGeomCell* prev_mcell1 = 0;
      SlimMergeGeomCell* prev_mcell2 = 0;
      SlimMergeGeomCell* mcell1 = 0;
      Point p1; //
      SlimMergeGeomCell* mcell2 = 0;
      Point p2;

      mcell1 = *(cluster1->get_time_cells_set_map().begin()->second.begin());
      p1 = mcell1->center();

      while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2) {
        prev_mcell1 = mcell1;
        prev_mcell2 = mcell2;

        // find the closest point and merged cell in cluster2
        std::pair<SlimMergeGeomCell*, Point> temp_results = cluster2->get_closest_point_mcell(p1);
        p2 = temp_results.second;
        mcell2 = temp_results.first;
        // find the closest point and merged cell in cluster1
        temp_results = cluster1->get_closest_point_mcell(p2);
        p1 = temp_results.second;
        mcell1 = temp_results.first;
      }
      double dis = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));

      if (dis < dis_save) { dis_save = dis; }

      prev_mcell1 = 0;
      prev_mcell2 = 0;

      mcell1 = *(cluster1->get_time_cells_set_map().rbegin()->second.begin());
      p1 = mcell1->center();

      while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2) {
        prev_mcell1 = mcell1;
        prev_mcell2 = mcell2;

        // find the closest point and merged cell in cluster2
        std::pair<SlimMergeGeomCell*, Point> temp_results = cluster2->get_closest_point_mcell(p1);
        p2 = temp_results.second;
        mcell2 = temp_results.first;
        // find the closest point and merged cell in cluster1
        temp_results = cluster1->get_closest_point_mcell(p2);
        p1 = temp_results.second;
        mcell1 = temp_results.first;
      }
      dis = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));

      if (dis < dis_save) { dis_save = dis; }
    }

    if (dis_save < dis_cut) { merge_clusters.push_back(temp_cluster); }
  }

  if (merge_clusters.size() > 0) {
    merge_clusters.push_back(main_cluster);
    Cluster* ncluster = new Cluster(main_cluster_id);
    for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
      Cluster* ocluster = *(it1);
      SMGCSelection& mcells = ocluster->get_mcells();
      for (auto it2 = mcells.begin(); it2 != mcells.end(); it2++) {
        SlimMergeGeomCell* mcell = (*it2);
        // std::cout << ocluster->get_cluster_id() << " " << mcell << std::endl;
        int time_slice = mcell->GetTimeSlice();
        ncluster->AddCell(mcell, time_slice);
      }
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
    main_cluster = ncluster;
  }
}

void TimingTPCBundle::add_bundle(TimingTPCBundle* bundle, Double_t* cos_pe_low, Double_t* cos_pe_mid)
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

  std::vector<double>& pes = bundle->get_pred_pmt_light();
  for (size_t i = 0; i != pred_pmt_light.size(); i++) {
    pred_pmt_light.at(i) += pes.at(i);
  }
  examine_bundle(cos_pe_low, cos_pe_mid);
}

bool TimingTPCBundle::examine_beam_bundle()
{
   // Store data in vectors instead of histograms
   std::vector<double> measured_dist(m_nchan);
   std::vector<double> predicted_dist(m_nchan);
   
   double pe[m_nchan], pe_err[m_nchan];
   double pred_pe[m_nchan];

   // Fill initial arrays
   for (int i = 0; i != m_nchan; i++) {
       pe[i] = flash->get_PE(i);
       pe_err[i] = flash->get_PE_err(i);
       pred_pe[i] = pred_pmt_light.at(i);
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
   double temp_ks_dis = 1.0;
   boost::math::kolmogorov_smirnov_test(norm_measured, norm_predicted, &temp_ks_dis);

   // Calculate chi-square statistics
   double temp_chi2 = 0;
   double temp_ndf = 0;
   double max_chi2 = 0;
   int max_bin = -1;

   for (int j = 0; j != m_nchan; j++) {
       double cur_chi2 = 0;
       if (flag_close_to_PMT) {
           if (pe[j] - pred_pe[j] > 350 &&
               pe[j] > pred_pe[j] * 1.3) {
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
   double temp_ks_dis1 = 1.0;
   boost::math::kolmogorov_smirnov_test(measured_dist2, predicted_dist2, &temp_ks_dis1);

   if ((temp_ks_dis < 0.1 || temp_ks_dis1 < 0.05) &&
       (temp_chi2 < temp_ndf * 12 || temp_chi2 - max_chi2 < (temp_ndf - 1) * 6)) {
       return true;
   }

   return false;
}

bool TimingTPCBundle::examine_bundle(Double_t* cos_pe_low, Double_t* cos_pe_mid)
{
    // Store data in vectors instead of histograms
    std::vector<double> measured_dist(m_nchan);
    std::vector<double> predicted_dist(m_nchan);
    
    double pe[m_nchan], pe_err[m_nchan];
    double pred_pe[m_nchan];

    // Fill initial arrays
    for (int i = 0; i != m_nchan; i++) {
        pe[i] = flash->get_PE(i);
        pe_err[i] = flash->get_PE_err(i);
        pred_pe[i] = pred_pmt_light.at(i);
    }

    // Process data and fill distributions
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

    // Calculate KS test using Boost
    if (total_predicted > 0) {
        boost::math::kolmogorov_smirnov_test(measured_dist, predicted_dist, &ks_dis);
    }

    chi2 = 0;
    ndf = 0;
    double max_chi2 = 0;
    int max_bin = -1;

    // Calculate chi-square statistics
    for (int j = 0; j != m_nchan; j++) {
        double cur_chi2 = 0;

        if (flag_close_to_PMT) {
            if (pe[j] - pred_pe[j] > 350 &&
                pe[j] > pred_pe[j] * 1.3) {
                cur_chi2 = pow(pred_pe[j] - pe[j], 2) / (pow(pe_err[j], 2) + pow(pe[j] * 0.5, 2));
            }
            else {
                cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
            }
        }
        else {
            cur_chi2 = pow(pred_pe[j] - pe[j], 2) / pow(pe_err[j], 2);
        }
        chi2 += cur_chi2;

        if (cur_chi2 > max_chi2) {
            max_chi2 = cur_chi2;
            max_bin = j;
        }

        if (pe[j] == 0 && pred_pe[j] == 0) {}
        else {
            ndf++;
        }
    }

    if (pe[max_bin] == 0 && pred_pe[max_bin] > 0) {
        chi2 -= max_chi2 - 1;
    }

    // Check multiple consistency conditions
    flag_high_consistent = false;
    if (ks_dis < 0.06 && ndf >= 3 && chi2 < ndf * 36) { 
        flag_high_consistent = true; 
    }
    else if (ks_dis < 0.05 && ndf >= 6 && chi2 < ndf * 45) {
        flag_high_consistent = true;
    }
    else if (ks_dis < 0.12 && ndf >= 3 && chi2 < ndf * 25) {
        flag_high_consistent = true;
    }
    else if (flag_at_x_boundary && ndf >= 2 && chi2 < 9 * ndf && ks_dis < 0.12) {
        flag_high_consistent = true;
    }
    else if (flag_at_x_boundary && ndf >= 1 && chi2 < 3 * ndf && ks_dis < 0.12) {
        flag_high_consistent = true;
    }
    else if (chi2 < 4 * ndf && ndf >= 3 && ks_dis < 0.15) {
        flag_high_consistent = true;
    }
    else if (chi2 < 1.5 * ndf && ks_dis < 0.2 && ndf >= 3) {
        flag_high_consistent = true;
    }
    else if (ks_dis < 0.12 && ndf >= 5 && chi2 < ndf * 55 && flag_close_to_PMT) {
        flag_high_consistent = true;
    }
    else if (ks_dis < 0.14 && ndf >= 3 && chi2 < ndf * 6) {
        flag_high_consistent = true;
    }

    if (flag_high_consistent) { 
        return true; 
    }
    else {
        // Check fired PMTs statistics
        double ntot = 0;
        double ntot1 = 0;
        double nfired = 0;
        double nfired1 = 0;
        for (int j = 0; j != m_nchan; j++) {
            if (pred_pe[j] > 0.33) {
                ntot++;
                if (pe[j] > 0.33 * pred_pe[j]) nfired++;
            }
            if (pred_pmt_light.at(j) > 1.0) {
                ntot1++;
                if (pe[j] > 0.33 * pred_pmt_light.at(j)) nfired1++;
            }
        }

        // Handle special cases
        if (nfired <= 1 && (nfired1 != 0 && nfired1 > 0.75 * ntot1)) return true;
        if (nfired <= 2 && ks_dis > 0.8 && chi2 > 60 * ndf && ndf >= 6) return false;

        if (flag_at_x_boundary && (!flag_close_to_PMT)) {
            if (nfired == 0) { flag_potential_bad_match = true; }
            if (nfired == 1 && ntot <= 2 && nfired1 < 0.2 * ntot1) { flag_potential_bad_match = true; }
            if (nfired < 0.5 * ntot && ntot - nfired >= 2) { flag_potential_bad_match = true; }
        }
        else {
            if (nfired == 0) {
                flag_potential_bad_match = true;
                return false;
            }
            if (nfired == 1 && ntot <= 2 && nfired1 < 0.2 * ntot1) {
                flag_potential_bad_match = true;
                return false;
            }
            if (nfired < 0.5 * ntot && ntot - nfired >= 2) {
                flag_potential_bad_match = true;
                return false;
            }
        }
    }
    return true;
}
