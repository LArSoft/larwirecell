#ifndef WIRECELLDEV_QLMATCH_UTIL
#define WIRECELLDEV_QLMATCH_UTIL

#include "WireCellIface/ITensorSet.h"
#include "WireCellUtil/PointTree.h"

#include "TimingTPCBundle.h"

namespace WireCell::QLMatch {
  void dump_bee_3d(const PointCloud::Tree::Points::node_t& root, const std::string& fn);
  void dump_bee_flash(const ITensorSet::pointer& ts, const std::string& fn);
  void dump_bee_bundle(const FlashBundlesMap& f2bundle,
                       const std::map<WireCell::PointCloud::Facade::Cluster*, int>& cluster_idx_map,
                       const std::string& fn);
}

#endif