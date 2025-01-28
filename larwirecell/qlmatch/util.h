#ifndef WIRECELLDEV_QLMATCH_UTIL
#define WIRECELLDEV_QLMATCH_UTIL

#include "WireCellIface/ITensorSet.h"
#include "WireCellUtil/PointTree.h"

namespace WireCell::QLMatch {
  void dump_bee_3d(const PointCloud::Tree::Points::node_t& root, const std::string& fn);
  void dump_bee_flash(const ITensorSet::pointer& ts, const std::string& fn);
}

#endif