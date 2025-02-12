#include "OpFlashSource.h"
#include "TTimeStamp.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Units.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include <boost/multi_array.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

WIRECELL_FACTORY(wclsOpFlashSource,
                 wcls::OpFlashSource,
                 wcls::IArtEventVisitor,
                 WireCell::ITensorSetSource)

using namespace wcls;
using namespace WireCell;
using WireCell::Aux::SimpleTensor;
using WireCell::Aux::SimpleTensorSet;
using WireCell::String::format;

OpFlashSource::OpFlashSource() : Aux::Logger("OpFlashSource", "op") {}

OpFlashSource::~OpFlashSource() {}

WireCell::Configuration OpFlashSource::default_configuration() const
{
  Configuration cfg;
  cfg["art_tag"] = ""; // how to look up the opflashes
  return cfg;
}

void OpFlashSource::configure(const WireCell::Configuration& cfg)
{
  const std::string art_tag = cfg["art_tag"].asString();
  if (art_tag.empty()) {
    THROW(ValueError() << errmsg{"WireCell::OpFlashSource requires a source_label"});
  }
  m_inputTag = cfg["art_tag"].asString();
}

void OpFlashSource::visit(art::Event& event)
{
  log->debug("OpFlashSource::visit {}", m_inputTag);
  art::Handle<std::vector<recob::OpFlash>> opflashes;
  event.getByLabel(m_inputTag, opflashes);
  if (!opflashes.isValid()) {
    THROW(ValueError() << errmsg{"WireCell::OpFlashSource failed to get opflashes"});
  }
  for (auto const& opflash : *opflashes) {
    log->debug("OpFlash time: {} PEs: {}", opflash.Time(), opflash.PEs().size());
  }

  const auto nflashes = opflashes->size();
  // Create a 2D boost::multi_array with shape 16 x nflashes
  typedef boost::multi_array<double, 2> MultiArray;
  MultiArray array(boost::extents[nflashes][m_npmts + 1]);

  for (size_t iflash = 0; iflash < nflashes; ++iflash) {
    const auto& opflash = opflashes->at(iflash);
    array[iflash][0] = opflash.Time() * units::microsecond;
    const auto& pes = opflash.PEs();
    if (pes.size() > m_npmts) {
      raise<ValueError>(
        format("WireCell::OpFlashSource got unexpected number of PMTs expecting %d got %d",
               m_npmts,
               pes.size()));
    }
    for (size_t ipmt = 0; ipmt < pes.size(); ++ipmt) {
      array[iflash][ipmt + 1] = pes[ipmt];
    }
  }
  std::vector<size_t> shape = {array.shape()[0], array.shape()[1]};
  Json::Value md = Json::objectValue;
  auto tensor = std::make_shared<SimpleTensor>(shape, array.data(), md);
  ITensor::vector* itv = new ITensor::vector;
  itv->push_back(tensor);
  Configuration set_md;
  auto tset = std::make_shared<SimpleTensorSet>(event.event(), set_md, ITensor::shared_vector(itv));
  m_tensorsets.push_back(tset);
}

bool OpFlashSource::operator()(WireCell::ITensorSet::pointer& tensorset)
{
  if (m_tensorsets.empty()) { return false; }
  tensorset = m_tensorsets.front();
  m_tensorsets.pop_front();
  return true;
}