#include "OpFlashSource.h"
#include "TTimeStamp.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/NamedFactory.h"
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

OpFlashSource::OpFlashSource() {}

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
  art::Handle<std::vector<recob::OpFlash>> opflashes;
  event.getByLabel(m_inputTag, opflashes);
  if (!opflashes.isValid()) {
    THROW(ValueError() << errmsg{"WireCell::OpFlashSource failed to get opflashes"});
  }
  for (auto const& opflash : *opflashes) {
    std::cout << "OpFlash time: " << opflash.Time() << " " << opflash.PEs().size() << std::endl;
  }

  //   const auto nflashes = opflashes->size();
  //   // Create a 2D boost::multi_array with shape 16 x nflashes
  //   typedef boost::multi_array<double, 2> MultiArray;
  //   MultiArray array(boost::extents[16][nflashes]);

  //   // ...

  //   for (auto const& opflash : *opflashes) {
  //     auto tensor = std::make_shared<SimpleTensor<double>>(
  //       opflash.Time(), opflash.TimeWidth(), opflash.TotalPE());
  //     m_tensorsets.push_back(tensor);

  //     // Store the values in the boost::multi_array
  //     // Assuming you want to store the values in the first dimension (16)
  //     for (int i = 0; i < 16; ++i) {
  //       array[i][m_tensorsets.size() - 1] = tensor->at(i);
  //     }
  //   }

  //   // ...

  //   for (auto const& opflash : *opflashes) {
  //     auto tensor = std::make_shared<SimpleTensor<double>>(
  //       opflash.Time(), opflash.TimeWidth(), opflash.TotalPE());
  //     m_tensorsets.push_back(tensor);
  //   }
}

bool OpFlashSource::operator()(WireCell::ITensorSet::pointer& tensorset)
{
  if (m_tensorsets.empty()) { return false; }
  tensorset = m_tensorsets.front();
  m_tensorsets.pop_front();
  return true;
}