/**
    This is to input dynamic list of channels to DBChannelSelector
    e.g. misconfigured channels or bad channels
*/

#include "art/Framework/Principal/Event.h"

#include "ChannelSelectorDB.h"

#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibProvider.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(wclsChannelSelectorDB, wcls::ChannelSelectorDB, wcls::ChannelNoiseDB)

using namespace WireCell;

wcls::ChannelSelectorDB::ChannelSelectorDB() {}

wcls::ChannelSelectorDB::~ChannelSelectorDB() {}

void wcls::ChannelSelectorDB::visit(art::Event& event)
{
  std::cerr << "ChannelSelectorDB art::Event visit!\n";

  // FIXME: the current assumption in this code is that LS channel
  // numbers are identified with WCT channel IDs.  For MicroBooNE
  // this holds but in general some translation is needed here.
  auto const& gc = *lar::providerFrom<geo::Geometry>();
  auto nchans = gc.Nchannels();

  if (m_type == "bad") {
    auto const& csvc = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    for (size_t ich = 0; ich < nchans; ++ich) {
      if (csvc.IsBad(event.time().value(), ich)) { m_bad_channels.push_back(ich); }
    }
  }

  if (m_type == "misconfigured") {
    const auto& esvc = art::ServiceHandle<lariov::ElectronicsCalibService const>()->GetProvider();

    std::vector<int> mc_channels;
    for (size_t ich = 0; ich < nchans; ++ich) {
      if (esvc.ExtraInfo(event.time().value(), ich).GetBoolData("is_misconfigured")) {
        m_miscfg_channels.push_back(ich);
      }
    }
  }
}

void wcls::ChannelSelectorDB::configure(const WireCell::Configuration& cfg)
{
  m_type = get<std::string>(cfg, "type", "misconfigured");
  if (m_type != "bad" && m_type != "misconfigured") {
    THROW(
      ValueError() << errmsg{String::format("Channel type \"%s\" cannot be identified.", m_type)});
  }
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
