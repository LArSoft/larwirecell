/**
   This channel noise database class inherits from OmniChannelNoiseDB
   and thus may be fully configured via WCT configuration.

   It additionally provides some select and configurable means to
   override this configuration with information taken from Art
   services.

   - If `bad_channel_policy` is set then take bad channels from the
   ChannelStatusService.  Example WCT config:

      bad_channel: {policy: "replace"},
      bad_channel: {policy: "union"},

   - If `misconfig_channel_policy` is set then take misconfigured
   channels from ElectronicsCalibService. Example WCT config:

      misconfig_channel: {
           policy: "replace",
           from: {gain:  4.7*wc.mV/wc.fC, shaping: 1.1*wc.us},
           to:   {gain: 14.0*wc.mV/wc.fC, shaping: 2.2*wc.us},
      }


   Fixme: original noise filter also used the following and needs to
   be handled somewhere:

    - DetPedestalService: to learn and restore a per-channel baseline.

*/

#include "ChannelNoiseDB.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibProvider.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(wclsChannelNoiseDB,
                 wcls::ChannelNoiseDB,
                 wcls::IArtEventVisitor,
                 WireCell::IChannelNoiseDatabase)

using namespace WireCell;

wcls::ChannelNoiseDB::ChannelNoiseDB()
  : OmniChannelNoiseDB(), m_bad_channel_policy(kNothing), m_misconfig_channel_policy(kNothing)
{}

wcls::ChannelNoiseDB::~ChannelNoiseDB() {}

void wcls::ChannelNoiseDB::visit(art::Event& event)
{
  if ((!m_bad_channel_policy) && (!m_misconfig_channel_policy)) {
    return; // no override
  }

  // FIXME: the current assumption in this code is that LS channel
  // numbers are identified with WCT channel IDs.  For MicroBooNE
  // this holds but in general some translation is needed here.
  auto nchans = art::ServiceHandle<geo::WireReadout const>()->Get().Nchannels();

  if (m_bad_channel_policy) {
    auto const& csvc = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    std::vector<int> bad_channels;
    for (size_t ich = 0; ich < nchans; ++ich) {
      if (csvc.IsBad(ich)) { bad_channels.push_back(ich); }
    }
    OmniChannelNoiseDB::set_bad_channels(bad_channels);
  }

  if (m_misconfig_channel_policy) {
    const auto& esvc = art::ServiceHandle<lariov::ElectronicsCalibService const>()->GetProvider();

    std::vector<int> mc_channels;
    for (size_t ich = 0; ich < nchans; ++ich) {
      if (esvc.ExtraInfo(ich).GetBoolData("is_misconfigured")) { mc_channels.push_back(ich); }
    }
    OmniChannelNoiseDB::set_misconfigured(
      mc_channels, m_fgstgs[0], m_fgstgs[1], m_fgstgs[2], m_fgstgs[3]);
  }
}

wcls::ChannelNoiseDB::OverridePolicy_t wcls::ChannelNoiseDB::parse_policy(const Configuration& jpol)
{
  if (jpol.empty()) {
    THROW(ValueError() << errmsg{"ChannelNoiseDB: empty override policy given"});
  }

  std::string pol = jpol.asString();

  if (pol == "union") { return kUnion; }

  if (pol == "replace") { return kReplace; }

  THROW(ValueError() << errmsg{"ChannelNoiseDB: unknown override policy given: " + pol});
}

void wcls::ChannelNoiseDB::configure(const WireCell::Configuration& cfg)
{
  // forward
  OmniChannelNoiseDB::configure(cfg);

  auto jbc = cfg["bad_channel"];
  if (!jbc.empty()) { m_bad_channel_policy = parse_policy(jbc["policy"]); }

  auto jmc = cfg["misconfig_channel"];
  if (!jmc.empty()) {
    m_misconfig_channel_policy = parse_policy(jmc["policy"]);

    // stash this for later
    m_fgstgs[0] = jmc["from"]["gain"].asDouble();
    m_fgstgs[1] = jmc["from"]["shaping"].asDouble();
    m_fgstgs[2] = jmc["to"]["gain"].asDouble();
    m_fgstgs[3] = jmc["to"]["shaping"].asDouble();
  }
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
