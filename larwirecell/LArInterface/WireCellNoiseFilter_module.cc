// Framework includes
#include "art/Framework/Principal/Event.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibProvider.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"

#include "lardataobj/RawData/RawDigit.h"

#include "WireCellUtil/Units.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellSigProc/Microboone.h"
#include "WireCellSigProc/OmnibusNoiseFilter.h"
#include "WireCellSigProc/SimpleChannelNoiseDB.h"

#include <numeric> // iota
#include <string>

using namespace WireCell;

namespace noisefilteralg {

  class WireCellNoiseFilter : public art::EDProducer {

  public:
    explicit WireCellNoiseFilter(fhicl::ParameterSet const& pset);
    virtual ~WireCellNoiseFilter();

    void produce(art::Event& evt);
    void reconfigure(fhicl::ParameterSet const& pset);

  private:
    void DoNoiseFilter(art::Event const& evt,
                       const std::vector<raw::RawDigit>&,
                       std::vector<raw::RawDigit>&) const;

    //******************************
    //Variables Taken from FHICL File
    std::string fDigitModuleLabel; // label for rawdigit module
    bool fDoNoiseFiltering;        // Allows for a "pass through" mode
    size_t fNumTicksToDropFront;   // If we are truncating then this is non-zero
    size_t fWindowSize;            // Number of ticks in the output RawDigit

    // services
  }; //end class Noise

  //-------------------------------------------------------------------
  WireCellNoiseFilter::WireCellNoiseFilter(fhicl::ParameterSet const& pset) : EDProducer(pset)
  {
    this->reconfigure(pset);
    produces<std::vector<raw::RawDigit>>();
  }

  //-------------------------------------------------------------------
  WireCellNoiseFilter::~WireCellNoiseFilter() {}

  //-------------------------------------------------------------------
  void
  WireCellNoiseFilter::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDigitModuleLabel = pset.get<std::string>("DigitModuleLabel", "daq");
    fDoNoiseFiltering = pset.get<bool>("DoNoiseFiltering", true);
    fNumTicksToDropFront = pset.get<size_t>("NumTicksToDropFront", 2400);
    fWindowSize = pset.get<size_t>("WindowSize", 6400);
  }

  //-------------------------------------------------------------------
  void
  WireCellNoiseFilter::produce(art::Event& evt)
  {
    // Recover services we will need
    const lariov::DetPedestalProvider& pedestalValues =
      art::ServiceHandle<lariov::DetPedestalService const>()->GetPedestalProvider();

    art::Handle<std::vector<raw::RawDigit>> rawDigitHandle;
    evt.getByLabel(fDigitModuleLabel, rawDigitHandle);

    // Define the output vector (in case we don't do anything)
    std::unique_ptr<std::vector<raw::RawDigit>> filteredRawDigit(new std::vector<raw::RawDigit>);

    if (rawDigitHandle.isValid() && rawDigitHandle->size() > 0) {
      const std::vector<raw::RawDigit>& rawDigitVector(*rawDigitHandle);

      // Make sure we have the correct window size (e.g. window size = 9600 but data is 9595)
      size_t windowSize(std::min(fWindowSize, rawDigitVector.at(0).NADC()));

      if (fNumTicksToDropFront + windowSize > rawDigitVector.at(0).NADC())
        throw cet::exception("WireCellNoiseFilter")
          << "Ticks to drop + windowsize larger than input buffer\n";

      if (fDoNoiseFiltering)
        DoNoiseFilter(evt, rawDigitVector, *filteredRawDigit);
      else {
        // Enable truncation
        size_t startBin(fNumTicksToDropFront);
        size_t stopBin(startBin + windowSize);

        raw::RawDigit::ADCvector_t outputVector(windowSize);

        for (const auto& rawDigit : rawDigitVector) {
          if (rawDigit.NADC() < windowSize) continue;

          const raw::RawDigit::ADCvector_t& rawAdcVec = rawDigit.ADCs();

          unsigned int channel = rawDigit.Channel();
          float pedestal = pedestalValues.PedMean(evt.time().value(), channel);

          std::copy(
            rawAdcVec.begin() + startBin, rawAdcVec.begin() + stopBin, outputVector.begin());

          filteredRawDigit->emplace_back(
            raw::RawDigit(channel, outputVector.size(), outputVector, raw::kNone));
          filteredRawDigit->back().SetPedestal(pedestal, 2.0);
        }
      }
    }

    //filtered raw digits
    evt.put(std::move(filteredRawDigit));
  }

  void
  WireCellNoiseFilter::DoNoiseFilter(art::Event const& e,
                                     const std::vector<raw::RawDigit>& inputWaveforms,
                                     std::vector<raw::RawDigit>& outputWaveforms) const
  {
    auto const runNum = e.run();

    // Recover services we will need
    const lariov::ChannelStatusProvider& channelStatus =
      art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();
    const lariov::DetPedestalProvider& pedestalValues =
      art::ServiceHandle<lariov::DetPedestalService const>()->GetPedestalProvider();
    const lariov::ElectronicsCalibProvider& elec_provider =
      art::ServiceHandle<lariov::ElectronicsCalibService const>()->GetProvider();
    const geo::GeometryCore& geometry = *lar::providerFrom<geo::Geometry>();
    auto const clock_data =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

    const unsigned int n_channels = inputWaveforms.size();

    // S&C microboone sampling parameter database
    const double tick = sampling_rate(clock_data); // 0.5 * units::microsecond;
    const size_t nsamples = inputWaveforms.at(0).NADC();
    const size_t windowSize = std::min(fWindowSize, nsamples);

    // Q&D microboone channel map
    std::vector<int> uchans(geometry.Nwires(0)), vchans(geometry.Nwires(1)),
      wchans(geometry.Nwires(2));
    const int nchans = uchans.size() + vchans.size() + wchans.size();
    std::iota(uchans.begin(), uchans.end(), 0);
    std::iota(vchans.begin(), vchans.end(), vchans.size());
    std::iota(wchans.begin(), wchans.end(), vchans.size() + uchans.size());

    // Q&D nominal baseline
    const double unombl = 2048.0, vnombl = 2048.0, wnombl = 400.0;

    // Q&D miss-configured channel database
    std::vector<int> miscfgchan;
    for (int channelIdx = 0; channelIdx < nchans; channelIdx++) {
      if (elec_provider.ExtraInfo(channelIdx).GetBoolData("is_misconfigured")) {
        miscfgchan.push_back(channelIdx);
      }
    }

    const double from_gain_mVfC = 4.7, to_gain_mVfC = 14.0, from_shaping = 1.0 * units::microsecond,
                 to_shaping = 2.0 * units::microsecond;

    // Recover bad channels from the database
    std::vector<int> bad_channels;
    for (int channelIdx = 0; channelIdx < nchans; channelIdx++)
      if (channelStatus.IsBad(e.time().value(), channelIdx)) bad_channels.push_back(channelIdx);

    // Q&D RC+RC time constant - all have same.
    const double rcrc = 1.0 * units::millisecond;
    std::vector<int> rcrcchans(nchans);
    std::iota(rcrcchans.begin(), rcrcchans.end(), 0);

    //harmonic noises
    std::vector<int> harmonicchans(uchans.size() + vchans.size());
    std::iota(harmonicchans.begin(), harmonicchans.end(), 0);

    std::vector<int> special_chans;
    special_chans.push_back(2240);

    WireCell::SigProc::SimpleChannelNoiseDB::mask_t h36kHz(0, 169, 173);
    WireCell::SigProc::SimpleChannelNoiseDB::mask_t h108kHz(0, 513, 516);
    WireCell::SigProc::SimpleChannelNoiseDB::mask_t hspkHz(0, 17, 19);
    WireCell::SigProc::SimpleChannelNoiseDB::multimask_t hharmonic;
    hharmonic.push_back(h36kHz);
    hharmonic.push_back(h108kHz);
    WireCell::SigProc::SimpleChannelNoiseDB::multimask_t hspecial;
    hspecial.push_back(h36kHz);
    hspecial.push_back(h108kHz);
    hspecial.push_back(hspkHz);

    float u_resp_array[120] = {
      0.0,          0.0,          0.0,          0.0,          0.0,          0.0,
      0.0,          0.0,          0.0,          0.0,          0.364382,     0.387949,
      0.411053,     0.433979,     0.456863,     0.479746,     0.502641,     0.52554,
      0.548441,     0.57134,      0.591765,     0.609448,     0.626848,     0.644094,
      0.661364,     0.678859,     0.695231,     0.710462,     0.726147,     0.742373,
      0.761332,     0.783313,     0.806325,     0.830412,     0.857676,     0.888412,
      0.920705,     0.954624,     0.990242,     1.02766,      1.06121,      1.09027,
      1.12037,      1.15157,      1.18392,      1.21748,      1.25229,      1.28824,
      1.32509,      1.36256,      1.40051,      1.43907,      1.47857,      1.51933,
      1.56134,      1.60404,      1.72665,      1.94005,      2.16994,      2.42041,
      2.69475,      3.07222,      3.67375,      4.60766,      5.91864,      7.30178,
      8.3715,       8.94736,      8.93705,      8.40339,      7.2212,       5.76382,
      3.8931,       1.07893,      -3.52481,     -11.4593,     -20.4011,     -29.1259,
      -34.9544,     -36.9358,     -35.3303,     -31.2068,     -25.8614,     -20.3613,
      -15.3794,     -11.2266,     -7.96091,     -5.50138,     -3.71143,     -2.44637,
      -1.57662,     -0.99733,     -0.62554,     -0.393562,    -0.249715,    -0.15914,
      -0.100771,    -0.062443,    -0.037283,    -0.0211508,   -0.0112448,   -0.00552085,
      -0.00245133,  -0.000957821, -0.000316912, -8.51679e-05, -2.21299e-05, -1.37496e-05,
      -1.49806e-05, -1.36935e-05, -9.66758e-06, -5.20773e-06, -7.4787e-07,  3.71199e-06,
      8.17184e-06,  1.26317e-05,  1.70916e-05,  2.15514e-05,  2.60113e-05,  3.04711e-05};
    float v_resp_array[120] = {
      0.0,          0.0,          0.0,          0.0,          0.0,          0.0,
      0.0,          0.0,          0.0,          0.0,          0.0865303,    0.0925559,
      0.0983619,    0.104068,     0.109739,     0.115403,     0.121068,     0.126735,
      0.132403,     0.138072,     0.143739,     0.149408,     0.155085,     0.160791,
      0.166565,     0.172454,     0.178514,     0.184795,     0.191341,     0.198192,
      0.205382,     0.212944,     0.220905,     0.229292,     0.238129,     0.247441,
      0.257256,     0.267601,     0.278502,     0.28999,      0.298745,     0.304378,
      0.310105,     0.315921,     0.321818,     0.327796,     0.333852,     0.339967,
      0.346098,     0.352169,     0.358103,     0.363859,     0.36945,      0.374915,
      0.380261,     0.385401,     0.39016,      0.394378,     0.39804,      0.401394,
      0.405145,     0.410714,     0.4205,       0.437951,     0.467841,     0.516042,
      0.587738,     0.694157,     0.840763,     1.01966,      1.22894,      1.5612,
      2.12348,      3.31455,      5.59355,      9.10709,      14.1756,      18.4603,
      19.9517,      17.4166,      10.6683,      1.40656,      -10.0638,     -19.034,
      -23.654,      -24.0558,     -21.4418,     -17.3229,     -12.9485,     -9.08912,
      -6.05941,     -3.86946,     -2.38669,     -1.43678,     -0.853335,    -0.503951,
      -0.296551,    -0.173029,    -0.0990099,   -0.0547172,   -0.0287882,   -0.0142758,
      -0.00661815,  -0.00284757,  -0.00115702,  -0.000456456, -0.000183439, -8.04214e-05,
      -4.20533e-05, -2.62903e-05, -1.64098e-05, -6.68039e-06, 3.04903e-06,  1.27784e-05,
      2.25079e-05,  3.22373e-05,  4.19667e-05,  5.16961e-05,  6.14255e-05,  7.11549e-05};
    WireCell::Waveform::realseq_t u_resp(nsamples);
    WireCell::Waveform::realseq_t v_resp(nsamples);
    for (int i = 0; i != 120; i++) {
      u_resp.at(i) = u_resp_array[i];
      v_resp.at(i) = v_resp_array[i];
    }
    WireCell::Waveform::compseq_t u_resp_freq = WireCell::Waveform::dft(u_resp);
    WireCell::Waveform::compseq_t v_resp_freq = WireCell::Waveform::dft(v_resp);

    int uplane_time_shift = 79;
    int vplane_time_shift = 82;

    // do the coherent subtraction
    std::vector<std::vector<int>> channel_groups;
    for (unsigned int i = 0; i != 172; i++) {
      //for (int i=150;i!=151;i++){
      std::vector<int> channel_group;
      for (int j = 0; j != 48; j++) {
        channel_group.push_back(i * 48 + j);
      }
      channel_groups.push_back(channel_group);
    }

    auto noise = new WireCell::SigProc::SimpleChannelNoiseDB;
    // initialize
    noise->set_sampling(tick, nsamples);
    // set nominal baseline
    noise->set_nominal_baseline(uchans, unombl);
    noise->set_nominal_baseline(vchans, vnombl);
    noise->set_nominal_baseline(wchans, wnombl);

    noise->set_response(uchans, u_resp_freq);
    noise->set_response(vchans, v_resp_freq);

    noise->set_response_offset(uchans, uplane_time_shift);
    noise->set_response_offset(vchans, vplane_time_shift);

    noise->set_pad_window_front(uchans, 20);
    noise->set_pad_window_back(uchans, 10);
    noise->set_pad_window_front(vchans, 10);
    noise->set_pad_window_back(vchans, 10);
    noise->set_pad_window_front(wchans, 10);
    noise->set_pad_window_back(wchans, 10);

    // set misconfigured channels
    noise->set_gains_shapings(miscfgchan, from_gain_mVfC, to_gain_mVfC, from_shaping, to_shaping);
    // do the RCRC
    noise->set_rcrc_constant(rcrcchans, rcrc);
    // set initial bad channels
    noise->set_bad_channels(bad_channels);
    // set the harmonic filter
    noise->set_filter(harmonicchans, hharmonic);
    noise->set_filter(special_chans, hspecial);
    noise->set_channel_groups(channel_groups);

    for (unsigned int i = 0; i != uchans.size(); i++) {
      if (uchans.at(i) < 100) {
        noise->set_min_rms_cut_one(uchans.at(i), 1);
        noise->set_max_rms_cut_one(uchans.at(i), 5);
      }
      else if (uchans.at(i) < 2000) {
        noise->set_min_rms_cut_one(uchans.at(i), 1.9);
        noise->set_max_rms_cut_one(uchans.at(i), 11);
      }
      else {
        noise->set_min_rms_cut_one(uchans.at(i), 0.9);
        noise->set_max_rms_cut_one(uchans.at(i), 5);
      }
    }
    for (unsigned int i = 0; i != vchans.size(); i++) {
      if (vchans.at(i) < 290 + (int)uchans.size()) {
        noise->set_min_rms_cut_one(vchans.at(i), 1);
        noise->set_max_rms_cut_one(vchans.at(i), 5);
      }
      else if (vchans.at(i) < 2200 + (int)uchans.size()) {
        noise->set_min_rms_cut_one(vchans.at(i), 1.9);
        noise->set_max_rms_cut_one(vchans.at(i), 11);
      }
      else {
        noise->set_min_rms_cut_one(vchans.at(i), 1);
        noise->set_max_rms_cut_one(vchans.at(i), 5);
      }
    }

    // these are the one after the Hardware Fix
    if (runNum > 8000) {
      for (unsigned int i = 0; i != uchans.size(); i++) {
        if (uchans.at(i) < 600) {
          noise->set_min_rms_cut_one(uchans.at(i), 1 + (1.7 - 1) / 600. * i);
          noise->set_max_rms_cut_one(uchans.at(i), 5);
        }
        else if (uchans.at(i) < 1800) {
          noise->set_min_rms_cut_one(uchans.at(i), 1.7);
          noise->set_max_rms_cut_one(uchans.at(i), 11);
        }
        else {
          noise->set_min_rms_cut_one(uchans.at(i), 1 + (1.7 - 1) / 600. * (2399 - i));
          noise->set_max_rms_cut_one(uchans.at(i), 5);
        }
      }
      for (unsigned int i = 0; i != vchans.size(); i++) {
        if (vchans.at(i) < 600 + (int)uchans.size()) {
          noise->set_min_rms_cut_one(vchans.at(i), 0.8 + (1.7 - 0.8) / 600. * i);
          noise->set_max_rms_cut_one(vchans.at(i), 5);
        }
        else if (vchans.at(i) < 1800 + (int)uchans.size()) {
          noise->set_min_rms_cut_one(vchans.at(i), 1.7);
          noise->set_max_rms_cut_one(vchans.at(i), 11);
        }
        else {
          noise->set_min_rms_cut_one(vchans.at(i), 0.8 + (1.7 - 0.8) / 600. * (2399 - i));
          noise->set_max_rms_cut_one(vchans.at(i), 5);
        }
      }
    }

    noise->set_min_rms_cut(wchans, 1.3);
    noise->set_max_rms_cut(wchans, 8.0);

    //Define database object
    std::shared_ptr<WireCell::IChannelNoiseDatabase> noise_sp(noise);

    auto one = new WireCell::SigProc::Microboone::OneChannelNoise;
    one->set_channel_noisedb(noise_sp);
    std::shared_ptr<WireCell::IChannelFilter> one_sp(one);
    auto many = new WireCell::SigProc::Microboone::CoherentNoiseSub;
    many->set_channel_noisedb(noise_sp);
    std::shared_ptr<WireCell::IChannelFilter> many_sp(many);

    //define noisefilter object
    WireCell::SigProc::OmnibusNoiseFilter bus;
    bus.set_channel_filters({one_sp});
    bus.set_grouped_filters({many_sp});
    bus.set_channel_noisedb(noise_sp);

    // Enable truncation
    size_t startBin(fNumTicksToDropFront);
    size_t stopBin(startBin + windowSize);

    //load waveforms into traces
    WireCell::ITrace::vector traces;
    for (unsigned int ich = 0; ich < n_channels; ich++) {
      if (inputWaveforms.at(ich).NADC() < windowSize) continue;

      const raw::RawDigit::ADCvector_t& rawAdcVec = inputWaveforms.at(ich).ADCs();

      WireCell::ITrace::ChargeSequence charges;

      charges.resize(nsamples);

      std::transform(rawAdcVec.begin(), rawAdcVec.end(), charges.begin(), [](auto& adcVal) {
        return float(adcVal);
      });

      unsigned int chan = inputWaveforms.at(ich).Channel();
      WireCell::SimpleTrace* st = new WireCell::SimpleTrace(chan, 0.0, charges);
      traces.push_back(WireCell::ITrace::pointer(st));
    }

    //Load traces into frame
    WireCell::SimpleFrame* sf = new WireCell::SimpleFrame(0, 0, traces);
    WireCell::IFrame::pointer frame = WireCell::IFrame::pointer(sf);
    WireCell::IFrame::pointer quiet = NULL;

    //Do filtering
    bus(frame, quiet);

    //std::cout << "HERE" << std::endl;
    //return;
    if (quiet == NULL) return;

    //Output results
    std::vector<short> waveform(windowSize);

    auto quiet_traces = quiet->traces();
    for (auto quiet_trace : *quiet_traces.get()) {
      //int tbin = quiet_trace->tbin();
      unsigned int channel = quiet_trace->channel();

      auto& quiet_charges = quiet_trace->charge();

      // Recover the database version of the pedestal, we'll offset the waveforms so it matches
      float pedestal = pedestalValues.PedMean(e.time().value(), channel);

      std::transform(quiet_charges.begin() + startBin,
                     quiet_charges.begin() + stopBin,
                     waveform.begin(),
                     [pedestal](auto charge) { return std::round(charge + pedestal); });

      outputWaveforms.emplace_back(raw::RawDigit(channel, waveform.size(), waveform, raw::kNone));
      outputWaveforms.back().SetPedestal(pedestal, 1.75);
    }

    return;
  }

  DEFINE_ART_MODULE(WireCellNoiseFilter)

} //end namespace WireCellNoiseFilteralg
