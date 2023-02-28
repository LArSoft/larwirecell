/** A WCT component which passes through IFrames untouched and saves
 * their content as:

 - waveform content as vector collections of either raw::RawDigit or
   recob::Wire.

 - summaries as vector<double>

 - channel mask maps as vector<int> holding channel numbers

 It can be configured to scale waveform or summary values by some constant.
*/

#ifndef LARWIRECELL_COMPONENTS_FRAMESAVER
#define LARWIRECELL_COMPONENTS_FRAMESAVER

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcls {

  class FrameSaver : public IArtEventVisitor,
                     public WireCell::IFrameFilter,
                     public WireCell::IConfigurable {
  public:
    FrameSaver();
    virtual ~FrameSaver();

    /// IArtEventVisitor
    virtual void produces(art::ProducesCollector& collector);
    virtual void visit(art::Event& event);

    /// IFrameFilter
    virtual bool operator()(const WireCell::IFrame::pointer& inframe,
                            WireCell::IFrame::pointer& outframe);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

    /// A WCT frame includes the concept of a "channel mask map"
    /// which is a named association between a set of channels and
    /// a number of bin ranges (the masks) for each channel.  The
    /// name is meant to connote some semantic meaning to the
    /// masks (eg, "bad" or "noisy").  The bin ranges are
    /// expressed as a begin and end index into the associated
    /// channel waveform (aka the "trace"). The end index is one
    /// past the intended coverage in the usual C++ iterator
    /// convention.  This information is stored in two
    /// vecotor<int> data products which are given instance names
    /// based on the channel mask map name.  The first vector<int>
    /// named "<name>channels" holds a simple list of channel
    /// numbers which have at least one mask.  The second
    /// vector<int> named "<name>masks" is a flattened 3xN array
    /// of entries like: [(channel, begin, end), ...].  This
    /// contrivied storage format is to avoid creating an data
    /// product class.  Both are produced so that the job may be
    /// configured to drop one or both.
    typedef std::vector<int> channel_list;
    typedef std::vector<int> channel_masks;

  private:
    // ordered.
    std::map<int, geo::View_t> m_chview;

    WireCell::IFrame::pointer m_frame;
    std::vector<std::string> m_frame_tags, m_summary_tags;
    std::vector<double> m_frame_scale, m_summary_scale;

    typedef std::function<float(const std::vector<float>& tsvals)> summarizer_function;
    std::unordered_map<std::string, summarizer_function> m_summary_operators;

    int m_nticks;
    bool m_digitize, m_sparse;
    Json::Value m_cmms, m_pedestal_mean;
    double m_pedestal_sigma;

    void save_as_raw(art::Event& event);
    void save_as_cooked(art::Event& event);
    void save_summaries(art::Event& event);
    void save_cmms(art::Event& event);
    void save_empty(art::Event& event);
  };
}

#endif
