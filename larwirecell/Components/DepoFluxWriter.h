#ifndef LARWIRECELL_COMPONENTS_DEPOFLUXWRITER
#define LARWIRECELL_COMPONENTS_DEPOFLUXWRITER

#include "canvas/Utilities/InputTag.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDepoSetFilter.h"
#include "WireCellUtil/Binning.h"

namespace wcls {

  class DepoFluxWriter : public IArtEventVisitor,
                         public WireCell::IDepoSetFilter,
                         public WireCell::IConfigurable {
  public:
    /// IArtEventVisitor
    virtual void produces(art::ProducesCollector& collector);
    virtual void visit(art::Event& event);

    /// IDepoSetFilter
    virtual bool operator()(const WireCell::IDepoSet::pointer& indepos,
                            WireCell::IDepoSet::pointer& outdepos);

    /// IConfigurable
    virtual WireCell::Configuration default_configuration() const;
    virtual void configure(const WireCell::Configuration& config);

  private:
    // Configuration:
    //
    // anodes - a string or array of string naming IAnodePlane
    // instances.
    //
    std::vector<WireCell::IAnodePlane::pointer> m_anodes;

    // field_response - name of an IFieldResponse.
    double m_speed{0}, m_origin{0};

    // simchan_label - name by which to save results to art::Event
    std::string m_simchan_label;
    // sparse - whether to save all channels or only those with
    // nonzero depo charge. Default is true.
    bool m_sparse;

    // tick - the sample time over which to integrate depo flux
    // into time bins.
    //
    // window_start - the start of an acceptance window in NOMINAL
    // TIME.
    //
    // window_duration - the duration of the acceptance window.
    WireCell::Binning m_tbins;

    // nsigma - number of sigma at which to truncate a depo Gaussian
    double m_nsigma{3.0};

    // referenced_time - An arbitrary time that is SUBTRACTED from
    // all NOMINAL TIMES prior to setting IDE tdc.
    //
    // time_offsets - An arbitrary time that is ADDED to NOMINAL
    // TIMES for each plane prior to setting IDE.
    std::vector<int> m_tick_offsets;

    // energy - The nominal, constant "energy" of the step
    // associated with the depo.  If zero, the depo will be
    // queried for the energy.  Default is zero.
    double m_energy{0};

    // smear_tran - Extra smearing applied to the depo in the
    // transverse direction of each plane.  This given in units of
    // pitch (ie, smear_tran=2.0 would additionally smear over 2
    // pitch distances) and it is added in quadrature with the
    // Gaussian sigma of the depo.  If zero (default) no extra
    // smearing is applied.
    std::vector<double> m_smear_tran{0};

    // smear_long - Extra smearing applied to the depo in the
    // longitudinal direction.  This given in units of tick (ie,
    // smear_long=2.0 would additionally smear across two ticks)
    // and added in quadrature with the Gaussian sigma of the
    // depo.  If zero (default) no extra smearing is applied.
    double m_smear_long{0};

    // sed_label - The art::Event label for a vector of
    // SimEnergyDeposit.  If given, the IDepo::id() is interpreted
    // as an index into this vector and both trackID and
    // origTrackId are transfered to the IDE.  If empty (default),
    // the IDepo::id() is set directly as the IDE trackID and no
    // origTrackID is set.
    std::string m_sed_label;

    std::string m_debug_file{""};

    // A queue of depos from WCT side
    std::vector<WireCell::IDepo::pointer> m_depos;

    WireCell::IAnodeFace::pointer find_face(const WireCell::IDepo::pointer& depo);
  };

}

#endif
// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
