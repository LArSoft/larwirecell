#ifndef WIRECELL_AIML_TRUTH2H5
#define WIRECELL_AIML_TRUTH2H5

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include <string>

#include "hdf5.h"

namespace WireCell::AIML {
  class Truth2h5 : public Aux::Logger,
                   public wcls::IArtEventVisitor,
                   public IFrameFilter,
                   public IConfigurable {
  public:
    Truth2h5();
    ~Truth2h5() override;

    struct NeutrinoInfo {
      bool valid{false};
      int pdg{0};
      int ccnc{0};
      int intType{0};
      double energy{0.0};
      double vx{0.0};
      double vy{0.0};
      double vz{0.0};
    };

    // IFrameFilter
    bool operator()(const input_pointer& in, output_pointer& out) override;

    // IConfigurable
    void configure(const WireCell::Configuration& config) override;
    WireCell::Configuration default_configuration() const override;

    // IArtEventVisitor
    void visit(art::Event& event) override;

  private:
    void reset();
    void ensure_file();
    void write_event(const IFrame& frame);

    std::string m_mctruth_label;
    std::string m_output_file;
    NeutrinoInfo m_info;
    hid_t m_file;
  };
}

#endif
