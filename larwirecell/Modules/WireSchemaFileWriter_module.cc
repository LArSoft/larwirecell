/** WireSchemaFileWriter (Option 1)
 *
 *  An art analyzer that writes a Wire-Cell Toolkit "wires file" (.json, .json.bz2
 *  or .json.gz) derived directly from the LArSoft geometry + wire-readout
 *  (channel map) services configured for the job.
 *
 *  Run it as a one-shot precursor over a single (even empty) event.  Because the
 *  wires come from the SAME channel map the "real" jobs use, the file is
 *  guaranteed consistent with that LArSoft channel numbering.  Record the
 *  dunecore/geometry provenance together with the emitted file.
 *
 *  FHiCL:
 *    physics.analyzers.wireswriter: {
 *      module_type: WireSchemaFileWriter
 *      output:     "dunevd10kt-wires-fromgeo.json.bz2"  # required
 *      grouping:   "tpcset"      # tpcset (default) | per-tpc | pair-tpc
 *      correction: 4             # 1=load 2=order 3=direction 4=pitch (default)
 *      validate:   true          # throw if the result violates wire-schema conventions
 *    }
 *  The LArSoft channel map + sorter are whatever services_dune.fcl selects for
 *  the geometry in use; this module simply reads the resulting service.
 */

#include "larwirecell/Components/WireSchemaBuilder.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "WireCellUtil/WireSchema.h"

#include <string>

namespace wcls {

    class WireSchemaFileWriter : public art::EDAnalyzer {
      public:
        explicit WireSchemaFileWriter(fhicl::ParameterSet const& pset);

      private:
        void beginJob() override;
        void analyze(const art::Event&) override {}  // nothing per-event

        std::string fOutput;
        std::string fGrouping;
        int fCorrection;
        bool fValidate;
        bool fDone{false};
    };

    WireSchemaFileWriter::WireSchemaFileWriter(fhicl::ParameterSet const& pset)
        : art::EDAnalyzer(pset)
        , fOutput(pset.get<std::string>("output"))
        , fGrouping(pset.get<std::string>("grouping", "tpcset"))
        , fCorrection(pset.get<int>("correction",
                          static_cast<int>(WireCell::WireSchema::Correction::pitch)))
        , fValidate(pset.get<bool>("validate", true))
    {
    }

    void WireSchemaFileWriter::beginJob()
    {
        if (fDone) return;
        fDone = true;

        auto const* geom = art::ServiceHandle<geo::Geometry const>().get();
        auto const& wires = art::ServiceHandle<geo::WireReadout const>()->Get();

        GroupingPolicy policy;
        policy.grouping = grouping_from_string(fGrouping);
        policy.correction = static_cast<WireCell::WireSchema::Correction>(fCorrection);
        policy.validate = fValidate;

        mf::LogInfo("WireSchemaFileWriter")
            << "building WCT wire schema from detector '" << geom->DetectorName()
            << "' (" << wires.Nchannels() << " channels), grouping=" << fGrouping;

        auto store = build_wire_schema(*geom, wires, policy);

        WireCell::WireSchema::dump(fOutput.c_str(), store);

        mf::LogInfo("WireSchemaFileWriter")
            << "wrote " << fOutput << ": " << store.anodes().size() << " anodes, "
            << store.faces().size() << " faces, " << store.planes().size()
            << " planes, " << store.wires().size() << " wires";
    }

}  // namespace wcls

DEFINE_ART_MODULE(wcls::WireSchemaFileWriter)
