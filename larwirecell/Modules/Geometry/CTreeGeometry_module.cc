// Dump TPC / Wire Geometry and mapping
// Chao Zhang (chao@bnl.gov) 2/7/2018
// adapted by Wenqiang Gu (wgu@bnl.gov) 8/30/2020

// LArSoft includes
#include "larcore/Geometry/AuxDetGeometry.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/fwd.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// C++ Includes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

namespace {

  class CTreeGeometry : public art::EDAnalyzer {
  public:
    explicit CTreeGeometry(fhicl::ParameterSet const& pset);

  private:
    void beginJob() override;
    void analyze(const art::Event& evt) override;

    void saveChannelWireMap();
    void printGeometry();

    // the parameters we'll read from the .fcl
    bool fSaveChannelWireMap;

    art::ServiceHandle<geo::Geometry> fGeom;
    geo::WireReadoutGeom const* fWireReadoutGeom{
      &art::ServiceHandle<geo::WireReadout const>()->Get()};
    geo::AuxDetGeometryCore const* fAuxDetGeom{
      art::ServiceHandle<geo::AuxDetGeometry>()->GetProviderPtr()};

    // Geometry Tree Leafs
    int fNchannels;
    vector<int> channel_starts; // vector os channel starts on each plane
    vector<int> channel_ends;   // vector os channel starts on each plane

  }; // class CTreeGeometry

  //-----------------------------------------------------------------------
  CTreeGeometry::CTreeGeometry(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset), fSaveChannelWireMap{pset.get<bool>("saveChannelWireMap")}
  {}

  //-----------------------------------------------------------------------
  void CTreeGeometry::beginJob()
  {
    fNchannels = fWireReadoutGeom->Nchannels();

    // Save Channel Map to text file.
    if (fSaveChannelWireMap) { saveChannelWireMap(); }

    printGeometry();
  }

  //-----------------------------------------------------------------------
  void CTreeGeometry::saveChannelWireMap()
  {
    ofstream out;
    out.open("ChannelWireGeometry.txt");
    double xyzStart[3];
    double xyzEnd[3];
    out << "# " << fGeom->GDMLFile() << "\n";
    out << "# channel\ttpc\tplane\twire\tsx\tsy\tsz\tex\tey\tez\n";
    int current_plane = 0;
    channel_starts.push_back(0);
    channel_ends.push_back(0);
    for (int i = 0; i < fNchannels; i++) {
      std::vector<geo::WireID> wireids = fWireReadoutGeom->ChannelToWire(i);
      int nWires = wireids.size();
      for (int j = 0; j < nWires; j++) {
        geo::WireID wid = wireids.at(j);
        int cstat = wid.Cryostat;
        int tpc = wid.TPC;
        int plane = wid.Plane;
        int wire = wid.Wire;

        int plane_id = plane + tpc * 10 + cstat * 100;
        if (plane_id != current_plane) {
          current_plane = plane_id;
          channel_starts.push_back(i);
          channel_ends.push_back(i);
        }
        else {
          channel_ends[channel_ends.size() - 1] = i;
        }

        fWireReadoutGeom->WireEndPoints(wid, xyzStart, xyzEnd);

        out << i << "\t" << cstat * 2 + tpc << "\t" << plane << "\t" << wire << "\t";
        for (int i = 0; i < 3; i++) {
          out << xyzStart[i] << "\t";
        }
        for (int i = 0; i < 3; i++) {
          out << xyzEnd[i] << "\t";
        }
        out << "\n";
      }
    }
    out.close();
  }

  //-----------------------------------------------------------------------
  void CTreeGeometry::printGeometry()
  {
    cout << "Detector Name: " << fGeom->DetectorName() << endl;
    cout << "GDML file: " << fGeom->GDMLFile() << endl;
    cout << "TPC (Active) Locations: " << endl;
    for (geo::TPCGeo const& TPC : fGeom->Iterate<geo::TPCGeo>()) {
      // get center in world coordinates
      auto const center = TPC.GetCenter();
      double tpcDim[3] = {TPC.ActiveHalfWidth(), TPC.ActiveHalfHeight(), 0.5 * TPC.ActiveLength()};
      double xmin = center.X() - tpcDim[0];
      double xmax = center.X() + tpcDim[0];
      double ymin = center.Y() - tpcDim[1];
      double ymax = center.Y() + tpcDim[1];
      double zmin = center.Z() - tpcDim[2];
      double zmax = center.Z() + tpcDim[2];
      cout << "\t[" << xmin << ", " << xmax << ", " << ymin << ", " << ymax << ", " << zmin << ", "
           << zmax << "]" << endl;
    } // for all TPC

    int size = channel_starts.size();
    cout << size << " planes: first channels: ";
    for (int i = 0; i < size; i++) {
      cout << channel_starts[i] << ", ";
    }
    cout << endl;
    size = channel_ends.size();
    cout << size << " planes: last channels: ";
    for (int i = 0; i < size; i++) {
      cout << channel_ends[i] << ", ";
    }
    cout << endl;

    cout << "fNchannels: " << fNchannels << endl;
    cout << "fNOpDet: " << fGeom->NOpDets() << endl;
    cout << "fAuxDetectors: " << fAuxDetGeom->NAuxDets() << endl;
    cout << endl;
  }

  //-----------------------------------------------------------------------
  void CTreeGeometry::analyze(const art::Event&) {}

  DEFINE_ART_MODULE(CTreeGeometry)

} // namespace
