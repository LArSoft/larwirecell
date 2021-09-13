// Dump TPC / Wire Geometry and mapping
// Chao Zhang (chao@bnl.gov) 2/7/2018
// adapted by Wenqiang Gu (wgu@bnl.gov) 8/30/2020

#ifndef CTreeGeometry_module
#define CTreeGeometry_module

// LArSoft includes
//#include "lardata/Utilities/DetectorProperties.h"
#include "lardata/Utilities/GeometryUtilities.h"
// #include "Utilities/LArProperties.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"
// #include "nusimdata/SimulationBase/MCParticle.h"
// #include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"


// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
// #include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"


// C++ Includes
#include <map>
#include <vector>
// #include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
// #include <string>
// #include <sstream>
// #include <cmath>

// #ifdef __MAKECINT__
// #pragma link C++ class vector<vector<int> >+;
// #pragma link C++ class vector<vector<float> >+;
// #endif

using namespace std;

namespace {

class CTreeGeometry : public art::EDAnalyzer {
public:

    explicit CTreeGeometry(fhicl::ParameterSet const& pset);
    virtual ~CTreeGeometry();

    void beginJob();
    void endJob();
    void analyze(const art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);
    void saveChannelWireMap();
    void printGeometry();


private:

    // the parameters we'll read from the .fcl
    bool fSaveChannelWireMap;

    art::ServiceHandle<geo::Geometry> fGeom;

    // Geometry Tree Leafs
    int fNcryostats;
    int fNTPC;
    vector<float> fTPC_x;  // TPC length in x
    vector<float> fTPC_y;  // TPC length in y
    vector<float> fTPC_z;  // TPC length in z
    // int fNplanes; // unused
    vector<int> fPlane_type;  // plane type: 0 == induction, 1 == collection
    vector<int> fPlane_view;  // wire orientation: 0 == U, 1 == V, 2 == X
    vector<double> fPlane_wirepitch;  // wire pitch of each plane
    vector<double> fPlane_wireangle;  // wire angle (to vertical) of each plane
    vector<int> fPlane_wires;  // number of wires in each plane
    int fNchannels;
    vector<int> channel_starts; // vector os channel starts on each plane
    vector<int> channel_ends; // vector os channel starts on each plane
    //int fNOpDets; // unused

    // Event Tree Leafs
    //int fEvent; // unused
    //int fRun; // unused
    //int fSubRun; // unused

   }; // class CTreeGeometry


//-----------------------------------------------------------------------
CTreeGeometry::CTreeGeometry(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    reconfigure(parameterSet);
}


//-----------------------------------------------------------------------
CTreeGeometry::~CTreeGeometry()
{
}


//-----------------------------------------------------------------------
void CTreeGeometry::reconfigure(fhicl::ParameterSet const& p){
    fSaveChannelWireMap = p.get< bool >("saveChannelWireMap");
}


//-----------------------------------------------------------------------
void CTreeGeometry::beginJob()
{

    fNcryostats = fGeom->Ncryostats();  // 1

    fNTPC = fGeom->NTPC();
    for (int i=0; i<fNTPC; i++) {
        fTPC_x.push_back(fGeom->DetHalfWidth(i)*2);
        fTPC_y.push_back(fGeom->DetHalfHeight(i)*2);
        fTPC_z.push_back(fGeom->DetLength(i));
    }

   //  fNplanes = fGeom->Nplanes();
   // cout << "#planes: " << fNplanes << endl;
    // for (int i=0; i<fNplanes; i++) {
    //     fPlane_type.push_back(fGeom->SignalType(geo::PlaneID(0, 0, i)));
    //     fPlane_view.push_back(fGeom->Plane(i).View());
    //     // fPlane_wirepitch[i] = fGeom->WirePitch(fPlane_view[i]);  // this doesn't seem to return the correct value!
    //     fPlane_wirepitch.push_back(fGeom->WirePitch(fPlane_view[i], 1, 0));  // this doesn't seem to return the correct value);
    //     fPlane_wireangle.push_back(fGeom->WireAngleToVertical(fGeom->Plane(i).View()));
    //     fPlane_wires.push_back(fGeom->Nwires(i));
    // }

    fNchannels = fGeom->Nchannels();
    // cout << "#channel: " << fNchannels << endl;

    // Save Channel Map to text file.
    if (fSaveChannelWireMap) {
        saveChannelWireMap();
    }

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
    for (int i=0; i<fNchannels; i++) {
        std::vector<geo::WireID> wireids = fGeom->ChannelToWire(i);
        int nWires = wireids.size();
        for (int j=0; j<nWires; j++) {
            geo::WireID wid = wireids.at(j);
            int cstat = wid.Cryostat;
            int tpc = wid.TPC;
            int plane = wid.Plane;
            int wire = wid.Wire;

            int plane_id = plane + tpc*10 + cstat*100;
            if (plane_id != current_plane) {
                current_plane = plane_id;
                channel_starts.push_back(i);
                channel_ends.push_back(i);
            }
            else {
                channel_ends[channel_ends.size()-1] = i;
            }

            fGeom->WireEndPoints(cstat, tpc, plane, wire, xyzStart, xyzEnd);

            out << i << "\t" << cstat*2+tpc << "\t" << plane << "\t" << wire << "\t";
            for (int i=0; i<3; i++) {
                out << xyzStart[i] << "\t";
            }
            for (int i=0; i<3; i++) {
                out << xyzEnd[i] << "\t";
            }
            out << "\n";
        }
    }
    out.close();

}


//-----------------------------------------------------------------------
void CTreeGeometry::endJob()
{
}


//-----------------------------------------------------------------------
void CTreeGeometry::printGeometry()
{
    cout << "Detector Name: " << fGeom->DetectorName() << endl;
    cout << "GDML file: " << fGeom->GDMLFile() << endl;
    // cout << "fNTPC: " << fNTPC << endl;
    // for (int i=0; i<fNTPC; i++) {
    //     cout << "\tTPC " << i << ": " << fTPC_x[i] << ", " << fTPC_y[i] << ", " << fTPC_z[i] << endl;
    // }
    cout << "TPC (Active) Locations: " << endl;
    for (geo::TPCGeo const& TPC: fGeom->IterateTPCs()) {
        // get center in world coordinates
        double origin[3] = {0.};
        double center[3] = {0.};
        TPC.LocalToWorld(origin, center);
        double tpcDim[3] = {TPC.ActiveHalfWidth(), TPC.ActiveHalfHeight(), 0.5*TPC.ActiveLength() };
        double xmin = center[0] - tpcDim[0];
        double xmax = center[0] + tpcDim[0];
        double ymin = center[1] - tpcDim[1];
        double ymax = center[1] + tpcDim[1];
        double zmin = center[2] - tpcDim[2];
        double zmax = center[2] + tpcDim[2];
        cout << "\t[" << xmin << ", " << xmax << ", " << ymin << ", " << ymax
             << ", " << zmin << ", " << zmax << "]" << endl;
    } // for all TPC

    int size = channel_starts.size();
    cout << size << " planes: first channels: ";
    for (int i=0; i<size ; i++) {
        cout << channel_starts[i] << ", ";
    }
    cout << endl;
    size = channel_ends.size();
    cout << size << " planes: last channels: ";
    for (int i=0; i<size ; i++) {
        cout << channel_ends[i] << ", ";
    }
    cout << endl;
    // for (geo::TPCGeo const& TPC: fGeom->IterateTPCs()) {
    //     for (geo::PlaneGeo const& plane: TPC.IteratePlanes()) {
    //         cout << plane.ID() <<
    //         " channels: " fGeom->PlaneWireToChannel(plane.FirstWire().node())<< endl;
    //      }

    // }

    // cout << "fNplanes: " << fNplanes << endl;
    // for (int i=0; i<fNplanes; i++) {
    //     cout
    //         << "\tplane " << i
    //         << "( type: " << fPlane_type[i]
    //         << ", view: " << fPlane_view[i]
    //         << ", wirepitch: " << fPlane_wirepitch[i]
    //         << ", wire angle: " << fPlane_wireangle[i]
    //         << ", wires: " << fPlane_wires[i]
    //         << ")" << endl;
    // }
    cout << "fNchannels: " << fGeom->Nchannels() << endl;
    cout << "fNOpDet: " << fGeom->NOpDets() << endl;
    cout << "fAuxDetectors: " << fGeom->NAuxDets() << endl;
    cout << endl;
}

//-----------------------------------------------------------------------
void CTreeGeometry::analyze( const art::Event& event )
{

    // fEvent  = event.id().event();
    // fRun    = event.run();
    // fSubRun = event.subRun();

    // printEvent();

}

DEFINE_ART_MODULE(CTreeGeometry)

} // namespace

#endif // CTreeGeometry_module
