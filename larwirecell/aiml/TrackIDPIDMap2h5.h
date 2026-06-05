#ifndef WIRECELL_AIML_TRACKIDPIDMAP2H5
#define WIRECELL_AIML_TRACKIDPIDMAP2H5

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

#include <array>
#include <hdf5.h>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace WireCell::AIML {
  class TrackIDPIDMap2h5 : public Aux::Logger,
                           public wcls::IArtEventVisitor,
                           public IFrameFilter,
                           public IConfigurable {
  public:
    TrackIDPIDMap2h5();
    ~TrackIDPIDMap2h5() override;

    // IFrameFilter
    bool operator()(const input_pointer& in, output_pointer& out) override;

    // IConfigurable
    void configure(const WireCell::Configuration& config) override;
    WireCell::Configuration default_configuration() const override;

    // IArtEventVisitor
    void visit(art::Event& event) override;

  private:
    void cache_simchannels(const std::vector<sim::SimChannel>& simchs);
    void populate_trackid_pid_map();
    void clear_cache();
    void ensure_file();
    void write_mapping(int frame_ident);
    void write_mapping_simchnl(int frame_ident);
    void write_mc_json(int frame_ident);

    // JSON helpers (CellTree-compatible)
    bool keep_mc(int tid) const;
    static std::string pdg_name(int pdg);
    bool dump_mc_json_node(int tid, std::ostream& out) const;
    void dump_mc_json(std::ostream& out) const;

    std::string m_simchannel_label;
    std::string m_particle_label;
    std::string m_output_file;
    bool m_save_mc_json{false};
    bool m_save_extended_mcpart{false};  // dump extra MCParticle fields to HDF5

    std::vector<sim::SimChannel> m_simchannels;

    // --- MCParticle-based maps (all stored G4 tracks, regardless of TPC ionization) ---
    // Basic (always dumped):
    std::unordered_map<int, int>   m_trackid_to_pid;        // PDG code
    std::unordered_map<int, int>   m_trackid_to_motherid;   // parent track ID
    std::unordered_map<int, int>   m_trackid_to_process;    // G4 creation process code
    std::unordered_map<int, int>   m_trackid_to_motherpid;  // parent's PDG code
    // Start/end positions [x, y, z, t] in cm from MCParticle trajectory
    std::unordered_map<int, std::array<float,4>> m_trackid_to_start;
    std::unordered_map<int, std::array<float,4>> m_trackid_to_end;
    // Start/end momentum [px, py, pz, E] in GeV from MCParticle
    std::unordered_map<int, std::array<float,4>> m_trackid_to_startmom;
    std::unordered_map<int, std::array<float,4>> m_trackid_to_endmom;
    // Daughter track IDs
    std::unordered_map<int, std::vector<int>> m_trackid_to_daughters;
    // Trajectory points [x, y, z] for each track
    std::unordered_map<int, std::vector<std::array<float,3>>> m_trackid_to_traj;

    // Extended (only when save_extended_mcpart is true):
    std::unordered_map<int, int>   m_trackid_to_endprocess;  // G4 termination process code
    std::unordered_map<int, int>   m_trackid_to_status;      // StatusCode (1=tracked)
    std::unordered_map<int, float> m_trackid_to_mass;        // PDG mass [GeV]
    std::unordered_map<int, int>   m_trackid_to_ndaughters;  // number of daughters
    std::unordered_map<int, int>   m_trackid_to_ntrajpts;    // number of trajectory points

    // --- SimChannel-based maps (only tracks that ionized in TPC, via ParticleInventoryService) ---
    std::unordered_map<int, int>   m_simchnl_trackid_to_pid;       // PDG code from pi_serv
    std::unordered_map<int, int>   m_simchnl_trackid_to_motherid;  // parent track ID from pi_serv
    std::unordered_map<int, int>   m_simchnl_trackid_to_process;   // G4 process code from pi_serv
    std::unordered_map<int, int>   m_simchnl_trackid_to_motherpid; // parent's PDG code from pi_serv
    std::unordered_map<int, float> m_simchnl_trackid_to_energy;    // summed deposited energy [MeV]

    // G4 process name → integer code (following CellTree convention)
    static const std::map<std::string, int> s_process_map;

    hid_t m_file;
  };
}

#endif
