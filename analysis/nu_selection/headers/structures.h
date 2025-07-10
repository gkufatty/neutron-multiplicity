#pragma once

#include <string>
#include <vector>
#include "TVector3.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this should be SRProxy.h, but there is an include error for that now. Alternatively, you can use SetBranchStatus function in TreeLoader, but it does not work for the common branch (to do)


// ======================
// CONSTANTS AND CONFIGURATION
// ======================


struct ParticleCode {
    // Masses in MeV/c^2
    static constexpr double neutron_mass = 939.5654;
    static constexpr double proton_mass  = 938.2702;
    // PDG Codes
    static constexpr int argon     = 1000180400;
    static constexpr int chlorine  = 1000170360;
    static constexpr int sulfur    = 1000160320;
    static constexpr int muon      = 13;
    static constexpr int neutron   = 2112;
    static constexpr int electron  = 11;
    static constexpr int pi0       = 111;
    static constexpr int pip       = 211;
    static constexpr int eta       = 221;
    static constexpr int proton    = 2212;
    static constexpr int numu      = 14;
    static constexpr int nue       = 12;
    static constexpr int photon    = 22;
};


struct TrueVertexSelection_np {
    std::string file_name;
    int tree_entry;
    int event;
    int vtx;
    bool insignal = false;
    float neutron_t0 = -1;
    float proton_t0 = -1;
    float proton_energy = -1;
    float proton_length = -1;
    float proton_distance = -1;
    // float best_match_length=-1;
    // float best_match_energy=-1;
    // float best_match_distance=-1;
    // TPC proton_tpc = TPC::Outside;
    // TPC vertex_tpc = TPC::Outside;
    // bool different_tpc = false;
};

namespace constants {
    constexpr double TPC_DIST = 3.07;
    constexpr double FV_DIST = 5;
    constexpr double X_BOUND = 63.931;
    constexpr double Y_BOUND = 62.076;
    constexpr double Z_BOUND = 64.3163;
    const TVector3 BEAM_DIR(0., 0., 1.);
}



// ======================
// DATA STRUCTURES
// ======================

//Store info of reco vertices
struct RecoSignal{
    std::string file_name; 
    int event;
    int tree_entry;
    int vtx_r;
    int vtx_t;
    bool insignal;
    float muon_t0;
};

//store best match info add length, energy and distance
struct PartBestMatch{
    int pdg;
    int type;
    int time;
    int interaction_idx;
    int particle_idx;
    TVector3 start;
    TVector3 end;
    int parent;
    bool coincidence;
    float length;
    float energy;
    float distance;
};

struct Counters {
    int nu_true = 0, nu_reco = 0, matched = 0;
    int prim_neutrons = 0; // Number of primary neutrons in true interactions
    int secp_true = 0; // Number of secondary protons in true interactions
    int np_true = 0; // Number of true neutron induced protons
    int n_fv = 0, matched_fv = 0;
    int n_2x2 = 0, matched_2x2 = 0;
    int n_mx2 = 0, matched_mx2 = 0;
    int n_1prim = 0, matched_1prim = 0;
    int n_1sec = 0, matched_1sec = 0;
    int n_less3sec = 0, matched_less3sec = 0;
    int QES = 0;
    int MEC = 0;
    int RES = 0;
    std::vector<RecoSignal> passed_events;  // Store event info that pass cuts
};




struct RecoVertexSelection {
    bool in_fv;
    bool has_muon_exiting;
    bool muon_crossed_mx2;
    int n_prim;
    int n_sec;
    bool is_signal;
    bool all_cuts_passed;
    float muon_t0;         // New
};
