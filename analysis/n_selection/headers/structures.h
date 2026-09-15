#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <limits>
#include "selection_input.h"
#include "TVector3.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this should be SRProxy.h, but there is an include error for that now. Alternatively, you can use SetBranchStatus function in TreeLoader, but it does not work for the common branch (to do)


// ======================
// CONSTANTS AND CONFIGURATION
// ======================

enum class TPC {
    TPC0 = 0,
    TPC1 = 1,
    TPC2 = 2,
    TPC3 = 3,
    TPC4 = 4,
    TPC5 = 5,
    TPC6 = 6,
    TPC7 = 7,
    Outside = -1
};


struct TPCMapping {
    int module;
    int io_group;
    TPC tpc;
    std::array<double, 2> x_bound;
    std::array<double, 2> y_bound;
    std::array<double, 2> z_bound;
};


constexpr std::array<TPCMapping, 8> tpc_map = {{
    {0, 1, TPC::TPC0, {63.931, 33.49 }, {-61.85, 61.85}, {2.68, 64.3163}},
    {0, 2, TPC::TPC1, {33.49 , 3.07}, {-61.85, 61.85}, {2.68, 64.3163}},
    {1, 3, TPC::TPC2, {63.931,33.49 }, {-61.85, 61.85}, {-64.3163, -2.68}},
    {1, 4, TPC::TPC3, {33.49 , 3.07}, {-61.85, 61.85}, {-64.3163, -2.68}},
    {2, 5, TPC::TPC4, {-3.07,-33.49}, {-61.85, 61.85}, {2.68, 64.3163}},
    {2, 6, TPC::TPC5, {-33.49 ,-63.931}, {-61.85, 61.85}, {2.68, 64.3163}},
    {3, 7, TPC::TPC6, {-3.07, -33.49 }, {-61.85, 61.85}, {-64.3163, -2.68}},
    {3, 8, TPC::TPC7, {-33.49,-63.931}, {-61.85, 61.85}, {-64.3163, -2.68}}
}}; 

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

struct PartBestMatch{
    bool valid = false;
    int interaction_idx = -1;
    int particle_idx = -1;
    int type = -1;
    int pdg = 0;
    float overlap = std::numeric_limits<float>::quiet_NaN();
    float energy = std::numeric_limits<float>::quiet_NaN();
    int parent = -1;
    int direct_parent_pdg = 0;
    float time = std::numeric_limits<float>::quiet_NaN();
    float length = std::numeric_limits<float>::quiet_NaN();
    float distance = std::numeric_limits<float>::quiet_NaN();
    caf::SRVector3D start;
    caf::SRVector3D end;
    bool neutron_induced = false;
    // CAF TrueParticleID::PartType for the direct neutron parent.
    int neutron_parent_type = -1;
    int neutron_parent_idx = -1;
    int neutron_parent_g4id = -1;
    bool in_signal = false;
};

struct RecoProtonInfo {
    int reco_pdg;
    int reco_type;
    float reco_energy;
    int reco_vtx;
    caf::SRVector3D reco_start;
    caf::SRVector3D reco_end;
    float reco_len;
    float reco_dist;
    bool reco_different_tpc = false;
    PartBestMatch bm;
    // Input metadata
    std::string input_file;
    std::int64_t input_event = -1;
    int input_vtx_t;
    bool input_matched;
    bool input_has_truth_match = false;
    int input_true_track_multiplicity = -1;
    int input_reco_track_multiplicity = -1;
    int input_true_primary_neutron_count = -1;
    int input_true_secondary_neutron_count = -1;
    float selected_truth_vertex_t0 = std::numeric_limits<float>::quiet_NaN();
    float matched_truth_vertex_t0 = std::numeric_limits<float>::quiet_NaN();
    float matched_particle_dt = std::numeric_limits<float>::quiet_NaN();
    bool same_truth_interaction = false;
    bool coincidence;
};

// One entry in the truth-centric denominator. Only true secondary protons
// with a direct neutron parent and a start point in the fiducial volume are
// included. The end point is deliberately recorded rather than selected on.
struct TrueNeutronProtonInfo {
    std::string input_file;
    std::int64_t input_event = -1;
    int truth_vtx = -1;
    int truth_particle_idx = -1;
    int truth_type = caf::TrueParticleID::kSecondary;
    int truth_pdg = ParticleCode::proton;
    float truth_energy = std::numeric_limits<float>::quiet_NaN();
    float truth_length = std::numeric_limits<float>::quiet_NaN();
    float truth_distance = std::numeric_limits<float>::quiet_NaN();
    float truth_time = std::numeric_limits<float>::quiet_NaN();
    float truth_vertex_time = std::numeric_limits<float>::quiet_NaN();
    float truth_particle_dt = std::numeric_limits<float>::quiet_NaN();
    caf::SRVector3D truth_start;
    caf::SRVector3D truth_end;
    bool truth_end_in_fv = false;
    int truth_parent_g4id = -1;
    int direct_parent_pdg = 0;
    int neutron_parent_type = -1;
    int neutron_parent_idx = -1;
    int neutron_parent_g4id = -1;

    // Event-wide reverse truth matching. "Any overlap" means that at least
    // one reco particle names this truth particle with finite overlap > 0.
    // A reco match additionally requires this truth particle to be that reco
    // particle's largest positive truth overlap.
    bool has_any_reco_overlap = false;
    float max_any_reco_overlap = std::numeric_limits<float>::quiet_NaN();
    bool has_reco_match = false;
    int n_reco_matches = 0;
    int best_reco_vtx = -1;
    int best_reco_particle_idx = -1;
    int best_reco_pdg = 0;
    bool best_reco_primary = false;
    float best_reco_overlap = std::numeric_limits<float>::quiet_NaN();
    bool best_reco_in_selected_interaction = false;
    bool best_reco_start_in_fv = false;
    bool best_reco_end_in_fv = false;
    bool best_reco_passes_select_n = false;
    bool has_selected_candidate = false;
};


struct Counters {
    int secp_reco = 0;  // all reco protons
    int matched_secp = 0; //store protons correctly reconstructed
    int matched_secp_pdg = 0;
    int matched_secp_type = 0;
    int np_reco = 0; // matched neutron-induced proton in the selected truth interaction
    int missed_sec = 0; // potential secondary protons
    int light_candidate = 0; // light candidate protons
};


// struct TrueVertexSelection_np {
//     std::string file_name;
//     int tree_entry;
//     int event;
//     int vtx;
//     bool insignal = false;
//     float neutron_t0 = -1;
//     float proton_t0 = -1;
//     float proton_energy = -1;
//     float proton_length = -1;
//     float proton_distance = -1;
//     // float best_match_length=-1;
//     // float best_match_energy=-1;
//     // float best_match_distance=-1;
//     TPC proton_tpc = TPC::Outside;
//     TPC vertex_tpc = TPC::Outside;
//     bool different_tpc = false;
// };


// struct ProtonsAnalysis {
//     std::vector<RecoProtonInfo> protons;
//     std::vector<TrueVertexSelection_np> true_protons;  // <-- Fixed: no longer a reference

//     // Reco histograms
//     TH1D* h_proton_energy = nullptr;
//     TH1D* h_proton_length = nullptr;
//     TH1D* h_proton_distance = nullptr;

//     // Best match histograms
//     TH1D* h_match_energy = nullptr;
//     TH1D* h_match_length = nullptr;
//     TH1D* h_match_distance = nullptr;

//     // // All True Protons histograms
//     TH1D* h_fulltrue_energy = nullptr;
//     TH1D* h_fulltrue_length = nullptr;
//     TH1D* h_fulltrue_distance = nullptr;

//     TH1D* h_allproton_energy = nullptr;
//     TH1D* h_allproton_length = nullptr;
//     TH1D* h_allproton_distance = nullptr;

//     // Best match histograms
//     TH1D* h_allmatch_energy = nullptr;
//     TH1D* h_allmatch_length = nullptr;
//     TH1D* h_allmatch_distance = nullptr;


//     // Constructor
//     ProtonsAnalysis() {
//         h_proton_energy = new TH1D("h_proton_energy", "", 20, 0, 100);
//         h_proton_length = new TH1D("h_proton_length", "", 80, 0, 100);
//         h_proton_distance = new TH1D("h_proton_distance", "", 80, 0, 100);

//         h_match_energy = new TH1D("h_match_energy", "", 20, 0, 100);
//         h_match_length = new TH1D("h_match_length", "", 80, 0, 100);
//         h_match_distance = new TH1D("h_match_distance", "", 80, 0, 100);

//         h_fulltrue_energy = new TH1D("h_fulltrue_energy", "", 20, 0, 100);
//         h_fulltrue_length = new TH1D("h_fulltrue_length", "", 80, 0, 100);
//         h_fulltrue_distance = new TH1D("h_fulltrue_distance", "", 80, 0, 100);

//         h_allproton_energy = new TH1D("h_allproton_energy", "", 20, 0, 100);
//         h_allproton_length = new TH1D("h_allproton_length", "", 80, 0, 100);
//         h_allproton_distance = new TH1D("h_allproton_distance", "", 80, 0, 100);

//         h_allmatch_energy = new TH1D("h_allmatch_energy", "", 20, 0, 100);
//         h_allmatch_length = new TH1D("h_allmatch_length", "", 80, 0, 100);
//         h_allmatch_distance = new TH1D("h_allmatch_distance", "", 80, 0, 100);
        
//     }

//     ~ProtonsAnalysis() {
//         delete h_proton_energy;
//         delete h_proton_length;
//         delete h_proton_distance;
//         delete h_match_energy;
//         delete h_match_length;
//         delete h_match_distance;
//         delete h_fulltrue_energy;
//         delete h_fulltrue_length;
//         delete h_fulltrue_distance;
//         //all reco protons
//         delete h_allproton_energy;
//         delete h_allproton_length;
//         delete h_allproton_distance;
//         delete h_allmatch_energy;
//         delete h_allmatch_length;
//         delete h_allmatch_distance;
 
//     }

//     void Fill(
//         const RecoProtonInfo& proton, 
//         const TrueVertexSelection_np& true_proton
//         ){
        
//         // Only store real reco protons (non-dummy)
//         if (proton.reco_pdg != -1 && proton.reco_type != -1) {
//             protons.push_back(proton);
//         }

//         // Only store real true protons
//         if (true_proton.proton_energy > 0 && 
//             true_proton.proton_length > 0 && 
//             true_proton.proton_distance > 0) 
//         {
//             true_protons.push_back(true_proton);
//             h_fulltrue_energy->Fill(true_proton.proton_energy);
//             h_fulltrue_length->Fill(true_proton.proton_length);
//             h_fulltrue_distance->Fill(true_proton.proton_distance);
//         }

//         // Histogram reco and match info only if valid and coincident neutron-induced
//         if (proton.is_neutron_induced && proton.is_coincident) {
//             h_proton_energy->Fill(proton.reco_energy);
//             h_proton_length->Fill(proton.reco_len);
//             h_proton_distance->Fill(proton.reco_dist);

//             h_match_energy->Fill(proton.best_match_energy);
//             h_match_length->Fill(proton.best_match_length);
//             h_match_distance->Fill(proton.best_match_distance);
//         }

//         if (proton.reco_pdg == ParticleCode::proton) {
//             h_allproton_energy->Fill(proton.reco_energy);
//             h_allproton_length->Fill(proton.reco_len);
//             h_allproton_distance->Fill(proton.reco_dist);

//             h_allmatch_energy->Fill(proton.best_match_energy);
//             h_allmatch_length->Fill(proton.best_match_length);
//             h_allmatch_distance->Fill(proton.best_match_distance);
//         }
//     }

//     void PrintUniqueMatchedPDGs() const {
//         std::map<int, int> pdg_counts;

//         for (const auto& proton : protons) {
//             if (proton.best_match_pdg !=-1) {  // Ignore dummy or unassigned
//                 pdg_counts[proton.best_match_pdg]++;
//             }
//         }

//         std::cout << "best_match_pdg counts for reco secondary protons:" << std::endl;
//         for (const auto& [pdg, count] : pdg_counts) {
//             std::cout << "PDG: " << pdg << " → Count: " << count << std::endl;
//         }

//     }

//     // Saves histograms to file
//     void Save(const std::string& output_path) {
//         TFile output(output_path.c_str(), "RECREATE");
        
//         // Save histograms
//         h_proton_energy->Write();
//         h_proton_length->Write();
//         h_proton_distance->Write();

        
//         // Create and save a tree with proton information
//         TTree proton_tree("proton_tree", "Proton Information");
        
//         // Variables to store in the tree
//         int pdg, type, vtx_idx, match_pdg, match_type;
//         float energy, length, distance;
//         bool neutron_induced, coincident, proton_t0;
        
//         // Set up branches
//         proton_tree.Branch("pdg", &pdg);
//         proton_tree.Branch("type", &type);
//         proton_tree.Branch("vtx_idx", &vtx_idx);
//         proton_tree.Branch("energy", &energy);
//         proton_tree.Branch("length", &length);
//         proton_tree.Branch("distance", &distance);
//         proton_tree.Branch("match_pdg", &match_pdg);
//         proton_tree.Branch("match_type", &match_type);
//         proton_tree.Branch("neutron_induced", &neutron_induced);
//         proton_tree.Branch("coincident", &coincident);
//         proton_tree.Branch("proton_t0", &proton_t0);
        
//         // Fill the tree
//         for (const auto& proton : protons) {
//             pdg = proton.reco_pdg;
//             type = proton.reco_type;
//             vtx_idx = proton.reco_vtx;
//             energy = proton.reco_energy;
//             length = proton.reco_len;
//             distance = proton.reco_dist;
//             match_pdg = proton.best_match_pdg;
//             match_type = proton.best_match_type;
//             neutron_induced = proton.is_neutron_induced;
//             coincident = proton.is_coincident;
            
//             proton_tree.Fill();
//         }
        
//         proton_tree.Write();
//         output.Close();
//     }

//     // Generates a PDF with the histograms
//     void SaveToPDF(const std::string& output_pdf) {
//         TCanvas* canvas = new TCanvas("canvas", "Reco vs Truth Proton Comparison", 800, 600);
//         canvas->Print((output_pdf + "[").c_str());  // Open multi-page PDF

//         int kRecoColor = TColor::GetColor("#67b8bd");  // Reco: teal
//         int kBMColor = TColor::GetColor("#ba75ff");    // Best match: purple
//         int kTruthColor = TColor::GetColor("#224a8e"); // Truth: dark blue

//         auto draw_comparison = [&](
//             TH1D* h_reco, TH1D* h_BM, TH1D* h_truth,
//             const std::string& title, const std::string& xaxis,
//             const std::string& filename,
//             double xmin = -1, double xmax = -1){
//             canvas->Clear();
//             canvas->cd();
//             gPad->SetRightMargin(0.08);
//             gPad->SetTopMargin(0.08);
//             gPad->SetLeftMargin(0.13);
//             gPad->SetBottomMargin(0.12);
//             gPad->SetLogy();

//             if (xmin >= 0 && xmax > xmin) {
//                 h_reco->GetXaxis()->SetRangeUser(xmin, xmax);
//                 h_BM->GetXaxis()->SetRangeUser(xmin, xmax);
//                 h_truth->GetXaxis()->SetRangeUser(xmin, xmax);
//             }

//             h_reco->SetLineColor(kRecoColor);
//             h_reco->SetLineWidth(2);
//             h_truth->SetLineColor(kTruthColor);
//             h_truth->SetLineWidth(2);
//             h_BM->SetLineColor(kBMColor);
//             h_BM->SetLineWidth(2);

//             h_reco->SetTitle((title + ";" + xaxis + ";Entries").c_str());
//             h_reco->SetTitleOffset(2, ""); 
//             h_reco->SetTitleSize(0.055, "");        // Title font size
//             h_reco->SetLabelSize(0.045, "XY");      // Tick label font size
//             h_reco->SetTitleSize(0.05, "XY");       // Axis title font size
//             h_reco->SetTitleOffset(1.6, "Y");       // Y-axis title spacing
//             h_reco->SetTitleOffset(1.2, "X");       // X-axis title spacing


//             double max_val = std::max({h_reco->GetMaximum(), h_truth->GetMaximum(), h_BM->GetMaximum()});
//             h_reco->SetMaximum(1.2 * max_val);

//             h_reco->Draw("HIST");
//             h_BM->Draw("HIST SAME");
//             h_truth->Draw("HIST SAME");

//             dunestyle::WIP(static_cast<ETextAlign>(kHAlignLeft + kVAlignTop), false);


//             // TLegend* leg = new TLegend(0.55, 0.72, 0.85, 0.88);
//             // leg->AddEntry(h_reco, "SPINE Reconstruction", "l");
//             // leg->AddEntry(h_BM, "Truth Best Match", "l");
//             // leg->AddEntry(h_truth, "Full Truth Sample", "l");
//             // leg->SetBorderSize(0);
//             // leg->SetFillStyle(0);
//             // leg->SetTextSize(0.035);
//             // leg->Draw();

//             canvas->Print(output_pdf.c_str());
//         };

        

//         // np Comparison plots
//         draw_comparison(h_proton_energy, h_match_energy, h_fulltrue_energy, "Neutron Induced Proton Energy - RHC #scale[1.0]{1 #times10^{19}} POT", "Energy [MeV]", "energy",0,100);
//         draw_comparison(h_proton_length, h_match_length, h_fulltrue_length, "Neutron Induced Proton Track Length - RHC #scale[1.0]{1 #times10^{19}} POT", "Length [cm]", "length", 0,40);
//         draw_comparison(h_proton_distance, h_match_distance, h_fulltrue_distance, "Neutron Induced Proton-Vtx Distance - RHC #scale[1.0]{1 #times10^{19}} POT", "Distance [cm]", "distance", 0,25);

//         // all reco protons comparison plots
//         draw_comparison(h_allproton_energy, h_allmatch_energy, h_fulltrue_energy, "Reconstructed Secondary Proton Energy - RHC #scale[1.0]{1#times10^{19}} POT", "Energy [MeV]", "energy",0,100);
//         draw_comparison(h_allproton_length, h_allmatch_length, h_fulltrue_length, "Reconstructed Secondary Proton Track Length - RHC #scale[1.0]{1#times10^{19}} POT", "Length [cm]", "length",0,30);
//         draw_comparison(h_allproton_distance, h_allmatch_distance, h_fulltrue_distance, "Reconstructed Secondary Proton-Vtx Distance - RHC #scale[1.0]{1#times10^{19}} POT", "Distance [cm]", "distance",0,35);

//         canvas->Print((output_pdf + "]").c_str());  // Close PDF

//         delete canvas;
//     }

// };





// ======================
