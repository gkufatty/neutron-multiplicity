// ======================
// HEADER SECTION
// ======================

// ROOT dependencies
#include "TTree.h"
#include "TSystem.h"  // For gSystem
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TVector3.h"

#include "TGraph.h"
#include "TLegend.h"
#include <TAxis.h>
#include "TGaxis.h"
#include "TSystem.h"  // For gSystem
#include <array>
#include <string>
#include <map> 
#include <stdexcept>
#include <cmath>  // For TMath functions
#include <utility> // For std::pair
#include <set>

// Standard dependencies
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sstream> 
#include <vector>

// CAF dependencies
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this should be SRProxy.h, but there is an include error for that now. Alternatively, you can use SetBranchStatus function in TreeLoader, but it does not work for the common branch (to do)

#include "DUNEStyle.h"
#include "cuts.h"  // Include the header file for cuts
#include "structures.h"  // Include the header file for structures
#include "utils.h"  // Include the header file for utils


// ======================
// FUNCTION DECLARATIONS
// ======================

// Main workflow
int select_2x2_neutrons(std::string input_file_list, const char* version ="6.2", const char* mode = "RHC");

// ======================
// MAIN FUNCTION
// ======================


int main(int argc, char** argv){

    char cwd[1024];
    if(getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cout << "Current working directory: " << cwd << std::endl;
    }

    if(argc != 3){
        std::cout << "\nUSAGE: " << argv[0] << " <input_caf_file_list> <mode>\n" << std::endl;
        return 1;
    }

    // std::string input_file_list = argv[1];
    
    std::string version = argv[1];
    std::string mode = argv[2];
    std::string input_file_list = "../nu_selection/output_files/MiniRun" + version + "_" + mode + "_QE_like_ixn.txt";
    std::cout << "Using MiniRun" << version << " in " << mode << " mode." << std::endl;
    select_2x2_neutrons(input_file_list, version.c_str(), mode.c_str());
    
    return 0;
}

// ======================
// FUNCTION DEFINITIONS
// ======================


int select_2x2_neutrons(std::string input_file_list, const char* version, const char* mode) {
    std::cout << "Running neutrons selector..." << std::endl;
    Counters stats;
    std::vector<RecoProtonInfo> all_protons;
    std::ifstream caf_list(input_file_list);

    if (!caf_list.is_open()) {
        std::cerr << "Error: File not found -> " << input_file_list << std::endl;
        return 1;
    }

    // Skip header
    std::string header;
    std::getline(caf_list, header);

    std::vector<InputCAFRow> input_rows;
    std::string line;
    int total_lines = 0;

    while (std::getline(caf_list, line)) {
        total_lines++;
        std::stringstream ss(line);
        std::string file_name, event_str, vtx_r_str, vtx_t_str, matched_str;

        std::getline(ss, file_name, ',');
        std::getline(ss, event_str, ',');
        std::getline(ss, vtx_r_str, ',');
        std::getline(ss, vtx_t_str, ',');
        std::getline(ss, matched_str, ',');

        if (file_name.empty() || vtx_r_str.empty()) continue;

        try {
            InputCAFRow row;
            row.file_name = file_name;
            row.event = std::stoi(event_str);
            row.vtx_r = std::stoi(vtx_r_str);
            row.vtx_t = std::stoi(vtx_t_str);
            row.matched_true_signal = (matched_str == "1" || matched_str == "true");
            input_rows.push_back(row);
        } catch (...) {
            std::cerr << "Error parsing line: " << line << std::endl;
        }
    }

    std::cout << "Total number of CAF rows (excluding header): " << total_lines << std::endl;
    std::cout << "Loaded " << input_rows.size() << " rows.\n";

    std::set<std::string> unique_files;
    for (const auto& row : input_rows)
        unique_files.insert(row.file_name);

    TChain* caf_tree = new TChain("cafTree");
    for (const auto& f : unique_files)
        caf_tree->Add(f.c_str());

    std::cout << "Loaded " << unique_files.size() << " unique CAF files into TChain.\n";

    auto sr = new caf::StandardRecord;
    caf_tree->SetBranchAddress("rec", &sr);
    int tree_entries = caf_tree->GetEntries();
    std::cout << "Number of entries in TChain: " << tree_entries << std::endl;

    for (int n = 0; n < tree_entries; n++) {
        caf_tree->GetEntry(n);
        std::string full_path = caf_tree->GetFile()->GetName();

        for (const auto& row : input_rows) {
            if (row.file_name == full_path && row.vtx_r >= 0 &&
                row.vtx_r < static_cast<int>(sr->common.ixn.dlp.size())) {

                const auto& reco = sr->common.ixn.dlp[row.vtx_r];
                std::vector<RecoProtonInfo> protons_vtx = process_protons(row.vtx_r, reco, sr, stats, mode);

                for (auto& p : protons_vtx) {
                    p.input_file = row.file_name;
                    p.input_event = row.event;
                    p.input_vtx_t = row.vtx_t;
                    p.input_matched = row.matched_true_signal;
                    all_protons.push_back(p);
                }
            }
        }
    }

    std::ofstream out("../n_selection/output_files/protons_MiniRun" + std::string(version) + "_" + std::string(mode) + ".txt");
    out << "file,event,rvtx,rtvtx,matched,"
        << "rpdg,rtype,rlen,rdist,rE,rstart,rend,"
        << "t_int_idx,tpart_idx,tpdg,ttype,tlen,tdist,tE,"
        << "tstart,tend,tparent,tt0,coincidence,insignal,ninduced\n";

    for (const auto& p : all_protons) {
        out << p.input_file << ","
            << p.input_event << ","
            << p.reco_vtx << ","
            << p.input_vtx_t << ","
            << p.input_matched << ","
            << p.reco_pdg << ","
            << p.reco_type << ","
            << p.reco_len << ","
            << p.reco_dist << ","
            << p.reco_energy << ","
            << p.reco_start << ","
            << p.reco_end << ","
            << p.bm.interaction_idx << ","
            << p.bm.particle_idx << ","
            << p.bm.pdg << ","
            << p.bm.type << ","
            << p.bm.length << ","
            << p.bm.distance << ","
            << p.bm.energy << ","
            << p.bm.start << ","
            << p.bm.end << ","
            << p.bm.parent << ","
            << p.bm.time << ","
            << p.coincidence << ","
            << p.bm.in_signal << ","
            << p.bm.neutron_induced << "\n";
    }

    out.close();
    std::cout << "[INFO] Saved " << all_protons.size() << " protons to output file.\n";
    std::cout << "[INFO] Output saved to /output_files/protons_" << version << "_" << mode << ".txt\n";
    std::cout << "[INFO] There are " << stats.light_candidate << " light candidate protons\n";
    std::cout << "[INFO] There are " << stats.missed_sec << " prim protons that match true secondary protons\n";
        // === ROOT Output ===
    std::string root_output_name = "../n_selection/output_files/protons_MiniRun" + std::string(version) + "_" + std::string(mode) + ".root";
    TFile* fout = new TFile(root_output_name.c_str(), "RECREATE");
    TTree* t = new TTree("protons", "Reconstructed Proton Info");

    // === Define variables for branches ===
    std::string file; int event, rvtx, rtvtx, matched;
    int rpdg, rtype; float rlen, rdist, rE;
    int t_int_idx, tpart_idx, tpdg, ttype;
    float tlen, tdist, tE;
    int tparent; float tt0;
    int coincidence, insignal, ninduced;

    // === Set up branches ===
    t->Branch("file", &file);
    t->Branch("event", &event);
    t->Branch("rvtx", &rvtx);
    t->Branch("rtvtx", &rtvtx);
    t->Branch("matched", &matched);
    t->Branch("rpdg", &rpdg);
    t->Branch("rtype", &rtype);
    t->Branch("rlen", &rlen);
    t->Branch("rdist", &rdist);
    t->Branch("rE", &rE);
    t->Branch("t_int_idx", &t_int_idx);
    t->Branch("tpart_idx", &tpart_idx);
    t->Branch("tpdg", &tpdg);
    t->Branch("ttype", &ttype);
    t->Branch("tlen", &tlen);
    t->Branch("tdist", &tdist);
    t->Branch("tE", &tE);
    t->Branch("tparent", &tparent);
    t->Branch("tt0", &tt0);
    t->Branch("coincidence", &coincidence);
    t->Branch("insignal", &insignal);
    t->Branch("ninduced", &ninduced);

    // === Fill tree ===
    for (const auto& p : all_protons) {
        file = p.input_file;
        event = p.input_event;
        rvtx = p.reco_vtx;
        rtvtx = p.input_vtx_t;
        matched = p.input_matched;

        rpdg = p.reco_pdg;
        rtype = p.reco_type;
        rlen = p.reco_len;
        rdist = p.reco_dist;
        rE = p.reco_energy;

        t_int_idx = p.bm.interaction_idx;
        tpart_idx = p.bm.particle_idx;
        tpdg = p.bm.pdg;
        ttype = p.bm.type;
        tlen = p.bm.length;
        tdist = p.bm.distance;
        tE = p.bm.energy;
        tparent = p.bm.parent;
        tt0 = p.bm.time;

        coincidence = p.coincidence;
        insignal = p.bm.in_signal;
        ninduced = p.bm.neutron_induced;

        t->Fill();
    }

    // === Write and close ===
    t->Write();
    fout->Close();
    std::cout << "[INFO] Wrote ROOT file to: " << root_output_name << std::endl;


    delete sr;
    delete caf_tree;
    return 0;
}

