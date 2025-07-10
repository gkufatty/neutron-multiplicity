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
#include "cuts.h"
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
            if (
                row.file_name == full_path && row.vtx_r >= 0 &&
                row.vtx_r < static_cast<int>(sr->common.ixn.dlp.size())
                ){
                const auto& target_ixn = sr->common.ixn.dlp[row.vtx_r];
                
                // === Step 1: find r_vtx tpc ===
                TPC rvtx_tpc = determineTPC(target_ixn.vtx);
                // find tpc of the vertex. 
                // === Step 2: Iterate over all reco vertices in the event ===
                for (size_t v = 0; v < sr->common.ixn.dlp.size(); ++v) {
                    const auto& vtx = sr->common.ixn.dlp[v];

                    // === Step 3: Loop through its particles ===
                    for (const auto& part : vtx.part.dlp) {
                        if (!(kIsVtxFV(part.start))|| !(kIsVtxFV(p.end))) continue;
                        if (part.pdg == 2212 && part.primary) { // proton
                            PartBestMatch bm = FindParticleBestMatch(part.truth, part.truthOverlap, sr, mode);
                            if (bm.neutron_induced){
                                stats.recovered_n_induced++;
                                TPC part_tpc = determineTPC(part.start);
                                //std::cout << "proton found in tpc " << tpcToString(part_tpc) << " vertex target in " << tpcToString(rvtx_tpc)<< std::endl;
                                if(part_tpc != rvtx_tpc) {
                                    stats.diff_tpc++;
                                }
                                else{
                                    stats.same_tpc++;
                                }
                            }

                            // Optionally: store into a RecoProtonInfo struct, etc.
                            // RecoProtonInfo pinfo{sr->hdr.run, sr->hdr.subrun, sr->hdr.evt, v, start};
                            // all_protons.push_back(pinfo);
                        }
                    }
                }
            }
        }
    }
    std::cout<< "Total n induced protons to recover: " << stats.recovered_n_induced << std::endl;
    std::cout<< "Total reco protons in different TPC than vertex: " << stats.diff_tpc << std::endl; 
    std::cout<< "Total reco protons in same TPC as vertex: " << stats.same_tpc << std::endl;

    delete sr;
    delete caf_tree;
    return 0;
}

