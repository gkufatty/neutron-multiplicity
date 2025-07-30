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
int select_2x2_nu(std::string input_file_list, const char* version ="6.2", const char* mode = "RHC");

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
    std::string input_file_list = "input_files/MiniRun" + version + "_"+ mode+".txt";
    std::cout << "Using MiniRun" << version << " in " << mode << " mode." << std::endl;
    select_2x2_nu(input_file_list, version.c_str(), mode.c_str());

    
    return 0;
}

// ======================
// FUNCTION DEFINITIONS
// ======================

// === True Analysis Functions ===




int select_2x2_nu(std::string input_file_list,const char* version, const char* mode){
    std::cout << "Running neutrino selector..." << std::endl;

    std::ifstream caf_list(input_file_list.c_str());
    if (!caf_list.is_open()) {
        std::cerr << Form("File %s not found", input_file_list.c_str()) << std::endl;
        return 1;
    }

    // Count CAF files
    int total_files = 0;
    std::string line;
    while (std::getline(caf_list, line)) total_files++;
    std::cout << "Total number of CAF trees = " << total_files << std::endl;

    caf_list.clear();
    caf_list.seekg(0, std::ios::beg);

    // Load into TChain
    TChain *caf_tree = new TChain("cafTree");
    std::string tmp;
    while (caf_list >> tmp) {
        caf_tree->Add(tmp.c_str());
    }

    // Link CAF
    caf::StandardRecord* sr = new caf::StandardRecord;
    caf_tree->SetBranchAddress("rec", &sr);
    int tree_entries = caf_tree->GetEntries();
    std::cout << "Number of entries: " << tree_entries << std::endl;


    caf::StandardRecord* out_sr = new caf::StandardRecord();

    std::string out_filename = "./output_files/MiniRun" + std::string(version) +"_"+std::string(mode) + "_eff_purity.csv";
    std::ofstream file_out(out_filename);

    std::vector<std::string> cuts;
    cuts.push_back("All");
    cuts.push_back("FV");
    cuts.push_back("2x2");
    cuts.push_back("Mx2");
    cuts.push_back("OnePrim");
    cuts.push_back("QE-like");

    for(int j = 0; j < cuts.size(); j++){
    
    // Stats
    Counters stats;
    std::set<std::string> unique_file_paths;
    std::vector<TrueVertexSelection_np> np_true_vertices;
    std::cout << "Getting results for cut: " << cuts.at(j) << std::endl;
            // Loop over entries
        for (int n = 0; n < tree_entries; n++) {
            caf_tree->GetEntry(n);
            std::string full_path = caf_tree->GetFile()->GetName();  // Use full path

            auto& all_vertices = np_true_vertices;
            for (size_t i = 0; i < sr->mc.nu.size(); i++) {
                if (QELikeSignal(sr, i, mode)) {
                    stats.nu_true++;
                    stats.prim_neutrons += sr->mc.nu[i].nneutron;
                    int event_id = sr->meta.nd_lar.event;
                    auto np_list = look_sec_particles(i, sr,stats);
                    all_vertices.insert(all_vertices.end(), np_list.begin(), np_list.end());
                }

            }

            for (size_t i = 0; i < sr->common.ixn.dlp.size(); ++i) {
                auto [signal, best_match_idx] = process_reco_interaction(sr, i,cuts.at(j));
                if (best_match_idx == -1) continue;
                if (signal){
                    stats.nu_reco++;
                    bool matched_true_signal = QELikeSignal(sr, best_match_idx, mode);
                    if (matched_true_signal) {
                        stats.matched++;
                    }
                }
            }        
        }
        float nu_efficiency = static_cast<double>(stats.matched) / stats.nu_true;
        float nu_purity = static_cast<double>(stats.matched) / stats.nu_reco;
        file_out << nu_efficiency << "," << nu_purity << "\n";
        print_stats(stats);



    }
    file_out.close();
   


  
    return 0;
}




