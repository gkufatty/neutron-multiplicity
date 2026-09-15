#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "structures.h"
#include "utils.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
// Detach addresses before the stack-owned StandardRecord goes out of scope,
// including when a read or validation throws.
struct BranchAddressGuard {
    TTree& tree;
    ~BranchAddressGuard() { tree.ResetBranchAddresses(); }
};

std::string VectorCsvField(const caf::SRVector3D& value) {
    std::ostringstream text;
    text << value;
    return CsvField(text.str());
}

void ProcessSelectedFile(const std::string& file_name,
                         const std::vector<const InputCAFRow*>& rows,
                         const char* mode, Counters& stats,
                         std::vector<RecoProtonInfo>& all_protons,
                         std::size_t& processed, std::size_t& zero_candidates) {
    try {
        std::unique_ptr<TFile> input(TFile::Open(file_name.c_str(), "READ"));
        if (!input || input->IsZombie()) throw std::runtime_error("Cannot open CAF");
        TTree* tree = nullptr;
        input->GetObject("cafTree", tree);
        if (!tree) throw std::runtime_error("Missing cafTree");
        auto* branch = tree->GetBranch("rec");
        if (!branch) throw std::runtime_error("Missing rec branch");
        caf::StandardRecord record;
        caf::StandardRecord* sr = &record;
        BranchAddressGuard guard{*tree};
        branch->SetAutoDelete(false);
        if (tree->SetBranchAddress("rec", &sr) < 0)
            throw std::runtime_error("Cannot bind rec branch; check CAF dictionaries");
        const auto event_index = BuildEventEntryIndex(*tree, sr);
        std::map<std::int64_t, std::vector<const InputCAFRow*>> selected_entries;
        for (const auto* row : rows)
            selected_entries[ResolveEventEntry(event_index, *row)].push_back(row);

        for (const auto& entry : selected_entries) {
            if (tree->GetEntry(entry.first) <= 0 || !sr)
                throw std::runtime_error("Cannot read selected entry " + std::to_string(entry.first));
            for (const auto* row : entry.second) {
                try {
                    if (sr->meta.nd_lar.event != row->event)
                        throw std::runtime_error("Event changed between index and selected read");
                    if (row->vtx_r < 0 || static_cast<std::size_t>(row->vtx_r) >= sr->common.ixn.dlp.size())
                        throw std::runtime_error("reco_ixn out of bounds");
                    if (row->has_truth_match && (row->vtx_t < 0 ||
                        static_cast<std::size_t>(row->vtx_t) >= sr->mc.nu.size()))
                        throw std::runtime_error("truth_ixn out of bounds");
                    const auto& reco = sr->common.ixn.dlp[row->vtx_r];
                    auto protons = process_protons(row->vtx_r, reco, sr, stats, mode);
                    ++processed;
                    if (protons.empty()) ++zero_candidates;
                    for (auto& proton : protons) {
                        AttachInputMetadata(proton, *row);
                        all_protons.push_back(std::move(proton));
                    }
                } catch (const std::exception& error) {
                    throw std::runtime_error(SelectionRowContext(*row) + ": " + error.what());
                }
            }
        }
    } catch (const std::exception& error) {
        throw std::runtime_error(file_name + ": " + error.what());
    }
}

void WriteOutputs(const std::vector<RecoProtonInfo>& all_protons,
                  const char* version, const char* mode) {
    std::filesystem::create_directories("output_files");
    const std::string stem = "output_files/protons_MiniRun" + std::string(version) + "_" + mode;
    std::ofstream out(stem + ".txt");
    if (!out) throw std::runtime_error("Cannot create " + stem + ".txt");
    out << "file,event,rvtx,rtvtx,matched,"
        << "rpdg,rtype,rlen,rdist,rE,rstart,rend,"
        << "t_int_idx,tpart_idx,tpdg,ttype,toverlap,tlen,tdist,tE,"
        << "tstart,tend,tparent,tt0,coincidence,insignal,ninduced,"
        << "has_truth_match,true_track_multiplicity,reco_track_multiplicity,"
        << "true_primary_neutron_count,true_secondary_neutron_count,has_particle_truth_match\n";

    for (const auto& p : all_protons) {
        out << CsvField(p.input_file) << ","
            << p.input_event << ","
            << p.reco_vtx << ","
            << p.input_vtx_t << ","
            << p.input_matched << ","
            << p.reco_pdg << ","
            << p.reco_type << ","
            << p.reco_len << ","
            << p.reco_dist << ","
            << p.reco_energy << ","
            << VectorCsvField(p.reco_start) << ","
            << VectorCsvField(p.reco_end) << ","
            << p.bm.interaction_idx << ","
            << p.bm.particle_idx << ","
            << p.bm.pdg << ","
            << p.bm.type << ","
            << p.bm.overlap << ","
            << p.bm.length << ","
            << p.bm.distance << ","
            << p.bm.energy << ","
            << VectorCsvField(p.bm.start) << ","
            << VectorCsvField(p.bm.end) << ","
            << p.bm.parent << ","
            << p.bm.time << ","
            << p.coincidence << ","
            << p.bm.in_signal << ","
            << p.bm.neutron_induced << ","
            << p.input_has_truth_match << ","
            << p.input_true_track_multiplicity << ","
            << p.input_reco_track_multiplicity << ","
            << p.input_true_primary_neutron_count << ","
            << p.input_true_secondary_neutron_count << ","
            << p.bm.valid << "\n";
    }

    out.close();
    if (!out) throw std::runtime_error("Failed writing " + stem + ".txt");

    const std::string root_output_name = stem + ".root";
    std::unique_ptr<TFile> fout(TFile::Open(root_output_name.c_str(), "RECREATE"));
    if (!fout || fout->IsZombie() || !fout->IsWritable())
        throw std::runtime_error("Cannot create " + root_output_name);
    TTree tree("protons", "Reconstructed Proton Info");
    tree.SetDirectory(nullptr);
    TTree* t = &tree;

    // === Define variables for branches ===
    std::string file; Long64_t event; int rvtx, rtvtx, matched;
    int rpdg, rtype; float rlen, rdist, rE;
    int t_int_idx, tpart_idx, tpdg, ttype;
    float toverlap, tlen, tdist, tE;
    int tparent; float tt0;
    int coincidence, insignal, ninduced;
    int has_truth_match, true_track_multiplicity, reco_track_multiplicity;
    int true_primary_neutron_count, true_secondary_neutron_count, has_particle_truth_match;

    // === Set up branches ===
    t->Branch("file", &file);
    t->Branch("event", &event, "event/L");
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
    t->Branch("toverlap", &toverlap);
    t->Branch("tlen", &tlen);
    t->Branch("tdist", &tdist);
    t->Branch("tE", &tE);
    t->Branch("tparent", &tparent);
    t->Branch("tt0", &tt0);
    t->Branch("coincidence", &coincidence);
    t->Branch("insignal", &insignal);
    t->Branch("ninduced", &ninduced);
    t->Branch("has_truth_match", &has_truth_match);
    t->Branch("true_track_multiplicity", &true_track_multiplicity);
    t->Branch("reco_track_multiplicity", &reco_track_multiplicity);
    t->Branch("true_primary_neutron_count", &true_primary_neutron_count);
    t->Branch("true_secondary_neutron_count", &true_secondary_neutron_count);
    t->Branch("has_particle_truth_match", &has_particle_truth_match);

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
        toverlap = p.bm.overlap;
        tlen = p.bm.length;
        tdist = p.bm.distance;
        tE = p.bm.energy;
        tparent = p.bm.parent;
        tt0 = p.bm.time;

        coincidence = p.coincidence;
        insignal = p.bm.in_signal;
        ninduced = p.bm.neutron_induced;

        has_truth_match = p.input_has_truth_match;
        true_track_multiplicity = p.input_true_track_multiplicity;
        reco_track_multiplicity = p.input_reco_track_multiplicity;
        true_primary_neutron_count = p.input_true_primary_neutron_count;
        true_secondary_neutron_count = p.input_true_secondary_neutron_count;
        has_particle_truth_match = p.bm.valid;
        if (t->Fill() < 0) throw std::runtime_error("Failed filling proton output tree");
    }

    // === Write and close ===
    fout->cd();
    if (t->Write() <= 0) throw std::runtime_error("Failed writing " + root_output_name);
    fout->Close();
    if (fout->TestBit(TFile::kWriteError))
        throw std::runtime_error("ROOT write error: " + root_output_name);
    std::cout << "[INFO] Saved " << all_protons.size() << " proton candidates to "
              << stem << ".txt and " << root_output_name << '\n';


}
} // namespace

int select_2x2_neutrons(const std::string& input_csv, const char* version, const char* mode) {
    const auto input_rows = ReadSelectionCsv(input_csv);
    std::map<std::string, std::vector<const InputCAFRow*>> files;
    std::size_t unmatched = 0, backgrounds = 0, skipped_events = 0;
    for (const auto& row : input_rows) {
        unmatched += !row.has_truth_match;
        backgrounds += row.has_truth_match && !row.matched_true_signal;
        // Temporary policy: unavailable event IDs are excluded before CAF
        // lookup. Other event IDs still require a unique matching entry.
        if (row.event == -1) {
            ++skipped_events;
            continue;
        }
        files[row.file_name].push_back(&row);
    }
    const auto eligible_rows = input_rows.size() - skipped_events;
    std::cout << "[INFO] Loaded " << input_rows.size()
              << " selected interactions; unmatched: " << unmatched
              << "; matched backgrounds: " << backgrounds << std::endl;
    std::cout << "[INFO] Skipped " << skipped_events << " CSV rows with event == -1; "
              << eligible_rows << " interactions eligible in " << files.size()
              << " CAF files" << std::endl;
    Counters stats;
    std::vector<RecoProtonInfo> all_protons;
    std::size_t processed = 0, zero_candidates = 0, completed_files = 0;
    for (const auto& file : files) {
        ProcessSelectedFile(file.first, file.second, mode, stats, all_protons,
                            processed, zero_candidates);
        if (++completed_files % 25 == 0 || completed_files == files.size())
            std::cout << "[INFO] Processed " << completed_files << '/' << files.size()
                      << " files, " << processed << " selected interactions" << std::endl;
    }
    if (processed != eligible_rows)
        throw std::runtime_error("Processed interaction count differs from eligible CSV row count");
    // Output starts only after every non-skipped interaction has resolved and run.
    WriteOutputs(all_protons, version, mode);
    std::cout << "[INFO] Selected interactions processed: " << processed
              << "; with zero proton candidates: " << zero_candidates
              << "; skipped (event == -1): " << skipped_events << '\n';
    std::cout << "[INFO] Secondary proton candidates: " << stats.secp_reco
              << "; matched neutron-induced secondary protons: " << stats.np_reco << '\n';
    return 0;
}

#ifndef N_SELECTION_NO_MAIN
int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "USAGE: " << argv[0] << " <selection.csv> <version> <RHC|FHC>\n";
        return 1;
    }
    try {
        const std::string version = argv[2];
        const std::string mode = argv[3];
        if (version.empty() || version.find_first_not_of("0123456789.") != std::string::npos)
            throw std::runtime_error("Version must contain only digits and dots, e.g. 6.5");
        if (mode != "RHC" && mode != "FHC")
            throw std::runtime_error("Mode must be RHC or FHC");
        return select_2x2_neutrons(argv[1], version.c_str(), mode.c_str());
    } catch (const std::exception& error) {
        std::cerr << "[ERROR] " << error.what() << '\n';
        return 1;
    }
}
#endif
