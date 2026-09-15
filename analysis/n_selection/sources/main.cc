#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "structures.h"
#include "utils.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
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
                         std::vector<TrueNeutronProtonInfo>& all_true_protons,
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
            std::set<int> selected_reco_vtxs;
            std::set<int> selected_truth_vtxs;
            for (const auto* row : entry.second) {
                try {
                    if (sr->meta.nd_lar.event != row->event)
                        throw std::runtime_error("Event changed between index and selected read");
                    if (row->vtx_r < 0 || static_cast<std::size_t>(row->vtx_r) >= sr->common.ixn.dlp.size())
                        throw std::runtime_error("reco_ixn out of bounds");
                    if (row->has_truth_match && (row->vtx_t < 0 ||
                        static_cast<std::size_t>(row->vtx_t) >= sr->mc.nu.size()))
                        throw std::runtime_error("truth_ixn out of bounds");
                    selected_reco_vtxs.insert(row->vtx_r);
                    if (row->has_truth_match)
                        selected_truth_vtxs.insert(row->vtx_t);
                    const auto& reco = sr->common.ixn.dlp[row->vtx_r];
                    auto protons = process_protons(row->vtx_r, reco, sr, stats, mode);
                    ++processed;
                    if (protons.empty()) ++zero_candidates;
                    for (auto& proton : protons) {
                        AttachInputMetadata(proton, *row);
                        if (row->has_truth_match)
                            proton.selected_truth_vertex_t0 = sr->mc.nu[row->vtx_t].time;
                        if (proton.bm.valid) {
                            proton.matched_truth_vertex_t0 =
                                sr->mc.nu[proton.bm.interaction_idx].time;
                            proton.matched_particle_dt = proton.bm.time -
                                proton.matched_truth_vertex_t0;
                        }
                        if (proton.coincidence && proton.bm.neutron_induced &&
                            proton.same_truth_interaction)
                            ++stats.np_reco;
                        all_protons.push_back(std::move(proton));
                    }
                } catch (const std::exception& error) {
                    throw std::runtime_error(SelectionRowContext(*row) + ": " + error.what());
                }
            }

            const std::vector<int> selected_reco_indices(
                selected_reco_vtxs.begin(), selected_reco_vtxs.end());
            for (const int truth_vtx : selected_truth_vtxs) {
                try {
                    auto true_protons = process_true_neutron_protons(
                        truth_vtx, selected_reco_indices, sr);
                    for (auto& proton : true_protons) {
                        proton.input_file = file_name;
                        proton.input_event = sr->meta.nd_lar.event;
                        all_true_protons.push_back(std::move(proton));
                    }
                } catch (const std::exception& error) {
                    throw std::runtime_error(
                        "event " + std::to_string(sr->meta.nd_lar.event) +
                        ", truth_ixn " + std::to_string(truth_vtx) + ": " +
                        error.what());
                }
            }
        }
    } catch (const std::exception& error) {
        throw std::runtime_error(file_name + ": " + error.what());
    }
}

void WriteTruthTree(TFile& output,
                    const std::vector<TrueNeutronProtonInfo>& true_protons,
                    const std::string& output_name) {
    TTree tree("true_neutron_protons", "True Neutron-Induced Protons");
    tree.SetDirectory(nullptr);

    std::string file;
    Long64_t event;
    int truth_ixn, true_part_idx, true_type, true_pdg;
    float true_energy, true_length, true_distance, true_time;
    float truth_vertex_time, truth_particle_dt;
    float true_start_x, true_start_y, true_start_z;
    float true_end_x, true_end_y, true_end_z;
    int true_end_in_fv;
    int truth_parent_g4id, direct_parent_pdg;
    int neutron_parent_type, neutron_parent_idx, neutron_parent_g4id;
    int has_any_reco_overlap;
    float max_any_reco_overlap;
    int has_reco_match, n_reco_matches;
    int best_reco_ixn, best_reco_part_idx, best_reco_pdg;
    int best_reco_primary;
    float best_reco_overlap;
    int best_reco_in_selected_interaction;
    int best_reco_start_in_fv, best_reco_end_in_fv;
    int best_reco_passes_select_n, has_selected_candidate;

    tree.Branch("file", &file);
    tree.Branch("event", &event, "event/L");
    tree.Branch("truth_ixn", &truth_ixn);
    tree.Branch("true_part_idx", &true_part_idx);
    tree.Branch("true_type", &true_type);
    tree.Branch("true_pdg", &true_pdg);
    tree.Branch("true_energy", &true_energy);
    tree.Branch("true_length", &true_length);
    tree.Branch("true_distance", &true_distance);
    tree.Branch("true_time", &true_time);
    tree.Branch("truth_vertex_time", &truth_vertex_time);
    tree.Branch("truth_particle_dt", &truth_particle_dt);
    tree.Branch("true_start_x", &true_start_x);
    tree.Branch("true_start_y", &true_start_y);
    tree.Branch("true_start_z", &true_start_z);
    tree.Branch("true_end_x", &true_end_x);
    tree.Branch("true_end_y", &true_end_y);
    tree.Branch("true_end_z", &true_end_z);
    tree.Branch("true_end_in_fv", &true_end_in_fv);
    tree.Branch("truth_parent_g4id", &truth_parent_g4id);
    tree.Branch("direct_parent_pdg", &direct_parent_pdg);
    tree.Branch("neutron_parent_type", &neutron_parent_type);
    tree.Branch("neutron_parent_idx", &neutron_parent_idx);
    tree.Branch("neutron_parent_g4id", &neutron_parent_g4id);
    tree.Branch("has_any_reco_overlap", &has_any_reco_overlap);
    tree.Branch("max_any_reco_overlap", &max_any_reco_overlap);
    tree.Branch("has_reco_match", &has_reco_match);
    tree.Branch("n_reco_matches", &n_reco_matches);
    tree.Branch("best_reco_ixn", &best_reco_ixn);
    tree.Branch("best_reco_part_idx", &best_reco_part_idx);
    tree.Branch("best_reco_pdg", &best_reco_pdg);
    tree.Branch("best_reco_primary", &best_reco_primary);
    tree.Branch("best_reco_overlap", &best_reco_overlap);
    tree.Branch("best_reco_in_selected_interaction",
                &best_reco_in_selected_interaction);
    tree.Branch("best_reco_start_in_fv", &best_reco_start_in_fv);
    tree.Branch("best_reco_end_in_fv", &best_reco_end_in_fv);
    tree.Branch("best_reco_passes_select_n", &best_reco_passes_select_n);
    tree.Branch("has_selected_candidate", &has_selected_candidate);

    for (const auto& proton : true_protons) {
        file = proton.input_file;
        event = proton.input_event;
        truth_ixn = proton.truth_vtx;
        true_part_idx = proton.truth_particle_idx;
        true_type = proton.truth_type;
        true_pdg = proton.truth_pdg;
        true_energy = proton.truth_energy;
        true_length = proton.truth_length;
        true_distance = proton.truth_distance;
        true_time = proton.truth_time;
        truth_vertex_time = proton.truth_vertex_time;
        truth_particle_dt = proton.truth_particle_dt;
        true_start_x = proton.truth_start.x;
        true_start_y = proton.truth_start.y;
        true_start_z = proton.truth_start.z;
        true_end_x = proton.truth_end.x;
        true_end_y = proton.truth_end.y;
        true_end_z = proton.truth_end.z;
        true_end_in_fv = proton.truth_end_in_fv;
        truth_parent_g4id = proton.truth_parent_g4id;
        direct_parent_pdg = proton.direct_parent_pdg;
        neutron_parent_type = proton.neutron_parent_type;
        neutron_parent_idx = proton.neutron_parent_idx;
        neutron_parent_g4id = proton.neutron_parent_g4id;
        has_any_reco_overlap = proton.has_any_reco_overlap;
        max_any_reco_overlap = proton.max_any_reco_overlap;
        has_reco_match = proton.has_reco_match;
        n_reco_matches = proton.n_reco_matches;
        best_reco_ixn = proton.best_reco_vtx;
        best_reco_part_idx = proton.best_reco_particle_idx;
        best_reco_pdg = proton.best_reco_pdg;
        best_reco_primary = proton.best_reco_primary;
        best_reco_overlap = proton.best_reco_overlap;
        best_reco_in_selected_interaction =
            proton.best_reco_in_selected_interaction;
        best_reco_start_in_fv = proton.best_reco_start_in_fv;
        best_reco_end_in_fv = proton.best_reco_end_in_fv;
        best_reco_passes_select_n = proton.best_reco_passes_select_n;
        has_selected_candidate = proton.has_selected_candidate;
        if (tree.Fill() < 0)
            throw std::runtime_error(
                "Failed filling true neutron-induced proton tree");
    }

    output.cd();
    if (tree.Write() <= 0)
        throw std::runtime_error("Failed writing truth tree to " + output_name);
}

void WriteOutputs(const std::vector<RecoProtonInfo>& all_protons,
                  const std::vector<TrueNeutronProtonInfo>& all_true_protons,
                  const char* version, const char* mode) {
    std::filesystem::create_directories("output_files");
    const std::string stem = "output_files/protons_MiniRun" + std::string(version) + "_" + mode;
    std::ofstream out(stem + ".txt");
    if (!out) throw std::runtime_error("Cannot create " + stem + ".txt");
    out << "file,event,rvtx,rtvtx,matched,"
        << "rpdg,rtype,rlen,rdist,rE,rstart,rend,"
        << "t_int_idx,tpart_idx,tpdg,ttype,toverlap,tlen,tdist,tE,"
        << "tstart,tend,tparent,direct_parent_pdg,tt0,selected_truth_vertex_t0,"
        << "matched_truth_vertex_t0,matched_particle_dt,"
        << "coincidence,insignal,ninduced,"
        << "neutron_parent_type,neutron_parent_idx,neutron_parent_g4id,"
        << "same_truth_interaction,"
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
            << p.bm.direct_parent_pdg << ","
            << p.bm.time << ","
            << p.selected_truth_vertex_t0 << ","
            << p.matched_truth_vertex_t0 << ","
            << p.matched_particle_dt << ","
            << p.coincidence << ","
            << p.bm.in_signal << ","
            << p.bm.neutron_induced << ","
            << p.bm.neutron_parent_type << ","
            << p.bm.neutron_parent_idx << ","
            << p.bm.neutron_parent_g4id << ","
            << p.same_truth_interaction << ","
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
    int tparent, direct_parent_pdg;
    float tt0, selected_truth_vertex_t0, matched_truth_vertex_t0;
    float matched_particle_dt;
    int coincidence, insignal, ninduced;
    int neutron_parent_type, neutron_parent_idx, neutron_parent_g4id;
    int same_truth_interaction;
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
    t->Branch("direct_parent_pdg", &direct_parent_pdg);
    t->Branch("tt0", &tt0);
    t->Branch("selected_truth_vertex_t0", &selected_truth_vertex_t0);
    t->Branch("matched_truth_vertex_t0", &matched_truth_vertex_t0);
    t->Branch("matched_particle_dt", &matched_particle_dt);
    t->Branch("coincidence", &coincidence);
    t->Branch("insignal", &insignal);
    t->Branch("ninduced", &ninduced);
    t->Branch("neutron_parent_type", &neutron_parent_type);
    t->Branch("neutron_parent_idx", &neutron_parent_idx);
    t->Branch("neutron_parent_g4id", &neutron_parent_g4id);
    t->Branch("same_truth_interaction", &same_truth_interaction);
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
        direct_parent_pdg = p.bm.direct_parent_pdg;
        tt0 = p.bm.time;
        selected_truth_vertex_t0 = p.selected_truth_vertex_t0;
        matched_truth_vertex_t0 = p.matched_truth_vertex_t0;
        matched_particle_dt = p.matched_particle_dt;

        coincidence = p.coincidence;
        insignal = p.bm.in_signal;
        ninduced = p.bm.neutron_induced;
        neutron_parent_type = p.bm.neutron_parent_type;
        neutron_parent_idx = p.bm.neutron_parent_idx;
        neutron_parent_g4id = p.bm.neutron_parent_g4id;
        same_truth_interaction = p.same_truth_interaction;

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
    WriteTruthTree(*fout, all_true_protons, root_output_name);
    fout->Close();
    if (fout->TestBit(TFile::kWriteError))
        throw std::runtime_error("ROOT write error: " + root_output_name);
    std::cout << "[INFO] Saved " << all_protons.size() << " proton candidates to "
              << stem << ".txt and " << root_output_name << "; saved "
              << all_true_protons.size()
              << " truth neutron-induced protons to the ROOT file\n";


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
    std::vector<TrueNeutronProtonInfo> all_true_protons;
    std::size_t processed = 0, zero_candidates = 0, completed_files = 0;
    for (const auto& file : files) {
        ProcessSelectedFile(file.first, file.second, mode, stats, all_protons,
                            all_true_protons, processed, zero_candidates);
        if (++completed_files % 25 == 0 || completed_files == files.size())
            std::cout << "[INFO] Processed " << completed_files << '/' << files.size()
                      << " files, " << processed << " selected interactions" << std::endl;
    }
    if (processed != eligible_rows)
        throw std::runtime_error("Processed interaction count differs from eligible CSV row count");
    // Output starts only after every non-skipped interaction has resolved and run.
    WriteOutputs(all_protons, all_true_protons, version, mode);
    std::cout << "[INFO] Selected interactions processed: " << processed
              << "; with zero proton candidates: " << zero_candidates
              << "; skipped (event == -1): " << skipped_events << '\n';
    std::cout << "[INFO] Secondary proton candidates: " << stats.secp_reco
              << "; matched neutron-induced secondary protons in the selected "
              << "truth interaction: " << stats.np_reco << '\n';
    const auto unreconstructed = std::count_if(
        all_true_protons.begin(), all_true_protons.end(),
        [](const TrueNeutronProtonInfo& proton) {
            return !proton.has_any_reco_overlap;
        });
    std::cout << "[INFO] Start-contained true neutron-induced protons: "
              << all_true_protons.size() << "; without any positive reco overlap: "
              << unreconstructed << '\n';
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
