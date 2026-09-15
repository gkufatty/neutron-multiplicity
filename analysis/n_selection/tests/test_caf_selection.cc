#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "structures.h"
#include "utils.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>

int select_2x2_neutrons(const std::string&, const char*, const char*);

namespace {
int checks = 0;
void Check(bool condition, const std::string& message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
void Throws(const std::function<void()>& action, const std::string& expected) {
    try { action(); }
    catch (const std::runtime_error& error) {
        Check(std::string(error.what()).find(expected) != std::string::npos, error.what());
        return;
    }
    throw std::runtime_error("Expected failure: " + expected);
}
void WriteCaf(const std::string& path, const std::vector<int>& events) {
    TFile file(path.c_str(), "RECREATE");
    caf::StandardRecord record;
    caf::StandardRecord* sr = &record;
    TTree tree("cafTree", "Test CAF");
    tree.Branch("rec", &sr);
    for (const int event : events) {
        record = caf::StandardRecord{};
        record.meta.nd_lar.event = event;
        record.common.ixn.dlp.resize(2);
        auto& reco = record.common.ixn.dlp[0];
        reco.vtx = caf::SRVector3D(10, 0, 10);
        caf::SRRecoParticle proton;
        proton.pdg = 2212;
        proton.primary = false;
        proton.E = 0.05;
        proton.start = caf::SRVector3D(15, 0, 10);
        proton.end = caf::SRVector3D(16, 0, 10);
        reco.part.dlp.push_back(proton); // Unmatched candidate must survive.
        if (event == 26) {
            record.mc.nu.resize(1);
            auto& nu = record.mc.nu[0];
            nu.iscc = false; // Legacy insignal must remain false.
            nu.prim.resize(1);
            nu.prim[0].pdg = 2112;
            nu.prim[0].G4ID = 123;
            nu.sec.resize(1);
            nu.sec[0].pdg = 2212;
            nu.sec[0].parent = 123;
            nu.sec[0].p.E = 1;
            nu.sec[0].start_pos = proton.start;
            nu.sec[0].end_pos = proton.end;
            caf::TrueParticleID id;
            id.ixn = 0; id.part = 0; id.type = caf::TrueParticleID::kSecondary;
            proton.truth = {id};
            proton.truthOverlap = {0.8f};
            reco.part.dlp.push_back(proton);
        }
        tree.Fill();
    }
    tree.Write();
    tree.ResetBranchAddresses();
}
const std::string header =
    "file_name,event,reco_ixn,truth_ixn,has_truth_match,matched_true_signal,"
    "true_track_multiplicity,reco_track_multiplicity,true_primary_neutron_count,"
    "true_secondary_neutron_count\n";
void WriteSelection(const std::string& contents) {
    std::ofstream csv("selection.csv");
    csv << header << contents;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    try {
        std::filesystem::current_path(argv[1]);
        const std::string caf = (std::filesystem::current_path() / "test,\"caf\".root").string();
        WriteCaf(caf, {99, 26, -1, -1, -2});
        const auto path = CsvField(caf);
        WriteSelection(path + ",26,0,0,1,1,7,2,5,0\n" +
                       path + ",26,1,-1,0,0,-1,0,-1,-1\n" +
                       path + ",-1,0,-1,0,0,-1,1,-1,-1\n" +
                       path + ",-2,0,-1,0,0,-1,1,-1,-1\n" +
                       "absent-skipped-file.root,-1,0,-1,0,0,-1,0,-1,-1\n");
        Check(select_2x2_neutrons("selection.csv", "6.5", "RHC") == 0, "Selection failed");

        TFile result("output_files/protons_MiniRun6.5_RHC.root", "READ");
        TTree* tree = nullptr;
        result.GetObject("protons", tree);
        Check(tree && tree->GetEntries() == 3, "Expected exactly three selected candidates");
        Check(std::string(tree->GetLeaf("event")->GetTypeName()) == "Long64_t", "64-bit event branch");
        int valid_count = 0;
        for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);
            auto value = [&](const char* name) { return tree->GetLeaf(name)->GetValue(); };
            Check(value("event") == 26 || value("event") == -2,
                  "Unselected or skipped event leaked into output");
            Check(value("rvtx") == 0, "Zero-candidate interaction emitted a proton");
            if (value("event") == 26) {
                Check(value("matched") == 1 && value("insignal") == 0 && value("has_truth_match") == 1,
                      "Input signal flag must be independent of legacy insignal");
                Check(value("true_track_multiplicity") == 7 && value("reco_track_multiplicity") == 2 &&
                      value("true_primary_neutron_count") == 5 && value("true_secondary_neutron_count") == 0,
                      "Metadata not propagated");
            } else {
                Check(value("has_truth_match") == 0 && value("rtvtx") == -1 &&
                      value("true_primary_neutron_count") == -1, "Unmatched interaction metadata");
            }
            if (value("has_particle_truth_match")) {
                ++valid_count;
                Check(value("coincidence") == 1 && value("ninduced") == 1, "Neutron parent match");
                Check(std::abs(value("toverlap") - 0.8) < 1e-6,
                      "Selected truth overlap");
            } else {
                Check(value("t_int_idx") == -1 && value("tpart_idx") == -1 &&
                      std::isnan(value("toverlap")) && std::isnan(value("tE")),
                      "Unavailable particle truth");
            }
        }
        Check(valid_count == 1, "Particle truth validity count");
        std::ifstream text("output_files/protons_MiniRun6.5_RHC.txt");
        std::vector<std::string> fields;
        Check(ReadCsvRecord(text, fields) && fields.size() == 33, "Text header width");
        int text_rows = 0;
        while (ReadCsvRecord(text, fields)) {
            Check(fields.size() == 33 && fields[0] == caf, "Text CSV quoting / width");
            ++text_rows;
        }
        Check(text_rows == 3, "Text/ROOT candidate count mismatch");

        WriteSelection(path + ",404,0,-1,0,0,-1,0,-1,-1\n");
        Throws([&] { select_2x2_neutrons("selection.csv", "6.6", "RHC"); }, "event not found");
        Check(!std::filesystem::exists("output_files/protons_MiniRun6.6_RHC.root"), "Failure created output");
        WriteSelection(path + ",26,2,-1,0,0,-1,0,-1,-1\n");
        Throws([&] { select_2x2_neutrons("selection.csv", "6.6", "RHC"); }, "reco_ixn out of bounds");
        WriteSelection(path + ",26,0,1,1,0,1,1,1,0\n");
        Throws([&] { select_2x2_neutrons("selection.csv", "6.6", "RHC"); }, "truth_ixn out of bounds");
        const auto duplicate = (std::filesystem::current_path() / "duplicate.root").string();
        WriteCaf(duplicate, {26, 26, -1, -1});
        WriteSelection(CsvField(duplicate) + ",26,0,-1,0,0,-1,0,-1,-1\n");
        Throws([&] { select_2x2_neutrons("selection.csv", "6.6", "RHC"); }, "ambiguous event ID");

        // All-skipped input is a successful empty selection; its CAFs need
        // not be opened, even if they do not exist or contain duplicate IDs.
        WriteSelection(CsvField(duplicate) + ",-1,0,-1,0,0,-1,0,-1,-1\n" +
                       "absent-skipped-file.root,-1,0,-1,0,0,-1,0,-1,-1\n");
        Check(select_2x2_neutrons("selection.csv", "6.7", "RHC") == 0,
              "All-skipped input must succeed");
        TFile empty_result("output_files/protons_MiniRun6.7_RHC.root", "READ");
        TTree* empty_tree = nullptr;
        empty_result.GetObject("protons", empty_tree);
        Check(empty_tree && empty_tree->GetEntries() == 0, "All-skipped ROOT output must be empty");
        std::ifstream empty_text("output_files/protons_MiniRun6.7_RHC.txt");
        Check(ReadCsvRecord(empty_text, fields) && fields.size() == 33 &&
              !ReadCsvRecord(empty_text, fields), "All-skipped text output must contain only the header");

        caf::StandardRecord sr;
        caf::TrueParticleID id;
        id.type = caf::TrueParticleID::kSecondary; id.ixn = 0; id.part = 0;
        Check(!FindParticleBestMatch({}, {}, &sr).valid, "Empty particle match");
        Throws([&] { FindParticleBestMatch({id}, {}, &sr); }, "lengths differ");
        Throws([&] { FindParticleBestMatch({id}, {1}, &sr); }, "interaction index out of bounds");
        sr.mc.nu.resize(1);
        Throws([&] { FindParticleBestMatch({id}, {1}, &sr); }, "particle index out of bounds");
        id.type = caf::TrueParticleID::kUnknown;
        Check(!FindParticleBestMatch({id}, {1}, &sr).valid, "Unsupported particle truth type");

        // Reproduce the three corrected truth-labeling cases. The first two
        // checks intentionally emit the corresponding avoidance diagnostics.
        auto& nu = sr.mc.nu[0];
        nu.prim.resize(1);
        nu.prim[0].pdg = ParticleCode::neutron;
        nu.prim[0].G4ID = -1;
        nu.sec.resize(1);
        nu.sec[0].pdg = ParticleCode::pip;
        nu.sec[0].parent = -1;
        nu.sec[0].p.E = 0.2f;
        id.type = caf::TrueParticleID::kSecondary;
        Check(!was_neutron_induced(0, 0, &sr),
              "Invalid neutron and parent IDs must not match");
        const auto zero_overlap = FindParticleBestMatch({id}, {0.0f}, &sr);
        Check(!zero_overlap.valid && std::isnan(zero_overlap.overlap),
              "Zero-overlap particle truth must be invalid");

        nu.prim[0].G4ID = 123;
        nu.sec[0].parent = 123;
        const auto pion_match = FindParticleBestMatch({id}, {0.75f}, &sr);
        Check(pion_match.valid && pion_match.neutron_induced &&
              std::abs(pion_match.overlap - 0.75f) < 1e-6f,
              "Positive-overlap pion truth match");
        Check(std::abs(pion_match.energy - 60.43f) < 0.02f,
              "Pion kinetic energy must use the pion mass");
        std::cout << "PASS: " << checks << " ROOT/CAF integration checks\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
