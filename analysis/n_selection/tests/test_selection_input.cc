#include "selection_input.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
int checks = 0;
void Check(bool condition, const std::string& message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}
void Throws(const std::function<void()>& action, const std::string& diagnostic) {
    try { action(); }
    catch (const std::runtime_error& error) {
        Check(std::string(error.what()).find(diagnostic) != std::string::npos,
              "Wrong diagnostic: " + std::string(error.what()));
        return;
    }
    throw std::runtime_error("Expected error containing: " + diagnostic);
}
const std::vector<std::string> header = {
    "file_name", "event", "reco_ixn", "truth_ixn", "has_truth_match",
    "matched_true_signal", "true_track_multiplicity", "reco_track_multiplicity",
    "true_primary_neutron_count", "true_secondary_neutron_count"
};
const std::vector<std::string> sample = {"a.root", "26", "2", "3", "1", "0", "7", "2", "5", "0"};
std::string Record(const std::vector<std::string>& fields) {
    std::string result;
    for (const auto& field : fields) {
        if (!result.empty()) result += ',';
        result += CsvField(field);
    }
    return result + "\r\n";
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) return 2;
    try {
        const auto temporary = std::filesystem::path(argv[1]) / "fixture.csv";
        auto parse = [&](const std::string& contents) {
            { std::ofstream out(temporary); out << contents; }
            return ReadSelectionCsv(temporary.string());
        };
        auto rows = parse(Record(header) + Record(sample));
        Check(rows.size() == 1 && rows[0].file_name == "a.root", "Unquoted path");
        Check(rows[0].has_truth_match && !rows[0].matched_true_signal, "Truth flags conflated");
        Check(rows[0].vtx_r == 2 && rows[0].vtx_t == 3 && rows[0].csv_record == 2, "Column mapping");
        Check(rows[0].true_track_multiplicity == 7 && rows[0].reco_track_multiplicity == 2 &&
              rows[0].true_primary_neutron_count == 5 && rows[0].true_secondary_neutron_count == 0,
              "Multiplicity mapping");

        auto changed = sample;
        changed[0] = "dir/a,\"quoted\"\nfile.root";
        changed[1] = "5000000000";
        rows = parse(Record(header) + Record(changed));
        Check(rows[0].file_name == changed[0] && rows[0].event == 5000000000LL,
              "Quoted path / multiline field / 64-bit event round trip");
        rows = parse(Record(header) + "\"a.root\",26,2,3,1,0,7,2,5,0");
        Check(rows[0].file_name == "a.root", "Quoted filename without final newline");

        auto reversed_header = header, reversed_row = sample;
        std::reverse(reversed_header.begin(), reversed_header.end());
        std::reverse(reversed_row.begin(), reversed_row.end());
        reversed_header.push_back("extra"); reversed_row.push_back("ignored");
        rows = parse(Record(reversed_header) + Record(reversed_row));
        Check(rows[0].event == 26 && rows[0].vtx_t == 3, "Header order / extra column");
        rows = parse(std::string("\xEF\xBB\xBF") + Record(header) + Record(sample));
        Check(rows.size() == 1, "UTF-8 BOM");

        changed = {"a.root", "-1", "0", "-1", "false", "false", "-1", "0", "-1", "-1"};
        rows = parse(Record(header) + Record(changed));
        Check(rows[0].event == -1 && !rows[0].has_truth_match && rows[0].vtx_t == -1 &&
              rows[0].true_primary_neutron_count == -1, "Unavailable truth sentinels");
        changed = sample; changed[2] = "3";
        rows = parse(Record(header) + Record(sample) + Record(changed));
        Check(rows.size() == 2, "Multiple selected interactions in one event");
        changed = sample; changed[0] = "b.root";
        Check(parse(Record(header) + Record(sample) + Record(changed)).size() == 2,
              "Same event ID in different files");

        Throws([&] { parse(Record(header) + Record(sample) + Record(sample)); }, "Duplicate");
        Throws([&] { parse(Record(header)); }, "No selected interactions");
        Throws([&] { parse(""); }, "Missing header");
        auto bad_header = header; bad_header.pop_back();
        Throws([&] { parse(Record(bad_header)); }, "Missing column");
        bad_header = header; bad_header[0] = "event";
        Throws([&] { parse(Record(bad_header)); }, "duplicate header");
        Throws([&] { parse(Record(header) + "a.root,1\n"); }, "Field count");
        Throws([&] { parse(Record(header) + "\"unterminated"); }, "Unterminated");
        Throws([&] { parse(Record(header) + "\"a.root\"junk,1\n"); }, "Invalid character");
        changed = sample; changed[1] = "26junk";
        Throws([&] { parse(Record(header) + Record(changed)); }, "Column 'event'");
        changed[1] = "9223372036854775808";
        Throws([&] { parse(Record(header) + Record(changed)); }, "invalid integer");
        changed = sample; changed[2] = "2147483648";
        Throws([&] { parse(Record(header) + Record(changed)); }, "Column 'reco_ixn'");
        changed[2] = "-1";
        Throws([&] { parse(Record(header) + Record(changed)); }, "must be nonnegative");
        changed = sample; changed[4] = "yes";
        Throws([&] { parse(Record(header) + Record(changed)); }, "invalid Boolean");
        changed = sample; changed[3] = "-1";
        Throws([&] { parse(Record(header) + Record(changed)); }, "inconsistent with has_truth_match");
        changed = {"a.root", "-1", "0", "-1", "0", "1", "-1", "0", "-1", "-1"};
        Throws([&] { parse(Record(header) + Record(changed)); }, "requires has_truth_match");
        changed[5] = "0"; changed[8] = "0";
        Throws([&] { parse(Record(header) + Record(changed)); }, "true_primary_neutron_count");
        changed = sample; changed[0] = "";
        Throws([&] { parse(Record(header) + ",26,2,3,1,0,7,2,5,0\n"); }, "empty path");
        Throws([&] { ReadSelectionCsv((temporary.parent_path() / "absent.csv").string()); }, "Cannot open");

        rows = parse(Record(header) + Record(sample));
        auto row = rows[0];
        EventEntryIndex index{{26, {8}}, {78, {0}}, {-1, {12}}};
        Check(ResolveEventEntry(index, row) == 8, "Metadata event must not be used as entry number");
        row.vtx_r = 3;
        Check(ResolveEventEntry(index, row) == 8, "Multiple interactions share the selected entry");
        row.event = -1;
        Check(ResolveEventEntry(index, row) == 12, "Unique negative event ID");
        index[-1].push_back(13);
        Throws([&] { ResolveEventEntry(index, row); }, "ambiguous event ID");
        row.event = 26; index[26].push_back(9);
        Throws([&] { ResolveEventEntry(index, row); }, "ambiguous event ID");
        row.event = 999;
        Throws([&] { ResolveEventEntry(index, row); }, "event not found");

        std::cout << "PASS: " << checks << " parser and event-resolution checks\n";
        if (argc == 3) {
            rows = ReadSelectionCsv(argv[2]);
            std::set<std::string> files;
            std::size_t unmatched = 0, backgrounds = 0, signals = 0, negative = 0;
            for (const auto& item : rows) {
                files.insert(item.file_name);
                unmatched += !item.has_truth_match;
                backgrounds += item.has_truth_match && !item.matched_true_signal;
                signals += item.matched_true_signal;
                negative += item.event < 0;
            }
            std::cout << "CSV: rows=" << rows.size() << " files=" << files.size()
                      << " unmatched=" << unmatched << " matched_backgrounds=" << backgrounds
                      << " signals=" << signals << " negative_events=" << negative << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
