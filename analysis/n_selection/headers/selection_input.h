#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <map>
#include <string>
#include <vector>

// Selection metadata is independent of ROOT so CSV and event resolution can
// be validated without loading CAF dictionaries.
struct InputCAFRow {
    std::string file_name;
    std::int64_t event = -1;
    int vtx_r = -1;
    int vtx_t = -1;
    bool has_truth_match = false;
    bool matched_true_signal = false;
    int true_track_multiplicity = -1;
    int reco_track_multiplicity = -1;
    int true_primary_neutron_count = -1;
    int true_secondary_neutron_count = -1;
    std::size_t csv_record = 0;
};

using EventEntryIndex = std::map<std::int64_t, std::vector<std::int64_t>>;

bool ReadCsvRecord(std::istream& input, std::vector<std::string>& fields);
std::vector<InputCAFRow> ReadSelectionCsv(const std::string& path);
std::string CsvField(const std::string& value);
std::string SelectionRowContext(const InputCAFRow& row);
std::int64_t ResolveEventEntry(const EventEntryIndex& index,
                               const InputCAFRow& row);
