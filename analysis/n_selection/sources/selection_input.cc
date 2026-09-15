#include "selection_input.h"

#include <array>
#include <charconv>
#include <fstream>
#include <set>
#include <stdexcept>
#include <tuple>

bool ReadCsvRecord(std::istream& input, std::vector<std::string>& fields) {
    fields.clear();
    std::string field;
    bool quoted = false;
    bool closed_quote = false;
    bool started = false;
    char ch;
    while (input.get(ch)) {
        started = true;
        if (quoted) {
            if (ch == '"') {
                if (input.peek() == '"') {
                    input.get(ch);
                    field += '"';
                } else {
                    quoted = false;
                    closed_quote = true;
                }
            } else {
                field += ch;
            }
        } else if (ch == ',') {
            fields.push_back(field);
            field.clear();
            closed_quote = false;
        } else if (ch == '\n' || ch == '\r') {
            if (ch == '\r' && input.peek() == '\n') input.get(ch);
            fields.push_back(field);
            return true;
        } else if (ch == '"' && field.empty() && !closed_quote) {
            quoted = true;
        } else {
            if (closed_quote || ch == '"')
                throw std::runtime_error("Invalid character in CSV field " +
                                         std::to_string(fields.size() + 1));
            field += ch;
        }
    }
    if (input.bad() || (!input.eof() && input.fail()))
        throw std::runtime_error("CSV read failure");
    if (quoted) throw std::runtime_error("Unterminated quoted CSV field " +
                                        std::to_string(fields.size() + 1));
    if (!started) return false;
    fields.push_back(field);
    return true;
}

namespace {
template <typename T>
T ParseInteger(const std::string& value, const std::string& column) {
    T result{};
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    if (value.empty() || parsed.ec != std::errc{} ||
        parsed.ptr != value.data() + value.size())
        throw std::runtime_error("Column '" + column + "': invalid integer '" + value + "'");
    return result;
}

bool ParseBoolean(const std::string& value, const std::string& column) {
    if (value == "1" || value == "true") return true;
    if (value == "0" || value == "false") return false;
    throw std::runtime_error("Column '" + column + "': invalid Boolean '" + value + "'");
}
} // namespace

std::vector<InputCAFRow> ReadSelectionCsv(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open selection CSV: " + path);
    std::vector<InputCAFRow> rows;
    std::vector<std::string> fields;
    std::size_t record = 1;
    try {
        if (!ReadCsvRecord(input, fields)) throw std::runtime_error("Missing header");
        // Accept a UTF-8 BOM emitted by spreadsheet applications.
        if (!fields.empty() && fields.front().compare(0, 3, "\xEF\xBB\xBF") == 0)
            fields.front().erase(0, 3);
        std::map<std::string, std::size_t> columns;
        const auto width = fields.size();
        for (std::size_t i = 0; i < width; ++i) {
            if (fields[i].empty() || !columns.emplace(fields[i], i).second)
                throw std::runtime_error("Empty or duplicate header: '" + fields[i] + "'");
        }
        const std::array<const char*, 10> required = {{
            "file_name", "event", "reco_ixn", "truth_ixn", "has_truth_match",
            "matched_true_signal", "true_track_multiplicity", "reco_track_multiplicity",
            "true_primary_neutron_count", "true_secondary_neutron_count"
        }};
        for (const auto name : required)
            if (!columns.count(name)) throw std::runtime_error("Missing column '" + std::string(name) + "'");

        std::set<std::tuple<std::string, std::int64_t, int>> keys;
        while (++record, ReadCsvRecord(input, fields)) {
            if (fields.size() != width) throw std::runtime_error("Field count does not match header");
            auto field = [&](const std::string& name) -> const std::string& {
                return fields.at(columns.at(name));
            };
            auto integer = [&](const std::string& name) {
                return ParseInteger<int>(field(name), name);
            };
            InputCAFRow row;
            row.csv_record = record;
            row.file_name = field("file_name");
            row.event = ParseInteger<std::int64_t>(field("event"), "event");
            row.vtx_r = integer("reco_ixn");
            row.vtx_t = integer("truth_ixn");
            row.has_truth_match = ParseBoolean(field("has_truth_match"), "has_truth_match");
            row.matched_true_signal = ParseBoolean(field("matched_true_signal"), "matched_true_signal");
            row.true_track_multiplicity = integer("true_track_multiplicity");
            row.reco_track_multiplicity = integer("reco_track_multiplicity");
            row.true_primary_neutron_count = integer("true_primary_neutron_count");
            row.true_secondary_neutron_count = integer("true_secondary_neutron_count");
            if (row.file_name.empty()) throw std::runtime_error("Column 'file_name': empty path");
            if (row.vtx_r < 0) throw std::runtime_error("Column 'reco_ixn': must be nonnegative");
            if (row.reco_track_multiplicity < 0)
                throw std::runtime_error("Column 'reco_track_multiplicity': must be nonnegative");
            for (const auto name : {"truth_ixn", "true_track_multiplicity",
                                   "true_primary_neutron_count", "true_secondary_neutron_count"}) {
                const int value = integer(name);
                if ((row.has_truth_match && value < 0) || (!row.has_truth_match && value != -1))
                    throw std::runtime_error("Column '" + std::string(name) + "': inconsistent with has_truth_match");
            }
            if (row.matched_true_signal && !row.has_truth_match)
                throw std::runtime_error("Column 'matched_true_signal': requires has_truth_match");
            if (!keys.emplace(row.file_name, row.event, row.vtx_r).second)
                throw std::runtime_error("Duplicate (file_name, event, reco_ixn)");
            rows.push_back(row);
        }
        if (rows.empty()) throw std::runtime_error("No selected interactions");
    } catch (const std::exception& error) {
        throw std::runtime_error(path + ": CSV record " + std::to_string(record) + ": " + error.what());
    }
    return rows;
}

std::string CsvField(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) return value;
    std::string encoded = "\"";
    for (const char ch : value) {
        if (ch == '"') encoded += '"';
        encoded += ch;
    }
    return encoded + '"';
}

std::string SelectionRowContext(const InputCAFRow& row) {
    return row.file_name + ": CSV record " + std::to_string(row.csv_record) +
           ", event " + std::to_string(row.event) + ", reco_ixn " + std::to_string(row.vtx_r);
}

std::int64_t ResolveEventEntry(const EventEntryIndex& index, const InputCAFRow& row) {
    const auto found = index.find(row.event);
    if (found == index.end() || found->second.empty())
        throw std::runtime_error(SelectionRowContext(row) + ": event not found in cafTree");
    if (found->second.size() != 1)
        throw std::runtime_error(SelectionRowContext(row) +
            ": ambiguous event ID; export a local CAF entry number from the CSV producer "
            "and extend the reader to use that identifier");
    return found->second.front();
}
