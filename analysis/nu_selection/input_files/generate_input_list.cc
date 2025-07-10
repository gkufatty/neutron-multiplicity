#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>  // Requires C++17
#include <algorithm>

namespace fs = std::filesystem;

std::string remove_quotes(const std::string& s) {
    std::string result = s;
    result.erase(std::remove(result.begin(), result.end(), '\"'), result.end());
    return result;
}

int main(int argc, char** argv) {
    
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <folder_path> <version> <mode>\n";
        std::cerr << "Example: " << argv[0] << " /exp/dune/data/users/noeroy/prod/MiniRun6.2_1E19_RHC/MiniRun6.2_1E19_RHC.caf/CAF/ 6.2 RHC\n";
        return 1;
    }

    std::string folder_path = argv[1];
    std::string version = argv[2];
    std::string mode = argv[3];
    std::cout << "Generating input file list from all the folders in " << folder_path << std::endl;

    std::string out_filename = "MiniRun" + version + "_" + mode + ".txt";
    std::ofstream outfile(out_filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file: " << out_filename << "\n";
        return 1;
    }

    int file_count = 0;
    for (const auto& entry : fs::recursive_directory_iterator(folder_path)) {
        if (entry.is_regular_file()) {
            std::string abs_path = fs::absolute(entry.path()).string();
            outfile << remove_quotes(abs_path) << "\n";
            file_count++;
        }
    }

    outfile.close();
    std::cout << "Wrote " << file_count << " file paths to " << out_filename << "\n";
    return 0;
}

