#pragma once

#include <functional>
#include <cmath>
#include <string>
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include <vector>
#include <utility>
#include "TVector3.h"

int FindVertexBestMatch(
    const caf::SRVector3D& vertex,
    const std::vector<std::size_t>& vtx_overlaps_idx,
    const std::vector<float>& vtx_overlaps,
    const caf::StandardRecord* sr
);

PartBestMatch FindParticleBestMatch(
    const int reco_pdg, const int reco_type,
    const std::vector<caf::TrueParticleID>& overlapping_particles,
    const std::vector<float>& overlaps,
    const caf::StandardRecord* sr, bool verbose
);

bool was_neutron_induced(int vtx_idx, int part_idx, const caf::StandardRecord* sr);



// true analysis function
void process_truth_interactions(
    const caf::StandardRecord* sr, 
    Counters& stats, 
    int tree_entry, 
    const char* file_name, 
    const char* mode = "RHC"
);

std::pair<bool,int> process_reco_interaction(
    caf::StandardRecord* sr, 
    int i
);

std::vector<TrueVertexSelection_np> look_sec_particles(
    const int vtx_idx, 
    caf::StandardRecord* sr,
    Counters& stats);