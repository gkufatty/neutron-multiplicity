#pragma once

#include <functional>
#include <cmath>
#include <string>
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include <vector>
#include <string>
#include <utility>
#include "structures.h"  // Include the header file for structures
#include "TVector3.h"


TPC determineTPC(const caf::SRVector3D& pos);
std::string tpcToString(TPC tpc); 

int FindVertexBestMatch(
    const caf::SRVector3D& vertex,
    const std::vector<std::size_t>& vtx_overlaps_idx,
    const std::vector<float>& vtx_overlaps,
    const caf::StandardRecord* sr
);


bool was_neutron_induced(int vtx_idx, int part_idx, const caf::StandardRecord* sr);


PartBestMatch FindParticleBestMatch(
    //this returns pdg, type, vtx, and index in true 
    const std::vector<caf::TrueParticleID>& overlaping_particles, 
    const std::vector<float>& overlaps,
    const caf::StandardRecord* sr,
    const char* mode = "RHC");

std::vector<RecoProtonInfo> process_protons(
    int vtx_r, 
    const caf::SRInteraction& reco, 
    caf::StandardRecord* sr,
    Counters& stats,
    const char* mode = "RHC"
);

bool LightCandidate();
bool isolated_part(
    int vtx_r, 
    const caf::SRInteraction& reco, 
    caf::StandardRecord* sr,
    const caf::SRVector3D reco_start,
    float radius= 5.0f
);
