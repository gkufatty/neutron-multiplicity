#include "cuts.h"
#include "structures.h" 
#include "utils.h"  // Include the header file for utils
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TTree.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <utility> // for std::pair
#include <stdexcept>

namespace {
float TruthKineticEnergyMeV(const caf::SRTrueParticle& part) {
    if (!std::isfinite(part.p.E))
        return std::numeric_limits<float>::quiet_NaN();

    const TParticlePDG* particle =
        TDatabasePDG::Instance()->GetParticle(part.pdg);
    if (!particle)
        return std::numeric_limits<float>::quiet_NaN();

    return static_cast<float>((part.p.E - particle->Mass()) * 1000.0);
}
} // namespace

EventEntryIndex BuildEventEntryIndex(TTree& tree, caf::StandardRecord*& sr) {
    EventEntryIndex index;
    const Long64_t entries = tree.GetEntries();
    if (entries < 0) throw std::runtime_error("Cannot count cafTree entries");
    for (Long64_t entry = 0; entry < entries; ++entry) {
        if (tree.GetEntry(entry) <= 0 || !sr)
            throw std::runtime_error("Cannot read cafTree entry " + std::to_string(entry));
        index[sr->meta.nd_lar.event].push_back(entry);
    }
    return index;
}

void AttachInputMetadata(RecoProtonInfo& proton, const InputCAFRow& row) {
    proton.input_file = row.file_name;
    proton.input_event = row.event;
    proton.input_vtx_t = row.vtx_t;
    proton.input_matched = row.matched_true_signal;
    proton.input_has_truth_match = row.has_truth_match;
    proton.input_true_track_multiplicity = row.true_track_multiplicity;
    proton.input_reco_track_multiplicity = row.reco_track_multiplicity;
    proton.input_true_primary_neutron_count = row.true_primary_neutron_count;
    proton.input_true_secondary_neutron_count = row.true_secondary_neutron_count;
}

// ======================
// FUNCTION IMPLEMENTATIONS
// =======================



int FindVertexBestMatch(
    const caf::SRVector3D& vertex,
    const std::vector<std::size_t>& vtx_overlaps_idx,
    const std::vector<float>& vtx_overlaps,
    const caf::StandardRecord* sr){
    int best_match_idx = -1;
    float max_overlap = -1;
    
    for(size_t i = 0; i < vtx_overlaps.size(); i++) {
        float current_overlap = vtx_overlaps[i];
        size_t current_idx = vtx_overlaps_idx[i];
        
        if(current_overlap > max_overlap) {
            max_overlap = current_overlap;
            best_match_idx = current_idx;
        }
        
        if(best_match_idx != -1) {
            auto best_match_vertex = sr->mc.nu[best_match_idx].vtx;
            float delta_x = abs(best_match_vertex.x - vertex.x);
            float delta_y = abs(best_match_vertex.y - vertex.y);
            float delta_z = abs(best_match_vertex.z - vertex.z);
            
            if(!(delta_x < 5 && delta_y < 5 && delta_z < 5)) {
                best_match_idx = -1;
            }
        }
    }
    return best_match_idx;
}




bool was_neutron_induced(int vtx_idx, int part_idx, const caf::StandardRecord* sr){
    const auto& nu = sr->mc.nu[vtx_idx];
    const int parent_id = nu.sec[part_idx].parent;

    // GEANT4 track IDs are positive. A parent ID of zero means no parent,
    // while negative IDs represent unavailable ancestry information.
    if (parent_id <= 0) {
        for (const auto& primary : nu.prim) {
            if (primary.G4ID == parent_id &&
                primary.pdg == ParticleCode::neutron) {
                std::cerr << "[DEBUG] Avoided bug 1: rejected invalid neutron-parent "
                          << "ID match (interaction=" << vtx_idx
                          << ", secondary=" << part_idx
                          << ", parent_id=" << parent_id << ")\n";
                break;
            }
        }
        return false;
    }

    for (const auto& primary : nu.prim) {
        if (primary.G4ID <= 0) continue;
        if (primary.G4ID == parent_id)
            return primary.pdg == ParticleCode::neutron;
    }
    return false;
}


PartBestMatch FindParticleBestMatch(
    //this returns pdg, type, vtx, and index in true 
    const std::vector<caf::TrueParticleID>& overlaping_particles, 
    const std::vector<float>& overlaps,
    const caf::StandardRecord* sr,
    const char* mode){
    PartBestMatch bm{};
    float max_overlap = 0.0f;
    int best = -1;
    int first_zero_overlap = -1;
    if (overlaping_particles.size() != overlaps.size())
        throw std::runtime_error("Particle truth and truthOverlap lengths differ");
    
    for(size_t n_part =0; n_part < overlaping_particles.size(); n_part++){
        float current_overlap = overlaps[n_part];
        if (std::isfinite(current_overlap) && current_overlap == 0.0f &&
            first_zero_overlap < 0)
            first_zero_overlap = static_cast<int>(n_part);
        if(std::isfinite(current_overlap) && current_overlap > max_overlap){
            max_overlap = current_overlap;
            best = static_cast<int>(n_part);
        }
    }

    if (best < 0) {
        // Report only a zero-overlap association that the old code would have
        // promoted all the way to a valid match.
        if (first_zero_overlap >= 0) {
            const auto& zero_id = overlaping_particles[first_zero_overlap];
            const bool supported = zero_id.type == 1 || zero_id.type == 3;
            const bool valid_interaction = zero_id.ixn >= 0 &&
                static_cast<std::size_t>(zero_id.ixn) < sr->mc.nu.size();
            bool valid_particle = false;
            if (supported && valid_interaction && zero_id.part >= 0) {
                const auto& zero_nu = sr->mc.nu[zero_id.ixn];
                const auto count = zero_id.type == 3 ? zero_nu.sec.size()
                                                     : zero_nu.prim.size();
                valid_particle = static_cast<std::size_t>(zero_id.part) < count;
            }
            if (supported && valid_interaction && valid_particle) {
                std::cerr << "[DEBUG] Avoided bug 2: rejected zero-overlap "
                          << "particle truth association (interaction="
                          << zero_id.ixn << ", particle=" << zero_id.part
                          << ", type=" << zero_id.type << ")\n";
            }
        }
        return bm;
    }
    const auto& id = overlaping_particles[best];
    // Other CAF truth types are unavailable to this neutrino-particle analysis.
    if (id.type != 1 && id.type != 3) return bm;
    if (static_cast<std::size_t>(id.ixn) >= sr->mc.nu.size())
        throw std::runtime_error("Particle truth interaction index out of bounds");
    const auto& nu = sr->mc.nu[id.ixn];
    const auto count = id.type == 3 ? nu.sec.size() : nu.prim.size();
    if (static_cast<std::size_t>(id.part) >= count)
        throw std::runtime_error("Particle truth particle index out of bounds");
    {
        bm.valid = true;
        bm.overlap = max_overlap;
        bm.interaction_idx = id.ixn;
        bm.particle_idx = id.part;
        bm.type = id.type;
        const auto& part = (bm.type == 3) ? nu.sec[bm.particle_idx]
            : nu.prim[bm.particle_idx];
        bm.pdg = part.pdg;
        bm.time = part.time;
        bm.start = part.start_pos;
        bm.end = part.end_pos;
        bm.parent = part.parent;
        bm.energy = TruthKineticEnergyMeV(part);
        bm.length = (bm.end - bm.start).Mag();
        bm.distance = (bm.start - nu.vtx).Mag();
        // bm.tpc = determineTPC(bm.start);
        // bm.vertex_tpc = determineTPC(nu.vtx);  
        // Conditionally set neutron_induced
        bm.neutron_induced = (bm.type == 3) ? was_neutron_induced(bm.interaction_idx, bm.particle_idx, sr) : false;
        // QELikeSignal's muon helpers access prim[0].
        bm.in_signal = !nu.prim.empty() && QELikeSignal(sr, bm.interaction_idx, mode);
    }
    return bm;
}


std::vector<RecoProtonInfo> process_protons(
    int vtx_r, 
    const caf::SRInteraction& reco, 
    caf::StandardRecord* sr,
    Counters& stats, 
    const char* mode) {
    caf::SRVector3D muon_start, muon_end;
    // Extract muon track from primary muon
    for (const auto& p : reco.part.dlp) {
        if (p.pdg == ParticleCode::muon && p.primary) {
            muon_start = p.start;
            muon_end = p.end;
            break;  // Use the first primary muon
        }
    }
    std::vector<RecoProtonInfo> reco_protons;
    for (const auto& p : reco.part.dlp) {
        if (p.pdg != ParticleCode::proton) continue; // only protons
        if(!(kIsVtxFV(p.start)) || !(kIsVtxFV(p.end))) continue; // Skip particles not fully contained in the FV   
        // if (!(around_vtx_region(p.start, reco.vtx, 0.5))) continue; // Skip particles within the vertex region
        // if (!around_track_region(p.start, muon_start, muon_end, 0.5)) continue;

        PartBestMatch bm_info = FindParticleBestMatch(p.truth, p.truthOverlap, sr, mode);
        // Keep reconstructed candidates without usable truth; bm.valid tells
        // downstream analyses whether truth quantities and flags are available.
        

        if(p.primary && bm_info.valid && bm_info.neutron_induced){
            stats.missed_sec++;
            bool LightCandidate = isolated_part(vtx_r, reco, sr,p.start,2.0f);
            // Check if the proton is a light candidate
            // if(LightCandidate){
            //     stats.light_candidate++;
            // }
            //do light check
        }
        if(!p.primary){    
            RecoProtonInfo proton{};
            stats.secp_reco++;        
            //reco information  
            proton.reco_pdg = p.pdg;
            proton.reco_type = 3; // secondary
            proton.reco_energy = p.E*1000;
            proton.reco_vtx = vtx_r;
            proton.reco_start = p.start;
            proton.reco_end = p.end;
            proton.reco_len = (p.start - p.end).Mag();
            proton.reco_dist = (p.start - reco.vtx).Mag();
            //bm information
            proton.bm = bm_info;  // direct assignment
            if(bm_info.valid && (proton.reco_pdg==bm_info.pdg) && (proton.reco_type==bm_info.type)){
                proton.coincidence = true;
                stats.matched_secp++;
            }
            else{
                proton.coincidence = false;
            }
            // Update statistics
            if(bm_info.pdg == proton.reco_pdg) stats.matched_secp_pdg++;
            if(bm_info.type == proton.reco_type) stats.matched_secp_type++;
            if(proton.coincidence && bm_info.neutron_induced){stats.np_reco++;}
            // Save proton info
            reco_protons.push_back(proton);  
        }
    }
    return reco_protons;
}


// bool LightCandidate(){
//     bool different_tpc = true;
//     if(isolated_part&& different_tpc){
//         return true;}
//     else{
//         return false;
//     }
// }
   

bool isolated_part(
    int vtx_r, 
    const caf::SRInteraction& reco, 
    caf::StandardRecord* sr,
    const caf::SRVector3D reco_start,
    float radius)
{
    for (const auto& p : reco.part.dlp){
        // Check if the start position is valid
        if (p.start.x != -999 && p.start.y != -999 && p.start.z != -999) {
            float dx = p.start.x - reco_start.x;
            float dy = p.start.y - reco_start.y;
            float dz = p.start.z - reco_start.z;
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (dist < radius) {
                return false;  // Found another particle too close
            }
        }
    }
    return true;  // No particle found within the radius
}
