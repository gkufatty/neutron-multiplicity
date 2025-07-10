#include "cuts.h"
#include "structures.h" 
#include "utils.h"  // Include the header file for utils
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this
#include <cassert>
#include <cstring>
#include <utility> // for std::pair

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
    bool neutron_induced = false;
    const auto& nu = sr->mc.nu[vtx_idx];
    const int parent_id = nu.sec[part_idx].parent;
    for(int i =0; i<nu.prim.size(); i++){
        //std::cout << "prim pdg" << nu.prim[i].pdg  << std::endl;
        if (nu.prim[i].G4ID == parent_id){
            //std::cout << "parent pdg" << nu.prim[i].pdg  << std::endl;
            if(nu.prim[i].pdg == 2112){
                neutron_induced = true;
            }
        }
    }
    return neutron_induced;
}


PartBestMatch FindParticleBestMatch(
    //this returns pdg, type, vtx, and index in true 
    const std::vector<caf::TrueParticleID>& overlaping_particles, 
    const std::vector<float>& overlaps,
    const caf::StandardRecord* sr,
    const char* mode){
    PartBestMatch bm{};
    float max_overlap = -1.0f;
    
    for(size_t n_part =0; n_part < overlaping_particles.size(); n_part++){
        float current_overlap = overlaps[n_part];
        int truth_ixn_idx = overlaping_particles[n_part].ixn;
        int truth_part_idx = overlaping_particles[n_part].part;
        auto truth_type = overlaping_particles[n_part].type;
        if(current_overlap>max_overlap){
            max_overlap = current_overlap;
            bm.interaction_idx = truth_ixn_idx;
            bm.particle_idx = truth_part_idx;
            bm.type = truth_type;
        }
    }

    if (bm.type == 1 || bm.type == 3) {
        const auto& nu = sr->mc.nu[bm.interaction_idx];
        const auto& part = (bm.type == 3) ? nu.sec[bm.particle_idx]
            : nu.prim[bm.particle_idx];
        bm.pdg = part.pdg;
        bm.time = part.time;
        bm.start = part.start_pos;
        bm.end = part.end_pos;
        bm.parent = part.parent;
        bm.energy = part.p.E*1000 - ParticleCode::proton_mass;;
        bm.length = (bm.end - bm.start).Mag();
        bm.distance = (bm.start - nu.vtx).Mag();
        // bm.tpc = determineTPC(bm.start);
        // bm.vertex_tpc = determineTPC(nu.vtx);  
        // Conditionally set neutron_induced
        bm.neutron_induced = (bm.type == 3) ? was_neutron_induced(bm.interaction_idx, bm.particle_idx, sr) : false;
        bm.in_signal = QELikeSignal(sr, bm.interaction_idx, mode);
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
        if (bm_info.interaction_idx == -1) continue; // Skip if no best match found 
        

        if(p.primary&& bm_info.neutron_induced){
            stats.missed_sec++;
            bool LightCandidate = isolated_part(vtx_r, reco, sr,p.start,2.0f);
            // Check if the proton is a light candidate
            // if(LightCandidate){
            //     stats.light_candidate++;
            // }
            //do light check
        }
        if(!p.primary){    
            RecoProtonInfo proton;
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
            if((proton.reco_pdg==bm_info.pdg) && (proton.reco_type==bm_info.type)){
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
