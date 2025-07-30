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


PartBestMatch FindParticleBestMatch(
    const int reco_pdg,
    const int reco_type,
    const std::vector<caf::TrueParticleID>& overlaping_particles, 
    const std::vector<float>& overlaps,
    const caf::StandardRecord* sr, 
    bool verbose){   
    PartBestMatch bm{};

    bm.coincidence = false;
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

    if (bm.type == 1) {  // Primary particle
        bm.pdg = sr->mc.nu[bm.interaction_idx].prim[bm.particle_idx].pdg;
        bm.time = sr->mc.nu[bm.interaction_idx].prim[bm.particle_idx].time;
        bm.start = sr->mc.nu[bm.interaction_idx].prim[bm.particle_idx].start_pos;
        bm.end = sr->mc.nu[bm.interaction_idx].prim[bm.particle_idx].end_pos;
        bm.parent = sr->mc.nu[bm.interaction_idx].prim[bm.particle_idx].parent;
        bm.energy= sr->mc.nu[bm.interaction_idx].prim[bm.particle_idx].p.E;
        bm.length = (bm.end-bm.start).Mag();
        bm.distance = (bm.start-sr->mc.nu[bm.interaction_idx].vtx).Mag();
    } 
    else if (bm.type == 3) {  // Secondary particle
        bm.pdg = sr->mc.nu[bm.interaction_idx].sec[bm.particle_idx].pdg;
        bm.time = sr->mc.nu[bm.interaction_idx].sec[bm.particle_idx].time;
        bm.start = sr->mc.nu[bm.interaction_idx].sec[bm.particle_idx].start_pos;
        bm.end = sr->mc.nu[bm.interaction_idx].sec[bm.particle_idx].end_pos;
        bm.parent = sr->mc.nu[bm.interaction_idx].sec[bm.particle_idx].parent;
        bm.energy= sr->mc.nu[bm.interaction_idx].sec[bm.particle_idx].p.E;
        bm.length = (bm.end-bm.start).Mag();
        bm.distance = (bm.start-sr->mc.nu[bm.interaction_idx].vtx).Mag();
    }
    //std::cout << bm.coincidence << std::endl;
    if((reco_pdg==bm.pdg) && (reco_type==bm.type)){
        bm.coincidence = true;
        //std::cout << bm.coincidence << std::endl;
    }
    return bm;
}

void process_truth_interactions(
    const caf::StandardRecord* sr, 
    Counters& stats, 
    int tree_entry, 
    const char* file_name, 
    const char* mode){
    for (size_t i = 0; i < sr->mc.nu.size(); i++) {
        if (QELikeSignal(sr, i, mode)) {
            std::cout << "Found a QELike signal in event " << sr->meta.nd_lar.event << std::endl;
            stats.nu_true++;
            int mode = sr->mc.nu[i].mode;
            switch(mode){
                case 1:stats.QES++; break;
                case 10:stats.MEC++; break;
                case 4:stats.RES++; break;
            }
        }
    }
}


// === Reco Analysis Functions ===

std::pair<bool,int> process_reco_interaction(
    caf::StandardRecord* sr, int i, std::string& reco_step){
    auto vertex = sr->common.ixn.dlp[i].vtx;
    auto vtx_overlaps = sr->common.ixn.dlp[i].truthOverlap;
    auto vtx_overlaps_idx = sr->common.ixn.dlp[i].truth; //index of the true interaction
    const auto& reco = sr->common.ixn.dlp[i];
    // Vertex matching
    int best_match_idx = FindVertexBestMatch(vertex, vtx_overlaps_idx, vtx_overlaps, sr);
    if (best_match_idx == -1) return {false, -1};
    if(reco_step=="All" && best_match_idx != -1){
        return {true,best_match_idx};
    }
    else if(reco_step=="FV"){
        return{kIsVtxFV(sr->common.ixn.dlp[i].vtx), best_match_idx};
    }
    else if(reco_step=="2x2"){
        return{RecoIsExitingMu(sr,i), best_match_idx};

    }
    else if(reco_step=="Mx2"){
        return{RecoIsBasicCCNuAr(sr,i), best_match_idx};
    }

    else if(reco_step=="OnePrim"){
        return{RecoHasOnePrimary(sr,i), best_match_idx};
    }

    else if(reco_step=="QE-like"){
        return {RecoIsQELikeSignal(sr, i), best_match_idx};
    }
}


std::vector<TrueVertexSelection_np> look_sec_particles(
    const int vtx_idx, 
    caf::StandardRecord* sr,
    Counters& stats) {
    const auto& nu = sr->mc.nu[vtx_idx];
    const auto& vtx_pos = nu.vtx;
    std::vector<TrueVertexSelection_np> result{}; 
    //iterate over the secondary particles of the vertex
    for (size_t j = 0; j < nu.sec.size(); ++j) {
        std::vector<int> n_induced_pdgs;

        bool n_induced = false;
        int pdg =nu.sec[j].pdg;
        const caf::SRVector3D start = nu.sec[j].start_pos;
        const caf::SRVector3D end = nu.sec[j].end_pos;
        //skip particles not fully contained in the FV
        if(!( kIsVtxFV(start)) || !(kIsVtxFV(end))) continue;

        if(was_neutron_induced(vtx_idx, j, sr)){
            n_induced_pdgs.push_back(pdg);
            n_induced = true;
        }


        if(pdg == ParticleCode::proton){
            stats.secp_true++;
            if(n_induced){
                stats.np_true++;
                TrueVertexSelection_np vtx_np{};
                vtx_np.vtx = vtx_idx;
                vtx_np.insignal = true;
                vtx_np.neutron_t0 = nu.sec[j].time;
                vtx_np.proton_t0 = nu.sec[j].time;
                vtx_np.proton_energy = nu.sec[j].p.E*1000 - ParticleCode::proton_mass;
                vtx_np.proton_length = (nu.sec[j].start_pos - nu.sec[j].end_pos).Mag();
                vtx_np.proton_distance = (nu.sec[j].start_pos - vtx_pos).Mag();
                result.push_back(vtx_np);
            }
        }
    }
    return result;
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


void print_stats(Counters& stats){
        // Print stats
    std::cout << "Total true QE-like interactions: " << stats.nu_true << std::endl;
    std::cout << "Total reconstructed QE-like interactions: " << stats.nu_reco << std::endl;
    std::cout << "Total matched QE-like interactions: " << stats.matched << std::endl;  
    std::cout << "Num of true neutron induced protons: " << stats.np_true << std::endl;
    std::cout << "Num of true secondary protons: " << stats.secp_true << std::endl;

    float nu_efficiency = static_cast<double>(stats.matched) / stats.nu_true;
    float nu_purity = static_cast<double>(stats.matched) / stats.nu_reco;
    std::cout << "Nu efficiency: " << nu_efficiency << std::endl;
    std::cout << "Nu purity: " << nu_purity << std::endl;
    std::cout << "There are " << stats.nu_true << " true QE-like interactions." << std::endl;
    std::cout << "There are " << stats.nu_reco << " reconstructed QE-like interactions with " << stats.matched << " matching true signal definition."<< std::endl;

}