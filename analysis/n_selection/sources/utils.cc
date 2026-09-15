#include "cuts.h"
#include "structures.h" 
#include "utils.h"  // Include the header file for utils
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TTree.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <utility> // for std::pair
#include <stdexcept>

namespace {
struct DirectParentMatch {
    int pdg = 0;
    bool induced = false;
    int type = -1;
    int particle_idx = -1;
    int g4id = -1;
};

float TruthKineticEnergyMeV(const caf::SRTrueParticle& part) {
    if (!std::isfinite(part.p.E))
        return std::numeric_limits<float>::quiet_NaN();

    const TParticlePDG* particle =
        TDatabasePDG::Instance()->GetParticle(part.pdg);
    if (!particle)
        return std::numeric_limits<float>::quiet_NaN();

    return static_cast<float>((part.p.E - particle->Mass()) * 1000.0);
}

DirectParentMatch FindDirectParent(
    int vtx_idx, int part_idx, const caf::StandardRecord* sr) {
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
        return {};
    }

    DirectParentMatch match;
    const auto inspect = [&](const auto& particles, int type) {
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const auto& candidate = particles[i];
            if (candidate.G4ID <= 0 || candidate.G4ID != parent_id) continue;

            // G4 track IDs identify one particle. A resolved non-neutron
            // parent therefore ends the direct-parent search.
            match.pdg = candidate.pdg;
            if (candidate.pdg == ParticleCode::neutron) {
                match.induced = true;
                match.type = type;
                match.particle_idx = static_cast<int>(i);
                match.g4id = candidate.G4ID;
            }
            return true;
        }
        return false;
    };

    if (inspect(nu.prim, caf::TrueParticleID::kPrimary)) return match;
    inspect(nu.sec, caf::TrueParticleID::kSecondary);
    return match;
}

bool IsTruthParticle(const caf::TrueParticleID& id,
                     int interaction_idx, int particle_idx) {
    return id.type == caf::TrueParticleID::kSecondary &&
        id.ixn == interaction_idx && id.part == particle_idx;
}

bool IsSelectedRecoInteraction(int reco_vtx,
                               const std::vector<int>& selected_reco_vtxs) {
    return std::find(selected_reco_vtxs.begin(), selected_reco_vtxs.end(),
                     reco_vtx) != selected_reco_vtxs.end();
}

bool PassesRecoProtonSelection(const caf::SRRecoParticle& particle) {
    return particle.pdg == ParticleCode::proton && !particle.primary &&
        kIsVtxFV(particle.start) && kIsVtxFV(particle.end);
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
    proton.same_truth_interaction = row.has_truth_match && proton.bm.valid &&
        proton.bm.interaction_idx == row.vtx_t;
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
    return FindDirectParent(vtx_idx, part_idx, sr).induced;
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
        // A neutron-induced proton has a neutron as its direct parent. Preserve
        // whether that neutron is stored as a true primary or true secondary.
        if (bm.type == caf::TrueParticleID::kSecondary) {
            const auto direct_parent = FindDirectParent(
                bm.interaction_idx, bm.particle_idx, sr);
            bm.direct_parent_pdg = direct_parent.pdg;
            bm.neutron_induced = direct_parent.induced;
            bm.neutron_parent_type = direct_parent.type;
            bm.neutron_parent_idx = direct_parent.particle_idx;
            bm.neutron_parent_g4id = direct_parent.g4id;
        }
        // QELikeSignal's muon helpers access prim[0].
        bm.in_signal = !nu.prim.empty() && QELikeSignal(sr, bm.interaction_idx, mode);
    }
    return bm;
}


std::vector<TrueNeutronProtonInfo> process_true_neutron_protons(
    int truth_vtx,
    const std::vector<int>& selected_reco_vtxs,
    const caf::StandardRecord* sr) {
    if (truth_vtx < 0 || static_cast<std::size_t>(truth_vtx) >= sr->mc.nu.size())
        throw std::runtime_error("Truth interaction index out of bounds");

    const auto& nu = sr->mc.nu[truth_vtx];
    std::vector<TrueNeutronProtonInfo> true_protons;
    for (std::size_t part_idx = 0; part_idx < nu.sec.size(); ++part_idx) {
        const auto& part = nu.sec[part_idx];
        if (part.pdg != ParticleCode::proton) continue;

        const auto parent = FindDirectParent(
            truth_vtx, static_cast<int>(part_idx), sr);
        if (!parent.induced) continue;
        // The production point defines truth eligibility. Keep the end point
        // only as an output feature so exiting protons remain measurable.
        if (!kIsVtxFV(part.start_pos)) continue;

        TrueNeutronProtonInfo proton;
        proton.truth_vtx = truth_vtx;
        proton.truth_particle_idx = static_cast<int>(part_idx);
        proton.truth_energy = TruthKineticEnergyMeV(part);
        proton.truth_length = (part.end_pos - part.start_pos).Mag();
        proton.truth_distance = (part.start_pos - nu.vtx).Mag();
        proton.truth_time = part.time;
        proton.truth_vertex_time = nu.time;
        proton.truth_particle_dt = part.time - nu.time;
        proton.truth_start = part.start_pos;
        proton.truth_end = part.end_pos;
        proton.truth_end_in_fv = kIsVtxFV(part.end_pos);
        proton.truth_parent_g4id = part.parent;
        proton.direct_parent_pdg = parent.pdg;
        proton.neutron_parent_type = parent.type;
        proton.neutron_parent_idx = parent.particle_idx;
        proton.neutron_parent_g4id = parent.g4id;

        // Invert the reco-to-truth associations across every reconstructed
        // interaction in the event. This separates any contribution to a reco
        // object from being that object's actual best truth match.
        for (std::size_t reco_vtx = 0; reco_vtx < sr->common.ixn.dlp.size(); ++reco_vtx) {
            const auto& reco = sr->common.ixn.dlp[reco_vtx];
            for (std::size_t reco_idx = 0; reco_idx < reco.part.dlp.size(); ++reco_idx) {
                const auto& reco_part = reco.part.dlp[reco_idx];
                if (reco_part.truth.size() != reco_part.truthOverlap.size())
                    throw std::runtime_error(
                        "Particle truth and truthOverlap lengths differ");

                float largest_overlap = 0.0f;
                int largest_idx = -1;
                for (std::size_t assoc_idx = 0;
                     assoc_idx < reco_part.truth.size(); ++assoc_idx) {
                    const float overlap = reco_part.truthOverlap[assoc_idx];
                    if (!std::isfinite(overlap) || overlap <= 0.0f) continue;
                    if (IsTruthParticle(reco_part.truth[assoc_idx], truth_vtx,
                                        static_cast<int>(part_idx))) {
                        proton.has_any_reco_overlap = true;
                        if (!std::isfinite(proton.max_any_reco_overlap) ||
                            overlap > proton.max_any_reco_overlap)
                            proton.max_any_reco_overlap = overlap;
                    }
                    if (overlap > largest_overlap) {
                        largest_overlap = overlap;
                        largest_idx = static_cast<int>(assoc_idx);
                    }
                }

                if (largest_idx < 0 ||
                    !IsTruthParticle(reco_part.truth[largest_idx], truth_vtx,
                                     static_cast<int>(part_idx)))
                    continue;

                proton.has_reco_match = true;
                ++proton.n_reco_matches;
                const bool in_selected_interaction = IsSelectedRecoInteraction(
                    static_cast<int>(reco_vtx), selected_reco_vtxs);
                const bool passes_select_n = in_selected_interaction &&
                    PassesRecoProtonSelection(reco_part);
                proton.has_selected_candidate |= passes_select_n;

                if (!std::isfinite(proton.best_reco_overlap) ||
                    largest_overlap > proton.best_reco_overlap) {
                    proton.best_reco_vtx = static_cast<int>(reco_vtx);
                    proton.best_reco_particle_idx = static_cast<int>(reco_idx);
                    proton.best_reco_pdg = reco_part.pdg;
                    proton.best_reco_primary = reco_part.primary;
                    proton.best_reco_overlap = largest_overlap;
                    proton.best_reco_in_selected_interaction =
                        in_selected_interaction;
                    proton.best_reco_start_in_fv = kIsVtxFV(reco_part.start);
                    proton.best_reco_end_in_fv = kIsVtxFV(reco_part.end);
                    proton.best_reco_passes_select_n = passes_select_n;
                }
            }
        }
        true_protons.push_back(std::move(proton));
    }
    return true_protons;
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
        if (!(kIsVtxFV(p.start)) || !(kIsVtxFV(p.end))) continue; // Skip particles not fully contained in the FV
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
