
#include "cuts.h"  // Include the header file for cuts
#include "structures.h" 
#include <cassert>
#include <cstring>
#include <utility> // for std::pair

bool kIsCC(const caf::StandardRecord* sr, int idx){
    return sr->mc.nu[idx].iscc;
}

bool kIsNumu(const caf::StandardRecord* sr, int idx){
    return (abs(sr->mc.nu[idx].pdg) == 14);
}

bool kIsAntiNumu(const caf::StandardRecord* sr,int idx){
    return (abs(sr->mc.nu[idx].pdg) == 14);
}

bool kIsVtxCont(const caf::StandardRecord* sr, int idx){
    return (sr->mc.nu[idx].isvtxcont);
}

bool kIsVtxFV(const caf::SRVector3D& vertex){
    //const auto& vertex = sr->mc.nu[0].vtx;  // Fix the missing vertex definition
    return (
        std::abs(vertex.x) < constants::X_BOUND - constants::FV_DIST &&
        std::abs(vertex.x) > constants::TPC_DIST &&
        std::abs(vertex.y) < constants::Y_BOUND - constants::FV_DIST &&
        std::abs(vertex.z) < constants::Z_BOUND - constants::FV_DIST &&
        std::abs(vertex.z) > constants::FV_DIST
    );
}

bool kIsArgon(const caf::StandardRecord* sr, int idx){
    return ((sr->mc.nu[idx].targetPDG) == ParticleCode::argon);
}

bool AngleCrossMx2(const caf::StandardRecord* sr, int idx) {
    const auto& nu = sr->mc.nu[idx];
    const auto& p = nu.prim[0];
    if(abs(p.pdg)!= ParticleCode::muon) return false;
    double muon_angle;
    double dz = p.start_pos.z - p.end_pos.z;
    double dx = p.start_pos.x - p.end_pos.x;
    double dy = p.start_pos.y - p.end_pos.y;
    double track_length = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (track_length > 0) muon_angle = std::abs(dz / track_length);
    else muon_angle = 1.0;
    return muon_angle > 0.9;
}

bool EnergyCrossMx2(const caf::StandardRecord* sr, int idx) {
    const auto& nu = sr->mc.nu[idx];
    const auto p = nu.prim[0];
    if (abs(p.pdg) != ParticleCode::muon) return false;  
   // muon_energy = p.p.E;
    return p.p.E > 1.0 ;
}


bool IsBasicCCNuAr(const caf::StandardRecord* sr, int idx, const char* mode) {
    const auto& nu = sr->mc.nu[idx];
    if (strcmp(mode, "RHC") == 0){
        return ( kIsCC(sr,idx) && kIsAntiNumu(sr,idx)  && kIsVtxFV(nu.vtx) && kIsArgon(sr,idx) && EnergyCrossMx2(sr,idx)&& AngleCrossMx2(sr,idx) );
    }
    if(strcmp(mode, "FHC") == 0){
        return ( kIsCC(sr,idx) && kIsNumu(sr,idx)  &kIsVtxFV(nu.vtx) && kIsArgon(sr,idx) && EnergyCrossMx2(sr,idx)&& AngleCrossMx2(sr,idx) );
    }
    return false;
}

bool ContLargeProtons(const caf::StandardRecord* sr, int idx) {
    const auto& nu = sr->mc.nu[idx];
    for (const auto& p : nu.prim) {
        if (p.pdg == ParticleCode::proton) {
            double dx = p.end_pos.x - p.start_pos.x;
            double dy = p.end_pos.y - p.start_pos.y;
            double dz = p.end_pos.z - p.start_pos.z;
            double track_length = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (track_length > 2.0) return true;
        }
    }
    return false;
}

// bool QELikeSignal(const caf::StandardRecord* sr, int idx, const char* mode) {
//     const auto& nu = sr->mc.nu[idx];
//     return (
//         IsBasicCCNuAr(sr, idx, mode) &&
//         !ContLargeProtons(sr, idx) &&
//         nu.nneutron > 0 &&
//         nu.npip == 0 &&
//         nu.npim == 0 &&
//         nu.npi0 == 0 &&
//         (nu.mode == 1 || nu.mode == 4 || nu.mode == 10)
//     );
// }


bool QELikeSignal(const caf::StandardRecord* sr, int idx, const char* /*mode*/) {
    const auto& nu = sr->mc.nu[idx];

    if (abs(nu.pdg) != ParticleCode::numu || 
        !nu.iscc || 
        nu.targetPDG != ParticleCode::argon ||
        (nu.mode != 1 && nu.mode != 4 && nu.mode != 10) ||
        nu.nneutron < 1 ||
        nu.npip > 0 || nu.npim > 0 || nu.npi0 > 0 ||
        !kIsVtxFV(nu.vtx)) {
        return false;
    }

    bool has_muon = false;
    bool short_proton = true;
    float muon_energy = 0;
    float muon_angle = 1;

    for (const auto& p : nu.prim) {
        if (abs(p.pdg) == ParticleCode::muon) {
            has_muon = true;
            muon_energy = p.p.E;
            double dx = p.start_pos.x - p.end_pos.x;
            double dy = p.start_pos.y - p.end_pos.y;
            double dz = p.start_pos.z - p.end_pos.z;
            double track_length = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (track_length > 0)
                muon_angle = std::abs(dz / track_length);
        }

        if (p.pdg == ParticleCode::proton) {
            double dx = p.end_pos.x - p.start_pos.x;
            double dy = p.end_pos.y - p.start_pos.y;
            double dz = p.end_pos.z - p.start_pos.z;
            double track_length = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (track_length > 2.0) short_proton = false;
        }
    }

    return has_muon && short_proton && muon_energy > 1.0 && muon_angle > 0.9;
}




bool track_minerva_match(int nixn, int npart, const caf::StandardRecord* sr){
    bool mcOnly = true;
    double mnvOffsetX=-10; double mnvOffsetY=5;
    if (mcOnly){ mnvOffsetX=0; mnvOffsetY=0;}
    int maxPartMinerva=-999; int maxTypeMinerva=-999;  
    int maxIxnMinerva=-999;  
    double dirZExiting=-999;


    // Calculate 2x2 reco track quantities
    auto start_pos=sr->common.ixn.dlp[nixn].part.dlp[npart].start;
    auto end_pos=sr->common.ixn.dlp[nixn].part.dlp[npart].end;

    double diffVertexdZ=abs(start_pos.z-sr->common.ixn.dlp[nixn].vtx.z); 
    double diffVertexdX=abs(start_pos.x-sr->common.ixn.dlp[nixn].vtx.x);
    double diffVertexdY=abs(start_pos.y-sr->common.ixn.dlp[nixn].vtx.y);
    double diffVertex=TMath::Sqrt(diffVertexdZ*diffVertexdZ+diffVertexdY*diffVertexdY+diffVertexdX*diffVertexdX);
    
    // If track away from vertex return false
    if(diffVertex>5) return false;
    
    double dX=(end_pos.x-start_pos.x);
    double dY=(end_pos.y-start_pos.y);
    double dZ=(end_pos.z-start_pos.z);
    double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
    double dirX=dX/length; double dirY=dY/length; double dirZ=dZ/length;
    // Revert direction if not properly reconstructed
    if (dirZ<0){ dirZ=-dirZ; dirX=-dirX; dirY=-dirY; auto temp=start_pos; end_pos=start_pos; end_pos=temp;}

    int minervaPass=0; double dotProductDS=-999;
    double deltaExtrapY=-999; double deltaExtrapX=-999; 

    // Loop over MINERvA ixn
    for(int i=0; i<sr->nd.minerva.ixn.size(); i++){

        // Loop over MINERvA tracks
        for (int j=0; j<sr->nd.minerva.ixn[i].ntracks; j++){
            double dir_z=sr->nd.minerva.ixn[i].tracks[j].dir.z;
            double end_z=sr->nd.minerva.ixn[i].tracks[j].end.z;
            double end_x=sr->nd.minerva.ixn[i].tracks[j].end.x;
            double start_x=sr->nd.minerva.ixn[i].tracks[j].start.x;

            double end_y=sr->nd.minerva.ixn[i].tracks[j].end.y;
            double start_y=sr->nd.minerva.ixn[i].tracks[j].start.y;
            double start_z=sr->nd.minerva.ixn[i].tracks[j].start.z;

            int truthPart=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
            double dXMnv=(sr->nd.minerva.ixn[i].tracks[j].end.x-sr->nd.minerva.ixn[i].tracks[j].start.x);
            double dYMnv=(sr->nd.minerva.ixn[i].tracks[j].end.y-sr->nd.minerva.ixn[i].tracks[j].start.y);
            double dZMnv=(sr->nd.minerva.ixn[i].tracks[j].end.z-sr->nd.minerva.ixn[i].tracks[j].start.z);
            double lengthMinerva=TMath::Sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv);
                
            if (lengthMinerva<10) return false;

            double dirXMinerva=dXMnv/lengthMinerva;
            double dirYMinerva=dYMnv/lengthMinerva;
            double dirZMinerva=dZMnv/lengthMinerva;
            double dotProduct=dirXMinerva*dirX+dirYMinerva*dirY+dirZ*dirZMinerva;
            double extrapdZ=start_z-end_pos.z;
            double extrapY=dirY/dirZ*(extrapdZ)+end_pos.y-start_y;
            double extrapX=dirX/dirZ*(extrapdZ)+end_pos.x-start_x;
            double diffExtrap=TMath::Sqrt(TMath::Power(extrapY-start_y,2));


            if (dotProductDS<dotProduct && abs(extrapY-mnvOffsetY)<15  && 
                abs(TMath::ATan(dirXMinerva/dirZMinerva)-TMath::ATan(dirX/dirZ))<0.06 && 
                abs(TMath::ATan(dirYMinerva/dirZMinerva)-TMath::ATan(dirY/dirZ))<0.06 && 
                abs(extrapX-mnvOffsetX)<15){ dotProductDS=dotProduct;
                                            deltaExtrapY=extrapY;
                                            deltaExtrapX=extrapX;
                                            dirZExiting=dirZ;

                if (mcOnly){
                maxPartMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                maxTypeMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
                maxIxnMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;}	


                if (end_z>300){ minervaPass=1;} if(dirZExiting<dirZ){ dirZExiting=dirZ;}

                return true;
            } 
        
            }   
        }
    if(dotProductDS==-999){
        return false;
    }

}



// Count number of primary and secondary reconstructed particles
std::pair<int, int> getNumRecoSecPrim(const caf::StandardRecord* sr, int idx) {
    int n_prim = 0;
    int n_sec = 0;
    const auto& reco = sr->common.ixn.dlp[idx];
    for (size_t j = 0; j < reco.part.dlp.size(); ++j) {
        const auto& p = reco.part.dlp[j];
        (p.primary ? n_prim++ : n_sec++);
    }
    return {n_prim, n_sec}; 
}

std::pair<bool, bool> has_valid_muon(const caf::StandardRecord* sr, int idx) {
    const auto& reco = sr->common.ixn.dlp[idx];
    bool has_exiting_muon = false;
    bool cross_mx2 = false;

    for (size_t j = 0; j < reco.part.dlp.size(); ++j) {
        const auto& p = reco.part.dlp[j];
        if (p.pdg == ParticleCode::muon && p.primary && p.end.z > constants::Z_BOUND) {
            has_exiting_muon = true;
            cross_mx2 = track_minerva_match(idx, j, sr); 
        }
    }

    return {has_exiting_muon, cross_mx2}; 
}

// Basic CC ν-Ar selection based on vertex cont
bool RecoIsBasicCCNuAr(const caf::StandardRecord* sr, int idx) {
    auto [has_exiting_muon, cross_mx2] = has_valid_muon(sr, idx);
    return(has_exiting_muon && cross_mx2 && kIsVtxFV(sr->common.ixn.dlp[idx].vtx));
}

// Final reco QE-like selector with topology limits
bool RecoIsQELikeSignal(const caf::StandardRecord* sr, int idx, int prim_limit, int sec_limit) {
    auto [n_prim, n_sec] = getNumRecoSecPrim(sr, idx);  // 
    //return (RecoIsBasicCCNuAr(sr, idx) && (n_prim < prim_limit) && (n_sec < sec_limit) && (n_sec > 0));
    return (
        RecoIsBasicCCNuAr(sr, idx) &&
        n_prim == 1 &&
        n_sec < 3
    );
}



