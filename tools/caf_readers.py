# Tools to extract data from flat CAFs

import numpy as np
import pandas as pd
import uproot
from utils import ParticleCode

particle_data = ParticleCode()

class CafReader():
    '''
    This class will contain
    helper functions to extract reco
    and truth data from flat CAFs

    Class arguments:
        df: pandas data frame
            data frame containing 
            flat caf branches 


    Class attributes:
        reco_branches: list of ml-reco reco branches
        primary_branches: list of true primary branches
        secondary_branches: list of true secondary branches


    Class methods:
        reco_backtrack: Backtracking from reco to true particle
        get_true_ixn_data: Returns a dictionary with true particle data 

    '''


    def __init__(self,df):
        self.nu_branches=[
                "rec.mc.nu.E",
                "rec.mc.nu.Q2",
                "rec.mc.nu.W",
                "rec.mc.nu.baseline",
                "rec.mc.nu.bjorkenX",
                "rec.mc.nu.generator",
                "rec.mc.nu.genieIdx",
                "rec.mc.nu.genweight",
                "rec.mc.nu.hitnuc",
                "rec.mc.nu.id",
                "rec.mc.nu.imp_weight",
                "rec.mc.nu.inelasticity",
                "rec.mc.nu.iscc",
                "rec.mc.nu.ischarm",
                "rec.mc.nu.isseaquark",
                "rec.mc.nu.isvtxcont",
                "rec.mc.nu.mode",
                "rec.mc.nu.modq",
                "rec.mc.nu.momentum.x",
                "rec.mc.nu.momentum.y",
                "rec.mc.nu.momentum.z",
                "rec.mc.nu.nneutron",
                "rec.mc.nu.npi0",
                "rec.mc.nu.npim",
                "rec.mc.nu.npip",
                "rec.mc.nu.nprefsi",
                "rec.mc.nu.nprim",
                "rec.mc.nu.nproton",
                "rec.mc.nu.nsec",
                "rec.mc.nu.parent_dcy_E",
                "rec.mc.nu.parent_dcy_mode",
                "rec.mc.nu.parent_dcy_mom.x",
                "rec.mc.nu.parent_dcy_mom.y",
                "rec.mc.nu.parent_dcy_mom.z",
                "rec.mc.nu.parent_pdg",
                "rec.mc.nu.pdg",
                "rec.mc.nu.pdgorig",
                "rec.mc.nu.prod_vtx.x",
                "rec.mc.nu.prod_vtx.y",
                "rec.mc.nu.prod_vtx.z",
                "rec.mc.nu.q0",
                "rec.mc.nu.removalE",
                "rec.mc.nu.resnum",
                "rec.mc.nu.t", 
                "rec.mc.nu.targetPDG", 
                "rec.mc.nu.time", 
                "rec.mc.nu.vtx.x", 
                "rec.mc.nu.vtx.y", 
                "rec.mc.nu.vtx.z", 
                "rec.mc.nu.xsec", 
                "rec.mc.nu.xsec_cvwgt" 
            ]




        self.reco_branches= [
                "rec.common.ixn.dlp.id",
                "rec.common.ixn.dlp.part.dlp..length",
                "rec.common.ixn.dlp.part.dlp..totarraysize",
                "rec.common.ixn.dlp.part.dlp.E",
                "rec.common.ixn.dlp.part.dlp.E_method",
                "rec.common.ixn.dlp.part.dlp.contained",
                "rec.common.ixn.dlp.part.dlp.end.x",
                "rec.common.ixn.dlp.part.dlp.end.y",
                "rec.common.ixn.dlp.part.dlp.end.z",
                "rec.common.ixn.dlp.part.dlp.p.x",
                "rec.common.ixn.dlp.part.dlp.p.y",
                "rec.common.ixn.dlp.part.dlp.p.z",
                "rec.common.ixn.dlp.part.dlp.pdg",
                "rec.common.ixn.dlp.part.dlp.primary",
                "rec.common.ixn.dlp.part.dlp.score",
                "rec.common.ixn.dlp.part.dlp.start.x",
                "rec.common.ixn.dlp.part.dlp.start.y",
                "rec.common.ixn.dlp.part.dlp.start.z",
                "rec.common.ixn.dlp.part.dlp.tgtA",
                "rec.common.ixn.dlp.part.dlp.truth..length",
                "rec.common.ixn.dlp.part.dlp.truth..totarraysize",
                "rec.common.ixn.dlp.part.dlp.truth.ixn",
                "rec.common.ixn.dlp.part.dlp.truth.part",
                "rec.common.ixn.dlp.part.dlp.truth.type",
                "rec.common.ixn.dlp.part.dlp.truth..idx",
                "rec.common.ixn.dlp.part.dlp.truthOverlap..length",
                "rec.common.ixn.dlp.part.dlp.truthOverlap..totarraysize",
                "rec.common.ixn.dlp.part.dlp.truthOverlap",
                "rec.common.ixn.dlp.part.dlp.truthOverlap..idx",
                "rec.common.ixn.dlp.part.dlp..idx",
            ]
        
        self.reco_nu_branches=[
                "rec.common.ixn.dlp.part.ndlp",
                "rec.common.ixn.dlp.vtx.x",
                "rec.common.ixn.dlp.vtx.y",
                "rec.common.ixn.dlp.vtx.z"
        ]
        
        self.primary_branches = [
                "rec.mc.nu.prim.G4ID",
                "rec.mc.nu.prim.end_pos.x",
                "rec.mc.nu.prim.end_pos.y",
                "rec.mc.nu.prim.end_pos.z",
                "rec.mc.nu.prim.end_process",
                "rec.mc.nu.prim.interaction_id",
                "rec.mc.nu.prim.p.E",
                "rec.mc.nu.prim.p.px",
                "rec.mc.nu.prim.p.py",
                "rec.mc.nu.prim.p.pz",
                "rec.mc.nu.prim.pdg",
                "rec.mc.nu.prim.start_pos.x",
                "rec.mc.nu.prim.start_pos.y",
                "rec.mc.nu.prim.start_pos.z",
                "rec.mc.nu.prim.start_process",
                "rec.mc.nu.prim.time",
                "rec.mc.nu.prim..idx"
            ]
        
        self.secondary_branches = [
                "rec.mc.nu.sec.G4ID",
                "rec.mc.nu.sec.ancestor_id.ixn",
                "rec.mc.nu.sec.ancestor_id.part",
                "rec.mc.nu.sec.ancestor_id.type",
                "rec.mc.nu.sec.end_pos.x",
                "rec.mc.nu.sec.end_pos.y",
                "rec.mc.nu.sec.end_pos.z",
                "rec.mc.nu.sec.end_process",
                "rec.mc.nu.sec.interaction_id",
                "rec.mc.nu.sec.p.E",
                "rec.mc.nu.sec.p.px",
                "rec.mc.nu.sec.p.py",
                "rec.mc.nu.sec.p.pz",
                "rec.mc.nu.sec.parent",
                "rec.mc.nu.sec.pdg",
                "rec.mc.nu.sec.start_pos.x",
                "rec.mc.nu.sec.start_pos.y",
                "rec.mc.nu.sec.start_pos.z",
                "rec.mc.nu.sec.start_process",
                "rec.mc.nu.sec.time",
                "rec.mc.nu.sec..idx"
            ]

        self.minerva_branches = [
                "rec.nd.minerva.ixn.tracks..length",
                "rec.nd.minerva.ixn.tracks..totarraysize",
                "rec.nd.minerva.ixn.tracks.E",
                "rec.nd.minerva.ixn.tracks.Evis",
                "rec.nd.minerva.ixn.tracks.dir.x",
                "rec.nd.minerva.ixn.tracks.dir.y",
                "rec.nd.minerva.ixn.tracks.dir.z",
                "rec.nd.minerva.ixn.tracks.end.x",
                "rec.nd.minerva.ixn.tracks.end.y",
                "rec.nd.minerva.ixn.tracks.end.z",
                "rec.nd.minerva.ixn.tracks.enddir.x",
                "rec.nd.minerva.ixn.tracks.enddir.y",
                "rec.nd.minerva.ixn.tracks.enddir.z",
                "rec.nd.minerva.ixn.tracks.len_cm",
                "rec.nd.minerva.ixn.tracks.len_gcm2",
                "rec.nd.minerva.ixn.tracks.qual",
                "rec.nd.minerva.ixn.tracks.start.x",
                "rec.nd.minerva.ixn.tracks.start.y",
                "rec.nd.minerva.ixn.tracks.start.z",
                "rec.nd.minerva.ixn.tracks.truth..length",
                "rec.nd.minerva.ixn.tracks.truth..totarraysize",
                "rec.nd.minerva.ixn.tracks.truth.ixn",
                "rec.nd.minerva.ixn.tracks.truth.part",
                "rec.nd.minerva.ixn.tracks.truth.type",
                "rec.nd.minerva.ixn.tracks.truth..idx",
                "rec.nd.minerva.ixn.tracks.truthOverlap..length",
                "rec.nd.minerva.ixn.tracks.truthOverlap..totarraysize",
                #"rec.nd.minerva.ixn.tracks.truthOverlap",
                "rec.nd.minerva.ixn.tracks.truthOverlap..idx",
                "rec.nd.minerva.ixn.tracks..idx",
                "rec.nd.minerva.nixn"
            ]


    # Methods     
    def reco_backtrack(self,ixn_index,verbose=False):
        # As it is it backtracks primary reco protons! 
        n_ixn = len(self.df['rec.common.ixn.dlp.id'])
        n_particles = self.df['rec.common.ixn.dlp.part.dlp..length'][ixn_index]
        if(ixn_index==0):
            n_pre = 0
        else: 
            n_pre = np.sum(self.df['rec.common.ixn.dlp.part.dlp..length'][:ixn_index]) 

        if(verbose):
            print('Number of interactions in this event ',n_ixn)
            print('Printing info for interaction ', ixn_index)
            print(f'This interaction has {n_particles} reco particles')

        # Scan over particles associated to this ixn 
        for ip in range(n_pre,n_pre + n_particles):
            #print(ip)
            # Check if particle is primary 
            is_primary = self.df['rec.common.ixn.dlp.part.dlp.primary'][ip]
            if(is_primary==False):
                continue 
            else:
                # Get reco PDG
                reco_pdg = self.df['rec.common.ixn.dlp.part.dlp.pdg'][ip]
                if(verbose): print(f"Particle {ip} is a {reco_pdg} primary")
                if(reco_pdg==particle_data.Proton()):
                    truth_length = self.df['rec.common.ixn.dlp.part.dlp.truth..length'][ip]
                    truth_pre = np.sum(self.df['rec.common.ixn.dlp.part.dlp.truth..length'][:ip])
                    if(verbose):print(f'With {truth_length} overlapping particles')
                    max_overlap = 0
                    best_match_idx = 0
                    best_match_type = 0
                    # Need to do the same trick as for scan over reco particles
                    for tp in range(truth_pre,truth_pre+truth_length):
                        temp_overlap = self.df['rec.common.ixn.dlp.part.dlp.truthOverlap'][tp]
                        temp_tp_match = self.df['rec.common.ixn.dlp.part.dlp.truth.part'][tp]
                        temp_tp_type = self.df['rec.common.ixn.dlp.part.dlp.truth.type'][tp] # 1 prim and 3 second
                        if(verbose):print(f'True ({temp_tp_type}) particle {temp_tp_match} has this amount of overlap {temp_overlap}')
                        if(temp_overlap>max_overlap):
                            max_overlap=temp_overlap
                            best_match_idx=temp_tp_match
                            best_match_type=temp_tp_type


                    if(best_match_type==3):
                        best_match_pdg = self.df['rec.mc.nu.sec.pdg'][best_match_idx]
                        if(verbose):print('true secondary, PDG of best match: ', best_match_pdg)
                    elif(best_match_type==1):
                        best_match_pdg = self.df['rec.mc.nu.prim.pdg'][best_match_idx]
                        if(verbose):print('true primary, PDG of best match: ', best_match_pdg)


    def get_true_ixn_data(self,my_event,ixn_id,verbose=False):
        nprim = my_event['rec.mc.nu.prim..length'][ixn_id]
        nprim_pre = np.sum(my_event['rec.mc.nu.prim..length'][:ixn_id])
        nsec = my_event['rec.mc.nu.sec..length'][ixn_id]
        nsec_pre = np.sum(my_event['rec.mc.nu.sec..length'][:ixn_id])
        out_dict = {}

        
        if(verbose):
            print(f"Extracting true data for ixn {ixn_id}")
            print(f"Number of primary particles {nprim}")
            print(f"Number of pre-prim {nprim_pre}")
            print(f"Number of secondary particles {nsec}")
            print(f"Number of pre-sec {nsec_pre}")

        # Get nu info 
        for branch in self.nu_branches:
            if(np.isscalar(my_event[branch])):
                out_dict[branch] = my_event[branch]
            else:
                out_dict[branch] = my_event[branch][ixn_id]

        # Get primaries
        for branch in self.primary_branches:
            out_dict[branch] = my_event[branch][nprim_pre:nprim_pre+nprim] 

        # Get secondaries 
        for branch in self.secondary_branches:
            out_dict[branch] = my_event[branch][nsec_pre:nsec_pre+nsec]

        return out_dict 

    def get_reco_ixn_data(self,my_event,ixn_id,verbose=False):
        out_dict={}
        nreco_part = my_event["rec.common.ixn.dlp.part.dlp..length"][ixn_id]
        nreco_part_pre = np.sum(my_event["rec.common.ixn.dlp.part.dlp..length"][:ixn_id])
        part_pdg = my_event["rec.common.ixn.dlp.part.dlp.pdg"][nreco_part_pre:nreco_part_pre+nreco_part]
        part_start_x = my_event["rec.common.ixn.dlp.part.dlp.start.x"][nreco_part_pre:nreco_part_pre+nreco_part]
        part_start_z = my_event["rec.common.ixn.dlp.part.dlp.start.z"][nreco_part_pre:nreco_part_pre+nreco_part]
        part_end_x = my_event["rec.common.ixn.dlp.part.dlp.end.x"][nreco_part_pre:nreco_part_pre+nreco_part]
        part_end_z = my_event["rec.common.ixn.dlp.part.dlp.end.z"][nreco_part_pre:nreco_part_pre+nreco_part]

        if(verbose):
            print(f"Number of reconstructed particles {nreco_part}")
            print(f"PDG codes of reco parts {part_pdg}")
            print(f"start points {part_start_z}, {part_start_x}")
            print(f"end points {part_end_z}, {part_end_x}")


        for branch in self.reco_nu_branches:
            if(np.isscalar(my_event[branch])):
                out_dict[branch] = my_event[branch]
            else:
                out_dict[branch] = my_event[branch][ixn_id]

        # Get reco part info 
        for branch in self.reco_branches:
            if(np.isscalar(my_event[branch])):
                out_dict[branch] = my_event[branch]
            else:
                out_dict[branch] = my_event[branch][nreco_part_pre:nreco_part_pre+nreco_part]
    
    def get_minerva_data(self, my_event):
        out_dict= {}
        # Get MINERvA info
        for branch in self.minerva_branches:
            out_dict[branch] = my_event[branch]
        return out_dict 
    
    
    
    def dump_branches(self,my_event,data_level):

        if(data_level=="truth" or data_level=="all"):
            print("Printing truth info...")
            # Get nu info 
            for branch in self.nu_branches:
                print(branch, my_event[branch])

            # Get primaries
            for branch in self.primary_branches:
                print(branch, my_event[branch])
            # Get secondaries 
            for branch in self.secondary_branches:
                print(branch, my_event[branch])
            print("==============================")

        
        if(data_level=="reco" or data_level=="all"):
            print("Printing reco info...")

            # Get spine/ML-reco 
            for branch in self.reco_branches:
                print(branch, my_event[branch])

            # Get MINERvA
            for branch in self.minerva_branches:
                print(branch, my_event[branch])
            print("==============================")

    
        