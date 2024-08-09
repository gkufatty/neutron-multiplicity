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


        # Get primaries
        for branch in self.primary_branches:
            out_dict[branch] = my_event[branch][nprim_pre:nprim_pre+nprim] 

        # Get secondaries 
        for branch in self.secondary_branches:
            out_dict[branch] = my_event[branch][nsec_pre:nsec_pre+nsec]

        return out_dict 
    
        