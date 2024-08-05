import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import uproot
import pandas as pd
import numpy as np
import sys



'''
Utility functions to run
the example notebooks
'''

def load_keys(keys_set):
    if(keys_set=="n-Ar"):
        input_keys = "./cfg/caf_keys.txt"
    keys_list = open(input_keys,'r')
    data = []
    for line in keys_list:
        line = line.strip()
        data.append(line)
    return data


def load_dataset(n_files,location):
    print("Openning MiniRun 5 beta 2.a CAFs")
    print("Reading ", n_files, " files")

    if(location=="nersc"):
        input_list = "./cfg/minirun_5_beta2a.txt"
    else:
        input_list = "./cfg/minirun_5_beta2a_fnal.txt"
    file_list = open(input_list, 'r')
    data = []
    df = pd.DataFrame(data)
    counter = 0
    counter_max = n_files 
    for line in file_list:
        if(counter > counter_max):
            break
        else:
            line = line.strip()
            print("Reading", line)
            caf_file = uproot.open(line)
            caf_tree = caf_file['cafTree']
            # Create an empty dictionary to store branch data
            caf_data_dict = {}
            # Import list of keys for n-Ar
            caf_keys = load_keys("n-Ar")
            #Iterate over the branch names in the TTree
            for branch_name in caf_keys:
                # Use branch_name as the key and fetch the data using .array()
                caf_data_dict[branch_name] = caf_tree[branch_name].array(library="np")
            # Create a Pandas DataFrame from the dictionary
            df_temp = pd.DataFrame(caf_data_dict)
            df = pd.concat([df, df_temp])
            counter+=1
    return df 



def reco_backtrack(df_test,ixn_index):
    n_ixn = len(df_test['rec.common.ixn.dlp.id'])
    n_particles = df_test['rec.common.ixn.dlp.part.dlp..length'][ixn_index]
    if(ixn_index==0):
        n_pre = 0
    else: 
        n_pre = np.sum(df_test['rec.common.ixn.dlp.part.dlp..length'][:ixn_index]) 
    print('Number of interactions in this event ',n_ixn)
    print('Printing info for interaction ', ixn_index)
    print(f'This interaction has {n_particles} reco particles')

    # Scan over particles associated to this ixn 
    for ip in range(n_pre,n_pre + n_particles):
        #print(ip)
        # Check if particle is primary 
        is_primary = df_test['rec.common.ixn.dlp.part.dlp.primary'][ip]
        if(is_primary==False):
            continue 
        else:
            # Get reco PDG
            reco_pdg = df_test['rec.common.ixn.dlp.part.dlp.pdg'][ip]
            print(f"Particle {ip} is a {reco_pdg} primary")
            if(reco_pdg==2212):
                truth_length = df_test['rec.common.ixn.dlp.part.dlp.truth..length'][ip]
                truth_pre = np.sum(df_test['rec.common.ixn.dlp.part.dlp.truth..length'][:ip])
                print(f'With {truth_length} overlapping particles')
                max_overlap = 0
                best_match_idx = 0
                best_match_type = 0
                # Need to do the same trick as for scan over reco particles
                for tp in range(truth_pre,truth_pre+truth_length):
                    temp_overlap = df_test['rec.common.ixn.dlp.part.dlp.truthOverlap'][tp]
                    temp_tp_match = df_test['rec.common.ixn.dlp.part.dlp.truth.part'][tp]
                    temp_tp_type = df_test['rec.common.ixn.dlp.part.dlp.truth.type'][tp] # 1 prim and 3 second
                    print(f'True ({temp_tp_type}) particle {temp_tp_match} has this amount of overlap {temp_overlap}')
                    if(temp_overlap>max_overlap):
                        max_overlap=temp_overlap
                        best_match_idx=temp_tp_match
                        best_match_type=temp_tp_type


                if(best_match_type==3):
                    best_match_pdg = df_test['rec.mc.nu.sec.pdg'][best_match_idx]
                    print('true secondary, PDG of best match: ', best_match_pdg)
                elif(best_match_type==1):
                    best_match_pdg = df_test['rec.mc.nu.prim.pdg'][best_match_idx]
                    print('true primary, PDG of best match: ', best_match_pdg)


class ParticleCode():
    # Class containing pdg codes
    def __init__(self):

        # Massess
        self.neutron_mass = 939.5654 # MeV/c2 
        self.proton_mass = 938.2702 # MeV/c2


        # PDG Codes 
        self.argon = 1000180400
        self.chlorine = 1000170360
        self.sulfur = 1000160320 
        self.muon = 13
        self.neutron = 2112
        self.pi0 = 111
        self.pip = 211
        self.eta = 221
        self.proton = 2212
        self.numu = 14
        self.nue = 12
        self.photon = 22
