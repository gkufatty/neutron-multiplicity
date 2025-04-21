import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import uproot
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys



'''
Utility functions to run
the example notebooks
'''

def load_keys(keys_set):
    if(keys_set=="mr6"):
        input_keys = "./cfg/mr6_keys.txt"

    else:
        input_keys = "./cfg/caf_keys.txt"

    keys_list = open(input_keys,'r')
    data = []
    for line in keys_list:
        line = line.strip()
        data.append(line)
    return data


def load_dataset(n_files,location,mr='mr6'):
    print("Reading ", n_files, " files")
    print(f"Location selected {location}")
    if(location=="nersc" and mr=='mr5_fix'):
        print("Openning MiniRun 5 beta 2.a CAFs")
        input_list = "./cfg/minirun5_noe_fix.txt"

    elif(location=="nersc" and mr=='mr6'):
        print("Openning MiniRun 6 CAFs")
        input_list = "./cfg/minirun6.txt"

    elif(location=="nersc" and mr=='mr5_beta1'):
        print("Openning MiniRun 5 beta 1 CAFs")
        input_list = "./cfg/minirun5_beta1.txt"
    else:
        print("Openning MiniRun 5 beta 2.a CAFs")
        input_list = "./cfg/minirun_5_beta2a_fnal.txt"

    with open(input_list, 'r') as file:
        # Convert the file object to a list
        lines = list(file)
    data = []
    df = pd.DataFrame(data)
    counter = 0
    counter_max = n_files 
    for ifile in tqdm(range(counter_max)):
        line =lines[ifile]
        line = line.strip()
        #print("Reading", line)
        caf_file = uproot.open(line)
        caf_tree = caf_file['cafTree']
        # Create an empty dictionary to store branch data
        caf_data_dict = {}
        # Import list of keys for n-Ar
        caf_keys = load_keys(mr)
        #Iterate over the branch names in the TTree
        for branch_name in caf_keys:
            # Use branch_name as the key and fetch the data using .array()
            caf_data_dict[branch_name] = caf_tree[branch_name].array(library="np")
        # Create a Pandas DataFrame from the dictionary
        df_temp = pd.DataFrame(caf_data_dict)
        df_temp['file_number'] = np.nan
        df_temp['file_entry'] = np.nan
        for i in range(len(df_temp)):
            df_temp.at[i,'file_number'] = ifile
            df_temp.at[i,'file_entry'] = i
        df = pd.concat([df, df_temp])
        counter+=1
    return df 

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
        self.electron = 11 
        self.pi0 = 111
        self.pip = 211
        self.eta = 221
        self.proton = 2212
        self.numu = 14
        self.nue = 12
        self.photon = 22