# neutrino_signal_analysis.py

import numpy as np
import pandas as pd
import h5py
import re
import uproot
import glob
import copy
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
from utils_analysis import ParticleCode, load_dataset
sys.path.insert(0, '../tools/')
from caf_readers import CafReader
from mx2_matching import Mx2_DS_Match


# Custom matplotlib settings
plt.style.use('./cfg/my_custom_plot.mplstyle')
import matplotlib.ticker

class MyLocator(matplotlib.ticker.AutoMinorLocator):
    def __init__(self, n=4):
        super().__init__(n=n)
        
matplotlib.ticker.AutoMinorLocator = MyLocator        
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True


#------------------------------------------------------------------------------------------------------
# neutrino_signal_analysis.py
# Analysis class to encapsulate methods
class NeutrinoAnalysis:
    def __init__(self, n_files=60, location="nersc", type="mr6", tpc_dist=5, xbound=63.931, ybound=62.076, zbound=64.3163):
        # Initialization and dataset loading
        self.df = load_dataset(n_files, location, type)
        self.caf_reader = CafReader(self.df)
        self.pdg_tab = ParticleCode()
        
        self.tpc_dist = tpc_dist
        self.xbound = xbound
        self.ybound = ybound
        self.zbound = zbound

        # Constants for interaction counting
        self.ScatteringMode = {'CC kQE': 1, 'CC kDIS': 3, 'CC kRes': 4, 'CC kCoh': 5, 'CC kMEC': 10}
        self.CurrentMode = {'NC': 0, 'RockMuons': 2}
        self.scattering_counters = {mode: 0 for mode in self.ScatteringMode.values()}
        self.current_counters = {mode: 0 for mode in self.CurrentMode.values()}
        self.insignal_counters = {mode: 0 for mode in self.ScatteringMode.values()}

    # Utility function to calculate cosL
    def calculate_cosL(self, ev_data, ip):
        dz = ev_data['rec.mc.nu.prim.start_pos.z'][ip] - ev_data['rec.mc.nu.prim.end_pos.z'][ip]
        start_pos = np.array([
            ev_data['rec.mc.nu.prim.start_pos.x'][ip], 
            ev_data['rec.mc.nu.prim.start_pos.y'][ip],
            ev_data['rec.mc.nu.prim.start_pos.z'][ip]
        ])
        end_pos = np.array([
            ev_data['rec.mc.nu.prim.end_pos.x'][ip], 
            ev_data['rec.mc.nu.prim.end_pos.y'][ip],
            ev_data['rec.mc.nu.prim.end_pos.z'][ip]
        ])
        track_length = np.linalg.norm(start_pos - end_pos)
        return abs(dz / track_length)

    # Step 1: Define the true neutrino signal
    def true_neutrino_signal(self, ev_data): 
        signal = np.zeros_like(ev_data['rec.mc.nu.vtx.x'], dtype=bool)
        mask_t = (
            (abs(ev_data['rec.mc.nu.vtx.x']) < self.xbound - self.tpc_dist) &
            (abs(ev_data['rec.mc.nu.vtx.x']) > self.tpc_dist) &
            (abs(ev_data['rec.mc.nu.vtx.y']) < self.ybound - self.tpc_dist) &
            (abs(ev_data['rec.mc.nu.vtx.z']) > self.tpc_dist) &
            (abs(ev_data['rec.mc.nu.vtx.z']) < self.zbound - self.tpc_dist) &
            (ev_data['rec.mc.nu.targetPDG'] == self.pdg_tab.argon) &
            (ev_data['rec.mc.nu.iscc'] == 1) &
            (abs(ev_data['rec.mc.nu.pdg']) == self.pdg_tab.numu) 
        #    (ev_data['rec.mc.nu.mode']==1)  
        )
        signal[mask_t] = True
        
        for ix in range(ev_data['rec.mc.nu..length']):
            if mask_t[ix]: 
                n_particles = ev_data['rec.mc.nu.prim..length'][ix]
                n_pre = np.sum(ev_data['rec.mc.nu.prim..length'][:ix]) if ix > 0 else 0
                for ip in range(n_pre, n_pre + n_particles):
                    pdg = ev_data['rec.mc.nu.prim.pdg'][ip]    
                    if pdg == self.pdg_tab.muon: 
                        Elep = ev_data['rec.mc.nu.prim.p.E'][ip]
                        cosL = self.calculate_cosL(ev_data, ip)
                        if cosL < 0.9 and Elep < 1:
                            signal[ix] = False
        return signal

    # Step 2: Define the reconstructed neutrino signal in FV
    def reco_neutrino_signal(self,ev_data, ev, mode, QE=False):
        mask_fv = (
            (abs(ev_data['rec.common.ixn.dlp.vtx.x']) < self.xbound - self.tpc_dist) &
            (abs(ev_data['rec.common.ixn.dlp.vtx.x']) > self.tpc_dist) &
            (abs(ev_data['rec.common.ixn.dlp.vtx.y']) < self.ybound - self.tpc_dist) &
            (abs(ev_data['rec.common.ixn.dlp.vtx.z']) > self.tpc_dist) &
            (abs(ev_data['rec.common.ixn.dlp.vtx.z']) < self.zbound - self.tpc_dist)
        )
        # If QE is True, apply QE-specific signal selection; otherwise, use fiducial volume mask
        vertices = mask_fv if not QE else self.QE_reco_signal(ev_data, mask_fv)
        if mode == 'FV': reco_vertices = vertices
        elif mode == '2x2':
            track_out_2x2 = self.track_selection_signal(ev_data)
            reco_vertices = vertices & track_out_2x2
        elif mode == 'mx2':
            match_mx2 = self.track_match_mx2(ev_data, ev)
            reco_vertices = vertices & match_mx2
        else: 
            raise ValueError("Invalid mode. Options: FV, 2x2, mx2")
        return reco_vertices

    # Step 2.5 select tracks that exit 2x2
    def track_selection_signal(self,ev_data):
        #select vertices that are primary and exit 2x2
        z_end_track = ev_data['rec.common.ixn.dlp.part.dlp.end.z']
        is_reco_primary = ev_data['rec.common.ixn.dlp.part.dlp.primary']
        #array to store signal values
        track_out = np.zeros_like(ev_data['rec.common.ixn.dlp.vtx.x'], dtype=bool)

        for vtx in range(len(ev_data['rec.common.ixn.dlp.vtx.x'])):
            n_reco_tracks = ev_data['rec.common.ixn.dlp.part.dlp..length'][vtx]
            n_pre = np.sum(ev_data['rec.common.ixn.dlp.part.dlp..length'][:vtx]) if vtx > 0 else 0
            for ip in range(n_pre, n_pre + n_reco_tracks):
                if (is_reco_primary[ip] and z_end_track[ip]>=self.zbound):
                        track_out[vtx] = True
                        break 
        return track_out

    def track_match_mx2(self,ev_data,ev):

        minerva_data = self.caf_reader.get_minerva_data(self.df.iloc[ev])
        is_reco_primary = ev_data['rec.common.ixn.dlp.part.dlp.primary']

        x_start_track  = ev_data['rec.common.ixn.dlp.part.dlp.start.x']
        y_start_track = ev_data['rec.common.ixn.dlp.part.dlp.start.y']
        z_start_track = ev_data['rec.common.ixn.dlp.part.dlp.start.z']

        x_end_track = ev_data['rec.common.ixn.dlp.part.dlp.end.x']
        y_end_track = ev_data['rec.common.ixn.dlp.part.dlp.end.y']
        z_end_track = ev_data['rec.common.ixn.dlp.part.dlp.end.z']

        track_Mx2 = np.zeros_like(ev_data['rec.common.ixn.dlp.vtx.x'], dtype=bool)

        for vtx in range(len(ev_data['rec.common.ixn.dlp.vtx.x'])):
            n_reco_tracks = ev_data['rec.common.ixn.dlp.part.dlp..length'][vtx]
            n_pre = np.sum(ev_data['rec.common.ixn.dlp.part.dlp..length'][:vtx]) if vtx > 0 else 0
            for ip in range(n_pre, n_pre + n_reco_tracks):
                    if (is_reco_primary[ip] and z_end_track[ip] >= self.zbound):
                        muon_track = [[x_start_track[ip], y_start_track[ip],z_start_track[ip]],[x_end_track[ip],y_end_track[ip],z_end_track[ip]]]
                        reco_vtx = [ev_data['rec.common.ixn.dlp.vtx.x'][vtx],ev_data['rec.common.ixn.dlp.vtx.y'][vtx], ev_data['rec.common.ixn.dlp.vtx.z'][vtx]]
                        m_index,ds_dp,deltaX,deltaY,exit_minerva = Mx2_DS_Match(muon_track, minerva_data, reco_vtx, mc=True,verbose=False)
                        if (exit_minerva == True and m_index!=None):
                            track_Mx2[vtx] = True
        return track_Mx2

    # This function selects only reco tracks that produced one primary particle
    # These are vertices that are potential QE interactions. This is simplest mode
    def QE_reco_signal(self,ev_data, mask_fv):
        QE_signal = mask_fv.copy()
        for vtx in range(ev_data['rec.common.ixn.dlp..length']):     
            n_vtx = ev_data['rec.common.ixn.dlp.part.dlp..length'][vtx]
            n_pre = np.sum(ev_data['rec.common.ixn.dlp.part.dlp..length'][:vtx]) if vtx > 0 else 0
            #iterate over particles in vtx
            count = 0
            for ip in range(n_pre, n_pre+n_vtx):
                count +=  ev_data['rec.common.ixn.dlp.part.dlp.primary'][ip]
                if count > 1: break 
            if count != 1: QE_signal[vtx] = 0  #discard anything else 
        return QE_signal

    # Step 3: Match reco vertices and true signal interactions
    def matching(self,ev_data, vtx_index, signal, reco_vertices_signal):

        n_vtx = ev_data['rec.common.ixn.dlp.truth..length'][vtx_index]
        n_pre = np.sum(ev_data['rec.common.ixn.dlp.truth..length'][:vtx_index]) if vtx_index > 0 else 0 
        max_overlap = 0
        best_match_idx = None

        for ip in range(n_pre, n_pre + n_vtx):
            temp_overlap = ev_data['rec.common.ixn.dlp.truthOverlap'][ip]
            temp_tp_ixn = ev_data['rec.common.ixn.dlp.truth'][ip]
            if temp_overlap > max_overlap:
                max_overlap = temp_overlap
                best_match_idx = temp_tp_ixn

        if best_match_idx is not None:
            delta_x = np.abs(ev_data['rec.mc.nu.vtx.x'][best_match_idx] - ev_data['rec.common.ixn.dlp.vtx.x'][vtx_index])
            delta_y = np.abs(ev_data['rec.mc.nu.vtx.y'][best_match_idx] - ev_data['rec.common.ixn.dlp.vtx.y'][vtx_index])
            delta_z = np.abs(ev_data['rec.mc.nu.vtx.z'][best_match_idx] - ev_data['rec.common.ixn.dlp.vtx.z'][vtx_index])
            if signal[best_match_idx] and delta_x < 5 and delta_y < 5 and delta_z < 5:
                reco_vertices_signal[vtx_index] = True
                return True, best_match_idx, reco_vertices_signal
        return False, best_match_idx, reco_vertices_signal

    #Step 4: get the type of interaction being reconstructed
    def type_interaction(self,ev_data, best_match_idx): 
        scattering_mode = ev_data['rec.mc.nu.mode'][best_match_idx]
        current_mode = ev_data['rec.mc.nu.iscc'][best_match_idx]
        if ev_data['rec.mc.nu.id'][best_match_idx] > 1E9:
            self.current_counters[self.CurrentMode['RockMuons']] += 1  # Increment rock_muons counter
        if scattering_mode in self.ScatteringMode.values() and current_mode ==1:
            self.scattering_counters[scattering_mode] += 1  # Increment the counter for the specific mode
        elif scattering_mode in self.ScatteringMode.values() and current_mode ==0:
            self.current_counters[self.CurrentMode['NC']] += 1  # Increment the counter for NC
        return self.scattering_counters, self.current_counters

    # Helper function to calculate efficiency and purity
    def efficiency_purity(self,count_match_FV, count_reco, num_neutrinos_signal_t):
        purity = count_match_FV / count_reco if count_reco > 0 else 0
        efficiency = count_match_FV / num_neutrinos_signal_t if num_neutrinos_signal_t > 0 else 0
        return efficiency, purity

        # Plotting functions for interaction distributions
    def plot_interactions(self, count_match_FV, pdf_pages=None):
        mode_labels = (
            [key for key, value in self.ScatteringMode.items() if value in self.scattering_counters] +
            [key for key, value in self.CurrentMode.items() if value in self.current_counters] +
            ['True Signal']
        )
        mode_counts = (
            [self.scattering_counters[self.ScatteringMode[label]] for label in self.ScatteringMode.keys()] +
            [self.current_counters[self.CurrentMode[label]] for label in self.CurrentMode.keys()] +
            [count_match_FV]
        )
        
        if all(value == 0 for value in mode_counts):  # Check if all values are zero
            print("No data to plot for interaction types.")
            return

        plt.figure(figsize=(7, 7))
        wedges, texts, autotexts = plt.pie(mode_counts, autopct='%1.1f%%', startangle=90, pctdistance=0.85)
        plt.legend(wedges, mode_labels, title="Interaction Modes", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))
        plt.title('Interaction Type of Reconstructed Vertices')
        plt.axis('equal')
        if pdf_pages:
            pdf_pages.savefig()
        plt.close()

    def plot_signal_interactions(self, true_counters, pdf_pages=None):
        labels = [label for label, code in self.ScatteringMode.items() if code in true_counters]
        sizes = [true_counters[self.ScatteringMode[label]] for label in labels]
        
        if all(value == 0 for value in sizes):  # Check if all values are zero
            print("No data to plot for scattering mode of vertices matching signal.")
            return

        plt.figure(figsize=(7, 7))
        plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        plt.title('Scattering Mode of Vertices Matching Signal')
        plt.axis('equal')
        if pdf_pages:
            pdf_pages.savefig()
        plt.close()

    

    def display_summary(self,efficiency, purity,pdf_pages=None):

        # Print and plot results    
        summary_text = f"Interaction Summary:\n"
        total_interactions = sum(self.scattering_counters.values()) + sum(self.current_counters.values())
        for mode, count in {**self.scattering_counters, **self.current_counters}.items():
            percent = (count / total_interactions) * 100 if total_interactions > 0 else 0
            summary_text += f"{mode}: {percent:.2f}%\n"
        summary_text += f"\nEfficiency: {efficiency:.2%}\nPurity: {purity:.2%}\n"  
        # Display in terminal
        print(summary_text)    
        # Save summary text as a page in the PDF
        if pdf_pages:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.text(0.1, 0.5, summary_text, fontsize=12, ha='left', wrap=True)
            ax.axis('off')
            pdf_pages.savefig(fig)
            plt.close(fig)

 
    def run_analysis(self, output_path, pdf_path):
        num_vertices, count_match_FV, num_match_signal = 0, 0, 0
        filtered_indices = []

        for ev in range(len(self.df)):
            ev_data = self.df.iloc[ev]
            true_interactions = self.true_neutrino_signal(ev_data)
            reco_vertices = self.reco_neutrino_signal(ev_data, ev, mode='mx2', QE=False)
            match_vertices_signal = np.zeros_like(reco_vertices, dtype=bool)
            num_vertices += np.sum(reco_vertices)
            for vtx_index, value in enumerate(reco_vertices):
                if value:
                    match, best_match_idx, match_vertices_signal = self.matching(ev_data, vtx_index, true_interactions, match_vertices_signal)
                    if match:
                        count_match_FV += 1
                        filtered_indices.append(ev)

        # Calculate efficiency and purity
        efficiency, purity = self.efficiency_purity(count_match_FV, num_vertices, num_match_signal)

        # Save filtered data
        filtered_data_df = self.df.iloc[filtered_indices]
        filtered_data_df.to_csv(output_path, index=False)

        # Save plots and summary to PDF
        with PdfPages(pdf_path) as pdf_pages:
            self.plot_interactions(count_match_FV, pdf_pages)
            self.plot_signal_interactions(self.insignal_counters, pdf_pages)
            self.display_summary(efficiency, purity, pdf_pages)

# Main function
if __name__ == "__main__":
    analysis = NeutrinoAnalysis()
    output_csv_path = '/global/homes/g/gkufatty/projects/neutron-multiplicity/analysis/neutrino_signal.csv'
    output_pdf_path = '/global/homes/g/gkufatty/projects/neutron-multiplicity/analysis/neutrino_signal_plots.pdf'
    analysis.run_analysis(output_csv_path, output_pdf_path)
