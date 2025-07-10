# How to run the analysis
There are three parts. Neutrino selection, neutron selection, and light study. If you are using a MiniRun version that is not in the folder, then there is a pre-step to generate the list of files you will use

## Analysis Setup

This does DUNE setup, set up dependencies and produces the input txt file with the list of CAFs to be use. 

To setup DUNE
```bash
/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer shell --shell=/bin/bash -B /cvmfs,/exp,/nashome,/pnfs/dune,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.co$
ipc --pid /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest 
```
Setup c++ dependencies
```bash
source neutron-multiplicity/analysis/setup_dependencies.sh
```
If you want to use a new MiniRun version, generate the files providing version and mode (RHC/FHC)
```bash
cd neutron-multiplicity/analysis/nu_selection/input_files/
source macro.sh
./generate_input_list <folder_path> <version> <MODE>
```
for example:
```bash
source macro.sh
./generate_input_list /exp/dune/data/users/noeroy/prod/MiniRun6.2_1E19_RHC/MiniRun6.2_1E19_RHC.caf/CAF/ 6.2 RHC
```

## Running the neutrino selection
Compile and then run the selection
```bash
cd neutron-multiplicity/analysis/nu_selection/
source macro.sh
./select_nus <version> <MODE>
```
for example:
```bash
source macro.sh
./select_nus 6.2 RHC
```
This generates two files in output_files, one with only the list of cafs files that passed the neutrino selection and other one with also the vertex informationn (event, file, vertex index, etc.). This will be the input to run neutron selection and light study.

## Running the neutron selection
Look up for neuton-induced protons. 
```bash
cd neutron-multiplicity/analysis/n_selection/
source macro.sh
./select_n <version> <MODE>
```
for example:
```bash
source macro.sh
./select_n 6.2 RHC
```
This generates the reconstructed protons list files with their info in root and txt formats. Generate the plots from root files. For rhc files:
```bash
cd analysis/n_selection/plots/rhc
root -l -b -q 'plot_particles.C("/exp/dune/app/users/gkufatty/dune_projects/neutrino_selection/analysis/n_selection/output_files/protons_MiniRun<version>_RHC.root")'
```
for example, 
```bash
cd analysis/n_selection/plots/rhc
root -l -b -q 'plot_particles.C("/exp/dune/app/users/gkufatty/dune_projects/neutrino_selection/analysis/n_selection/output_files/protons_MiniRun6.2_RHC.root")'
```

## Running the light study
