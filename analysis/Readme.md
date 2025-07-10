# How to run the analysis
There are three parts. Neutrino selection, neutron selection, and light study. If you are using a MiniRun version that is not in the folder, then there is a pre-step to generate the list of files you will use

##Analysis Setup

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
If you want to use a new MiniRun version, generate the files
```bash
source neutron-multiplicity/analysis/nu_selection/input_files/macro.sh
<folder_path> <version> <mode>
```
for example 
```bash
./generate_input_list /exp/dune/data/users/noeroy/prod/MiniRun6.2_1E19_RHC/MiniRun6.2_1E19_RHC.caf/CAF/ 6.2 RHC
```

```bash

```

```bash

```

```bash

```

## Running the neutrino selection
source compile_macro in analysis/nu_selection folder 
it generates a select_nus file. 
./select_nus -i input_file -o output_file -c config_file
for example
```bash
./select_nus 6.3 RHC    
```
Then run the neutrons selection:
```bash
./select_neutrons -i input_file -o output_file -c config_file


## Running the neutron selection


## Running the light study
