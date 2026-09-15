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
Look for neutron-induced proton candidates using a neutrino-selection CSV.
The CSV must contain `file_name,event,reco_ixn,truth_ixn,has_truth_match,matched_true_signal,true_track_multiplicity,reco_track_multiplicity,true_primary_neutron_count,true_secondary_neutron_count`.
Columns are read by name; quoted and unquoted fields are supported.
The older five-column QE-like interaction files need conversion or regeneration.
```bash
cd neutron-multiplicity/analysis/n_selection/
source macro.sh
./select_n <selection.csv> <version> <MODE>
```
for example:
```bash
source macro.sh
./select_n /global/homes/l/lmlepin/2x2_trackMultStudies/spineAna/minirun6.5_neutrino_selection_09-14-2026.csv 6.5 RHC
```
The selector matches `event` to `rec.meta.nd_lar.event` within each CAF file.
Rows with `event == -1` are temporarily skipped before CAF lookup. The terminal
summary reports skipped rows separately from processed interactions and those
with zero candidates. A file referenced only by skipped rows is not opened;
an all-skipped input produces empty outputs successfully. Other missing or
ambiguous event IDs still cause an error. Event IDs are never treated as entry
numbers.

The text output and the ROOT `protons` tree contain one row per reconstructed
proton candidate. The input interaction's truth flag and four multiplicities
are preserved. The ROOT `event` branch is a signed 64-bit integer.
`has_particle_truth_match` indicates whether particle truth quantities are
available; unmatched candidates remain in the output with invalid truth
indices and unavailable scalar quantities.

The same ROOT file also contains a separate `true_neutron_protons` tree. It has
one row per direct neutron-induced true secondary proton whose start position
passes the fiducial-volume cut. The true end position is recorded, but it is not
used to reject truth rows. Event-wide reverse-association and `select_n`
selection-status branches distinguish protons with no reco overlap, protons
reconstructed but rejected by the candidate cuts, and selected candidates.
`matched` retains the input signal flag; `insignal` still uses the existing
QE-like definition for the particle's truth interaction. Interaction counts are
repeated per candidate and must not be summed over candidate rows. Selected
interactions yielding zero candidates are counted in the terminal summary.

Run the ROOT-independent input tests from the repository root:
```bash
bash analysis/n_selection/tests/run_tests.sh [selection.csv]
```
With ROOT and duneanaobj configured, also run the synthetic CAF integration tests:
```bash
bash analysis/n_selection/tests/run_root_tests.sh
```
The inspected MiniRun6.5 CSV contains 25 rows with `event == -1`. These are
excluded, leaving 8,964 interactions eligible for processing. To recover skipped
rows in a future analysis, export an unambiguous identifier and extend the reader
to use it. See `n_selection/CSV_COMPATIBILITY_REPORT.txt` for validation results
and the exact scope of the change.

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
