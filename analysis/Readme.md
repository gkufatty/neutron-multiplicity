### How to run the analysis
There are two parts. Neutrino selection and the neutrons selection. 
First run the neutrino selection: 
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