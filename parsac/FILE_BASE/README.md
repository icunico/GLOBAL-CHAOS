### Run Model Calibration 
- file needed: 
   1. BOUSSOLE_calibration.xml        -> configure calibration
   2. TotalCHLA.obs                   -> observation file 
   3. run_calibration_template.sbatch -> run on HPC machines
- run:
```
sbatch run_calibration_template.sbatch
```


### Run Model Sensitivity 
- file needed: 
   1. BOUSSOLE_sensitivity.xml        -> configure sensitivity
   2. sensitivity_sample.sh           -> launcher of the sensitivity sampling phase
   3. run_sensitivity_template.sbatch -> run sensitivity on HPC machines
   4. sensitivity_run.sh              -> launcher of the sensitivity run phase

- run:
  1. sampling step: `bash sensitivity_sample.sh `
  2. ensemble model run step: `bash sensitivity_run.sh`
  3. analysis step: `bash sensitivity_analyze.sh`

- results are stored in BOUSSOLE_sensitivity_CV.txt


