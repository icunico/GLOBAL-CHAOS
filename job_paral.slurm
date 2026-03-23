#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mem=30gb
#SBATCH --job-name=lat_Trends
#SBATCH --account OGS23_PRACE_IT_0
#SBATCH --partition g100_usr_prod
##SBATCH --qos=g100_qos_dbg

eval "$(/g100/home/userexternal/icunico0/programmi/ls/bin/conda shell.bash hook)"
#conda activate copernicus
conda activate conda-sto

module load intel

date
# Con 5 processi 
mpirun -np 5 ./tau_areasfolta /g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/regional_dataset_CHL_10_20_-30_-15.nc CHL
#python lyapunov_Guido.py regional_dataset_CHL_10_20_-30_-15.nc CHL
#mpirun -np 4 python tempoWolf.py
date

