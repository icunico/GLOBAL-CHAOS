#!/bin/bash
#SBATCH --job-name=CHAOS
#SBATCH --time=24:00:00
#SBATCH --account=OGS23_PRACE_IT_0
#SBATCH -p boost_usr_prod
#SBATCH --mail-type=ALL
##SBATCH --mail-user=plazzari@ogs.it
#SBATCH --mem=30G

ulimit -s unlimited

date

export RANKS_PER_NODE=1

cd /leonardo/home/userexternal/icunico0/seamless-notebooks/setups/GLOBAL_CHOAS/GLOBAL-CHAOS

# Run the script with the specific year
python Lorenz_embedding.py 
date
