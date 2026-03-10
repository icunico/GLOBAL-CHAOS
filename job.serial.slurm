#!/bin/bash
#SBATCH --job-name=CHAOS
#SBATCH --time=24:00:00
#SBATCH --account=OGS_test2528_0
#SBATCH -p dcgp_usr_prod
#SBATCH --mail-type=ALL
##SBATCH --mail-user=plazzari@ogs.it
#SBATCH --mem=30G

ulimit -s unlimited

date

export RANKS_PER_NODE=1

cd /leonardo_work/OGS_test2528_0/plazzari/GLOBAL-CHAOS/

# Run the script with the specific year
python Lorenz_embedding.py 
date
