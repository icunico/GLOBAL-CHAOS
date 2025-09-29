#!/bin/bash
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_prod
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 48
#SBATCH --time=04:00:00
#SBATCH --mem=300gb
#SBATCH --output=%j.out
#SBATCH --error=%j.err
eval "$(/g100/home/userexternal/gocchipi/my_libraries/conda-2023/bin/conda shell.bash hook)"
conda activate seamless-1d
export PYTHONUNBUFFERED=TRUE

# Generate secret for Parallel Python authentication between server and client
secret=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 32`

# Pick IP port for Parallel Python communication (IANA: use 49152 to 65535 for dynamic/private ports)
port=`python -c 'import random; print(random.randint(49152,65535))'`

# run the servers
echo Starting ppserver on $SLURM_JOB_NODELIST:$port
srun --output %j-%3t.out --error %j-%3t.err ppserver -w $SLURM_CPUS_PER_TASK -t 3600 -s $secret -p $port -r &

hostnames=`scontrol show hostnames $SLURM_JOB_NODELIST | sed -z 's/\n/:'"$port"',/g'`
echo Host names provided to parsac: $hostnames

# launch the parsac job
echo Starting parsac
WRKDIR=$CINECA_SCRATCH/BOUSSOLE_CALIBRATION/
xml_file=BOUSSOLE_calibration.xml
parsac calibration run --ppservers $hostnames --ncpus 0 --secret=$secret --tempdir $WRKDIR --maxiter 5 $xml_file
echo Done
