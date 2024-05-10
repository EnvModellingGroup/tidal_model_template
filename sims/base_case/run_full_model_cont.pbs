#!/bin/bash
#SBATCH --job-name=modern_uk                 # Job name
#SBATCH --mail-type=FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jon.hill@york.ac.uk      # Where to send mail
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1                    # Number of CPU cores per task
#SBATCH --time=4:00:00                       # Time limit hrs:min:sec
#SBATCH --output=thetis_%j.log               # Standard output and error log
#SBATCH --account=ENV-TSUNAMI-2019           # Project account

module load firedrake

unset PYTHONPATH


. /mnt/scratch/projects/env-tsunami-2019/firedrake_dec23/bin/activate
# we need to add the library path so the various .so files are picked up
# make sure the path matches the one above!
export LD_LIBRARY_PATH=/mnt/scratch/projects/env-tsunami-2019/firedrake_dec23/src/petsc/default/lib/:$LD_LIBRARY_PATH


export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1


mpiexec -n 32 python tidal_model_cont.py
