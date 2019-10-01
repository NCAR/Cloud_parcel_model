#!/bin/bash
#SBATCH --time=00:5:00		#walltime
#SBATCH --job-name=parcel_err	
#SBATCH --output=stdout_graham
#SBATCH --error=output_error
#SBATCH --account=def-yaumanko
#SBATCH --ntasks=1 		#number of MPI processes
#SBATCH --mem-per-cpu=1024M 	#memory per unit

srun ./cloud_parcel
