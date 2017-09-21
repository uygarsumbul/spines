#!/bin/sh
#
#SBATCH --account=stats          # The account name for the job.
#SBATCH --job-name=heatMap1       # The job name.
#SBATCH -c 1                     # The number of cpu cores to use.
#SBATCH --time=299:00            # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=5gb        # The memory the job will use per cpu core.
#SBATCH -e arrayjob1-%a.err       # Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
#SBATCH -o arrayjob1-%a.out       # Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID

module load matlab/2015a
matlab -nodisplay -r "allHeatMapsAndSpecialNodes_function1, exit"

date
# End of script