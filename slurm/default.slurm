#!/bin/bash
#SBATCH --job-name=rare
#SBATCH --output=logs/cluster_logs/%x-%j.out
#SBATCH --mail-user=pschloss@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --account=pschloss1
#SBATCH --partition=standard
#SBATCH --time=110:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g

### Set Up Environment
### Need to access the internet
source /etc/profile.d/http_proxy.sh

### Conda environment:
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rare

### This chunk changes directory to the directory from which it was submitted
echo "Changing directory to $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR

### Run script
snakemake --profile config_files/

echo "done."
