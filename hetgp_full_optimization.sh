#!/bin/bash
#SBATCH --job-name=OMFittingMultiSeed

#SBATCH --mem=1024G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=11

#SBATCH --time=368:00:00
#SBATCH --qos=infinite

#SBATCH --output=out.txt
#SBATCH --error=err.txt

ml R/3.6.3-foss-2018b

R CMD BATCH  /scicore/home/smith/reiker/GitRepos/om_fitting/hetgp_full_optim_small.r
