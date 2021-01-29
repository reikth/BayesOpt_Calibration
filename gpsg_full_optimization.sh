#!/bin/bash
#SBATCH --job-name=OMFittingGPSGValidation

#SBATCH --mem=900G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=11

#SBATCH --time=168:00:00
#SBATCH --qos=1week

#SBATCH --output=out.gpsg.txt
#SBATCH --error=err.gpsg.txt

ml R/3.6.3-foss-2018b

R CMD BATCH  /scicore/home/smith/reiker/GitRepos/om_fitting/2020_01_Full_BayesOpt_GPSG_mars_brnn_ranger.R
