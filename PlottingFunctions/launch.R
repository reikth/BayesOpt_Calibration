############################################################
# LAUNCH
#
# Launch analyses and plots for Malaria Modelling Consortium
# benchmarking exercises.
#
# Written by A.J.Shattock - andrewjames.shattock@unibas.ch
############################################################

# Clear global environment
rm(list = ls())

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

if(!require(pacman)){install.packages("pacman")};require(pacman)
source("myRfunctions.R")
source("options.R")
source("scenarios.R")
source("results.R")

# Tidy up
if (interactive()) clc()  # Clear console
if (interactive()) clf()  # Close figures

# Set options (see options.R)
o = set_options(do_step =  1:2)

# 1) Run OpenMalaria for defined scenarios
run_scenarios(o)  # See scenarios.R

# 2) Produce results
run_results(o)  # See results.R

# Finish up
message("* Finished!")

