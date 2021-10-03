############################################################
# OPTIONS
#
# Define simulation settings and assumptions.
#
# Written by A.J.Shattock - andrewjames.shattock@unibas.ch
############################################################

source("directories.R")

pacman::p_load(plyr,dplyr)


# ---------------------------------------------------------
# Set a series of options for the analysis.
# ---------------------------------------------------------
set_options = function(do_step = c(1,2), quiet = FALSE) {
  
  if (quiet == FALSE) message("* Setting options")
  
  # Reset R's most annoying default options
  options(stringsAsFactors = FALSE, scipen = 999,
          dplyr.summarise.inform = FALSE)
  
  # Initiate options list
  o = list(do_step = do_step)
  
  # Detect user - required for checking cluster jobs
  o$user = Sys.info()[["user"]]
  
  # Set directories
  o = set_dirs(o)  # See directories.R
  
  # ---- Time settings ----
  
  # Years of analysis
  o$time$start_year = 2000
  o$time$end_year   = 2020
  
  # Number of years for 'secondary' burn in
  burn_in_years = 30
  
  # Set the burn in relative to start year
  o$time$burn_in = o$time$start_year - burn_in_years
  
  # ---- Scenario and seasonality settings ----
  
  # Vector of EIR values to run
  o$sim$eir = c(0.25,0.5,0.75,1,1.1,1.25,1.35,1.5,1.75,2,2.5,3,4,5,6,7,8,9,10,12,14,16,18,20,22,25,30,35,40,45,50,64,73,80,100, 128, 150,200,256, 512) #2 ^ c(-1 : 10)  # 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024
  
  # Access to care
  o$sim$access = 0.44
  
  # Run different model variants
  o$sim$models = c("0000GA", "0000GP", "0000GPSG")# , "0068")#, "0131", "0132", "0133", "0670")
  
  # Define functional form of non-perennial, seasonal setting
  season_daily = 1 + sin(2 * pi * ((1 : 365) / 365))
  season_idx   = round(1 + seq(0, 365, length.out = 13))[-13] # First day of month index
  
  # Index for first of month and normalise
  season_month = season_daily[season_idx]
  season_month = season_month / max(season_month)
  
  # Append to a dataframe with a perennial setting
  o$sim$seasons = data.frame(z1 = 1, z2 = season_month)
  
  # ---- Simulation settings ----
  
  # OpenMalaria population size
  o$sim$pop_size = 50000
  
  # Age group upper bounds (first bin bounded below by 0)
  o$sim$age_groups = c(0.5,1,2,5,10,15,20,100)
  
  # Number of seeds
  o$sim$n_seeds = 10
  
  # Proportion of mosquito biting outdoors
  outdoor_biting = 0.2
  
  # Paramter values for outdoors/indoors mosquito biting
  o$sim$outdoor = outdoor_biting
  o$sim$indoor  = 1 - outdoor_biting
  
  # ---- Post-processing settings ----
  
  # Rerun only postprocessing for simulation - often useful when debugging/experimenting
  o$post$postprocess_only = TRUE  # For normal functionality this should be set to false
  
  # Define which indicators to extract from simulations
  o$post$metrics = c("prevalence_2to10", "inputEIR", "simulatedEIR","nHost","nPatent","nUncomp", "nSevere","nTreatments1","expectedSevere", "expectedDirectDeaths", "expectedIndirectDeaths")
  
  # ---- Plotting settings ----
  
  # Set plotting flags
  o$plot$eir_error = TRUE
  o$plot$eir2prev  = TRUE
  
  return(o)
}

