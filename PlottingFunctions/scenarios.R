############################################################
# SCENARIOS
#
# Create, simulate, and post-process all OM scenarios.
#
# Written by A.J.Shattock - andrewjames.shattock@unibas.ch
############################################################

source("convert_access.R")
source("openmalaria_setup.R")
source("openmalaria_run.R")

pacman::p_load(stringr,utils)


# ---------------------------------------------------------
# Parent function for creating and simulating all scenarios.
# ---------------------------------------------------------
run_scenarios = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  # Generate full list of scenarios
  sim_list <<- generate_scenario_list(o)
  
  # Number of these scenarios
  n_sim <<- thou_sep(nrow(sim_list$sim_df))
  
  message("* Simulating ", n_sim, " scenarios")
  
  # When debugging/experimenting we may want to only postprocess
  if (!o$post$postprocess_only) {
    
    # Generate sed pattern replacements for full simulation xml files
    generate_sed_patterns(o, sim_list$param_df)  # See openmalaria_setup.R

    # Generate scenario xml files from the base xml and sed patterns
    generate_xml_files(o, sim_list$param_df,
                       sed_path = o$pth$sim_sed,
                       xml_path = o$pth$sim_xml)  # See openmalaria_setup.R
    
    # Run OpenMalaria for the calibration EIR samples
    run_model(o, sim_list$sim_df,
              xml_path = o$pth$sim_xml,
              sim_path = o$pth$sim_out)  # See openmalaria_run.R
  }
  
  # Post process all OpenMalaria simulations and summarise over seeds
  postprocess(o, sim_list$sim_df)  # See openmalaria_run.R
}

# ---------------------------------------------------------
# Generate list of scenarios to run in dataframe form.
# ---------------------------------------------------------
generate_scenario_list = function(o) {
  
  # Initiate parameter set df with all EIR and access combinations
  param_df = merge(o$sim$eir, o$sim$access)
  
  # Number of seasons and seeds
  n_season = ncol(o$sim$seasons)
  n_seeds  = o$sim$n_seeds
  
  # Repeat parameter sets for each season and seed
  param_df = merge(param_df, as.data.frame(1 : n_season))
  param_df = merge(as.data.frame(1 : n_seeds), param_df)
  
  # ---- Dataframe for unique parameter sets ----
  
  # Rename columns of dataframe
  names(param_df) = c("seed", "eir", "access", "season")
  
  # Format sim IDs for all parameter sets
  param_df$sim_id = get_sim_ids(o, param_df)
  
  # Reorder columns of parameter set dataframe
  param_df = param_df[, c("sim_id", "eir", "access", "season", "seed")]
  
  # ---- Dataframe of all simulations (inc. model variant) ----
  
  # Repeat unique parameter set dataframe for each model variant
  sim_df = merge(data.frame(model = o$sim$models), param_df)
  
  # Update sim IDs to include model ID for all simulations
  sim_df$sim_id = get_sim_ids(o, sim_df)
  
  # Reorder columns of simulation dataframe
  sim_df = sim_df[, c("sim_id", "eir", "access", "season", "seed", "model")]
  
  # Append these two dataframes to an output list
  sim_list = list(sim_df   = sim_df, 
                  param_df = param_df)
  
  return(sim_list)
}

# ---------------------------------------------------------
# Convert scenario details to sim IDs.
# ---------------------------------------------------------
get_sim_ids = function(o, df) {
  
  # Convert EIR value to EIR index
  eir_idx = match(df$eir, o$sim$eir)
  
  # Format sim ID components with sufficient padding
  eir    = paste0("e", sprintf("%02i", eir_idx))
  access = paste0("a", sprintf("%02i", df$access * 100))
  season = paste0("z", sprintf("%01i", df$season))
  seed   = paste0("s", sprintf("%02i", df$seed))
  
  # Bind columns into a dataframe
  id_df = data.frame(eir, access, season, seed)
  
  # Include model variant column if necessary
  if (any(names(df) == "model"))
    id_df$model = paste0("m", df$model)

  # Collapse all columns to create vector of sim IDs
  sim_ids = unite(id_df, "x", sep = "")$x
  
  return(sim_ids)
}

