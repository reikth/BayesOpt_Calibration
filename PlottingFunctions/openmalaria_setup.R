############################################################
# OPENMALARIA SETUP
#
# Prepare for OpenMalaria simulations.
#
# Written by A.J.Shattock - andrewjames.shattock@unibas.ch
############################################################

source("convert_access.R")

pacman::p_load(pracma,stringr,tgp,qdap,roperators,tidyr,xlsx)


# ---------------------------------------------------------
# Generate sed replacement patterns for OpenMalaria xml files.
# ---------------------------------------------------------
generate_sed_patterns = function(o, param_df) {
  
  message("  - Generating scenarios")
  
  # ---- Format fixed parameters for this country ----
  
  # Define list of constant parameters to be set in all xml files
  constant_params = c("pop_size", "outdoor", "indoor")
  time_params     = c("start_year", "end_year", "burn_in")
  
  # Extract constant and efficacy parameters from o list
  constant_values = extract_param_values(o$sim,  constant_params)
  time_values     = extract_param_values(o$time, time_params)
  
  # End year needs to be modified for the xml file
  #
  # NOTE: +1 for end of year, then +1 as last time point is not reported
  time_values[time_params == "end_year"] %+=% 2
  
  # Convert access to care to 5 day probabilties for use in xml files
  param_df$access = pmax(convert_access(param_df$access * 100), 0.04)
  
  # Parameter names for seasonality values 
  season_params = paste0("seasonality", 1 : 12)
  
  # ---- Loop through scenarios ----
  
  # Number of scenarios
  n_scenarios = nrow(param_df)
  
  # Initiate progress bar
  pb = start_progress_bar(n_scenarios)
  
  # Loop through scenarios
  for (i in 1 : n_scenarios) {
    s = param_df[i, ]
    
    # Extract seasonality values from options dataframe
    season_values = o$sim$seasons[[paste0("z", s$season)]]
    
    # Order initial parameter values: EIR then seed
    param_names  = c("eir", "access", "seed", season_params)
    param_values = c(s$eir, s$access, s$seed, season_values)
    
    # Then append constants and time parameters
    param_names  = c(param_names,  constant_params, time_params)
    param_values = c(param_values, constant_values, time_values)
    
    # Format parameter values to not use scientific notation
    param_format = format(param_values, scientific = FALSE, 
                          trim = TRUE, drop0trailing = TRUE)
    
    # Create sed command replacement pattern for each parameter
    sed_patterns = paste0("s$@", param_names ,"@$", param_format, "$g;")
    
    # ---- Save file to be used by bash script ----
    
    # Save these sed patterns to a text file
    write.table(sed_patterns, file = paste0(o$pth$sim_sed, s$sim_id, ".txt"), 
                quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
}

# ---------------------------------------------------------
# Extract parameter values from o sub-list.
# ---------------------------------------------------------
extract_param_values = function(o_list, params) {
  
  # Preallocate vector
  values = rep(NA, length(params))
  
  # Loop through parameters to extract values
  for (i in 1 : length(params)) {
    
    # Extract value from o sub list
    values[i] = o_list[[params[i]]]
  }
  
  return(values)
}

# ---------------------------------------------------------
# Generate scenario xml files from base file using sed patterns.
# ---------------------------------------------------------
generate_xml_files = function(o, param_df, sed_path = NA, xml_path = NA) {
  
  message("  - Generating xml files")
  
  # Number of scenarios
  n_param_sets  = nrow(param_df)
  n_models      = length(o$sim$models)
  n_simulations = n_param_sets * n_models
  
  # Initiate progress bar
  pb = start_progress_bar(n_simulations)
  
  # Loop through parameter sets
  for (i in 1 : n_param_sets) {
    param_id = param_df[i, ]$sim_id
    
    # Loop through model variants
    for (j in 1 : n_models) {
      model = o$sim$models[j]
      
      # Construct full sim ID
      sim_id = paste0(param_id, "m", model)
      
      # Concatenate file names for files representing this simulation
      sed_file = paste0(sed_path, param_id, ".txt")
      xml_file = paste0(xml_path, sim_id, ".xml")
      
      # Point to xml base file for this model variant
      base_file = paste0(o$pth$input, "R", model, ".xml")
      
      # Replace all sed patterns and save to a new xml file
      sed_command = paste("sed -f", sed_file, base_file, ">", xml_file)
      
      # Invoke this command
      system(sed_command)
      
      # Update progress bar
      setTxtProgressBar(pb, n_models * (i - 1) + j)
    }
  }
  
  # Close progress bar
  close(pb)
}

