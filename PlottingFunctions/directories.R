############################################################
# SET DIRECTORIES
#
# Set and get directories in one place in the name of consistency
# and ease. Also creates any directories that do not currently exist.
#
# Outputs a list of relevant directories (within o$pth) which 
# can be referenced elsewhere.
#
# Written by A.J.Shattock - andrewjames.shattock@unibas.ch
############################################################

source("myRfunctions.R")

# ---------------------------------------------------------
# Define paths for project inputs and outputs.
# ---------------------------------------------------------
set_dirs = function(o) {
  
  # Initiate file path lists
  pth_in = pth_out = list()
  
  # ---- Code and resource locations ----
  
  # Base path to code repositories
  base_stm = file.path("~/GitRepos")
  
  # Parent paths to all input files relating to this project
  pth_in$code  = file.path(base_stm, "mmc_age_plots")
  pth_in$input = file.path(pth_in$code, "input")
  
  # Path to OpenMalaria resource files
  pth_in$om_files = file.path(base_stm, "OM_schema38")
  
  # ---- Output directories and files ----
  
  # Parent path to all output files relating to this project
  stem_ouput = file.path(pth_in$code, "output")
  
  # Paths to simulation files
  pth_out$sim_sed = file.path(stem_ouput, "1_sed_files")
  pth_out$sim_xml = file.path(stem_ouput, "2_xml_files")
  pth_out$sim_out = file.path(stem_ouput, "3_om_output")
  
  # Paths to results and figures
  pth_out$results = file.path(stem_ouput, "4_results")
  pth_out$figures = file.path(stem_ouput, "5_figures")
  
  # ---- Create directories if needed ----
  
  # Make all output directories
  make_output_dirs(o, pth_out)
  
  # Initiate path field
  o$pth = list()
  
  # Append paths to o list
  o = append_dirs(o, pth_in)
  o = append_dirs(o, pth_out)
  
  return(o)
}

# ---------------------------------------------------------
# Make all output directories if they do not already exist.
# ---------------------------------------------------------
make_output_dirs = function(o, pth_out) {
  
  # Extract all path names in list
  pth_names = names(pth_out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = pth_out[[pth_name]]
    
    # If it does not already exist, create it
    if (!dir.exists(this_pth)) {
      dir.create(this_pth, recursive = TRUE)
    }
  }
}

# ---------------------------------------------------------
# Concatenate separators and append directories to o list.
# ---------------------------------------------------------
append_dirs = function(o, pth) {
  
  # Loop through these path names
  for (pth_name in names(pth)) {
    this_pth = pth[[pth_name]]
    
    # Add a file separator to end of output paths
    o$pth[[pth_name]] = paste0(this_pth, file_sep())
  }
  
  return(o)
}

