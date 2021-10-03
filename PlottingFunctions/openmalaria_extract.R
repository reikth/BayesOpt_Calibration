################################################################################
# OPENMALARIA EXTRACT
#
# Extract OpenMalaria metrics from series of simulation output files. The default
# behaviour is to (sum) aggregate metrics by time and age, resulting in one single
# output per metric per simulation. This behaviour can be changed (see input 
# descriptions below), resulting in a long-format dataframe.
#
# Aside from any OpenMalaria metric*, this function can also compute the
# following 'special' metrics:
#  prevalence        - Overall population prevalence
#  prevalence_2to10  - Prevalence among 2-10 year olds
#  daly_direct       - DALYs from direct malaria mortality and morbidity
#  daly_all          - DALYs from all malaria-related mortality and morbidity
#
# *NOTE: See https://github.com/SwissTPH/openmalaria/wiki/XmlMonitoring
#
# Required inputs:
#  files       - Vector of full paths to OpenMalaria simulation outputs
#  metrics     - Vector of OM indicators to extract
#  age_groups  - Vector of age group upper bounds (first group bounded below by 0)
#  
# Optional inputs:
#  time_points      - Vector of time points to extract data for [Default: all time points]
#  sum_age          - Bool for summing results by age [TRUE]
#  sum_time         - Bool for summing results by time [FALSE]
#  quiet            - Bool for turning off messaging to console (aside from fatal errors) [TRUE]
#  discounting*     - Discounting rate per year (value between 0 and 1) [0]
#  age_weighting*   - Bool for turning on age weighting for DALY calculations [FALSE]
#  expected_deaths* - Bool for 'expected' mortality in DALY calculations [FALSE]
#
# *NOTE: Inputs only used to calculate DALYs
#
# Example usage:
#  df = run_extract(files = paste0(sim_path, sim_file_names, ".txt"),
#                   metrics = c("prevalence_2to10", "daly_direct", "nSevere", "nDirDeaths"),
#                   age_groups = c(1, 2, 5, 10, 40, 70, 100),
#                   discounting = 0.03, expected_deaths = TRUE)
#
# Adapted from code witten by K.Galactionova (e.galactionova@unibas.ch)
# which itself was based on work by M.Penny (melissa.penny@unibas.ch)
#
# Written by A.J.Shattock - andrewjames.shattock@unibas.ch
################################################################################

source("openmalaria_metrics.R")

pacman::p_load(plyr,dplyr,data.table,stringr)

# ---------------------------------------------------------
# Extract outputs from a load of OpenMalaria simulations.
# ---------------------------------------------------------
run_extract = function(files, metrics, age_groups, ...) {
  
  # Turn off dplyr grouping messages 
  options(dplyr.summarise.inform = FALSE)
  
  # Initiate list to store sim outputs
  sims_list = list()
  
  # Initiate a progress bar
  pb = start_progress_bar(length(files))
  
  # Loop through the simulation files
  for (i in 1 : length(files)) {
    file = files[i]
    
    # Throw an error if this file does not exist
    if (!file.exists(file))
      stop("Simulation '", file, "' not found")
    
    # Extract the metrics of interest (nested function)
    sims_list[[i]] = sim_extract(file, metrics, age_groups,sum_age = FALSE) #, ... ## TR edit after Error:"'...' used in an incorrect context" + added sum_age=FALSE
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
  
  # Bind list of dataframes into single dataframe
  sims_df = as.data.frame(rbindlist(sims_list))
  
  return(sims_df)
}

# ---------------------------------------------------------
# Extract outputs from a single OpenMalaria simulation.
# ---------------------------------------------------------
sim_extract = function(file, metrics, age_groups, time_points = NULL, 
                       sum_age = TRUE, sum_time = FALSE, 
                       discounting = 0, age_weighting = FALSE, 
                       expected_deaths = FALSE, quiet = TRUE) {
  
  # Extract simulation name/id from file path
  sim_name = gsub(".txt", "", basename(file))
  
  if (!quiet) message("Extracting from", sim_name)
  
  # ---- Which metrics to extract ----
  
  # Special metrics we can calculate
  daly_metrics = c("daly_direct", "daly_all")
  prev_metrics = c("prevalence", "prevalence_2to10")
  
  # Flags for running DALY / prevalence calculations
  do_daly = any(daly_metrics %in% metrics)
  do_prev = any(prev_metrics %in% metrics)
  
  # Which raw metrics do we need to extract
  raw_metrics = setdiff(metrics, c(daly_metrics, prev_metrics))
  
  # Raw metrics needed for prevalence calculations
  if (do_prev) raw_metrics = c(raw_metrics, "nHost", "nPatent")
  
  # Raw metrics needed for DALY calculations
  if (do_daly) {
    raw_metrics = c(raw_metrics, "nUncomp", "nSevere", "nSeq")
    
    # Two approaches possible for counting deaths
    if (expected_deaths == FALSE) death_metrics = c("nDirDeaths", "nIndDeaths")
    else death_metrics = c("expectedDirectDeaths", "expectedIndirectDeaths")
    
    # Append the chosen death indicators into metric vector
    raw_metrics = c(raw_metrics, death_metrics)
  }
  
  # Remove any duplication
  raw_metrics = unique(raw_metrics)
  
  # ---- DALY equation parameter values ----
  
  # NOTE: All DALY parameter values taken from Katya's Stata code - I'll try to 
  #       find references to the original papers
  
  if (do_daly) {
    
    # Initiate lists to store parameters and variables for DALY calculations
    y = u = n = p = list() # y := years; u := utility weight; n := number;
    
    # Years of life associated with type of disability
    y$mild     = 0.057534247 # Uncomplicated cases
    y$moderate = 0.057534247 # Sequelae cases
    y$severe   = 0.05        # Severe cases
    
    # DALY utility weight for each type of disability
    u$mild     = 0.005 # Uncomplicated cases
    u$moderate = 0.053 # Sequelae cases
    u$severe   = 0.210 # Severe cases
    
    # Proportional split of uncomplicated into mild and moderate
    p$mild     = 0.663 / (0.663 + 0.332)
    p$moderate = 0.332 / (0.663 + 0.332)
    
    # Age weighting parameters
    b = 0.04   # Age-weighting function parameter
    c = 0.1658 # Adjustment constant for age weights
    
    # Constants common to YLL and YLD
    r = discounting       # Discount rate
    k = age_weighting * 1 # Age-weighting toggle
  }
  
  # Difference between bin edges for all age groups
  age_diff = diff(c(0, age_groups))
  
  # Use the upper bound vector to determine lower bound thresholds
  age_lower = c(0, age_groups[1 : (length(age_groups) - 1)])
  
  # Convert into mean of the age bins
  age_mid = age_lower + age_diff / 2
  
  # ---- Load simulation results and extract metrics of interest ----
  
  # Load simulation result
  sim = read.table(file, sep = "\t")
  colnames(sim) = c("time", "age_group", "metric", "value")
  
  # Load full list of possible OpenMalaria metrics (aka 'events' or 'outputs')
  all_metrics = openmalaria_metrics()  # See openmalaria_metrics.R
  
  # If no time input given, use all timepoints by default
  if (is.null(time_points)){time_points = unique(sim$time)}
  
  # TODO: Throw error if user-defined time indices outside of those available
  
  # Initiate output value list
  out_list = list()
  
  # Iterate through raw metrics (aka 'measures' or 'events')
  for (metric_name in raw_metrics) {
    metric_values = all_metrics[[metric_name]]
    
    # Row indices of OM results table associated with this metric
    #
    # NOTE : Use is_element here as metric index could actually be a vector of several indicies
    metric_idx = is.element(sim$metric, metric_values)
    
    # Extract time, age group, and value of this metric
    metric_table = sim[metric_idx, c("time", "age_group", "value")]
    
    # Reduce table to only time points of interest (between current_year and analysis_year)
    metric_table = metric_table[is.element(metric_table$time, time_points), ]
    
    # If mutliple metrics (as with allDeaths) we need some groupby action
    if (length(metric_values) > 1) {
      
      # Group by multiple columns (age_group, time)
      metric_table = metric_table %>% 
        dplyr::group_by(age_group, time) %>% 
        dplyr::summarise(value = sum(value)) %>% 
        as.data.frame()
    }
    
    # If metric does not have age-segregation
    if (all(metric_table$age_group == 0)) {
      out_list[[metric_name]] = as.matrix(metric_table$value)
      
    } else {
      
      # Reshape into a matrix (of size n_time_pts x n_age)
      age_time_matrix = unstack(metric_table, value~age_group)
      colnames(age_time_matrix) = age_mid
      
      # Store matrix in output array
      out_list[[metric_name]] = as.matrix(age_time_matrix)
    }
  }
  
  # Additional formatting required for DALY calculations
  if (do_daly) {
    
    # Rename death outcomes for generalisability
    out_list$death_direct = out_list[[death_metrics[1]]]
    out_list$death_all    = out_list[[death_metrics[2]]] + out_list$death_direct
    
    # Rename indicators used to quantify mild, moderate, and severe cases
    n$mild     = out_list$nUncomp * p$mild
    n$moderate = out_list$nUncomp * p$moderate
    n$severe   = out_list$nSevere - out_list$nSeq - out_list$death_direct
  }
  
  # ---- DALY calculations ----
  
  # Initiate output dataframe
  out_df = NULL
  
  # Check flag for DALY calculations
  if (do_daly) {
    
    # Age matrix (of size n_timepoints x n_age_groups)
    a = t(matrix(age_mid, nrow = length(age_groups), ncol = length(time_points)))
    
    # Life expetency for each age
    l = life_expectancy(a)
    
    # Years of life with disability (YLD)
    years_life_disability = 
      calculate_YLD(n, y, u, a, r, k, c, b, "mild") + 
      calculate_YLD(n, y, u, a, r, k, c, b, "moderate") +
      calculate_YLD(n, y, u, a, r, k, c, b, "severe")
    
    # Years of life lost (YLL) ('direct' and 'all' malaria-related mortality)
    years_life_lost_all    = calculate_YLL(a, l, r, k, c, b, out_list$death_all)
    years_life_lost_direct = calculate_YLL(a, l, r, k, c, b, out_list$death_direct)
    
    # Finally, combine to calculate 'all' and 'direct' DALYs
    daly_all    = years_life_disability + years_life_lost_all
    daly_direct = years_life_disability + years_life_lost_direct
    
    # Only output metrics explictly requested - DALY all
    if ("daly_all" %in% metrics)
      out_df = format_df(out_df, daly_all, "daly_all")
    
    # Only output metrics explictly requested - DALY direct
    if ("daly_direct" %in% metrics)
      out_df = format_df(out_df, daly_direct, "daly_direct")
  }
  
  # ---- Prevalence calculations ----
  
  # Check flag for DALY calculations
  if (do_prev) {
    
    c("prevalence", "prevalence_2to10")
    
    # Do we need to calculate all age prevalence?
    if ("prevalence" %in% metrics) {
      
      # Number of people and number with malaria
      n_malaria = rowSums(out_list$nPatent)
      n_people  = rowSums(out_list$nHost)
      
      # Overall prevalece (all ages)
      prevalence = as.matrix(n_malaria / n_people)
      
      # Format and append into output dataframe
      out_df = format_df(out_df, prevalence, "prevalence")
    }
    
    # Do we need to calculate PfPR2-10?
    if ("prevalence_2to10" %in% metrics) {
      
      # Indices of people aged between 2 and 10
      age_2to10 = (o$sim$age_groups >= 2 & o$sim$age_groups <= 10)
      
      # TODO: Throw an error/warning if poorly defined age groups for PfPR2-10
      
      # Number of people and number with malaria (aged 2-10)
      n_malaria = rowSums(out_list$nPatent[, age_2to10])
      n_people  = rowSums(out_list$nHost[, age_2to10])
      
      # Prevalece among 2-10 year olds
      prevalence_2to10 = as.matrix(n_malaria / n_people)
      
      # Format and append into output dataframe
      out_df = format_df(out_df, prevalence_2to10, "prevalence_2to10")
    }
  }
  
  # ---- Append metrics and format output ----
  
  # Set of matrix indicators to additionally append
  append_metrics = setdiff(metrics, unique(out_df$metric))
  
  # Iterate through and append each indicator
  for (append_metric in append_metrics)
    out_df = format_df(out_df, out_list[[append_metric]], append_metric)
  
  # Rename dataframe columns
  setnames(out_df, "Var1", "time")
  setnames(out_df, "Var2", "age")
  
  # Add column explaining simulation name
  out_df$sim = sim_name
  
  # Reorder dataframe columns
  out_df = out_df[, c("sim", "metric", "time", "age", "value")]
  
  # Check flag for summing by age
  if (sum_age == TRUE) {
    
    # Group and sum value by all but age
    out_df = out_df %>% 
      dplyr::group_by_at(vars(-age, -value)) %>% 
      dplyr::summarise(value = sum(value)) %>% 
      as.data.frame()
  }
  
  # Check flag for summing by time
  if (sum_time == TRUE) {
    
    # Group and sum value by all but time
    out_df = out_df %>% 
      dplyr::group_by_at(vars(-time, -value)) %>% 
      dplyr::summarise(value = sum(value)) %>% 
      as.data.frame()
  }
  
  return(out_df)
}

# ---------------------------------------------------------
# Function to compute life expectancy.
# ---------------------------------------------------------
life_expectancy <- function(a) 52.993 - 0.6665*a + 0.0012*a^2

# ---------------------------------------------------------
# Function to compute years of life with disability (YLD).
# ---------------------------------------------------------
calculate_YLD = function(n, y, u, a, r, k, c, b, level) {
  
  # Extract values for this level
  n_cases = n[[level]]
  utility = u[[level]]
  y       = y[[level]]
  
  # I've split this into terms for the sake of readability and debugging
  term_1 = k*c*exp(r*a) / (r+b)^2
  term_2 = ((r+b) *  a    + 1) * exp(-(r+b) *  a)
  term_3 = ((r+b) * (a+y) + 1) * exp(-(r+b) * (a+y))
  term_4 = ((1-k) / r) * (1 - exp(-r*y))
  
  if (r == 0) term_4 = y
  
  # Combine the terms together for years of life with disability
  yld = n_cases * utility * (term_1 * (term_2 - term_3) + term_4)
  
  return(yld)
}

# ---------------------------------------------------------
# Function to compute years of life lost (YLL).
# ---------------------------------------------------------
calculate_YLL = function(a, l, r, k, c, b, n_deaths) {
  
  # Again, split this into terms for the sake of readability
  #
  # NOTE: The equation is exactly the same form as above so
  #       should be generalised to a single function
  term_1 = k*c*exp(r*a) / (r+b)^2
  term_2 = ((r+b) *  a    + 1) * exp(-(r+b) *  a)
  term_3 = ((r+b) * (a+l) + 1) * exp(-(r+b) * (a+l))
  term_4 = ((1-k) / r) * (1 - exp(-r*l))
  
  if (r == 0) term_4 = l
  
  # Combine the terms together for years of life lost
  yll = n_deaths * (term_1 * (term_2 - term_3) + term_4)
  
  return(yll)
}

# ---------------------------------------------------------
# Format metric and append into output dataframe
# ---------------------------------------------------------
format_df = function(out_df, metric_val, metric_name) {
  
  # Shape the (possibly trivial) matrix into dataframe
  metric_df = reshape2::melt(metric_val)
  
  # Trivialise age values if appropriate
  if (all(metric_df$Var2 == 1))
    metric_df$Var2 = NA
  
  # Append metric name
  metric_df$metric = metric_name
  
  # Bind to other results
  out_df = rbind(out_df, metric_df)
  
  return(out_df)
}

