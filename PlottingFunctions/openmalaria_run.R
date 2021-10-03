############################################################
# OPENMALARIA RUN
#
# Run scenario xml files and perform necessary post-processing.
#
# Written by A.J.Shattock - andrewjames.shattock@swisstph.ch
############################################################

source("openmalaria_extract.R")

pacman::p_load(readr,stringr,tidyr,matrixStats)

# ---------------------------------------------------------
# Run OpenMalaria for all generated xml files.
# ---------------------------------------------------------
run_model = function(o, sim_df, xml_path = NA, sim_path = NA) {
  
  message("  - Running OpenMalaria")
  
  # Number of scenarios
  n_scenarios = nrow(sim_df)
  
  # Create a new log file for the cluster jobs
  log_file = create_bash_log(o$pth$code, "scicore_log.txt")
  
  # Construct slurm array command for running in parallel
  slurm_array = paste0("--array=1-", n_scenarios)
  
  # Concatenate system command
  sys_command = paste("sbatch", slurm_array, "bash_openmalaria.sh", 
                      xml_path, sim_path, o$pth$om_files, log_file)
  
  # Invoke this command
  system(sys_command)
  
  # Wait for all cluster jobs to complete
  wait_for_jobs(o, log_file, n_scenarios)
}

# ---------------------------------------------------------
# Perform full OpenMalaria post processing procedures.
# ---------------------------------------------------------
postprocess = function(o, sim_df) {
  
  message("  - Performing post-processing")
  
  # Concatenate full file paths for all simulations
  sim_files = paste0(o$pth$sim_out, sim_df$sim_id, ".txt")
  
  # Extract all outcomes for this scenario (see openmalaria_extract.R)
  metrics_df = run_extract(files       = sim_files, 
                           metrics     = o$post$metrics, 
                           age_groups  = o$sim$age_groups)
  

  
  #drop first  survey
  metrics_df=metrics_df[!(metrics_df$time==1),]
  
  sev = metrics_df[metrics_df$metric=="expectedSevere","value"]
  host = metrics_df[metrics_df$metric=="nHost","value"] #these values are ok...
  # Save full post-processed output
  save_path = paste0(o$pth$results, "postprocessed.txt") #this is ok
  write.table(metrics_df, save_path, row.names = FALSE)
  
  # Convert to wide format using only last time point
  wide_df = metrics_df %>% 
    dplyr::rename(sim_id = sim) %>%
    #filter(time == max(time)) %>% 
    tidyr::spread(metric, value) %>% 
    dplyr::select(one_of("sim_id", "age", "time",o$post$metrics)) #TR added age #also still ok
  
  wide_df$age_sev_inc = wide_df$nSevere/wide_df$nHost #max=120
  
  wide_df$age_prev = wide_df$nPatent/wide_df$nHost
 
  # Write result dataframe to file
  save_path = paste0(o$pth$results, "all_out_wide.txt")
  write.table(wide_df, save_path, row.names = FALSE)
  # ---- EIR-PfPR2-10 results ----
  
  # Subset the wide dataframe with metrics of interest
  eir_prev_df = wide_df %>% 
    dplyr::select(sim_id, inputEIR, simulatedEIR, prevalence_2to10) %>%
    mutate(inputEIR     = (365 / 5) * inputEIR, 
           simulatedEIR = (365 / 5) * simulatedEIR)
  
  # Then join results to the sim ID dataframe
  eir_prev_df = left_join(sim_df, eir_prev_df, by = "sim_id")
  
  # Write result dataframe to file
  save_path = paste0(o$pth$results, "eir_prevalence.txt")
  write.table(eir_prev_df, save_path, row.names = FALSE)
  
  
  # ---- prevalence-incidence results ----
  
  # Subset the wide dataframe with metrics of interest
  prev_inc_df = wide_df %>% 
    dplyr::select(sim_id, time,age, nTreatments1,nHost,nPatent,nUncomp, nSevere,expectedDirectDeaths,expectedSevere,expectedIndirectDeaths) %>%
    group_by(sim_id,age) %>%
    summarise_at(c("nTreatments1","nHost","nPatent","nUncomp","nSevere","expectedDirectDeaths","expectedSevere","expectedIndirectDeaths"),sum) %>%
    na.omit() %>%
    mutate(
      age_prev = nPatent / nHost,
      nHost = nHost / length(levels(as.factor(wide_df$time))) #divide by number of surveys per year 
    )

  tmp = wide_df %>% 
    dplyr::select(sim_id, time,age, prevalence_2to10) %>%
    group_by(sim_id,age) %>%
    summarise_at(c("prevalence_2to10"),mean) %>%
    dplyr::select(sim_id,prevalence_2to10) %>%
    na.omit()
  prev_inc_df =  merge(prev_inc_df,tmp,by=c("sim_id"))
  
  
  
  # Then join results to the sim ID dataframe
  prev_inc_df = left_join(sim_df, prev_inc_df, by = "sim_id")
  
  # Write result dataframe to file
  save_path = paste0(o$pth$results, "prev_incidence.txt")
  write.table(prev_inc_df, save_path, row.names = FALSE)
  

  

}

