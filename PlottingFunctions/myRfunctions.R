############################################################
# MY R FUNCTIONS
#
# A series of helpful R functions.
#
# Written by A.J.Shattock - andrewjames.shattock@swisstph.ch
############################################################

pacman::p_load(rlist,stringr,fastmatch,zoo)

# ---------------------------------------------------------
# Assert that a condition is true. Throw error if not.
# ---------------------------------------------------------
assert = function(condition, errMsg) {
  if (!condition) {
    stop(errMsg)
  }
}

# ---------------------------------------------------------
# Clear the console.
# ---------------------------------------------------------
clc = function() {
  cat("\014")
}

# ---------------------------------------------------------
# Clear all figures.
# ---------------------------------------------------------
clf = function() {
  graphics.off()
}

# ---------------------------------------------------------
# Clear environment (and also console and figures by default).
# ---------------------------------------------------------
clear = function(clearConsole=TRUE, clearFigs=TRUE) {
  if (clearConsole == TRUE) {clc()}
  if (clearFigs == TRUE) {clf()}
  rm(list=ls(envir=globalenv()), envir=globalenv())
}

# ---------------------------------------------------------
# Is an element conatined within a vector? Returns boolean.
# ---------------------------------------------------------
contains = function(vec, el) {
  bool = is.element(el, vec)  # Simply for readability
  return(bool)
}

# ---------------------------------------------------------
# Create a log file (see function wait_for_jobs)
# ---------------------------------------------------------
create_bash_log = function(pth, file) {
  log_file = file.path(pth, file)
  if (file.exists(log_file)) file.remove(log_file)
  Sys.sleep(0.1)
  file.create(log_file)
  Sys.sleep(0.1)
  return(log_file)
}

# ---------------------------------------------------------
# Create dataframe with custom column headings in one line.
# ---------------------------------------------------------
dataFrame = function(x, colNames=NULL, rowNames=NULL, ...) {
  df = data.frame(x, row.names=rowNames, ...)  # Create dataframe
  if (!is.null(colNames)) 
    colnames(df) = colNames  # Rename column headers if need be
  return(df)
}

# ---------------------------------------------------------
# Drop one or more columns from a dataframe - improved readability.
# ---------------------------------------------------------
df_drop = function(df, drop_cols) {
  drop_idx = (names(df) %in% drop_cols)
  df = df[, !drop_idx]
  return(df)
}

# ---------------------------------------------------------
# Convert dataframe to a list using a column as list fields.
# ---------------------------------------------------------
df2list = function(df, col_idx, append_list = NULL) {
  if (is.null(append_list)) append_list = list()  # Initiate list if needed
  for (i in 1 : nrow(df)) {
    vec = unlist(unname(df[i, -col_idx]))  # Convert to vector
    vec = na.trim(vec, sides = "right")  # Remove trailing NAs
    append_list[[df[i, col_idx]]] = vec  # Append vector to list
  }
  return(append_list)
}

# ---------------------------------------------------------
# Create figure window and display plot.
# ---------------------------------------------------------
figure = function(figHandle) {
  operatingSystem = Sys.info()['sysname']
  if (operatingSystem == "Linux") plot(figHandle)
  if (operatingSystem == "Windows") {
    windows()  # Other techniques available
    print(figHandle)
  }
}

# ---------------------------------------------------------
# Platform specific file separator - for readability.
# ---------------------------------------------------------
file_sep = function() {
  platform_file_sep = .Platform$file.sep
  return(platform_file_sep)
}

# ---------------------------------------------------------
# Convert a field of a list into a vector.
# ---------------------------------------------------------
list2vec = function(l, idx=NA) {
  vals = unlist(lapply(l, '[[', idx))
  return(vals)
}

# ---------------------------------------------------------
# Number of running and pending jobs on the cluster.
# ---------------------------------------------------------
n_slurm_jobs = function(user = "shatto0000") {
  
  # Base sq command for user
  sq = paste("squeue -u", user)
  
  # Concatenate full commands
  slurm_running = paste(sq, "-t running | wc -l")
  slurm_pending = paste(sq, "-t pending | wc -l")
  
  # System call to determine number of slurm processes
  n_running = system(slurm_running, intern = TRUE)
  n_pending = system(slurm_pending, intern = TRUE)
  
  # Convert to numeric and minus 1 to get number of jobs
  n_running = as.numeric(n_running) - 1
  n_pending = as.numeric(n_pending) - 1
  
  # Compile into list
  n_jobs = list(running = n_running, pending = n_pending)
  
  return(n_jobs)
}

# ---------------------------------------------------------
# Sum whilst ignoring NAs - better readability.
# ---------------------------------------------------------
nasum = function(x) {
  y = sum(x, na.rm=TRUE)
  return(y)
}

# ---------------------------------------------------------
# Suppress output from a function call.
# ---------------------------------------------------------
quiet = function(x) { 
  sinkCon = file("sink.txt")
  sink(sinkCon, type="output")
  sink(sinkCon, type="message")
  on.exit(sink(type="output"))
  on.exit(sink(type="message"), add=TRUE)
  on.exit(file.remove("sink.txt"), add=TRUE)
  invisible(force(x)) 
} 

# ---------------------------------------------------------
# Initiate progress bar with normal-use options.
# ---------------------------------------------------------
start_progress_bar = function(n_tasks) {
  pb = txtProgressBar(min = 0, max = n_tasks,
                      initial = 0, width = 100, style = 3)
  return(pb)
}

# ---------------------------------------------------------
# Sumbit an array job to the SciCORE cluster.
# ---------------------------------------------------------
submit_cluster = function(bash_file, n_jobs, ...) {
  
  # Create a new log file for the cluster jobs
  log_file = create_bash_log("~", "scicore_log.txt")
  
  # Construct slurm array command for running in parallel
  slurm_array = paste0("--array=1-", n_jobs)
  
  # Concatenate system command
  sys_command = paste("sbatch", slurm_array, 
                      bash_file, log_file, ...)
  
  # Invoke this command
  system(sys_command)
  
  browser()  # Need to rejig wait_for_jobs a little...
  
  # Wait for all cluster jobs to complete
  wait_for_jobs(o, log_file, n_jobs)
}

# ---------------------------------------------------------
# Summing function for higher dimensional arrays.
# ---------------------------------------------------------
sumMat = function(x, sumDim) {
  dims  = dim(x)
  nDims = length(dims)
  if (nDims == 2) {  # 2D matrices included for readability
    if (sumDim == 1) {
      y = rowSums(x)  # Simply apply rowSums
    } else if (sumDim == 2) {
      y = colSums(x)  # Simply apply colSums
    }
  }
  if (nDims > 2) {  # Higher dimensional arrays are summed and squeezed
    applyDim = setdiff(1 : nDims, sumDim)
    y = apply(x, applyDim, sum)  # Use apply method
  }
  return(y)
}

# ---------------------------------------------------------
# Format a number with thousand mark separators.
# ---------------------------------------------------------
thou_sep = function(val) {
  format_val = format(val, scientific = FALSE,
                      trim = TRUE, 
                      drop0trailing = TRUE, 
                      big.mark = ",")
  return(format_val)
}

# ---------------------------------------------------------
# Repeat vector into a 2D matrix.
# ---------------------------------------------------------
vec2mat = function(v, nRow=1, nCol=1) {
  m = matrix(v, nrow=nRow, ncol=nCol*length(v), byrow=TRUE)
  return(m)
}

# ---------------------------------------------------------
# Waits until all jobs have finished (marked in a log file).
#
# Adapted from code written by Monica Golumbeanu.
# ---------------------------------------------------------
wait_for_jobs = function(o, log_file, n_lines, 
                         wait_time = 1, pause_time = 5) {
  
  # Wait for log file to be created
  while (!file.exists(log_file)) Sys.sleep(wait_time)
  
  # Initiate a progress bar
  pb = txtProgressBar(min = 0, max = n_lines, 
                      width = 100, style = 3)
  
  # Wait for at least one line to be written
  while (file.info(log_file)$size == 0) Sys.sleep(wait_time)
  
  # ---- Continuously update progress bar ----
  
  # Wait for all jobs to write to log file
  while (nrow(read.table(log_file)) < n_lines) {
    
    # Update progress bar
    setTxtProgressBar(pb, nrow(read.table(log_file)))
    
    # Check number of running and pending jobs
    n_jobs = n_slurm_jobs(user = o$user)
    
    # Have all jobs already finished?
    if (sum(unlist(n_jobs)) == 0) {
      
      # Pause - let any recently completed jobs write to log
      Sys.sleep(pause_time)
      
      break # Break out of while loop
    }
    
    # Wait before testing again
    Sys.sleep(wait_time)
  }
  
  # Finalise progress bar
  setTxtProgressBar(pb, n_lines)
  close(pb)
  
  # ---- Report any batch job failures ----
  
  # Number of jobs not reported in log file
  n_missing = n_lines - nrow(read.table(log_file))
  
  # Tell user how many jobs did not successfully complete
  if (n_missing > 0)
    message("!! Batch job finished with ", n_missing, " errors !!")
}

# ---------------------------------------------------------
# Writes a line of text to file (handling file connection).
# ---------------------------------------------------------
write_lines = function(txt, file_name) {
  file_conn = file(file_name, blocking = TRUE)
  writeLines(txt, file_conn)
  close(file_conn)
}

