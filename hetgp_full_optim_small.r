u = "" #HPC cluster user name
GitDir = "" # Local version of the Git repository
ExperimentDir = "" #Parent folder to experiments

ExperimentName = "Tests"
DataFile= "fieldData.txt"  #when not doing synthetic data validation, this should be set to "fieldData.txt"

Diagnostic_plots = FALSE
collection.of.full.gps.ON = FALSE
SCICORE = TRUE
INIT = TRUE
eval.ID = "param.eval.hetGP" #needed as a job identifier in EvalOM. Different IDs allow for running multiple fitting runs at the same time
runInitial_gp = TRUE # Allows for continuing
STEP1= TRUE # This is step is required for noise estimation.
noise=FALSE #for loading an existing noise esimation object: noise=TRUE
STEP2=TRUE
CONT1=FALSE
CONT2=FALSE
Continued_run=FALSE
Continue_optimization=FALSE
T_Process = FALSE #Set to true to implemnt Student t process: Student-t processes generalize GPs, keeping most of their benefits at almost no extra cost,offering an improved robustness to outliers and larger tail noise. Several choices exist in theliterature; see, for example, the work ofWang, Shi, and Lee(2017). [Binois 2016]


drop_point =FALSE #to drop contentious point from objective 7
############################################################
#
# Load Packages
#
############################################################


if(!require(pacman)){install.packages("pacman")}; require(pacman) 
pacman::p_load(Matrix, MASS, reshape2, tmvtnorm, gridExtra, doParallel, mlr, hetGP,dplyr, scater, tictoc, rlist,parallel,lubridate, randtoolbox,tmvtnorm) 

############################################################
#
# set up folders
#
############################################################

#set local.directory
local.directory = paste0(ExperimentDir,ExperimentName,'/') #working directory
om.path = paste0(local.directory,'openmalaria') #OpenMalaria directory. 


if(!Continue_optimization){
  #set up folder structure
  system(paste0("mkdir -p ",local.directory,"{PriorDraws/{LikelihoodComp,Output,Plots/{main/{updated,best},init},scenarios/{scenario_buffers,scenario_templates,scenarios_complete},Simulation},script,out}"))
  #copy necessary files from Git
  system(paste0("cp -R ",GitDir,"openmalaria ",local.directory))
  system(paste0("cp -R ",GitDir,"scenario_templates/* ",local.directory,"PriorDraws/scenarios/scenario_templates"))
  system(paste0("cp -R ",GitDir,"LikelihoodComp/* ",local.directory,"PriorDraws/LikelihoodComp/"))
  if (DataFile != "fieldData.txt"){
  system(paste0("cp ",GitDir,"simulatedFieldData/",DataFile," ",local.directory,"PriorDraws/LikelihoodComp/fieldData.txt"))
  } 
}
############################################################
#
# functions
#
############################################################
invlogit <- function(x) {1/(1+exp(-x))} # define this here for later use

############################################################
#
# load functions
#
############################################################
## Load functions for running individual simulation models: Eval.model.large() and Eval.model.small(): must return a list with $Objectives (a W dimensional vector of sub-objectives, failed objectives being NAs) and $Runtime (a single number giving the computational runtime of the model evaluation)


D = 23 #dimensions
lockBinding("D", globalenv())
W = 11 #output objectives
lockBinding("W", globalenv())
nseeds = 2 #number of seeds per point; previously = 2; in multi-seed run this was =5

setwd(local.directory)
source(paste0(GitDir,"EvalOM.parallel.R"))
source(paste0(GitDir,"UCB_acquisition_function.r"))

if(drop_point){source(paste0(GitDir,"LikelihoodComp/likelihood.scicore.drop_data.R"))} # for drop_point run: To overwrite the likelihood calculation file that is called in EvalOM
############################################################################
#############################################################################
#
###                     INITIALISE THE PROBLEM
#
############################################################################
############################################################################
### Initialise problem
## Initial design points may be selected deterministically or randomly; a minimum of ~50 (better = 500) initial function evaluations are required to ensure stable behaviour of this algorithm in the initial stages of exploration
## Since an estimate of the iid noise term on the small and large models is required one might choose to incorporate this into the initial design as below  
## The function evaluations should be stored in a list 'evaluations' with members $small and $large
## The example code below implements a Sobol sampling design with initial draws at the prior median to estimate iid noise terms 

total.time.init=0
if(runInitial_gp) { 
  tic("Initialisation")
  # Melissa: For testing
  #test values
  #Nrep.noise.small <- 100 #number of initial points (multiplied with nseeds for number of sim setups; *61 for number of actually simulated scenarios)
  #Nrep.noise.large <- 10 #10	#20 = apporx. 480 min = 8h
  Nrep.noise.small= 30 #= 150
  Nrep.noise.large=1
  #total time to evaluate test values = approx 21h20min
  if(STEP1){
    evaluations <- list()
    evaluations$small <- list()
    evaluations$large <- list()
    
    if(CONT1){
      load("eval.state.noise.logLH")
      start.rep.small = length(evaluations$small) + 1
      start.rep.large = length(evaluations$large) + 1
    } else{
      start.rep.small = 1
      start.rep.large = 1
    }
    
    
    if(start.rep.small <= (Nrep.noise.small*nseeds)){
      
      for (i in seq(start.rep.small,Nrep.noise.small*nseeds,nseeds)) {
        evaluations$small[[i]] <- list()
        evaluations$small[[i]]$parameters <- rep(0.5,D)
        
        for (j in 1:nseeds){
          evaluations$small[[i+j-1]] <- evaluations$small[[i]] 
          evaluations$small[[i+j-1]]$parameters <- evaluations$small[[i]]$parameters          
        }
      }
      
      #set up parallelisation
      n = detectCores()
      
      cl = makeCluster(n)
      registerDoParallel(cl)
      
      # do this every time before writing .bat files.Not the most elegant fix, but circumvents current problems with global vs local variable assignment...
      size <<- "small"
      foreach (i=1:length(evaluations$small)) %dopar% {  
        print(i)
        current.index=paste0("init_",sprintf("%04i",i))
        write.xmls.small(evaluations$small[[i]]$parameters)
      } #write XMLs
      
      stopCluster(cl)
      registerDoSEQ()
      
      init1 = toc()
      total.time.init = total.time.init +(init1$toc - init1$tic)
      
      #Run simulations
      size <<- "small"
      Eval.model.small()
      
      #timekeeping
      total.time.init = total.time.init +  period_to_seconds(hms(Runtime))
      tic("Initialisation_pt2")
      
      for (i in 1:length(evaluations$small)) {  
        #Runtime = NA #Remove later
        current.index=paste0("init_",sprintf("%04i",i))
        size <<- "small"
        
        if (as.numeric(system(paste0('ls ',local.directory,'PriorDraws/Output/',size,'/output_',current.index,' | wc -l'),intern=TRUE))<61){ #if the number of output files is less than 61, some runs have crashed and the location is invalid
          cat("Error: Bad location at index ",current.index,". Moving on.\n ", sep="")
          evaluations$small[[i]]$Objectives <- rep(NA,W)
        }else{
        evaluations$small[[i]] = list.merge(evaluations$small[[i]],obj.calc.small(evaluations$small[[i]]$parameters, dir="init/"))
        }# LF calculations times 61
      } 
      
      setwd(local.directory)
      save(evaluations,file="eval.state.noise.logLH") 
      init2 = toc()
      total.time.init = total.time.init + (init2$toc-init2$tic)
    }
    
    if(start.rep.large <= Nrep.noise.large*nseeds){
      tic("Initialisation continued...")
      
      for (i in seq(start.rep.large,Nrep.noise.large*nseeds,nseeds)) {
        evaluations$large[[i]] <- list()
        evaluations$large[[i]]$parameters <- rep(0.5,D)
        
        for (j in 1:nseeds){
          evaluations$large[[i+j-1]] <- evaluations$large[[i]] 
          evaluations$large[[i+j-1]]$parameters <- evaluations$large[[i]]$parameters 
          
        }  
      }
      
      for (i in 1:length(evaluations$large)) {  
        current.index=paste0("init_",sprintf("%04i",i))
        write.xmls.large(evaluations$large[[i]]$parameters)} 
      
      init2 = toc()
      total.time.init = total.time.init + (init2$toc - init2$tic)
      
      size <<- "large"
      Eval.model.large()
      #timekeeping
      total.time.init = total.time.init + period_to_seconds(hms(Runtime))
      tic("Initialisation_pt4")
      
      
      for (i in 1:length(evaluations$large)) {  
        #Runtime = NA #Remove later
        current.index=paste0("init_",sprintf("%04i",i))
        if (as.numeric(system(paste0('ls ',local.directory,'PriorDraws/Output/',size,'/output_',current.index,' | wc -l'),intern=TRUE))<61){ #if the number of output files is less than 61, some runs have crashed and the location is invalid
          cat("Error: Bad location at index ",current.index,". Moving on.\n ", sep="")
          evaluations$large[[i]]$Objectives <- rep(NA,W)
        }else{
          evaluations$large[[i]] = list.merge(evaluations$large[[i]],obj.calc.large(evaluations$large[[i]]$parameters, dir="init/"))  # LF calculations times 61
        }
        
      } 
      
      
      setwd(local.directory)
      save(evaluations,file="eval.state.noise.logLH")
      
      init3=toc()
      total.time.init = total.time.init + (init3$toc-init3$tic)
    }
  } else if (Continued_run) {load("eval.state.newer") # TR added for additional simulations to previous runs 13.3.2019
  } else if (noise){load("eval.state.noise.logLH")} ## <<added 20190218 MP>>##
  
  
  #Note: Ewan had included a noise estimation step here. - The "nugs" object from predict.hetGP/hetTP provides automated noise estimation 
  
  if(STEP2) { ## <<added 20190218 MP>>## ; ##TR: changed condition 20190320
    tic("Initialisation step2 part 1")
    
    evaluations <- list()
    evaluations$small <- list()
    evaluations$large <- list()
    

    Ninit.small=1000 #large number required to ensure stable behaviour of this algorithm in the initial stages of exploration
    Ninit.large=1 #not required in this version of the code
    if(CONT2){
      load("eval.state.init")
      start.it.small =length(evaluations$small) + 1
      start.it.large = length(evaluations$large) + 1
    } else{
      start.it.small = 1
      start.it.large = 1
    }
    
    sobol_sample_points = randtoolbox::sobol(Ninit.small, D)
    sobol_row=1
    for (i in seq(start.it.small,Ninit.small*nseeds,nseeds)) {
      evaluations$small[[i]] <- list()
      #evaluations$small[[i]]$parameters <-  runif(D,0,1) #values must be from hypercube for quant.to.prior function
      evaluations$small[[i]]$parameters <-  sobol_sample_points[sobol_row,] # Sobol sequence sampling for better coverage of the space
      
      for (j in 1:nseeds){
        evaluations$small[[i+j-1]] <- evaluations$small[[i]] 
        evaluations$small[[i+j-1]]$parameters <- evaluations$small[[i]]$parameters          
      }   
      sobol_row = sobol_row+1
    }    
    rm(sobol_sample_points,sobol_row)
    save(evaluations, file="eval.state.init") #save in between
    
    #set up parallelisation
    n = detectCores()
    
    cl = makeCluster(n)
    registerDoParallel(cl)
    
    foreach (i=1:length(evaluations$small)) %dopar% {  
      
      current.index=paste0("init_",sprintf("%04i",i))
      size <<-"small"
      write.xmls.small(evaluations$small[[i]]$parameters)} 
    
    
    #stop cluster
    parallel::stopCluster(cl)
    registerDoSEQ()
    
    init_s2_pt1 = toc()
    total.time.init = total.time.init + (init_s2_pt1$toc - init_s2_pt1$tic)
    size <<-"small"
    Eval.model.small() # run everything
    
    #timekeeping
    total.time.init = total.time.init + period_to_seconds(hms(Runtime))
    tic("Initialisation step2 part 2")
    
    for (i in 1:length(evaluations$small)) {  
      
      current.index=paste0("init_",sprintf("%04i",i))
      size<<-"small"
      
      if (as.numeric(system(paste0('ls ',local.directory,'PriorDraws/Output/',size,'/output_',current.index,' | wc -l'),intern=TRUE))<61){ #if the number of output files is less than 61, some runs have crashed and the location is invalid
        cat("Error: Bad location at index ",current.index,". Moving on.\n ", sep="")
        evaluations$small[[i]]$Objectives <- rep(NA,W)
      }else{
        current.index=paste0("init_",sprintf("%04i",i))
        evaluations$small[[i]] = list.merge(evaluations$small[[i]],obj.calc.small(evaluations$small[[i]]$parameters, dir="init/"))  # LF calculations times 61
      }
    }
    
    
    setwd(local.directory)
    save(evaluations,file="eval.state.init")
    init_s2_pt2 = toc()
    total.time.init = total.time.init + (init_s2_pt2$toc - init_s2_pt2$tic)
  }
  
  
  tic("Initialisation step2 part 3")
  
  sobol_sample_points = randtoolbox::sobol(Ninit.large, D)
  sobol_row=1
  for (i in seq(start.it.large,Ninit.large*nseeds,nseeds)) {
    evaluations$large[[i]] <- list()
    #evaluations$large[[i]]$parameters <-  runif(D,0,1) #values must be from hypercube for quant.to.prior function
    evaluations$large[[i]]$parameters <-  sobol_sample_points[sobol_row,] # Sobol sequence sampling for better coverage of the space
    
    for (j in 1:nseeds){
      evaluations$large[[i+j-1]] <- evaluations$large[[i]] 
      evaluations$large[[i+j-1]]$parameters <- evaluations$large[[i]]$parameters          
    }      
    sobol_row = sobol_row + 1
  }    
  rm(sobol_sample_points,sobol_row)
  
  #save updated evaluations object
  save(evaluations,file = "eval.state.init")
  
  cl = makeCluster(n)
  registerDoParallel(cl)
  
  foreach (i=1:length(evaluations$large)) %dopar% {  
    #size ="large"
    current.index=paste0("init_",sprintf("%04i",i))
    size <<-"large"
    write.xmls.large(evaluations$large[[i]]$parameters)} # write times 61
  
  
  
  #stop cluster
  parallel::stopCluster(cl)
  registerDoSEQ()
  
  init_s2_pt3 = toc()
  total.time.init = total.time.init + (init_s2_pt3$toc - init_s2_pt3$tic)
  
  #Run
  size <<-"large"
  Eval.model.large() 
  
  total.time.init = total.time.init + period_to_seconds(hms(Runtime))
  tic("Finalising initialisation")
  
  for (i in 1:length(evaluations$large)) {  
    
    current.index=paste0("init_",sprintf("%04i",i))
    size <<-"large"
    if (as.numeric(system(paste0('ls ',local.directory,'/PriorDraws/Output/',size,'/output_',current.index,' | wc -l'),intern=TRUE))<61){ #if the number of output files is less than 61, some runs have crashed and the location is invalid
      cat("Error: Bad location at index ",current.index,". Moving on.\n ", sep="")
      evaluations$large[[i]]$Objectives <- rep(NA,W)
    }else{
      current.index=paste0("init_",sprintf("%04i",i))
      evaluations$large[[i]] = list.merge(evaluations$large[[i]],obj.calc.large(evaluations$large[[i]]$parameters, dir="init/"))  # LF calculations times 61
      
    }
  } 
  
  setwd(local.directory)
  save(evaluations,file="eval.state.init")
  
  final.init = toc()
  total.time.init =total.time.init + (final.init$toc - final.init$tic)
  save(total.time.init, file = "total.time.init")
}else {load("eval.state.init")} ## TR: changed file name 20190320

#} else {load("eval.state.om.newer")}  ## <<added 20190218 MP>>##


############################################################
#
#  CHECKPOINT 2
#
############################################################
tic("Preparing main loop")
# Initialise proposal distributions for test particles
  proposal.means <- rep(0.5,D) #rep(0,D)
  proposal.covs <- diag(1,D)



if("eval.state.noise.logLH" %in% system("ls", intern=TRUE)){
  load("eval.state.noise.logLH")
  Nrep.noise.small = length(evaluations$small)
  Nrep.noise.large = length(evaluations$large)
}
small.repetitions <- matrix(-999999,nrow=Nrep.noise.small,ncol=W)
for (w in 1:W) {
  for (i in 1:Nrep.noise.small) {small.repetitions[i,w] <- evaluations$small[[i]]$Objectives[w]}
}
small.sum.of.squares <- colVars(small.repetitions,na.rm=T)
small.n.obs <- colSums(!is.na(small.repetitions))
large.repetitions <- matrix(-999999,nrow=Nrep.noise.large,ncol=W)
for (w in 1:W) {
  for (i in 1:Nrep.noise.large) {large.repetitions[i,w] <- evaluations$large[[i]]$Objectives[w]}
}
large.sum.of.squares <-colVars(large.repetitions,na.rm=T)
large.n.obs <- colSums(!is.na(large.repetitions))



proposal.means <- rep(0.5,D) 
proposal.covs <- diag(1,D)

prep = toc()
total.time.init = total.time.init + (prep$toc-prep$tic)


############################################################################
############################################################################
#
###                     DATA PREPARATION
#
############################################################################
############################################################################
#	
############################################################################
### (a) Prepare data: Generate or load training and test set
############################################################################

if(!Continue_optimization){
  time=list() #time excludes initalization phase
  
  tic("Preparing data")
  
  #Load evaluations object
  load("eval.state.init")
  master.data.small <- matrix(0,nrow=length(evaluations$small),ncol=(D*2+W))
  for (i in 1:length(evaluations$small)) {
    master.data.small[i,1:D] <- evaluations$small[[i]]$parameters-0.5
    master.data.small[i,(D+1):(2*D)] <- (evaluations$small[[i]]$parameters-0.5)^2
    master.data.small[i,(2*D+1):(2*D+W)] <- evaluations$small[[i]]$Objectives
  }
  colnames.small <- list()
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("PARAM",sprintf("%02i",i))}
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("SQ_PARAM",sprintf("%02i",i))}
  for (i in 1:W) {colnames.small[[length(colnames.small)+1]] <- paste0("OBJECTIVE",sprintf("%02i",i))}
  colnames.small <- unlist(colnames.small)
  colnames(master.data.small) <- colnames.small
  master.data.small <- as.data.frame(master.data.small)
  
  master.data.large <- matrix(0,nrow=length(evaluations$large),ncol=(D*2+W))
  for (i in 1:length(evaluations$large)) {
    master.data.large[i,1:D] <- evaluations$large[[i]]$parameters-0.5
    master.data.large[i,(D+1):(2*D)] <- (evaluations$large[[i]]$parameters-0.5)^2
    master.data.large[i,(2*D+1):(2*D+W)] <- evaluations$large[[i]]$Objectives
  }
  colnames.large <- list()
  for (i in 1:D) {colnames.large[[length(colnames.large)+1]] <- paste0("PARAM",sprintf("%02i",i))}
  for (i in 1:D) {colnames.large[[length(colnames.large)+1]] <- paste0("SQ_PARAM",sprintf("%02i",i))}
  for (i in 1:W) {colnames.large[[length(colnames.large)+1]] <- paste0("OBJECTIVE",sprintf("%02i",i))}
  colnames.large <- unlist(colnames.large)
  colnames(master.data.large) <- colnames.large
  master.data.large <- as.data.frame(master.data.large)
  
  #Log transform objectives for better emulator performance (minimzation still possible)
  master.data.small[,c((2*D+1):(2*D+W))] = log(master.data.small[,c((2*D+1):(2*D+W))]) 
  master.data.large[,c((2*D+1):(2*D+W))] = log(master.data.large[,c((2*D+1):(2*D+W))]) 
  
  
  #divide into test and training sets
  h1 = sample(seq(1,dim(master.data.small)[1],nseeds),dim(master.data.small)[1]/nseeds*0.1) %>%sort() ##sample indeces of unique design points (10% of total unique design points) 
  h2=h1
  if(nseeds>1){
    for(s in 2:nseeds-1){h2 = sort(c(h2, h1+s))} #add replicate indeces
  }
  
  master.data.small.test <- master.data.small[h2,] 
  master.data.small.train <- master.data.small[-h2,] #complement to the holdout
  
  
  
  #Repeat procedure with large simulation data set
  h1 = sample(seq(1,dim(master.data.large)[1],nseeds),dim(master.data.large)[1]/nseeds*0.1) %>%sort() ##sample numbers divisible by nseeds 1 and dim(master.data.small)[1] 800 rows = 10%
  h2=h1
  if(nseeds>1){
    for(s in 2:nseeds-1){h2 = sort(c(h2, h1+s))}
  }
  master.data.large.test <- master.data.large[h2,] 
  master.data.large.train <- master.data.large[-h2,] #complement to the holdout
  
  
  
  #aggregated test set by seeds for truth vs predicted plots
  objective.colnames = colnames.small[c((2*D+1):(2*D+W))]
  param.colnames = colnames.small[c(1:(2*D))]
  #divide data into predictors and response variable and convert to matrix. Necessary for hetGP functions
  X = master.data.small.test %>% dplyr::select(param.colnames) %>% as.matrix() 
  prdata.Z0=list()
  for(w in 1:W){
    obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
    Z = master.data.small.test %>% dplyr::select(obj.w) %>% as.matrix()
    prdata =find_reps(X=X, Z=Z, rescale = FALSE, normalize = FALSE)
    prdata.Z0[[w]] = prdata$Z0      
  }
  master.data.small.test.aggregated = cbind(prdata$X0, prdata.Z0[[1]])
  for(w in 2:W){master.data.small.test.aggregated = cbind(master.data.small.test.aggregated, prdata.Z0[[w]])}
  colnames(master.data.small.test.aggregated)[c((2*D+1):(2*D+W))] = objective.colnames
  
  
  #bind results from eval.state.noise.logLH to training frame for initial noise estimation
  
  #Load evaluations object
  load("eval.state.noise.logLH")
  noise.data.small <- matrix(0,nrow=length(evaluations$small),ncol=(D*2+W))
  for (i in 1:length(evaluations$small)) {
    noise.data.small[i,1:D] <- evaluations$small[[i]]$parameters-0.5
    noise.data.small[i,(D+1):(2*D)] <- (evaluations$small[[i]]$parameters-0.5)^2
    noise.data.small[i,(2*D+1):(2*D+W)] <- evaluations$small[[i]]$Objectives
  }
  
  colnames(noise.data.small) <- colnames.small
  noise.data.small <- as.data.frame(noise.data.small)
  
  #Log transform objectives in noise data
  noise.data.small[,c((2*D+1):(2*D+W))] = log(noise.data.small[,c((2*D+1):(2*D+W))]) 
  
  master.data.small.train = rbind(noise.data.small, master.data.small.train)
  
  
  noise.data.large <- matrix(0,nrow=length(evaluations$large),ncol=(D*2+W))
  for (i in 1:length(evaluations$large)) {
    noise.data.large[i,1:D] <- evaluations$large[[i]]$parameters-0.5
    noise.data.large[i,(D+1):(2*D)] <- (evaluations$large[[i]]$parameters-0.5)^2
    noise.data.large[i,(2*D+1):(2*D+W)] <- evaluations$large[[i]]$Objectives
  }
  colnames(noise.data.large) <- colnames.large
  noise.data.large <- as.data.frame(noise.data.large)
  noise.data.large[,c((2*D+1):(2*D+W))] = log(noise.data.large[,c((2*D+1):(2*D+W))]) 
  
  master.data.large.train = rbind(noise.data.large, master.data.large.train)
  
  
  #save test and training frames
  write.table(master.data.small.train,"master.data.small.train")
  write.table(master.data.small.test,"master.data.small.test")
  write.table(master.data.large.train,"master.data.large.train")
  write.table(master.data.large.test,"master.data.large.test")
  write.table(master.data.small.test.aggregated,"master.data.small.test.aggregated")
  
  n.init.sims.small = dim(unique(master.data.small.train[,c(1:2*D)]))[1] #save the number of unique design points for later
  
  load("eval.state.init")
  start.iteration=1
}else{
  colnames.small <- list()
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("PARAM",sprintf("%02i",i))}
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("SQ_PARAM",sprintf("%02i",i))}
  for (i in 1:W) {colnames.small[[length(colnames.small)+1]] <- paste0("OBJECTIVE",sprintf("%02i",i))}
  colnames.small <- unlist(colnames.small)
  colnames.large <- list()
  for (i in 1:D) {colnames.large[[length(colnames.large)+1]] <- paste0("PARAM",sprintf("%02i",i))}
  for (i in 1:D) {colnames.large[[length(colnames.large)+1]] <- paste0("SQ_PARAM",sprintf("%02i",i))}
  for (i in 1:W) {colnames.large[[length(colnames.large)+1]] <- paste0("OBJECTIVE",sprintf("%02i",i))}
  colnames.large <- unlist(colnames.large)
  
  
  master.data.small.test = read.table("master.data.small.test")
  master.data.small.test.aggregated = read.table("master.data.small.test.aggregated")
  master.data.small.train = read.table("master.data.small.train.updated")
  
  master.data.large.train = read.table("master.data.large.train")
  master.data.large.test = read.table("master.data.large.test")
  
  all_selected_points=read.table("all_selected_points")
  names(all_selected_points) = c(colnames.small[1:D], "lcb", "mpe" )
  load("eval.state.newer")
  load("time.hetgp")
  start.iteration=length(time)+1
  N.main.iterations = 200 
  n.init.sims.small = dim(unique(master.data.small.train[,c(1:2*D)]))[1]
}
#############################################################################
#
### (1) Fit hetGP model: 
#
############################################################################
### (1) Fit hetGP model: 
#			(a) Prepare data: Generate or load training and test sets
#     (b) Fit hetGP to training data: Hyperparameter learnin; Matern 5_2 Kernel



############################################################################
### (b) Fit hetGP to training data: Hyperparameter learnin; Matern 5_2 Kernel
############################################################################
start.time <- proc.time() #initialise system timekeeping

if(!Continue_optimization){
  total.time = 0 #initialise "true" CPU runtime timekeeping
  
  N.main.iterations = 200 
  
  all_selected_points <- data.frame(Date=as.Date(character()),
                   File=character(), 
                   User=character(), 
                   stringsAsFactors=FALSE)
  evaluations$best =list()
}else if(Continue_optimization){
  
  total.time=time[[start.iteration-1]]$time["elapsed"]
}

for (main.iterations in start.iteration:N.main.iterations) {
  #main.iterations=1 #for testing
  tic("Main program: main")
  cat("Main iterations =",main.iterations,",: Fitting hetGP model ...")
  
  
  n = detectCores()
  cl = makeCluster(min(n,W))
  registerDoParallel(cl)
  
  
  hetGP.model = foreach(w =1:W, .export=ls(globalenv()), .packages=c("hetGP", "dplyr")) %dopar% {
    output=list()
    obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
    setwd(local.directory)
    
    ##################################
    # hetGP
    
    #prepare data to train on
    hetgp.parnames = colnames.small[1:D] # linear predictors only
    col.names = c(obj.w, hetgp.parnames)
    obj.w.data = master.data.small.train %>%
      dplyr::select(col.names) %>% na.omit()
    
    nvar = length(col.names)-1
    lower <- rep(0.001, nvar) #These are the bounds for theta, a hyperparameter
    upper <- rep(10, nvar)
    
    #divide data into predictors and response variable and convert to matrix. Necessary for hetGP.
    X = obj.w.data %>% dplyr::select(hetgp.parnames) %>% as.matrix() 
    Z = obj.w.data %>% dplyr::select(obj.w) %>% as.matrix()
    
    #account for repetitions
    prdata =find_reps(X=X, Z=Z, rescale = FALSE, normalize = FALSE)
    
    #full model fitting, incl. hyperparameter tuning
    if(T_Process){
    #Student-t Process
        output$hetgp.model$hetgp.model = mleHetTP(X=list(X0 = prdata$X0, Z0=prdata$Z0, mult = prdata$mult), Z = prdata$Z, lower = rep(0.01, nvar), upper = rep(10, nvar), covtype = "Matern5_2")
    }else{
    #Gaussian Process
        output$hetgp.model$hetgp.model = mleHetGP(X=list(X0 = prdata$X0, Z0=prdata$Z0, mult = prdata$mult), Z = prdata$Z, lower = rep(0.01, nvar), upper = rep(10, nvar), covtype = "Matern5_2")
    }
    #hold out predictions
    X = as.matrix(master.data.small.test.aggregated[,1:D])
    output$hetgp.model$hetgp.test.pred = predict(output$hetgp.model$hetgp.model, x=X)$mean		  
    
    output      
  }
  
  #Stop parallelisation. Not sure we need to put registerDoSEQ twice, but it works. And it's required at least once according to stackoverflow on troubleshooting with doParallel
  registerDoSEQ() 
  stopCluster(cl) #end cluster
  registerDoSEQ() #register seqiential coding
  
  
  cat("Done!\n")
  
  
  #Plot hold out predictions:
  pdf(file=paste0("hetGP.holdout.",main.iterations,".blue.pdf"))
  p=list()
  for(w in 1:W){
    dat = cbind(master.data.small.test.aggregated[,46+w],hetGP.model[[w]]$hetgp.model$hetgp.test.pred) 
    colnames(dat) = c("truth","response")
    dat = dat %>% as.data.frame() %>% na.omit()
    r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
    
    p[[w]]= ggplot(dat, aes(truth, response)) + 
      geom_point(size = 1, alpha=0.35, colour = "deepskyblue4")+
      geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
      labs(title = paste0("Objective ",w, ", R2=",round(r2,3)),x = "Truth", y="hetGP Response", size= 8) +
      theme_minimal() +
      coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response), max(dat$truth, dat$response)), ylim = c(min(dat$truth, dat$response), max(dat$truth, dat$response))) + 
      theme(axis.line = element_line(colour = "black"),
            plot.title=element_text(size=7.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
  }
  
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
  
  weighted.truths = numeric(dim(master.data.small.test.aggregated)[1])
  for (w in 1:W) {weighted.truths <- weighted.truths + (weights.w[w]*exp(master.data.small.test.aggregated[,(2*D+w)]))} # for each parameter vector, compute
  
  weighted.hetGP = numeric(dim(master.data.small.test.aggregated)[1]) #nrows in test frame
  for (w in 1:W) {weighted.hetGP <- weighted.hetGP + (weights.w[w]*exp((hetGP.model[[w]]$hetgp.model$hetgp.test.pred)))} 
  
  dat = cbind(weighted.hetGP,weighted.truths) 
  colnames(dat) = c("truth","response")
  dat = dat %>% as.data.frame() %>% na.omit()
  
  r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
  p[[W+1]]= ggplot(dat, aes(truth, response)) + 
    geom_point(size = 1, alpha=0.35, colour = "deepskyblue4")+
    geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
    labs(title = paste0("w. sum, R2=",round(r2,4)),x = "Truth", y="hetGP Response", size= 8) +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
    theme(axis.line = element_line(colour = "black"),
          plot.title=element_text(size=7.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  multiplot(plotlist = p, cols = 4)
  dev.off()
  
  
  ############################################################################
  #
  ### (1) Update model N.Updates times: 
  #
  ############################################################################
  ### (1) Update model N.Updates times: 
  #			(a) Make Fine grid predictions
  #			(b) Choose min (UCB - SD) as next simulation point. These are loss terms! Minimize!
  #			(c) Simulate
  #			(d) Update model
  #     (e) Plot progress
  
  if(main.iterations==1){load("eval.state.init")}else{load("eval.state.newer")}
  Ninit.small = n.init.sims.small + dim(all_selected_points)[1]
  
  N.Updates = 1 #updating doesn't seem to improve estimation much
  n.select = 100
  n.seeds = 1 #how many seeds should be included for each new point? Since the nugget effect is already estimated, we're leaving this at one for now, but we could include more seeds here, in multi-seed run this was 5
  
  #data frame to store new simulation outputs
  update.data.small <- matrix(NA,N.Updates*n.select*n.seeds,ncol=(2*D+W))
  colnames.small <- list()
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("PARAM",sprintf("%02i",i))}
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("SQ_PARAM",sprintf("%02i",i))}
  for (i in 1:W) {colnames.small[[length(colnames.small)+1]] <- paste0("OBJECTIVE",sprintf("%02i",i))}
  colnames.small <- unlist(colnames.small)
  colnames(update.data.small) <- colnames.small
  update.data.small <- as.data.frame(update.data.small)
  
  cat("Main iterations =",main.iterations,",: Start updating hetGP ... \n")
  
  #Updating is a way of adding data to an existing hetGP model without re-estimating theta
  
  main = toc()
  main.elapsed = main$toc-main$tic
  
  total.time = total.time +main.elapsed #update total.time object
  
  
  #Updating is a way of adding data to an existing hetGP model without re-estimating theta
  for (update in 1:N.Updates){
    
    tic("Updating (pre simulation)") #start timekeeping up to simulation
    ############################################################################
    ### (a) Fine grid predictions
    ############################################################################ 
    
    #set up grid 
    current.quantile <- 0.5
    
    N.test <- 200000 # Number of random points at which to estimate UCB: 200,000
    
    #test.points <- matrix(runif(N.test*D),nrow=N.test,ncol=D)-0.5 #Draw points from uniform dist. Don't Sobol sample here because Sobol is deterministic
    
    #Alternative: mvnorm sampling:
    if(main.iterations==1){
    proposal.means <- rep(0,D)
    proposal.covs <- diag(1,D)
    }else{
    tmp = exp(master.data.small.train[,(2*D+1)])
    for(w in 2:W){tmp = tmp+ exp(master.data.small.train[,(2*D+w)])}
  
    proposal.means =as.numeric(as.vector(master.data.small.train[which.min(tmp),c(1:D)]))
    proposal.covs <- cov(data.matrix(master.data.small.train[c(Ninit.small:dim(master.data.small.train)[1]),c(1:D)]),use="complete.obs")
    rm(tmp)}
    
    test.points <- rtmvnorm(N.test,proposal.means,proposal.covs, lower = rep(-0.5,D), upper=rep(0.5,D),algorithm="gibbs") 
   
    
    test.data <- cbind(test.points,test.points^2,matrix(0,nrow=N.test,ncol=W))
    colnames.test <- list()
    for (d in 1:D) {colnames.test[[length(colnames.test)+1]] <- paste0("PARAM",sprintf("%02i",d))}
    for (d in 1:D) {colnames.test[[length(colnames.test)+1]] <- paste0("SQ_PARAM",sprintf("%02i",d))}
    for (w in 1:W) {colnames.test[[length(colnames.test)+1]] <- paste0("OBJECTIVE",sprintf("%02i",w))}
    colnames.test <- unlist(colnames.test)
    colnames(test.data) <- colnames.test
    test.data = as.data.frame(test.data)
    
    cat("Main iterations =",main.iterations,", update = ",update,": Computing posterior predictives ...")
    
    
    # hetGP predict
    registerDoSEQ()  #run this sequentially because it keeps crashing otherwise...
    
    hetGP.posterior.predictives = foreach(w =1:W, .export=ls(globalenv()), .packages=c("hetGP")) %dopar% {
      obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
      setwd(local.directory)		
      X = as.matrix(test.data[,1:D])
      posterior.predictives = predict(hetGP.model[[w]]$hetgp.model$hetgp.model, x=X)
      posterior.predictives    		
    }
    
    cat("Done!\n")
    
    ############################################################################
    ### (b) choose max UCB + SD (UCB acquisition function like in Bayesian Optimization)
    ############################################################################ 
    
    lcb = numeric(N.test)
    mpe = numeric(N.test) #minimum point estimate
    
    k = kappa_calc(t=Ninit.small,d=D,delta=0.1) #for computing LCB (see Brochu 20210 and Srinivas 20XX)
    #Weighted sum of upper confidence bounds of objective estimates
    weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
    #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
    
    for (w in 1:W) {lcb <- lcb + (weights.w[w]*exp((hetGP.posterior.predictives[[w]]$mean-k*(hetGP.posterior.predictives[[w]]$sd2+hetGP.posterior.predictives[[w]]$nugs))))} # lower confidence bound
    for (w in 1:W) {mpe <- mpe + (weights.w[w]*exp(hetGP.posterior.predictives[[w]]$mean))} # lower confidence bound
    

    prediction_data = cbind(test.points, as.data.frame(lcb), as.data.frame(mpe))
    names(prediction_data) = c(colnames.small[1:D], "lcb", "mpe" )
    prediction_data = prediction_data[order(prediction_data$lcb, decreasing = FALSE),]
    
    #every 10 iterations force pure exploitation
    if ((main.iterations %% 10) == 0) {
      prediction_data = prediction_data[order(prediction_data$mpe, decreasing = FALSE),]
    }
    
    selected_points = prediction_data[c(1:n.select),]
    selected_points = selected_points[rep(seq_len(nrow(selected_points)), each = n.seeds), ]

    all_selected_points = rbind(all_selected_points,selected_points)
    write.table(all_selected_points, file="all_selected_points")
    
    for(addition in 1:dim(selected_points)[1]){
      
      evaluations$small[[length(evaluations$small)+1]] <- list()
      evaluations$small[[length(evaluations$small)]]$parameters <- as.numeric(paste0(selected_points[addition,c(1:D)]+0.5))
    }
    
    
    cat("Main iterations =",main.iterations,", update = ",update,": Simulating ...")
    
    ############################################################################
    ### (c) simulate
    ############################################################################ 
    
    n = detectCores()
    
    cl = makeCluster(n)
    registerDoParallel(cl)
    
    size <<-"small"    # do this every time before writing .bat files.Not the most elegant fix, but circumvents current problems with global vs local variable assignment...
    
    foreach (a=1:dim(selected_points)[1]) %dopar% {  
      current.index<-paste0("main_",sprintf("%04i",main.iterations),"_",sprintf("%04i",a))
      setwd(local.directory)
      system(paste0('rm -f ./PriorDraws/Output/',size,'/output_',current.index,'/*'),intern=TRUE)
      write.xmls.small(evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$parameters)
    } #write XMLs
    
    stopCluster(cl)
    registerDoSEQ()
    
    
    
       
    time.pre.sim = toc()
    total.time = total.time + (time.pre.sim$toc-time.pre.sim$tic) #add elapsed time up to simulation to total.time
    
    
    #Run simulations
    size <<-"small"
    Eval.model.small()
    
    total.time = total.time + period_to_seconds(hms(Runtime))
    
    tic("Updating (post simulation)")
    observations  = list()
    setwd(local.directory)
    
    #Calculate LF for new points
    errors = 0
    for (a in 1:dim(selected_points)[1]) {  
      current.index=paste0("main_",sprintf("%04i",main.iterations),"_",sprintf("%04i",a))
      
      #if the number of output files is less than 61, some runs have crashed and the location is invalid
      if (system(paste0('ls ./PriorDraws/Output/',size,'/output_',current.index,'/ | wc -l'),intern=TRUE)!=61){ 
        cat("Error: Bad location... Try elsewhere.\n ")
        errors = errors+1
      }else{
        size <<- "small"
        evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]  = list.merge(evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]],obj.calc.small(evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$parameters, dir="main/updated/"))  # LF calculations
        
        save(evaluations, file="eval.state.newer") # note that Objectives in the evaluations object are not log transformed
        
        update.data.small[((update-1)*dim(selected_points)[1])+a,1:D] <- evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$parameters-0.5
        update.data.small[((update-1)*dim(selected_points)[1])+a,(D+1):(2*D)] <- (evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$parameters-0.5)^2
        update.data.small[((update-1)*dim(selected_points)[1])+a,(2*D+1):(2*D+W)] <- evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$Objectives %>% log() #include log transform for easier merging
        
        
        #create updated training object
        master.data.small.train <- rbind(master.data.small.train,update.data.small[((update-1)*dim(selected_points)[1])+a,])
        
      }
    }
    
    #As long as the new simulations didn't all fail
    if(errors != dim(selected_points)[1]){
      
      cat("Done!\n")
      cat("Main iterations =",main.iterations,", update = ",update,": Updating hetGP model ...")
      
      ############################################################################
      ### (d) update hetGP model with new data
      ############################################################################  
      
      n = detectCores()
      cl = makeCluster(min(n,W))
      registerDoParallel(cl)
      
      #update model without refitting theta
      hetGP.model = foreach(w =1:W, .export=ls(globalenv()), .packages=c("hetGP", "dplyr")) %dopar% {
        
        output=list()
        obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
        setwd(local.directory)
        
        #prepeare data
        hetgp.parnames = colnames.small[1:D]
        col.names = c(obj.w, hetgp.parnames)
        obj.w.data = update.data.small[c(((update-1)*(dim(selected_points)[1])+1):(update*dim(selected_points)[1])),] %>% #only add latest row each iteration because hetGP.model gets updated and should include previous rows
          dplyr::select(col.names) %>% na.omit()
        
        nvar = length(col.names)-1
        lower <- rep(0.001, nvar) #These are the bounds for theta, a hyperparameter
        upper <- rep(10, nvar)
        
        #predictors and response variables into matrix format for input into hetGP
        Xnew = obj.w.data %>% dplyr::select(hetgp.parnames) %>% as.matrix() 
        Znew = obj.w.data %>% dplyr::select(obj.w) %>% as.matrix()
        
        newdata =find_reps(X=Xnew, Z=Znew, rescale = FALSE, normalize = FALSE)
        
        #model updating
        output$hetgp.model$hetgp.model = update(object=hetGP.model[[w]]$hetgp.model$hetgp.model, X=newdata$X0, Z = newdata$Z0, maxit =0, maxiter = 0, ginit = 0.01, method = "quick")      
        
        #output$hetgp.model$hetgp.model = update(object=hetGP.model[[w]]$hetgp.model$hetgp.model, Xnew, Znew, maxit =0, maxiter = 0, ginit = 0.01, method = "quick")      
        
        #if(i==N.Updates){
        X = as.matrix(master.data.small.test.aggregated[,1:D])
        #output$hetgp.model$hetgp.test.pred = predict(hetGP.model[[1]]$hetgp.model$hetgp.model, x=X)$mean	  #updated prediction on the holdout set
        
        output$hetgp.model$hetgp.test.pred = predict(output$hetgp.model$hetgp.model, x=X)$mean	  #updated prediction on the holdout set
        #}
        
        output
      }	
      #close cluster
      registerDoSEQ() 
      stopCluster(cl)
      registerDoSEQ() 
      
      cat("Done!\n")
      
      #store all realizations beyond the initialization set in this data frame 
      write.table(update.data.small, file=paste0("update.data.small",main.iterations))  
    }else{cat("All new points in bad locations... Moving on. \n")}
    
    
    
    cat("Main iterations =",main.iterations,", update = ",update,": Plotting ...")
    #whenever the model is updated including re-estimating hyper-parameters, predict on the hold out:
    
    weighted.truths = numeric(dim(master.data.small.test.aggregated)[1])
    for (w in 1:W) {weighted.truths <- weighted.truths + (weights.w[w]*exp(master.data.small.test.aggregated[,2*D+w]))} # for each parameter vector, compute
    
    pdf(file=paste0("hetGP.holdout.",main.iterations,".updated.pdf"))
    p=list()
    r2 = list()
    for(w in 1:W){
      dat = cbind(master.data.small.test.aggregated[,(2*D+w)],hetGP.model[[w]]$hetgp.model$hetgp.test.pred) 
      colnames(dat) = c("truth","response")
      dat = dat %>% as.data.frame() %>% na.omit()
      r2[[w]] =cor(dat$response, dat$truth,use = "complete.obs")^2
      p[[w]]= ggplot(dat, aes(truth, response)) + 
        geom_point(size = 1, alpha=0.35, colour = "deepskyblue4")+
        geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
        labs(title = paste0("Objective ",w, ", R2=",round(r2[[w]],3)),x = "Truth", y="hetGP Response", size= 8) +
        theme_minimal() +
        coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response), max(dat$truth, dat$response)), ylim = c(min(dat$truth, dat$response), max(dat$truth, dat$response))) + 
        theme(axis.line = element_line(colour = "black"),
              plot.title=element_text(size=7.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) 
    }
    
    #weighted
    weighted.hetgp = numeric(dim(master.data.small.test.aggregated)[1]) #nrows in test frame
    weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
    #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
    
    for (w in 1:W) {weighted.hetgp <- weighted.hetgp + (weights.w[w]*exp((hetGP.model[[w]]$hetgp.model$hetgp.test.pred)))} 
    
    dat = cbind(weighted.hetgp,weighted.truths) 
    colnames(dat) = c("truth","response")
    dat = dat %>% as.data.frame() %>% na.omit()
    
    r2[[W+1]] =cor(dat$truth, dat$response,use = "complete.obs")^2
    p[[W+1]]= ggplot(dat, aes(truth, response)) + 
      geom_point(size = 1, alpha=0.35, colour = "deepskyblue4")+
      geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
      labs(title = paste0("w. sum, R2=",round(r2[[W+1]],3)),x = "Truth", y="hetGP Response", size= 8) +
      theme_minimal() +
      coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
      theme(axis.line = element_line(colour = "black"),
            plot.title=element_text(size=7.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
    
    
    
    
    multiplot(plotlist = p, cols = 4)
    dev.off()
    
    ############################################################################
    ### (e) plot current performance
    ############################################################################     
    
    
    time.post.sim = toc()
    total.time = total.time + (time.post.sim$toc-time.post.sim$tic)
    tic("Housekeeping")
    #time keeping
    time[[(main.iterations-1)*N.Updates+update]] = list()  
    time[[(main.iterations-1)*N.Updates+update]]$time= proc.time() - start.time 
    time[[(main.iterations-1)*N.Updates+update]]$cpu.time = total.time
    time[[(main.iterations-1)*N.Updates+update]]$r2 = r2
    time[[(main.iterations-1)*N.Updates+update]]$obj$min_lcb = min(lcb)
    time[[(main.iterations-1)*N.Updates+update]]$params_lcb = test.points[which.min(lcb),]
    
    time[[(main.iterations-1)*N.Updates+update]]$obj$mpe = min(mpe)
    time[[(main.iterations-1)*N.Updates+update]]$params_mpe = test.points[which.min(mpe),]    
    
    
    save(time, file="time.hetgp")
    time.frame <- matrix(0,nrow=length(time),(W+4))
    for (g in 1:length(time)) {
      time.frame[g,1] <- as.numeric(time[[g]]$time["elapsed"])
      time.frame[g,c(2:(2+W))] <- unlist(time[[g]]$r2)
      time.frame[g,(W+3)] <- time[[g]]$obj$min_lcb
      time.frame[g,(W+4)] <- time[[g]]$obj$mpe
    }
    
    
    colnames <- list()
    colnames[[length(colnames)+1]] <- paste0("elapsed")
    for (w in 1:W) {colnames[[length(colnames)+1]] <- paste0("r2_",sprintf("%02i",w))}
    colnames[[length(colnames)+1]] <- paste0("r2_w")
    colnames[[length(colnames)+1]] <- paste0("min_lcb")
    colnames[[length(colnames)+1]] <- paste0("mpe")
    
    colnames <- unlist(colnames)
    colnames(time.frame) <- colnames
    time.frame <- as.data.frame(time.frame)
    
    d <- melt(time.frame[,1:(W+2)], id.vars="elapsed")
    
    #Plot current learner performance
    pdf("current.learner.performance.hetgp.pdf")
    # Everything on the same plot
    print(ggplot(d, aes(elapsed,value, col=variable)) + 
            geom_line()) 
    dev.off()
    
    #Convergence plot LCB
    pdf("min_lcb.progress.hetgp.pdf")
    
    plot(time.frame$min_lcb~seq(1,length(time.frame$min_lcb)), type="l", col="red",xlab="iteration", ylab="current minimum", lty=2,ylim=c(min(time.frame$min_lcb,time.frame$mpe),max(time.frame$min_lcb,time.frame$mpe)))
    lines(x=seq(1,length(time.frame$min_lcb)), y=time.frame$mpe, col="black")
    
    plot(time.frame$min_lcb~time.frame$elapsed, type="l", col="red",xlab="time", ylab="current minimum", lty=2,ylim=c(min(time.frame$min_lcb,time.frame$mpe),max(time.frame$min_lcb,time.frame$mpe)))
    lines(x=time.frame$elapsed, y=time.frame$mpe, col="black")
    dev.off()
    
    #Plot Euclidian distance between current and previous best
    if(length(time)>1){
      pdf("distance.pdf")
      
      edist_mpe = vector(length=length(time))
      edist_mpe[1]=NA
      for(e in 2:length(time)){edist_mpe[e] = dist(rbind(time[[e]]$params_mpe, time[[e-1]]$params_mpe))}#calculates the Euclidian distance between current and previous mpe
      plot(edist_mpe~seq(1,length(time),1), type="l", ylab="Distance", xlab="Iteration", main="Euclidian distance between current and previous MPE", col="steelblue",lwd=2)
      
      
      edist_lcb = vector(length=length(time))
      edist_lcb[1]=NA
      for(e in 2:length(time)){edist_lcb[e] = dist(rbind(time[[e]]$params_lcb, time[[e-1]]$params_lcb))}#calculates the Euclidian distance between current and previous lcb
      plot(edist_lcb~seq(1,length(time),1), type="l", ylab="Distance", xlab="Iteration", main="Euclidian distance between current and previous LCB",col="tan1",lwd=2)
      dev.off()
      
    }
    
    #end time keeping
    housekeeping = toc()
    total.time = total.time + (housekeeping$toc-housekeeping$tic)
    
    
    
  }
  
  tic("Finalising main iterations after updates (post current best simulation)")
  
  #update training data for next main
  master.data.small.train <- rbind(master.data.small.train,update.data.small)
  write.table(master.data.small.train,"master.data.small.train.updated")
  
  #after N updates, compute allPlots with best parameter set.
  cat("Done!\n")
  cat("Main iterations =",main.iterations,": Evaluating current best ...")
  
  current.quantile <- 0.5
  
  N.test <- 500000 # Number of random points at which to estimate UCB: 500,000
  if ((main.iterations %% 15) == 0) {N.test = 5000000} #every 10 main iterations, do larger evaluation (5Million) of the space
  #mvnorm sampling???
  
   #test.points <- matrix(runif(N.test*D),nrow=N.test,ncol=D)-0.5 #Draw points from uniform dist. Don't Sobol sample here because Sobol is deterministic
    
    #Alternative: mvnorm sampling:
    if(main.iterations==1){
    proposal.means <- rep(0,D)
    proposal.covs <- diag(1,D)
    }else{
    tmp = exp(master.data.small.train[,(2*D+1)])
    for(w in 2:W){tmp = tmp+ exp(master.data.small.train[,(2*D+w)])}
  
    proposal.means =as.numeric(as.vector(master.data.small.train[which.min(tmp),c(1:D)]))
    proposal.covs <- cov(data.matrix(master.data.small.train[c(Ninit.small:dim(master.data.small.train)[1]),c(1:D)]),use="complete.obs")
    rm(tmp)}
    
    test.points <- rtmvnorm(N.test,proposal.means,proposal.covs, lower = rep(-0.5,D), upper=rep(0.5,D),algorithm="gibbs") 

  test.data <- cbind(test.points,test.points^2,matrix(0,nrow=N.test,ncol=W))
  colnames.test <- list()
  for (d in 1:D) {colnames.test[[length(colnames.test)+1]] <- paste0("PARAM",sprintf("%02i",d))}
  for (d in 1:D) {colnames.test[[length(colnames.test)+1]] <- paste0("SQ_PARAM",sprintf("%02i",d))}
  for (w in 1:W) {colnames.test[[length(colnames.test)+1]] <- paste0("OBJECTIVE",sprintf("%02i",w))}
  colnames.test <- unlist(colnames.test)
  colnames(test.data) <- colnames.test
  test.data = as.data.frame(test.data)
  
  cat("Done!\n")
  cat("Main iterations =",main.iterations,": Computing posterior predictives ...")
  
  
  # hetGP predict
  registerDoSEQ() 
  cat("Starting...\n")
  hetGP.posterior.predictives = foreach(w =1:W, .export=ls(globalenv()), .packages=c("hetGP")) %dopar% {
    obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
    setwd(local.directory)
    
    X = as.matrix(test.data[,1:D])
    posterior.predictives = predict(hetGP.model[[w]]$hetgp.model$hetgp.model, x=X)
    posterior.predictives    		
  }
  cat("registering doseq...\n")
  registerDoSEQ() 
  #stopCluster(cl)
  #registerDoSEQ() 
  
  cat("Done!\n")
  
  lcb = numeric(N.test)
  sds = numeric(N.test)
  
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
  for (w in 1:W) {lcb <- lcb + weights.w[w]*exp(hetGP.posterior.predictives[[w]]$mean)}
  for (w in 1:W) {sds <- sds + weights.w[w]*exp(hetGP.posterior.predictives[[w]]$sd2 + hetGP.posterior.predictives[[w]]$nugs)}
  
  lcb.point <- test.points[which.min(lcb),] #TR: min!
  lcb.index <- which.min(lcb) #TR: min!
  lcb.sd <- sds[lcb.index]
  
  evaluations$best[[length(evaluations$best)+1]] <- list()
  evaluations$best[[length(evaluations$best)]]$predicted =list()
  evaluations$best[[length(evaluations$best)]]$predicted$parameters <- lcb.point+0.5
  
  evaluations$best[[length(evaluations$best)]]$known =list()
  #evaluate best known
   minimum.known = numeric(dim(master.data.small.train)[1])
  for (w in 1:W) {minimum.known <- minimum.known + weights.w[w]*exp(master.data.small.train[,(2*D+w)])}
   
  evaluations$best[[length(evaluations$best)]]$known$parameters <- as.numeric(master.data.small.train[which.min(minimum.known),c(1:D)]+0.5)
  
  #if we're in the first iteration or we have a new best point, simulate
  if(main.iterations==1 || evaluations$best[[length(evaluations$best)]]$known$parameters!=evaluations$best[[length(evaluations$best)-1]]$known$parameters){
    
    cat("Main iterations =",main.iterations,": New minimum, simulating best...")
    #simulating best known, not best predicted	
    current.index<-paste0("main_",sprintf("%04i",main.iterations))
    size <- "small"
    write.xmls.small(evaluations$best[[length(evaluations$best)]]$known$parameters) # write times 61
    
    setwd(local.directory)
    system(paste0('rm -f ./PriorDraws/Output/',size,'/output_',current.index,'/*'),intern=TRUE)
    
    main.end1 = toc()
    total.time = total.time + (main.end1$toc - main.end1$tic)
    
    Eval.model.small()
    
    #Runtime = NA #Remove later
    total.time = total.time + period_to_seconds(hms(Runtime))
    
    tic("Finalising main iteration (post current best simulation)")
    observations  = list()
    setwd(local.directory)
    
    if (system(paste0('ls ./PriorDraws/Output/',size,'/output_',current.index,'/ | wc -l'),intern=TRUE)!=61){ #if the number of output files is less than 61, some runs have crashed and the location is invalid
      cat("Error: Bad location...\n ")
    }else{
      
      evaluations$best[[length(evaluations$best)]]$known = list.merge(evaluations$best[[length(evaluations$best)]]$known,obj.calc.small(evaluations$best[[length(evaluations$best)]]$known$parameters, dir="main/best/"))
      
      time[[(main.iterations-1)*N.Updates+update]]$min_LF_known = min(minimum.known, na.rm=T)
      time[[(main.iterations-1)*N.Updates+update]]$params_min_LF_known = master.data.small.train[which.min(minimum.known),c(1:D)]+0.5
      
    }
    
  }
  save(evaluations, file="eval.state.newer")
  bp=quant.to.prior(evaluations$best[[length(evaluations$best)]]$known$parameters) #best point estimate
  
  save(bp, file="best.params")
  
  main.end2 = toc()
  total.time = total.time + (main.end2$toc - main.end2$tic)
  
} 

save(hetgp.model, file="hetGP.model")
#Add stochastic gradient descent


# Repeat. (time to convergence defined manually by plots)


