u = "reiker" #Scicore user name
GitDir = "/scicore/home/smith/reiker/GitRepos/om_fitting/" # Local version of the Git repository
ExperimentDir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/" #Parent folder to experiments

ExperimentName = "2020_07_GPSG_synth_data_validation" #Name of this Experiment  (folder name)
hetgpName = NA #If using initialisation from hetGP
DataFile= "simulatedFieldData_2020_05_full_fit_small_GPSG_norm_samplingmain_0010.txt" #when not doing synthetic data validation, this should be set to "fieldData.txt"

# Fitting Options
Diagnostic_plots = FALSE
collection.of.full.gps.ON = FALSE
SCICORE = TRUE
INIT = TRUE
eval.ID = "param.evaluation.gpsg.valid" #needed as a job identifier in EvalOM. Different IDs allow for running multiple fitting runs at the same time
runInitial_gp = FALSE # Here, we're using the same setf initial simulations as in the hetGP run 
STEP1= TRUE # This is the noise estimation - currently not included but can be added back in. 
noise=FALSE #for loading an existing noise esimation object: noise=TRUE
STEP2=TRUE
CONT1=FALSE
CONT2=FALSE
Continued_run=FALSE
use_hetgp_init=FALSE #should the initialization be run or do we use initial simulations from the hetGP fitting run? Should be TRUE when comparing performance
Continue_optimization=FALSE

############################################################
#
# Load Packages
#
############################################################

if(!require(pacman)){install.packages("pacman",repos = "http://cran.us.r-project.org")}; require(pacman) 
pacman::p_load(Matrix, MASS, reshape2, tmvtnorm, gridExtra, doParallel, mlr, hetGP,dplyr, scater, tictoc, rlist,parallel, lubridate,tictoc,randtoolbox,tmvtnorm) 

install.packages("foreach",repos = "http://cran.us.r-project.org")#foreach needs to be updated from what is loaded with doParallel
if(!require(doFuture)){install.packages("doFuture",repos = "http://cran.us.r-project.org")}; require(doFuture) 

options(future.globals.maxSize= 6291456000) # Needed to allow >500Mb per thread when entering future. This number is equivalent to 6GB=6000Mb
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
nseeds = 3 #number of seeds per (initial) point


setwd(local.directory)
source(paste0(GitDir,"EvalOM.parallel.R"))
source(paste0(GitDir,"UCB_acquisition_function.r"))


############################################################################
#############################################################################
#
###                     INITIALISE THE PROBLEM
#
############################################################################
############################################################################
### Initialise problem
## Initial design points may be selected deterministically or randomly; a minimum of ~50 (better = 500) initial function evaluations are required to ensure stable behaviour of this algorithm in the initial stages of exploration; likewise draws must be split between small and large versions of model with a split no worse than ~1:10
## Since an estimate of the iid noise term on the small and large models is required one might choose to incorporate this into the initial design as below  
## The function evaluations should be stored in a list 'evaluations' with members $small and $large
## The example code below implements a simple uniform sampling design with initial draws at the prior median to estimate iid noise terms 

total.time.init=0
if(runInitial_gp) { ## TR: Above required regardless
  tic("Initialisation")
  # Melissa: For testing
  #test values
  #Nrep.noise.small <- 100 #number of initial points (multiplied with nseeds for number of sim setups; *61 for number of actually simulated scenarios)
  #Nrep.noise.large <- 10 #10	#20 = apporx. 480 min = 8h
  Nrep.noise.small= 100
  Nrep.noise.large=2
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
        evaluations$small[[i]] = list.merge(evaluations$small[[i]],obj.calc.small(evaluations$small[[i]]$parameters, dir="init/"))  # LF calculations times 61
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
        evaluations$large[[i]] = list.merge(evaluations$large[[i]],obj.calc.large(evaluations$large[[i]]$parameters, dir="init/"))  # LF calculations times 61
      } 
      
      
      setwd(local.directory)
      save(evaluations,file="eval.state.noise.logLH")
      
      init3=toc()
      total.time.init = total.time.init + (init3$toc-init3$tic)
    }
  } else if (Continued_run) {load("eval.state.newer") # TR added for additional simulations to previous runs 13.3.2019
  } else if (noise){load("eval.state.noise.logLH")} ## <<added 20190218 MP>>##
  
  
  #Note: Ewan had included a noise estimation step here. - The "nugs" object from predict.hetGP provides automated noise estimation 
  
  if(STEP2) { ## <<added 20190218 MP>>## ; ##TR: changed condition 20190320
    tic("Initialisation step2 part 1")
    
    evaluations <- list()
    evaluations$small <- list()
    evaluations$large <- list()
    
    #Ninit.small <- 8000#75 # 8000 * 5min = 40,000 min; 450 in parallel => 90 min
    #Ninit.large <- 100#20 # Probably between 10 and 50; 20 * 24min = 480 min = 8h; thus 13h for initial design
    Ninit.small=1000 #for testing. should take 6.5-10h ish (depending on the queue)
    Ninit.large=10
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
  
  total.time.init = total.time.init +  period_to_seconds(hms(Runtime))
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
}else if(use_hetgp_init){
  
  hetgp.directory = paste0(ExperimentDir,hetgpName,'/') #working directory
  setwd(hetgp.directory)
  ready = FALSE
  while(!ready){
    if("master.data.small.train" %in% system("ls",intern=TRUE)){ready=TRUE}
    Sys.sleep(10)
  }
  
  load("eval.state.init")} ## TR: changed file name 20190320

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


time=list() #time excludes initalization phase

tic("Enter main program loop")

setwd(local.directory)


start.time <- proc.time() #initialise system timekeeping
total.time = 0 #initialise "true" CPU runtime timekeeping

N.main.iterations = 200 

all_selected_points <- data.frame(Date=as.Date(character()),
                                  File=character(), 
                                  User=character(), 
                                  stringsAsFactors=FALSE)


for (main.iterations in 1:N.main.iterations) {
  
  tic("main")
  #main.iterations =1 				  
  #############################################################################
  #
  ### (1) Fit stacked generalisation model: 
  #
  ############################################################################
  ### (1) Fit stacked generalisation model: 
  #			(a)	Prepare data: Generate of load training and test data
  #			(b) Fit an ensemble of machine learning models to a ten-fold cross-validation of the currently available Eval.model.small outputs, producing an out-of-sample prediction for each observation; 
  #			(c) Fit a GP model to learn the optimal weights (on the members of the machine learning model ensemble) for predicting the Eval.model.small outputs given their out-of-sample machine learning predictions; 
  #			(d) Use the ensemble weights in fitting a joint model for the Eval.model.small and Eval.model.large outputs; 
  #			(e) Estimate runtimes and failure rates for the .small and .large models
  
  
  ############################################################################
  ### (a) Prepare data: Generate or load training and test set
  ############################################################################
  
  if(use_hetgp_init && main.iterations==1 && !Continue_optimization){
    
    setwd(hetgp.directory)
    load("eval.state.init")
    master.data.small.train <- read.table("master.data.small.train")
    master.data.small.test <- read.table("master.data.small.test")
    master.data.small.test <- read.table("master.data.small.test")
    master.data.large.train <- read.table("master.data.large.train")
    master.data.small.test.aggregated <- read.table("master.data.small.test.aggregated")	
    
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
    
    setwd(local.directory)
    write.table(master.data.small.train, file="master.data.small.train")
    write.table(master.data.small.train, file="master.data.small.test")
    write.table(master.data.small.train, file="master.data.large.train")
    write.table(master.data.small.train, file="master.data.large.test")
    write.table(master.data.small.test.aggregated, file="master.data.small.test.aggregated")
    
    n.init.sims.small = dim(unique(master.data.small.train[,c(1:2*D)]))[1] #save the number of unique design points for later
    
    save(evaluations, file="eval.state.init")
    evaluations$best = list()
    
  }else if(main.iterations==1 && !Continue_optimization){
    
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
    h1 = sample(seq(1,dim(master.data.small)[1],nseeds),dim(master.data.small)[1]*0.1) %>%sort() ##sample indeces of unique design points (10% of total unique design points) 
    h2=h1
    if(nseeds>1){
      for(s in 2:nseeds-1){h2 = sort(c(h2, h1+s))} #add replicate indeces
    }
    
    master.data.small.test <- master.data.small[h2,] 
    master.data.small.train <- master.data.small[-h2,] #complement to the holdout
    
    
    #Repeat procedure with large simulation data set
    h1 = sample(seq(1,dim(master.data.large)[1],nseeds),dim(master.data.large)[1]*0.1) %>%sort() ##sample numbers divisible by nseeds 1 and dim(master.data.small)[1] 800 rows = 10%
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
    master.data.small.test.aggregated = as.data.frame(master.data.small.test.aggregated)
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
    evaluations$best = list()
    
  }else if(Continue_optimization){
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
    #names(all_selected_points) = c(colnames.small[1:D], "lcb" )
    load("eval.state.newer")
    load("time.gpsg")
    total.time = time[[length(time)]]$time["elapsed"]
    main.iterations=length(time)+1
    N.main.iterations = 200 
    n.init.sims.small = dim(unique(master.data.small.train[,c(1:2*D)]))[1]
    
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
  }
  

  #############################################################################
  #
  #	(b) Fit an ensemble of machine learning models to a ten-fold cross-validation of the currently available Eval.model.small outputs, producing an out-of-sample prediction for each observation
  #
  ############################################################################
  
  cat("Main iterations = ",main.iterations,": Fitting Machine learning models ...", sep="")
  
  # set up parallelization
  n = detectCores()
  
  cl = parallel::makeCluster(min(n,W))
  registerDoFuture()
  future::plan(cluster, workers = cl) #tested different methods ("cluster", multicore, multiprocess, and sequential. "cluster" was the fastest. )
  
  ensemble.of.machine.learning.models =  foreach(w=1:W) %dopar% { #this takes ~4 minutes if parallelized
    output = list()
    output$mars.model = list()
    output$randomforest.model = list()
    output$brnn.model = list()
    
    
    #PROBLEM SET UP
    ### These descrobe things common to all learning procedures
    obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
    setwd(local.directory)
    
    #data set:
    parnames = colnames.small[1:(2*D)]# linear and squared predictors
    col.names = c(obj.w, parnames)
    obj.w.data =  master.data.small.train %>%
      select(col.names)
    master.data.large.sub =  master.data.large.train %>% select(col.names)
    #Create learning Task (makeRegrTask can't handle missing values)
    tsk = makeRegrTask(data = na.omit(obj.w.data), target = paste0("OBJECTIVE",sprintf("%02i",w))) #checkdata False for speed 	
    #Prepare for resampling:
    ### Resampling evaluates a learner on a given machine learning taslk using the selected resampling strategy
    rdesc = makeResampleDesc("CV", iters = 10) #defines resampling strategy: 10fold- Cross Validation
    rinst = makeResampleInstance(rdesc,tsk) #fixes the split: using risnt rather than rdesc makes this repeatable; see set.seed at the beginning 
    
    
    ##################################
    #MARS
    #note that MARS doesn't have hyperparameters to tune! [ getHyperPars(makeLearner("regr.mars")) = list() ]
    
    # Train learner:
    # Fit MARS to all available Eval.model.small outputs in training set 
    # note: add makeLearner command first if we want to change any default values (e.g. hyperparameter values)
    output$mars.model$mars.model = train("regr.mars", tsk) # tsk defined above
    # Save within sample estimates for Eval.model.small outputs
    output$mars.model$Ypred = predict(output$mars.model$mars.model, newdata=obj.w.data)$data$response	
    # Save out-of-sample estimates for Eval.model.large outputs		
    output$mars.model$Ypred.largepop = predict(output$mars.model$mars.model,newdata= master.data.large.sub)$data$response
    # Assess performance and get 10-fold CV predictions
    # The defaut performance measure for regression learners is MSE (but we don't actually need this, we just need the CV preds)
    mars.res = resample("regr.mars", tsk, rinst, show.info= FALSE) # carry out the resampling and store under res
    #save each out-of-sample prediction from the cross-validation on Eval.model.small outputs
    
    obj.w.data$index = NA
    obj.w.data[!is.na(obj.w.data[,obj.w]),"index"] = mars.res$pred$data[order(mars.res$pred$data$id),]$response
    output$mars.model$Ypred_CV  = obj.w.data$index
    
    ##################################
    #RANDOM FOREST	(No tuning)
    
    # Train learner:
    # Fit Ranger RF to all available Eval.model.small outputs in training set 
    # note: add makeLearner command first if we want to change any default values (e.g. hyperparameter values)
    output$randomforest.model$rf.model = train("regr.ranger", tsk) 
    # Save within sample estimates for Eval.model.small outputs
    output$randomforest.model$Ypred = predict(output$randomforest.model$rf.model, newdata=obj.w.data)$data$response    # Save best out-of-sample estimates for Eval.model.large outputs		
    # Save out-of-sample estimates for Eval.model.large outputs		
    output$randomforest.model$Ypred.largepop = predict(output$randomforest.model$rf.model,newdata= master.data.large.sub)$data$response
    # Assess performance and get 10-fold CV predictions
    # The defaut performance measure for regression learners is MSE (but we don't actually need this, we just need the CV preds)
    rf.res = resample("regr.ranger", tsk, rinst, show.info= FALSE) # carry out the resampling and store under res
    #save each out-of-sample prediction from the cross-validation on Eval.model.small outputs
    
    obj.w.data$index = NA
    obj.w.data[!is.na(obj.w.data[,obj.w]),"index"] = rf.res$pred$data[order(rf.res$pred$data$id),]$response
    output$randomforest.model$Ypred_CV    = obj.w.data$index
    
    
    ##################################
    # BRNN
    
    # Train learner:
    # Fit Ranger RF to all available Eval.model.small outputs in training set 
    # note: add makeLearner command first if we want to change any default values (e.g. hyperparameter values)
    output$brnn.model$brnn.model = train("regr.brnn", tsk) # the model already has the hyperparams stored; lambdamin under output$gbm.model$learner.model$lambda.min
    # Save within sample estimates for Eval.model.small outputs
    output$brnn.model$Ypred = predict(output$brnn.model$brnn.model, newdata=obj.w.data)$data$response	
    # Save out-of-sample estimates for Eval.model.large outputs		
    output$brnn.model$Ypred.largepop = predict(output$brnn.model$brnn.model,newdata= master.data.large.sub)$data$response
    # Assess performance and get 10-fold CV predictions
    # The defaut performance measure for regression learners is MSE (but we don't actually need this, we just need the CV preds)
    brnn.res = resample("regr.brnn", tsk, rinst, show.info= FALSE) # carry out the resampling and store under res
    #save each out-of-sample prediction from the cross-validation on Eval.model.small outputs
    
    obj.w.data$index = NA
    obj.w.data[!is.na(obj.w.data[,obj.w]),"index"] = brnn.res$pred$data[order(brnn.res$pred$data$id),]$response
    output$brnn.model$Ypred_CV = obj.w.data$index
    
    
    output	
    
  }
  
  cat("Done!\n")
  
  ## CLEANUP
  parallel::stopCluster(cl)
  plan(sequential)
  
  
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  
  weighted.truths = numeric(dim(master.data.small.test.aggregated)[1])
  for (w in 1:W) {weighted.truths <- weighted.truths + (weights.w[w]*exp(master.data.small.test.aggregated[,(2*D+w)]))} # for each parameter vector, compute
  
  #Plot hold out predictions:
  pdf(file=paste0("3ML.holdout.",main.iterations,".blue.pdf"))
  
  
  p.mars=list()	
  mars.pred=list()
  
  for(w in 1:W){
    mars.pred[[w]] = predict(ensemble.of.machine.learning.models[[w]]$mars.model$mars.model, newdata = master.data.small.test.aggregated[,1:46])
    dat = cbind(master.data.small.test.aggregated[,46+w],mars.pred[[w]]$data[,1]) 
    colnames(dat) = c("truth","response")
    dat = dat %>% as.data.frame() %>% na.omit()
    r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
    
    p.mars[[w]]= ggplot(dat, aes(truth, response)) + 
      geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
      geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
      labs(title = paste0("Objective ",w, ", R2=",round(r2,3)),x = "Truth", y="MARS Response", size= 8) +
      theme_minimal() +
      coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response), max(dat$truth, dat$response)), ylim = c(min(dat$truth, dat$response), max(dat$truth, dat$response))) + 
      theme(axis.line = element_line(colour = "black"),
            plot.title=element_text(size=7.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
  }
  
  # weighted predictions
  weighted.mars = numeric(dim(master.data.small.test.aggregated)[1]) #nrows in test frame
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
  
  for (w in 1:W) {weighted.mars <- weighted.mars + (weights.w[w]*exp((mars.pred[[w]]$data)))} 
  
  dat = cbind(weighted.mars,weighted.truths) 
  colnames(dat) = c("truth","response")
  dat = dat %>% as.data.frame() %>% na.omit()
  
  r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
  p.mars[[W+1]]= ggplot(dat, aes(truth, response)) + 
    geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
    geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
    labs(title = paste0("w. sum, R2=",round(r2,4)),x = "Truth", y="MARS Response", size= 8) +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
    theme(axis.line = element_line(colour = "black"),
          plot.title=element_text(size=7.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  multiplot(plotlist = p.mars, cols = 4)
  
  
  
  
  p.rf=list()	
  rf.pred=list()
  for(w in 1:W){
    rf.pred[[w]] = predict(ensemble.of.machine.learning.models[[w]]$randomforest.model$rf.model, newdata = master.data.small.test.aggregated[,1:46])
    dat = cbind(master.data.small.test.aggregated[,46+w],rf.pred[[w]]$data[,1]) 
    colnames(dat) = c("truth","response")
    dat = dat %>% as.data.frame() %>% na.omit()
    r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
    
    p.rf[[w]]= ggplot(dat, aes(truth, response)) + 
      geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
      geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
      labs(title = paste0("Objective ",w, ", R2=",round(r2,3)),x = "Truth", y="RF Response", size= 8) +
      theme_minimal() +
      coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response), max(dat$truth, dat$response)), ylim = c(min(dat$truth, dat$response), max(dat$truth, dat$response))) + 
      theme(axis.line = element_line(colour = "black"),
            plot.title=element_text(size=7.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
  }
  # weighted predictions
  weighted.rf = numeric(dim(master.data.small.test.aggregated)[1]) #nrows in test frame
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
  
  for (w in 1:W) {weighted.rf <- weighted.rf + (weights.w[w]*exp((rf.pred[[w]]$data)))} 
  
  dat = cbind(weighted.rf,weighted.truths) 
  colnames(dat) = c("truth","response")
  dat = dat %>% as.data.frame() %>% na.omit()
  
  r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
  p.rf[[W+1]]= ggplot(dat, aes(truth, response)) + 
    geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
    geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
    labs(title = paste0("w. sum, R2=",round(r2,4)),x = "Truth", y="RF Response", size= 8) +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
    theme(axis.line = element_line(colour = "black"),
          plot.title=element_text(size=7.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  multiplot(plotlist = p.rf, cols = 4)
  
  
  
  p.brnn=list()	
  brnn.pred=list()
  for(w in 1:W){
    brnn.pred[[w]] = predict(ensemble.of.machine.learning.models[[w]]$brnn.model$brnn.model, newdata = master.data.small.test.aggregated[,1:46])
    dat = cbind(master.data.small.test.aggregated[,46+w],brnn.pred[[w]]$data[,1]) 
    colnames(dat) = c("truth","response")
    dat = dat %>% as.data.frame() %>% na.omit()
    r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
    
    p.brnn[[w]]= ggplot(dat, aes(truth, response)) + 
      geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
      geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
      labs(title = paste0("Objective ",w, ", R2=",round(r2,3)),x = "Truth", y="BRNN Response", size= 8) +
      theme_minimal() +
      coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response), max(dat$truth, dat$response)), ylim = c(min(dat$truth, dat$response), max(dat$truth, dat$response))) + 
      theme(axis.line = element_line(colour = "black"),
            plot.title=element_text(size=7.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
  }
  
  
  # weighted predictions
  weighted.brnn = numeric(dim(master.data.small.test.aggregated)[1]) #nrows in test frame
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
  
  for (w in 1:W) {weighted.brnn <- weighted.brnn + (weights.w[w]*exp((brnn.pred[[w]]$data)))} 
  
  dat = cbind(weighted.brnn,weighted.truths) 
  colnames(dat) = c("truth","response")
  dat = dat %>% as.data.frame() %>% na.omit()
  
  r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
  p.brnn[[W+1]]= ggplot(dat, aes(truth, response)) + 
    geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
    geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
    labs(title = paste0("w. sum, R2=",round(r2,4)),x = "Truth", y="BRNN Response", size= 8) +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
    theme(axis.line = element_line(colour = "black"),
          plot.title=element_text(size=7.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  multiplot(plotlist = p.brnn, cols = 4)
  dev.off()	
  
  
  
  ############################################################################
  #  (b) Fit a hetGP model to learn the optimal weights
  ############################################################################
  ### Fit het GP on out of sample (cross-validation) predictions of each of the level 0 learners
  
  #start parallelization
  cl = parallel::makeCluster(min(n,W)) #cl is defined above but needs to be reset, otherwise "Error in serialize(data, node$con) : connection is not open" 
  
  registerDoFuture()
  future::plan(cluster, workers = cl) 
  
  collection.of.stacking.weights <- foreach(w=1:W,.export="master.data.small.train") %dopar% {
    opt = list()
    obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
    setwd(local.directory)
    
    #Prepare data; remember this order for later!
    Y.brnn <- ensemble.of.machine.learning.models[[w]]$brnn.model$Ypred_CV
    Y.mars <- ensemble.of.machine.learning.models[[w]]$mars.model$Ypred_CV
    Y.rf <- ensemble.of.machine.learning.models[[w]]$randomforest.model$Ypred_CV
    
    #Extract simulation results from training data frame
    eval(parse(text=paste0("Y <- master.data.small.train$",obj.w)))
    
    Y.all = cbind(Y.brnn, Y.mars, Y.rf,Y) %>% na.omit() # - note na.omit!
    X = subset(Y.all, select = -c(Y))
    Y = na.omit(Y)
    
    #standardise input data
    X =scale(X)
    opt$scaling=list()
    #important: store the scaling parameters for later use during prediction!
    opt$scaling$center = t(as.matrix(attributes(X)["scaled:center"][["scaled:center"]])) 
    opt$scaling$scale = t(as.matrix(attributes(X)["scaled:scale"][["scaled:scale"]]))	
    
    reps =find_reps(X=X, Z=Y, rescale = FALSE, normalize = FALSE) #rescaled manually above
    nvar = dim(X)[2]
    
    #fit hetGP
    opt$model = mleHetGP(X=list(X0 = reps$X0, Z0=reps$Z0, mult = reps$mult), Z = reps$Z, lower = rep(0.1, nvar), upper = rep(10, nvar), covtype = "Matern5_2")
    
    opt
    
  } 
  ## CLEANUP
  parallel::stopCluster(cl)
  plan(sequential)
  
  
  #Generate hold-out predictions for model assessment.
  
  #cl = makeCluster(min(n,W)) #
  #registerDoFuture()
  #plan(cluster, workers = cl)
  # not sure why this crashes in parallel mode... keep sequential, it's fast anyway. - Because the predict function spikes the CPU usage!
  
  
  holdout.preds = foreach(w =1:W, .export="master.data.small.test.aggregated") %dopar% {
    obj.w = paste0("OBJECTIVE",sprintf("%02i",w))
    
    out = list()
    out$brnn=list()
    out$mars=list()
    out$rf=list()
    out$gpsg=list()
    
    out$brnn = predict(ensemble.of.machine.learning.models[[w]]$brnn.model$brnn.model, newdata = master.data.small.test.aggregated[,1:46])
    out$mars = predict(ensemble.of.machine.learning.models[[w]]$mars.model$mars.model, newdata = master.data.small.test.aggregated[,1:46])
    out$rf = predict(ensemble.of.machine.learning.models[[w]]$randomforest.model$rf.model, newdata =master.data.small.test.aggregated[,1:46])
    
    #scale 
    Y.brnn.test.scaled = scale(as.data.frame(out$brnn), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.brnn"], 
                               scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.brnn"])
    Y.mars.test.scaled = scale(as.data.frame(out$mars), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.mars"], 
                               scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.mars"])
    Y.rf.test.scaled = scale(as.data.frame(out$rf), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.rf"], 
                             scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.rf"])
    
    out$gpsg= predict(collection.of.stacking.weights[[w]]$model, x=as.matrix(cbind(Y.brnn.test.scaled, Y.mars.test.scaled,Y.rf.test.scaled)))
    out
  }
  
  
  #Truth vs predicted plots
  p=list()
  for(w in 1:W){ 
    dat = cbind(master.data.small.test.aggregated[,(2*D+w)],as.data.frame(holdout.preds[[w]]$gpsg$mean)) 
    colnames(dat) = c("truth","response")
    dat = dat %>% as.data.frame() %>% na.omit()
    
    r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
    p[[w]]= ggplot(dat, aes(truth, response)) + 
      geom_point(size = 1, alpha=0.35, colour = "deepskyblue4")+
      geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
      labs(title = paste0("Objective ",w, ", R2=",round(r2,4)),x = "Truth", y="GPSG Response", size= 8) +
      theme_minimal() +
      coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
      theme(axis.line = element_line(colour = "black"),
            plot.title=element_text(size=7.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
  }
  
  
  #add weighted predictions
  weighted.truths = numeric(dim(master.data.small.test.aggregated)[1])
  weighted.preds = numeric(dim(master.data.small.test.aggregated)[1]) #nrows in test frame
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
  
  for (w in 1:W) {weighted.truths <- weighted.truths + (weights.w[w]*exp(master.data.small.test.aggregated[,(2*D+w)]))} # for each parameter vector, compute
  for (w in 1:W) {weighted.preds <- weighted.preds + (weights.w[w]*exp((holdout.preds[[w]]$gpsg$mean)))} 
  
  dat = cbind(weighted.preds,weighted.truths) 
  colnames(dat) = c("truth","response")
  dat = dat %>% as.data.frame() %>% na.omit()
  
  r2 =cor(dat$truth, dat$response,use = "complete.obs")^2
  p[[W+1]]= ggplot(dat, aes(truth, response)) + 
    geom_point(size = 1, alpha=0.35, colour = "deepskyblue4")+
    geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
    labs(title = paste0("w. sum, R2=",round(r2,4)),x = "Truth", y="Response", size= 8) +
    theme_minimal() +
    coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
    theme(axis.line = element_line(colour = "black"),
          plot.title=element_text(size=7.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  
  pdf(file=paste0("3ml.gpsg.holdout.",main.iterations,".blue.pdf"))
  multiplot(plotlist = p, cols = 4)
  dev.off()
  
  
  
  ############################################################################
  #
  ### (2) Update model N.Updates times: 
  #
  ############################################################################
  ### (1) Update model N.Updates times: 
  #			(a) Make Fine grid predictions
  #			(b) Choose min (UCB - SD) as next simulation point
  #			(c) Simulate
  #			(d) Update model
  
  
  
  if(main.iterations==1 && use_hetgp_init){setwd(hetgp.directory); load("eval.state.init"); setwd(local.directory)}else if (Continue_optimization){load("eval.state.newer")}

  Ninit.small = n.init.sims.small + dim(all_selected_points)[1]
  
  N.Updates = 1 #updating doesn't seem to improve estimation much
  n.select = 250
  n.seeds = 1 #how many seeds should be included for each new point? Since the nugget effect is already estimated, we're leaving this at one for now, but we could include more seeds here
  
  
  #data frame to store new simulation outputs
  update.data.small <- matrix(NA,N.Updates*n.select*n.seeds,ncol=(2*D+W))
  colnames.small <- list()
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("PARAM",sprintf("%02i",i))}
  for (i in 1:D) {colnames.small[[length(colnames.small)+1]] <- paste0("SQ_PARAM",sprintf("%02i",i))}
  for (i in 1:W) {colnames.small[[length(colnames.small)+1]] <- paste0("OBJECTIVE",sprintf("%02i",i))}
  colnames.small <- unlist(colnames.small)
  colnames(update.data.small) <- colnames.small
  update.data.small <- as.data.frame(update.data.small)
  
  cat("Main iterations =",main.iterations,",: Start updating algorithm ... \n")
  
  #Updating is a way of adding data to an existing hetGP model without re-estimating theta
  main = toc()
  main.elapsed = main$toc-main$tic
  
  total.time = total.time +main.elapsed #update total.time object
  

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
    
    test.points <- rtmvnorm(N.test,proposal.means,proposal.covs, lower = rep(-0.5,D), upper=rep(0.5,D),algorithm="gibbs") #+0.5
  
    test.data <- cbind(test.points,test.points^2,matrix(0,nrow=N.test,ncol=W))
    colnames.test <- list()
    for (d in 1:D) {colnames.test[[length(colnames.test)+1]] <- paste0("PARAM",sprintf("%02i",d))}
    for (d in 1:D) {colnames.test[[length(colnames.test)+1]] <- paste0("SQ_PARAM",sprintf("%02i",d))}
    for (w in 1:W) {colnames.test[[length(colnames.test)+1]] <- paste0("OBJECTIVE",sprintf("%02i",w))}
    colnames.test <- unlist(colnames.test)
    colnames(test.data) <- colnames.test
    test.data = as.data.frame(test.data)
    
    cat("Main iterations =",main.iterations,", update = ",update,": Computing posterior predictives ...")
    
    
    
    # GPSG predict
    
    cl = makeCluster(min(n,W)) #
    registerDoFuture()
    #plan(cluster, workers = cl)
    plan(sequential) #runs sequentially, unsure why it crashes in parallelization: Predict function takes too much memory. Add parallelization for submitted job, run sequentially for testing.
    
    gpsg.posterior.predictives = foreach(w =1:W, .export=c("ensemble.of.machine.learning.models","test.data","collection.of.stacking.weights")) %dopar% {
      obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
      setwd(local.directory)		
      
      posterior.predictions =list()
      
      posterior.predictions$brnn=list()
      posterior.predictions$mars=list()
      posterior.predictions$rf=list()
      posterior.predictions$gpsg=list()
      
      #ML holdout predictions
      posterior.predictions$brnn = predict(ensemble.of.machine.learning.models[[w]]$brnn.model$brnn.model, newdata = test.data[,1:(2*D)])
      posterior.predictions$mars = predict(ensemble.of.machine.learning.models[[w]]$mars.model$mars.model, newdata = test.data[,1:(2*D)])
      posterior.predictions$rf = predict(ensemble.of.machine.learning.models[[w]]$randomforest.model$rf.model, newdata = test.data[,1:(2*D)])
      
      #scale
      Y.brnn.test.scaled = scale(as.data.frame(posterior.predictions$brnn), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.brnn"], 
                                 scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.brnn"])
      Y.mars.test.scaled = scale(as.data.frame(posterior.predictions$mars), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.mars"], 
                                 scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.mars"])
      Y.rf.test.scaled = scale(as.data.frame(posterior.predictions$rf), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.rf"], 
                               scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.rf"])
      
      #stacked prediction
      posterior.predictions$gpsg= predict(collection.of.stacking.weights[[w]]$model, x=cbind(Y.brnn.test.scaled, Y.mars.test.scaled,Y.rf.test.scaled))
      
      
      posterior.predictions    		
    }
    parallel::stopCluster(cl)
    plan(sequential)
    
    
    cat("Done!\n")
    
    
  
    ############################################################################
    ### (b) choose max UCB + SD (UCB acquisition function like in Bayesian Optimization)
    ############################################################################ 
    
    lcb = numeric(N.test)
    mpe = numeric(N.test) #minimum point estimate
    
    k = kappa_calc(t=Ninit.small,d=3,delta=0.1) #for computing LCB (see Brochu 20210 and Srinivas 20XX)
    #Weighted sum of upper confidence bounds of objective estimates
    weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
    #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
    
    for (w in 1:W) {lcb <- lcb + (weights.w[w]*exp((gpsg.posterior.predictives[[w]]$gpsg$mean-k*(gpsg.posterior.predictives[[w]]$gpsg$sd2+gpsg.posterior.predictives[[w]]$gpsg$nugs))))} # lower confidence bound
    for (w in 1:W) {mpe <- mpe + (weights.w[w]*exp(gpsg.posterior.predictives[[w]]$gpsg$mean))} # lower confidence bound
    
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
    
    
    #add only one realisation to evaluations object
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
    
    
    setwd(local.directory)

    time.pre.sim = toc()
    total.time = total.time + (time.pre.sim$toc-time.pre.sim$tic) #add elapsed time up to simulation to total.time
    print(total.time)
    #Run simulations
    size <<- "small"
    Eval.model.small()
    
    total.time.init = total.time.init + period_to_seconds(hms(Runtime))
    print(total.time)
    
    tic("Updating (post simulation)")
    observations  = list()
    setwd(local.directory)
    
    
    #Calculate LF for new points
    errors = 0
    for (a in 1:dim(selected_points)[1]) {  
      current.index=paste0("main_",sprintf("%04i",main.iterations),"_",sprintf("%04i",a))
      update.data.small[((update-1)*dim(selected_points)[1])+a,(1:D)] <- evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$parameters-0.5
      update.data.small[((update-1)*dim(selected_points)[1])+a,(D+1):(2*D)] <- (evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$parameters-0.5)^2
      
      #if the number of output files is less than 61, some runs have crashed and the location is invalid
      if (system(paste0('ls ./PriorDraws/Output/',size,'/output_',current.index,'/ | wc -l'),intern=TRUE)!=61){ 
        cat("Error: Bad location... Try elsewhere.\n ")
        errors = errors+1
        cat("Error count = ",errors,".\n ")
        
      }else{
        size <<- "small"
        evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]  = list.merge(evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]],obj.calc.small(evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$parameters, dir="main/updated/"))  # LF calculations
        
        save(evaluations, file="eval.state.newer") # note that Objectives in the evaluations object are not log transformed
        

        update.data.small[((update-1)*dim(selected_points)[1])+a,(2*D+1):(2*D+W)] <- evaluations$small[[length(evaluations$small)-dim(selected_points)[1]+a]]$Objectives %>% log() #include log transform for easier merging
        
        
        #create updated training object
        master.data.small.train <- rbind(master.data.small.train,update.data.small[((update-1)*dim(selected_points)[1])+a,])
        
      }
    }
    
    #As long as the new simulations didn't all fail
    if(errors != n.select){
      cat("Done!\n")
      cat("Main iterations =",main.iterations,", update = ",update,": Updating Stacked model ... \n")
      
      ############################################################################
      ### (d) Refit ML algorithms but only UPDATE hetGP model with new data
      ############################################################################  
      cat("Updating ensemble of machine learning models... \n")
      
      
      #don't completely refit ML models, otherwise we have to completely retrain hetGP
      #just predict on new point with MLs and add new point (LF ~ MLpreds) to update hetGP
      
      #ML preds on new point; Xnew  = evaluations$small[[length(evaluations$small)]]$parameters;  => MLpreds
      plan(sequential) #so fast, doesn't matter if sequentialor parallel... can make this parallel later if we add more points per iteration
      MLpreds.new = foreach(w = 1:W, .export="n.select") %dopar% {
        
        output=list()
        obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
        setwd(local.directory)
        
        preds.new = list()
        preds.new$brnn = list()
        preds.new$mars = list()
        preds.new$rf = list()
        
        #prepeare data
        parnames = colnames.small[1:(2*D)]
        col.names = c(obj.w, parnames)
        obj.w.data = update.data.small[c(((update-1)*(dim(selected_points)[1])+1):(update*dim(selected_points)[1])),] %>% #only add latest row each iteration because hetGP.model gets updated and should include previous rows
          dplyr::select(parnames) %>% na.omit()
        

        #ML holdout predictions
        preds.new$brnn = predict(ensemble.of.machine.learning.models[[w]]$brnn.model$brnn.model, newdata = obj.w.data)$data$response
        preds.new$mars = predict(ensemble.of.machine.learning.models[[w]]$mars.model$mars.model, newdata = obj.w.data)$data$response
        preds.new$rf = predict(ensemble.of.machine.learning.models[[w]]$randomforest.model$rf.model, newdata = obj.w.data)$data$response
        
        preds.new
        
      }
      
      
      #for hetGP: Xnew= MLpreds ; Znew = evaluations$small[[length(evaluations$small)]]$Objectives
      
      cl = parallel::makeCluster(min(n,W)) 
      
      registerDoFuture()
      #future::plan(cluster, workers = cl) 
      future::plan(sequential) #for testing
      
      collection.of.stacking.weights <- foreach(w=1:W,.export=c("update.data.small","update", "n.select")) %dopar% {
        opt = list()
        obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
        setwd(local.directory)
        
        #Extract simulation results from update data
        eval(parse(text=paste0("Y <- update.data.small[c(((update-1)*(dim(selected_points)[1])+1):(update*dim(selected_points)[1])),'",obj.w,"']")))
        
        
        #standardise input data based on previous
        #scale
        Y.brnn.test.scaled = scale(as.data.frame(MLpreds.new[[w]]$brnn), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.brnn"], 
                                   scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.brnn"])
        Y.mars.test.scaled = scale(as.data.frame(MLpreds.new[[w]]$mars), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.mars"], 
                                   scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.mars"])
        Y.rf.test.scaled = scale(as.data.frame(MLpreds.new[[w]]$rf), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.rf"], 
                                 scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.rf"])
        
        #keep scaling while updating collection of stacking weights object
        opt$scaling=list()
        
        opt$scaling$center=collection.of.stacking.weights[[w]]$scaling$center
        opt$scaling$scale=collection.of.stacking.weights[[w]]$scaling$scale
        
        #prepare update data
        Y.all = cbind(Y.brnn.test.scaled, Y.mars.test.scaled, Y.rf.test.scaled,Y) %>% na.omit() # - note na.omit!
        colnames(Y.all) = c("Y.brnn", "Y.mars", "Y.rf", "Y")
        
        Xnew = subset(Y.all, select = -c(Y))
        Znew = na.omit(Y)
        
        
        #newdata =find_reps(X=Xnew, Z=Znew, rescale = FALSE, normalize = FALSE)
        #opt$model = update(object=collection.of.stacking.weights[[w]]$model, X=newdata$X0, Z = newdata$Z0, maxit =0, maxiter = 0, ginit = 0.01, method = "quick")      
        
        #update hetGP
        opt$model = update(object=collection.of.stacking.weights[[w]]$model, X=Xnew, Z = Znew, maxit =0, maxiter = 0, ginit = 0.01, method = "quick")      
        
        
        
        opt
        
      } 
      ## CLEANUP
      parallel::stopCluster(cl)
      plan(sequential)
      
      write.table(update.data.small, file=paste0("update.data.small",main.iterations))  
    }else{cat("All new points in bad locations... Moving on. \n")}
    
    
    
    cat("Main iterations =",main.iterations,", update = ",update,": Plotting ...")
    
    holdout.preds = foreach(w =1:W, .export=c("master.data.small.test.aggregated", "ensemble.of.machine.learning.models", "collection.of.stacking.weights")) %dopar% {
      obj.w = paste0("OBJECTIVE",sprintf("%02i",w))
      
      out = list()
      out$brnn=list()
      out$mars=list()
      out$rf=list()
      out$gpsg=list()
      
      out$brnn = predict(ensemble.of.machine.learning.models[[w]]$brnn.model$brnn.model, newdata = master.data.small.test.aggregated[,1:(2*D)])
      out$mars = predict(ensemble.of.machine.learning.models[[w]]$mars.model$mars.model, newdata = master.data.small.test.aggregated[,1:(2*D)])
      out$rf = predict(ensemble.of.machine.learning.models[[w]]$randomforest.model$rf.model, newdata =master.data.small.test.aggregated[,1:(2*D)])
      
      #scale 
      Y.brnn.test.scaled = scale(as.data.frame(out$brnn), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.brnn"], 
                                 scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.brnn"])
      Y.mars.test.scaled = scale(as.data.frame(out$mars), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.mars"], 
                                 scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.mars"])
      Y.rf.test.scaled = scale(as.data.frame(out$rf), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.rf"], 
                               scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.rf"])
      
      out$gpsg= predict(collection.of.stacking.weights[[w]]$model, x=as.matrix(cbind(Y.brnn.test.scaled, Y.mars.test.scaled,Y.rf.test.scaled)))
      out
    }
    
    
    #Truth vs predicted plots
    p=list()
    r2 = list()
    for(w in 1:W){ 
      dat = cbind(master.data.small.test.aggregated[,46+w],as.data.frame(holdout.preds[[w]]$gpsg$mean)) 
      colnames(dat) = c("truth","response")
      dat = dat %>% as.data.frame() %>% na.omit()
      
      r2[[w]] =cor(dat$truth, dat$response,use = "complete.obs")^2
      p[[w]]= ggplot(dat, aes(truth, response)) + 
        geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
        geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
        labs(title = paste0("Objective ",w, ", R2=",round(r2[[w]],3)),x = "Truth", y="GPSG Response", size= 8) +
        theme_minimal() +
        coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
        theme(axis.line = element_line(colour = "black"),
              plot.title=element_text(size=7.5),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()) 
    }
    
    
    #add weighted predictions
    weighted.truths = numeric(dim(master.data.small.test.aggregated)[1])
    weighted.preds = numeric(dim(master.data.small.test.aggregated)[1]) #nrows in test frame
    weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
    #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
    
    for (w in 1:W) {weighted.truths <- weighted.truths + (weights.w[w]*exp(master.data.small.test.aggregated[,2*D+w]))} # for each parameter vector, compute
    for (w in 1:W) {weighted.preds <- weighted.preds + (weights.w[w]*exp((holdout.preds[[w]]$gpsg$mean)))} 
    
    dat = cbind(weighted.preds,weighted.truths) 
    colnames(dat) = c("truth","response")
    dat = dat %>% as.data.frame() %>% na.omit()
    
    r2[[W+1]] =cor(dat$truth, dat$response,use = "complete.obs")^2
    p[[W+1]]= ggplot(dat, aes(truth, response)) + 
      geom_point(size = 1, alpha=0.2, colour = "deepskyblue4")+
      geom_abline(slope = 1, intercept = 0, colour = "black", size = .7,linetype = "dashed") + 
      labs(title = paste0("w. sum, R2=",round(r2[[W+1]],3)),x = "Truth", y="Response", size= 8) +
      theme_minimal() +
      coord_fixed(ratio = 1, xlim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE)), ylim = c(min(dat$truth, dat$response,na.rm=TRUE), max(dat$truth, dat$response,na.rm=TRUE))) + 
      theme(axis.line = element_line(colour = "black"),
            plot.title=element_text(size=7.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) 
    
    pdf(file=paste0("3ml.gpsg.holdout.",main.iterations,".updated.pdf"))
    multiplot(plotlist = p, cols = 4)
    dev.off()
    
    
    ############################################################################
    ### (e) plot current performance
    ############################################################################  
    
    
    time.post.sim = toc()
    total.time = total.time + (time.post.sim$toc-time.post.sim$tic)
    print(total.time)
    tic("Housekeeping")
    #time keeping
    time[[(main.iterations-1)*N.Updates+update]] = list()  
    time[[(main.iterations-1)*N.Updates+update]]$time= proc.time() - start.time 
    time[[(main.iterations-1)*N.Updates+update]]$cpu.time = total.time["elapsed"]
    time[[(main.iterations-1)*N.Updates+update]]$r2 = r2 
    time[[(main.iterations-1)*N.Updates+update]]$obj$min_lcb = min(lcb)
    time[[(main.iterations-1)*N.Updates+update]]$params_lcb = test.points[which.min(lcb),]
    
    time[[(main.iterations-1)*N.Updates+update]]$obj$mpe = min(mpe)
    time[[(main.iterations-1)*N.Updates+update]]$params_mpe = test.points[which.min(mpe),]    
    
    
    save(time, file="time.gpsg")
    time.frame <- matrix(0,nrow=length(time),(W+4))
    for (g in 1:length(time)) {
      time.frame[g,1] <- as.numeric(time[[g]]$cpu.time["elapsed"])
      time.frame[g,c(2:(2+W))] <- unlist(time[[g]]$r2)
      time.frame[g,(W+3)] <- time[[g]]$obj$min_lcb
      time.frame[g,(W+4)] <- time[[g]]$obj$mpe
    }
    
    
    colnames <- list()
    colnames[[length(colnames)+1]] <- paste0("cpu.time")
    for (w in 1:W) {colnames[[length(colnames)+1]] <- paste0("r2_",sprintf("%02i",w))}
    colnames[[length(colnames)+1]] <- paste0("r2_w")
    colnames[[length(colnames)+1]] <- paste0("min_lcb")
    colnames[[length(colnames)+1]] <- paste0("mpe")
    
    colnames <- unlist(colnames)
    colnames(time.frame) <- colnames
    time.frame <- as.data.frame(time.frame)
    
    d <- melt(time.frame[,1:(W+2)], id.vars="cpu.time")
    
    #Plot current learner performance
    pdf("current.learner.performance.gpsg.pdf")
    # Everything on the same plot
    print(ggplot(d, aes(cpu.time,value, col=variable)) + 
            geom_line()) 
    dev.off()
    
    #Convergence plot for current minimum
    pdf("min_lcb.progress.gpsg.pdf")
    
    plot(time.frame$min_lcb~seq(1,length(time.frame$min_lcb)), type="l", col="red",xlab="iteration", ylab="current minimum", lty=2,ylim=c(min(time.frame$min_lcb,time.frame$mpe),max(time.frame$min_lcb,time.frame$mpe)))
    lines(x=seq(1,length(time.frame$min_lcb)), y=time.frame$mpe, col="black")
    
    plot(time.frame$min_lcb~time.frame$cpu.time, type="l", col="red",xlab="time", ylab="current minimum", lty=2,ylim=c(min(time.frame$min_lcb,time.frame$mpe),max(time.frame$min_lcb,time.frame$mpe)))
    lines(x=time.frame$cpu.time, y=time.frame$mpe, col="black")
    
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
  
  N.test <- 500000 # Number of random points at which to estimate LCB: 500,000
 #mvnorm sampling???
  
  #test.points <- matrix(runif(N.test*D),nrow=N.test,ncol=D)-0.5 #Draw points from uniform dist. Don't Sobol sample here because Sobol is deterministic
  if(main.iterations==1){
    proposal.means <- rep(0,D)
    proposal.covs <- diag(1,D)
    }else{
    tmp = exp(master.data.small.train[,(2*D+1)])
    for(w in 2:W){tmp = tmp+ exp(master.data.small.train[,(2*D+w)])}
  
    proposal.means =as.numeric(as.vector(master.data.small.train[which.min(tmp),c(1:D)]))
    proposal.covs <- cov(data.matrix(master.data.small.train[c(Ninit.small:dim(master.data.small.train)[1]),c(1:D)]),use="complete.obs")
    rm(tmp)}
    
  test.points <- rtmvnorm(N.test,proposal.means,proposal.covs, lower = rep(-0.5,D), upper=rep(0.5,D),algorithm="gibbs") #+0.5
  
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
  # GPSG predict
  
  cl = makeCluster(min(n,W)) #
  registerDoFuture()
  #plan(cluster, workers = cl)
  plan(sequential) #runs sequentially, unsure why it crashes in parallelization: Predict function takes too much memory. Add parallelization for submitted job, run sequentially for testing.
  
  gpsg.posterior.predictives = foreach(w =1:W, .export=ls(globalenv())) %dopar% {
    obj.w = paste0("OBJECTIVE",sprintf("%02i",w))	
    setwd(local.directory)		
    
    posterior.predictions =list()
    
    posterior.predictions$brnn=list()
    posterior.predictions$mars=list()
    posterior.predictions$rf=list()
    posterior.predictions$gpsg=list()
    
    #ML holdout predictions
    posterior.predictions$brnn = predict(ensemble.of.machine.learning.models[[w]]$brnn.model$brnn.model, newdata = test.data[,c(1:(2*D))])
    posterior.predictions$mars = predict(ensemble.of.machine.learning.models[[w]]$mars.model$mars.model, newdata = test.data[,c(1:(2*D))])
    posterior.predictions$rf = predict(ensemble.of.machine.learning.models[[w]]$randomforest.model$rf.model, newdata = test.data[,c(1:(2*D))])
    
    #scale
    Y.brnn.test.scaled = scale(as.data.frame(posterior.predictions$brnn), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.brnn"], 
                               scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.brnn"])
    Y.mars.test.scaled = scale(as.data.frame(posterior.predictions$mars), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.mars"], 
                               scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.mars"])
    Y.rf.test.scaled = scale(as.data.frame(posterior.predictions$rf), center =collection.of.stacking.weights[[w]]$scaling$center[,"Y.rf"], 
                             scale =collection.of.stacking.weights[[w]]$scaling$scale[,"Y.rf"])
    
    #stacked prediction
    posterior.predictions$gpsg= predict(collection.of.stacking.weights[[w]]$model, x=cbind(Y.brnn.test.scaled, Y.mars.test.scaled,Y.rf.test.scaled))
    
    
    posterior.predictions    		
  }
  parallel::stopCluster(cl)
  plan(sequential)
  
  
  cat("Done!\n")
  
  lcb = numeric(N.test)
  sds = numeric(N.test)
  
  weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
  #Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)
  for (w in 1:W) {lcb <- lcb + weights.w[w]*exp(gpsg.posterior.predictives[[w]]$gpsg$mean)}
  for (w in 1:W) {sds <- sds + weights.w[w]*exp(gpsg.posterior.predictives[[w]]$gpsg$sd2+gpsg.posterior.predictives[[w]]$gpsg$nugs)}
  
  lcb.point <- test.points[which.min(lcb),] #TR: Min!
  lcb.index <- which.min(lcb) #TR: Min!
  lcb.sd <- sds[lcb.index]
  
  evaluations$best[[length(evaluations$best)+1]] <- list()
  evaluations$best[[length(evaluations$best)]]$predicted =list()
  evaluations$best[[length(evaluations$best)]]$predicted$parameters <- lcb.point+0.5
  
  evaluations$best[[length(evaluations$best)]]$known =list()
  #evaluate best known
  minimum.known = exp(numeric(dim(master.data.small.train)[1]))
  for (w in 1:W) {minimum.known <- minimum.known + weights.w[w]*exp(master.data.small.train[,(2*D+w)])}
  
  evaluations$best[[length(evaluations$best)]]$known$parameters <- as.numeric(master.data.small.train[which.min(minimum.known),c(1:D)]+0.5)
  
  
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
      time[[(main.iterations-1)*N.Updates+update]]$params_min_LF_known = min(minimum.known, na.rm=T)
      
    }
    
  }
  save(evaluations, file="eval.state.newer")
  bp=quant.to.prior(evaluations$best[[length(evaluations$best)]]$known$parameters) #best point estimate
  
  save(bp, file="best.params")
  
  main.end2 = toc()
  total.time = total.time + (main.end2$toc - main.end2$tic)
  
} 

save(ensemble.of.machine.learning.models, file="ensemble.of.machine.learning.models")
save(collection.of.stacking.weights, file="collection.of.stacking.weights")
