### R script to perform Bayesian optimisation on an openMalaria-style joint objective function using the "multi-fidelity, stacked generalisation with GPs" algorithm
#		Bayesian optimisation to select the best combination of parameters. To fit ot data. Iterative process, sequential design stratefy for global optimisation
#			# We want to minimise the loss function
#		"multi-fidelity, stacked generalisation with GPs" algorithm = some alorithm for data sampling and optimisation I think

## Noisy objective function, Obj.large(theta) of W sub-objectives in a D dimensional parameter space given by the D-hypercube with probability Fail.large[w] of simulation failure in the w-th objective and model cost Cost.large [to evaluate all sub-objectives simultaneously]
#		W objectives: I think we have 10 or 11 as outlined in Smith et al 2008 - I think this is what the objective function(s?) is (are?) based on
#		D dimensional parameter space

############################################################
#
# TO BEGIN
#
############################################################
# need to install all packages, see below

############################################################
#
# global variables
#
############################################################
User = "Theresa.scicore.Sep"  #"Theresa.local.Rosenbrock.test"  #choose from c("Theresa.local.Rosenbrock.test",Theresa.local", Melissa", "Theresa.scicore.test","Theresa.scicore.full","Theresa.scicore.molineaux","Theresa.scicore.full.prev.inc","Theresa.scicore.molineaux.prev.inc","Ewan") for filepaths


runInitial_gp = TRUE
Continued_run = FALSE
CONT = FALSE
STEP1 = FALSE
STEP2 = TRUE
CONT1 = FALSE
CONT2 = TRUE
SCICORE = TRUE
INIT  = TRUE
options(scipen=999)

############################################################
# D is number of parameters (23 for OM and 2 or more for Rosenbrook)
# W is the number of objectivs (11 for OM, 1 or more for RB)

# Melissa: For testing
D = 23 #no less than 2 for RB
lockBinding("D", globalenv())
W = 11 #1-2 to run on laptop
lockBinding("W", globalenv())

############################################################
#
# Load packages
#
############################################################
if(!require(pacman)){install.packages("pacman")}; require(pacman) #pacman::p_load(h2o)
pacman::p_load(h2o, Matrix, MASS, tmvtnorm, gridExtra, doParallel) # TODO ,doParallel)
if(!require(rmngb)){
  url <- "https://cran.r-project.org/src/contrib/Archive/rmngb/rmngb_0.6-1.tar.gz"
  pkgFile <- "rmngb_0.6-1.tar.gz"
  download.file(url = url, destfile = pkgFile)
  # Install package
  install.packages(pkgs=pkgFile, type="source", repos=NULL)
  # Delete package tarball
  unlink(pkgFile)
}; require(rmngb)

#cl <- makeCluster(W)
#registerDoParallel(cl)

############################################################
#
# set local directory
#
############################################################

if(User=="Theresa.scicore.Sep"){
  local.directory = '/scicore/home/smith/reiker/Paper_3_Model_Fitting/2019_09_16_Complete_redo/'
  om.path = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2019_09_16_Complete_redo/openmalaria/"
  u = "reiker"
}

setwd(local.directory)
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
## Load functions for running individual simulation models: Eval.model.large(theta) and Eval.model.small(theta): must return a list with $Objectives (a W dimensional vector of sub-objectives, failed objectives being NAs) and $Runtime (a single number giving the computational runtime of the model evaluation)
## By default the line below loads the toy Rosenbrock function model from the RosenMod.R script supplied: this example requires D>=2 and imagines the problem space to have structure in only Deff effective dimensions (specified in RosenMod.R)

#source("RosenMod.R")
source("./script/2019_09_EvalOM.parallel.R", echo = TRUE)

############################################################
#
# Initialise problem
#
############################################################
### Initialise problem
## Initial design points may be selected deterministically or randomly; a minimum of ~50 (better = 500) initial function evaluations are required to ensure stable behaviour of this algorithm in the initial stages of exploration; likewise draws must be split between small and large versions of model with a split no worse than ~1:10
## Since an estimate of the iid noise term on the small and large models is required one might choose to incorporate this into the initial design as below  
## The function evaluations should be stored in a list 'evaluations' with members $small and $large
## The example code below implements a simple uniform sampling design with initial draws at the prior median to estimate iid noise terms 


if(runInitial_gp) { ## TR: Above required regardless
# Melissa: For testing
#test values
Nrep.noise.small <- 200 #25  #200 = approx. 800 min = 13h20min
Nrep.noise.large <- 10 #10	#20 = apporx. 480 min = 8h
#total time to evaluate test values = approx 21h20min
if(STEP1){
evaluations <- list()
evaluations$small <- list()
evaluations$large <- list()
INIT = TRUE
if(CONT1){
	load("eval.state.noise.logLH")
	start.rep.small = length(evaluations$small) + 1
	start.rep.large = length(evaluations$large) + 1
} else{
	start.rep.small = 1
	start.rep.large = 1
}
if(start.rep.small <= Nrep.noise.small){
for (i in seq(start.rep.small,Nrep.noise.small,2)) {
  evaluations$small[[i]] <- list()
  evaluations$small[[i+1]] <- evaluations$small[[i]] 
  evaluations$small[[i]]$parameters <- rep(0.5,D)
  evaluations$small[[i+1]]$parameters <- evaluations$small[[i]]$parameters}


for (i in 1:Nrep.noise.small) {  
  current.index=paste0("init_",sprintf("%04i",i))
  write.xmls.small(evaluations$small[[i]]$parameters)} # write times 61

Eval.model.small() # run everything
observations  = list()
for (i in 1:Nrep.noise.small) {  
  Runtime = NA #Remove later
  current.index=paste0("init_",sprintf("%04i",i))
  observations[[i]] = obj.calc.small(evaluations$small[[i]]$parameters)  # LF calculations times 61
  evaluations$small[[i]]$Objectives <- observations[[i]]$Objectives
  } 
  
  setwd(local.directory)
  save(evaluations,file="eval.state.noise.logLH") 
 }
#OK until here
 
if(start.rep.large <= Nrep.noise.large){
for (i in seq(start.rep.large,Nrep.noise.large,2)) {
  evaluations$large[[i]] <- list()
  evaluations$large[[i+1]] <- evaluations$large[[i]]  
  evaluations$large[[i]]$parameters <- rep(0.5,D)
  evaluations$large[[i+1]]$parameters <- evaluations$large[[i]]$parameters}

for (i in 1:Nrep.noise.large) {  
  current.index=paste0("init_",sprintf("%04i",i))
  write.xmls.large(evaluations$large[[i]]$parameters)} 
  
Eval.model.large()
observations  = list()
for (i in 1:Nrep.noise.large) {  
  Runtime = NA #Remove later
  current.index=paste0("init_",sprintf("%04i",i))
  observations[[i]] = obj.calc.large(evaluations$large[[i]]$parameters)  # LF calculations times 61
  evaluations$large[[i]]$Objectives <- observations[[i]]$Objectives
  } 
  setwd(local.directory)
  save(evaluations,file="eval.state.noise.logLH") ## <<added 20190218 MP>>##
}
} else if (Continued_run) {load("eval.state.newer") # TR added for additional simulations to previous runs 13.3.2019
} else {load("eval.state.noise.logLH")} ## <<added 20190218 MP>>##

#ALL GOOD
############################################################
#
# estimate noise  - CHECKPOINT 1
#
############################################################
## Estimate iid noise terms from repetition design points using hierarhical model
#small.repetitions <- matrix(-999999,nrow=Nrep.noise.small,ncol=W)
#for (w in 1:W) {
#  for (i in 1:Nrep.noise.small) {small.repetitions[i,w] <- evaluations$small[[i]]$Objectives[w]}
#}
#small.sum.of.squares <- colVars(small.repetitions,na.rm=T)
#small.n.obs <- colSums(!is.na(small.repetitions))
#large.repetitions <- matrix(-999999,nrow=Nrep.noise.large,ncol=W)
#for (w in 1:W) {
#  for (i in 1:Nrep.noise.large) {large.repetitions[i,w] <- evaluations$large[[i]]$Objectives[w]}
#}
#large.sum.of.squares <-colVars(large.repetitions,na.rm=T)
#large.n.obs <- colSums(!is.na(large.repetitions))

#if (length(which(small.n.obs<2))!=0) {stop("Too many failed init sims for iid noise estimation!")}
#if (length(which(large.n.obs<2))!=0) {stop("Too many failed init sims for iid noise estimation!")}
###########################################################################
# COMPILE
###########################################################################
#compile("./script/noise_est.cpp") # TMB code for hierarchical noise estimation model
#dyn.load(dynlib("./script/noise_est"))

###########################################################################
###########################################################################
#obj <- MakeADFun(
#  data = list('W'=W,
#              'Nobs_small'=small.n.obs,
#              'Nobs_large'=large.n.obs,
#              'SS_small'=small.sum.of.squares,
#              'SS_large'=large.sum.of.squares
#  ),
#  
#  parameters = list(prior_mean_log_iid_noise_small=0,
#                    prior_sd_log_iid_noise_small=0,
#                    prior_mean_logit_iid_noise_factor_large=0,
#                    prior_sd_logit_iid_noise_factor_large=0,
#                    log_iid_noise_small=rep(0,W),
#                    logit_iid_noise_factor_large=rep(0,W)
#  ),
#  DLL = "noise_est",
#  random=c('log_iid_noise_small','logit_iid_noise_factor_large')
#)
#opt <- nlminb(obj$par,obj$fn,obj$gr,max.iter=2000)
#rep <- sdreport(obj)
#iid.noise.ests <- obj$report(c(rep$par.fixed,rep$par.random))

#if(runInitial_gp) { ## <<added 20190218 MP>>##
if(STEP2) { ## <<added 20190218 MP>>## ; ##TR: changed condition 20190320
evaluations <- list()
evaluations$small <- list()
evaluations$large <- list()

Ninit.small <- 8000#75 # 8000 * 5min = 40,000 min; 450 in parallel => 90 min
Ninit.large <- 100#20 # Probably between 10 and 50; 20 * 24min = 480 min = 8h; thus 13h for initial design
if(CONT2){
	load("eval.state.init")
	start.it.small =length(evaluations$small) + 1
	start.it.large = length(evaluations$large) + 1
} else{
	start.it.small = 1
	start.it.large = 1
}





if(start.it.small <= Ninit.small){
for (i in seq(start.it.small,Ninit.small,2)) {
  evaluations$small[[i]] <- list()
  evaluations$small[[i+1]] <- evaluations$small[[i]] 
  evaluations$small[[i]]$parameters <-  runif(D,0,1)
  evaluations$small[[i+1]]$parameters <- evaluations$small[[i]]$parameters} 
  #start.writing.small = as.numeric(system("ls ./PriorDraws/scenarios/scenarios_complete/ | wc -l", intern = TRUE))
  
  N = detectCores() -1
  cl = makeCluster(N)
  registerDoParallel(cl)
  foreach (i=start.writing.small:Ninit.small) %dopar% {  
	current.index=paste0("init_",sprintf("%04i",i))
	write.xmls.small(evaluations$small[[i]]$parameters)} # write times 61
  
  Eval.model.small() # run everything
  } #REMOVE THIS!!!!
  observations  = list()
  
  for (i in 1:Ninit.small) {  
	Runtime = NA #Remove later
	current.index=paste0("init_",sprintf("%04i",i))
	size<-"small"
	
	if (as.numeric(system(paste0('ls ',local.directory,'/PriorDraws/Output/',size,'/output_',current.index,' | wc -l'),intern=TRUE))<61){ #if the number of output files is less than 61, some runs have crashed and the location is invalid
		cat("Error: Bad location at index",current.index,". Moving on.\n ", sep="")
		evaluations$small[[i]]$Objectives <- rep(NA,W)
	}else{
	observations[[i]] = obj.calc.small(evaluations$small[[i]]$parameters)  # LF calculations times 61
	evaluations$small[[i]]$Objectives <- observations[[i]]$Objectives}
  } 
  
  setwd(local.directory)
  save(evaluations,file="eval.state.init2")
#} #PUT THIS BACL IN!!!!


if(start.it.large<=Ninit.large){
	for (i in seq(start.it.large,Ninit.large,2)) {
	evaluations$large[[i]] <- list()
	evaluations$large[[i+1]] <- evaluations$large[[i]] 
	evaluations$large[[i]]$parameters <-  runif(D,0,1)
	evaluations$large[[i+1]]$parameters <- evaluations$large[[i]]$parameters} 
	save(evaluations,file = "eval.state.init3")
	N = detectCores() -1
	cl = makeCluster(N)
	registerDoParallel(cl)
	
	start.writing.large = 1 # can be changed for testing or continued writing
  foreach (i=start.writing.large:Ninit.large) %dopar% {  
    #size ="large"
	current.index=paste0("init_",sprintf("%04i",i))
	write.xmls.large(evaluations$large[[i]]$parameters)} # write times 61

  Eval.model.large() # run everything
  observations  = list()
  
   
  for (i in 1:Ninit.large) {  
    size <- "large"
	Runtime = NA #Remove later
	current.index=paste0("init_",sprintf("%04i",i))
	
	if (as.numeric(system(paste0('ls ',local.directory,'/PriorDraws/Output/',size,'/output_',current.index,' | wc -l'),intern=TRUE))<61){ #if the number of output files is less than 61, some runs have crashed and the location is invalid
		cat("Error: Bad location at index",current.index,". Moving on.\n ", sep="")
		evaluations$large[[i]]$Objectives <- rep(NA,W)
	}else{
	observations[[i]] = obj.calc.large(evaluations$large[[i]]$parameters)  # LF calculations times 61
	evaluations$large[[i]]$Objectives <- observations[[i]]$Objectives}
  } 
  
  setwd(local.directory)
    save(evaluations,file="eval.state.init2")

	}
}
#} else {load("eval.state.om.newer")}  ## <<added 20190218 MP>>##
} else {load("eval.state.init")} ## TR: changed file name 20190320
}
 
############################################################
#
#  CHECKPOINT 2
#
############################################################

# Initialise proposal distributions for test particles
if(runInitial_gp) { ## << MP added >>
	proposal.means <- rep(0.5,D) #rep(0,D)
	proposal.covs <- diag(1,D)
}else{
	valTest = 50 #213 ## <<added 20190218 MP>>##
	proposal.means <- evaluations$small[[valTest]]$proposal.means #rep(0,D)
	proposal.covs <- diag(evaluations$small[[valTest]]$proposal.sds^2) #diag(1,D)
}
