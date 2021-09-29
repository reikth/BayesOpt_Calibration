# ===========================================================================================
# phaseApred.r
#
# Takes output from 61 scenarios and plots phase A predictions to a pdf file
# initial plotting from PhaseA by AR: 6.7.07
# New file (with errors fixed) MP: 14.05.2012 very very substantial changes to structure, data, 
#						obtaining predictions, output, naming conventions, plots, fix errors, 
#						include LF calculations
#
# Requrie : outputXXX.txt files (where XXX is the scenario number)
#			LF_#  - empty file where # is the current LF value
#			param_##.txt - file with parameter values and ## indicates param_id
#			fieldData.txt - field data file
# ===========================================================================================
#set this to your working dir. 
library(RMySQL,dplyr)
setwd("/scicore/home/smith/reiker/Paper_3_Model_Fitting/")


#Script for generating convergence plots 

u = "reiker" #Scicore user name
GitDir = "/scicore/home/smith/reiker/GitRepos/om_fitting/" # Local version of the Git repository
ExperimentDir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/" #Parent folder to experiments
outputDir = "/PriorDraws/Output/small/" #relative to the parent folder, followed by ID
xmlDir = "/PriorDraws/scenarios/scenarios_complete/small/"#relative to the parent folder, followed by ID
  
ExperimentNameHetGP = "2020_05_full_fit_small_hetGP_norm_sampling"
ExperimentNameGPSG = "2020_05_full_fit_small_GPSG_norm_sampling"


############################################################
#
# Load Packages
#
############################################################


if(!require(pacman)){install.packages("pacman")}; require(pacman) 
pacman::p_load(mlr,sensitivity,ggplot2,dplyr,tidyverse,hetGP,doParallel,viridis,lattice,colorspace,tgp,rlist,latex2exp)
############################################################
#
# set up folders
#
############################################################

#set local.directory
local.directory = paste0(ExperimentDir,ExperimentNameHetGP,'/') #working directory

###################

# ===========================================================================================
# CODES/IDS FIELD DATA AND MODEL RESULTS codes/ids for for field data and model predictions
# ===========================================================================================
# function fieldData_nameCodes -  In the field data, the codes used for the different measures are as follows:
fieldData_nameCodes <- function(){
  fieldDataCode <- list()
  fieldDataCode$nHost <- 0;  #number of hosts
  fieldDataCode$nInfect <- 1;  #number of infected hosts
  fieldDataCode$nExpectd <- 2;  #expected number of infected hosts
  fieldDataCode$npatent <- 3;  #number of patent hosts
  fieldDataCode$sumlogdens <- 5;  #sum of the logarithm of the density
  fieldDataCode$pyr <- 10;  #person-years
  fieldDataCode$totalpatentInf <- 13;  #total patent infections
  fieldDataCode$malMort <- 16;  #Malaria mortality rate
  fieldDataCode$recInc <- 18;  #recorded incidence disease
  fieldDataCode$pyrogenThres <- 19;  #pyrogenic threshold as para/leuc ratio
  fieldDataCode$rate <- 20;  #relative risk of severe compared to agegp 1
  fieldDataCode$prop <- 21;  #proportion (unspecified)
  fieldDataCode$allCauseIMR <- 22;  #all-cause infant mortality rate
  return(fieldDataCode)
}


# function modelResults_nameCodes -  for model results from openMalaria, the codes used for the different measures are as follows:
modelResults_nameCodes <- function(){
  modelResultsCode <- list()		
  modelResultsCode$nHost <- 0;   #Total number of humans
  modelResultsCode$nInfect <- 1;   #number of infected hosts
  modelResultsCode$nExpectd <- 2;   #expected number of infected hosts
  modelResultsCode$nPatent <- 3;   #number of patent hosts
  modelResultsCode$sumLogPyrogenThres <- 4;   #Sum of the log of the pyrogen threshold
  modelResultsCode$sumlogDens <- 5;   #Sum of the logarithm of the parasite density
  modelResultsCode$totalInfs <- 6;   #Total infections
  modelResultsCode$nTransmit <- 7;   # Infectiousness of human population to mosquitoes: sum(p(transmit_i)) across humans i, weighted by availability to mosquitoes. Single value, not per age-group.
  modelResultsCode$totalPatentInf <- 8;   #Total patent infections
  modelResultsCode$sumPyrogenThresh <- 10;   #Sum of the pyrogenic threshold
  modelResultsCode$nTreatments1 <- 11;   #number of treatments (1st line) (added to 1-day model in 24.1)
  modelResultsCode$nTreatments2 <- 12;   #number of treatments (2nd line) (added to 1-day model in 24.1)
  modelResultsCode$nTreatments3 <- 13;   #number of treatments (inpatient) (added to 1-day model in 24.1)
  modelResultsCode$nUncomp <- 14;   #number of episodes (uncomplicated)
  modelResultsCode$nSevere <- 15;   #number of episodes (severe)
  modelResultsCode$nSeq <- 16;   #recovered cases with sequelae
  modelResultsCode$nHospitalDeaths <- 17;   #deaths in hospital
  modelResultsCode$nIndDeaths <- 18;   #number of deaths (indirect)
  modelResultsCode$nDirDeaths <- 19;   #number of deaths (direct)
  modelResultsCode$nEPIVaccinations <- 20;   #number of EPI vaccine doses given
  modelResultsCode$allCauseIMR <- 21;   #all cause infant mortality rate (returned as a single number over whole intervention period, instead of from a survey interval)
  modelResultsCode$nMassVaccinations <- 22;   #number of Mass / Campaign vaccine doses given
  modelResultsCode$nHospitalRecovs <- 23;   #recoveries in hospital without sequelae
  modelResultsCode$nHospitalSeqs <- 24;   #recoveries in hospital with sequelae
  modelResultsCode$nIPTDoses <- 25;   #number of IPT Doses
  modelResultsCode$annAvgK <- 26;   #Annual Average Kappa. Calculated once a year as sum of human infectiousness weighted by initial EIR for that time of year.
  modelResultsCode$nNMFever <- 27;   #Number of episodes of non-malaria fever
  modelResultsCode$expectedDirectDeaths <- 74;   #Expected number of deaths (direct)
  modelResultsCode$expectedHospitalDeaths <- 75;   #Expected number of deaths in hospital
  modelResultsCode$expectedIndirectDeaths <- 76;   #Expected number of deaths (indirect)
  modelResultsCode$expectedSevere <- 78;   #Expected number of episodes (severe)
  return(modelResultsCode)
}

# 
# ===========================================================================================
# GET DENS BIAS VALUES FROM FILE
# ===========================================================================================

getDensBiasVals <- function(readDensBiasFromFile, runDir, param_file){
  if (readDensBiasFromFile ==0){
    # values from base model
    densBias <- list()
    densBias$garki <- 4.8
    densBias$nonGarki <- 0.18
  }
  if (readDensBiasFromFile ==1){
    densBias <- list()	
    densBias$garki <- system(paste0('grep "Density bias (Garki)" ',param_file, '| grep -Eo "[0-9]+\\.[0-9]+" '),intern=T) %>% as.numeric() #extracts decimal numbers from line that contains specified string
    densBias$nonGarki <- system(paste0('grep "Density bias (non Garki)" ',param_file, '| grep -Eo "[0-9]+\\.[0-9]+" '),intern=T) %>% as.numeric() 
  }
  return(densBias)
}


# ===========================================================================================
# OBECTIVES  - functions to get field data and model predictions and 
#				produce equivalent simulated data
# ===========================================================================================

#---------------------------------------------------------------------------------------------
# Age Pattern of Incidence after intervention
#---------------------------------------------------------------------------------------------
OBJ_AgePatternIncidenceAfterIntervention<-function(modelResultID, fieldDataAll, fieldDataID, densBias, oDirAbs, runID,scenarioNums){
  # scenarios 30
  
  # read in predictions for particular scenario number and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    # new field data
    # only need outcome 3 (nPatent) not 5 (densities)
    fieldNew[[scenario]] <- outputTemp[[scenario]][outputTemp[[1]]$V3==modelResultID$nHost | outputTemp[[1]]$V3==modelResultID$nPatent,]
    
    # check dimensions
    if (dim(fieldNew[[scenario]])[1]!= dim(fieldTemp[[scenario]])[1]){
      stop(paste("   ERROR in obj 1 for scenario", scenarioNums[scenario], sep=" "))
    }
  }
  return(fieldNew)
}

#---------------------------------------------------------------------------------------------
# Age Pattern of Prevalence And Parasite Density
#---------------------------------------------------------------------------------------------
OBJ_AgePatternPrevalenceAndDensity<-function(modelResultID, fieldDataAll, fieldDataID, densBias, oDirAbs, runID,scenarioNums){
  # scenarios 24,28,29,35,34 and 31
  
  # density bias vector, one for each scenario
  densBiasVec<-c(densBias$garki, densBias$garki, densBias$garki, densBias$nonGarki, densBias$nonGarki, densBias$nonGarki)
  
  # read in predictions for particular scenario number and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    
    # new field data
    fieldNew[[scenario]]<- outputTemp[[scenario]]
    # replace field sumlogDens with transformation based on loss function
    fieldNew[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$sumlogDens] <-outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$sumlogDens] - outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nPatent]*log(densBiasVec[scenario])
    # check dimensions
    if (scenario == 5){# ie 34
      if (dim(fieldNew[[scenario]][fieldNew[[scenario]]$V3!=modelResultID$totalInfs & fieldNew[[scenario]]$V3!=modelResultID$totalPatentInf,])[1] != dim(fieldTemp[[scenario]][fieldTemp[[scenario]]$V4!=fieldDataID$totalpatentInf,])[1]){
        stop(paste("   ERROR in obj 2&3 for scenario", scenarioNums[scenario], sep=" "))
      }
    }else{
      if (dim(fieldNew[[scenario]])[1]!= dim(fieldTemp[[scenario]])[1]){
        stop(paste("   ERROR in obj 2&3 for scenario", scenarioNums[scenario], sep=" "))
      }
    }
    # 
  }
  return(fieldNew)
  
}

#---------------------------------------------------------------------------------------------
# Age Pattern of Number of concurrent infections
#---------------------------------------------------------------------------------------------

OBJ_AgePatternMOI <- function(modelResultID, fieldDataAll, fieldDataID, oDirAbs, runID, scenarioNums){	
  # scenario 34
  
  # read in predictions for particular scenario number and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    
    # new field data		
    fieldNew[[scenario]]<- outputTemp[[scenario]]
    
    # check dimensions
    if (dim(fieldNew[[scenario]][fieldNew[[scenario]]$V3==modelResultID$totalPatentInf ,])[1] != dim(fieldTemp[[scenario]][fieldTemp[[scenario]]$V4==fieldDataID$totalpatentInf,])[1]){
      stop(paste("   ERROR in obj 4 for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
}

#---------------------------------------------------------------------------------------------
# Age Pattern of Incidence of Clinical Malaria in Senegal
#---------------------------------------------------------------------------------------------
OBJ_AgePatternIncidenceClinicalMalaria_Senegal <- function(modelResultID, fieldDataAll, fieldDataID,  oDirAbs, runID,  scenarioNums){
  # field data (scenarios 232 and 233)
  
  # The rate multiplier is the duration in years for which episodes are collected
  # which is 5 years
  rateMultiplier<-5
  
  # read in predictions for particular scenario number and transform if needed for new field data and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    
    # new field data		
    fieldNew[[scenario]]<- outputTemp[[scenario]]
    # transform
    # replace field new (but nHost, replace number 18)
    fieldNew[[scenario]]$V3[outputTemp[[scenario]]$V3==modelResultID$nHost] <- 18
    #For measuresd severe choose:
    #fieldNew[[scenario]]$V4[fieldNew[[scenario]]$V3==18] <-(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nUncomp] + outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nSevere])/(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nHost])/rateMultiplier
    #For expected:
    fieldNew[[scenario]]$V4[fieldNew[[scenario]]$V3==18] <-(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nUncomp] + outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$expectedSevere])/(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nHost])/rateMultiplier
    
    
    # check dimensions
    if (dim(fieldNew[[scenario]][fieldNew[[scenario]]$V3==18,])[1] != dim(fieldTemp[[scenario]][fieldTemp[[scenario]]$V4==fieldDataID$recInc,])[1]){
      stop(paste("   ERROR in obj 5a for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
}

#---------------------------------------------------------------------------------------------
# Age Pattern of Incidence of Clinical Malaria in Idete 
#---------------------------------------------------------------------------------------------

OBJ_AgePatternIncidenceClinicalMalaria_Idete <- function(modelResultID, fieldDataAll, fieldDataID,  oDirAbs, runID,  scenarioNums){
  # (scenario 49)
  
  # The rate multiplier is 1/access, access = 36%
  accessConst<-2.7975236259999985
  rateMultiplier<-accessConst
  
  # read in predictions for particular scenario number and transform if needed for new field data and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    
    # new field data		
    fieldNew[[scenario]]<- outputTemp[[scenario]]
    # transform
    # replace field new (but nHost, replace number 18)
    fieldNew[[scenario]]$V3[outputTemp[[scenario]]$V3==modelResultID$nHost] <- 18
    fieldNew[[scenario]]$V4[fieldNew[[scenario]]$V3==18] <-(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nUncomp])/(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nHost])/rateMultiplier
    
    # check dimensions
    if (dim(fieldNew[[scenario]][fieldNew[[scenario]]$V3==18,])[1] != dim(fieldTemp[[scenario]][fieldTemp[[scenario]]$V4==fieldDataID$recInc,])[1]){
      stop(paste("   ERROR in obj 5b for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
  
}

#---------------------------------------------------------------------------------------------
# Age Pattern of parasite density threshold for clinical attack
#---------------------------------------------------------------------------------------------
OBJ_AgePatternThresholdClinicalAttack <- function(modelResultID, fieldDataAll, fieldDataID, densBias,  oDirAbs, runID,  scenarioNums){
  # (scenario 234)
  
  # The rate multiplier associated with log parasite/leucocyt ratio and thus μ = 1/(8000densbias) is the non-garki density bias.
  biasPy <-1/(8000*densBias$nonGarki) 
  
  # read in predictions for particular scenario number and transform if needed for new field data and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    
    # new field data		
    # survey one only
    fieldNew[[scenario]]<- outputTemp[[scenario]][outputTemp[[scenario]]$V1==1&outputTemp[[scenario]]$V3==0,]
    # transform
    # replace field new (but nHost, replace number 19)# survey one only!!!
    fieldNew[[scenario]]$V3[fieldNew[[scenario]]$V3==modelResultID$nHost] <- 19	
    fieldNew[[scenario]]$V4[fieldNew[[scenario]]$V3==19] <- exp( log(biasPy) + outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$sumLogPyrogenThres & outputTemp[[scenario]]$V1==1]/outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3==modelResultID$nHost& outputTemp[[scenario]]$V1==1])
    
    # check dimensions
    if (dim(fieldNew[[scenario]])[1] != dim(fieldTemp[[scenario]])[1]){
      stop(paste("   ERROR in obj 6 for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
}

#---------------------------------------------------------------------------------------------
# Severe episodes vs prevalence
#---------------------------------------------------------------------------------------------
OBJ_SevereEpisodesVsPrevalence <- function(modelResultID, fieldDataAll, fieldDataID,  oDirAbs, runID,  scenarioNums){
  # field data (scenarios 501,502,503,504,505,506,507,508,509,510,511,512,514,515,516,517,518,519,520,521,522,523,524,525,526,527)
  
  AccessForSevereCases <- 0.48
  # read in predictions for particular scenario number and transform if needed for new field data and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  output <- list()
  output$Prev <- array(NA,dim=c(length(scenarioNums)))
  output$SevRate <- array(NA,dim=c(length(scenarioNums)))
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]

    outputTempTemp <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    # first age-group only (0-9yrs) 
    outputTemp[[scenario]] <- outputTempTemp[outputTempTemp$V2==1,]
    rm(outputTempTemp)
    pred_nhost<-mean(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3== modelResultID$nHost]) # mean num hosts
    pred_npatent<-mean(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3== modelResultID$nPatent]) # mean num patent
    pred_prevalence<-(pred_npatent/pred_nhost) 
    # multiply by 100 to get to percentage
    pred_severe<-sum(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3== modelResultID$expectedSevere]) 
    #pred_severe<-sum(outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V3== modelResultID$nSevere]) 
    #require rate episodes per 1000 person years (so divide pred$severe/2)
    pred_sevEpPer1000personYr<-AccessForSevereCases*(pred_severe/(2*pred_nhost))*1000 
    
    
    output$Prev[scenario] <- pred_prevalence
    output$SevRate[scenario] <-pred_sevEpPer1000personYr
    rm(pred_nhost, pred_npatent, pred_prevalence, pred_sevEpPer1000personYr)
  }
  # use predicted to create relations ship between prevalence and severe rate
  trialDat <- array(NA, dim = c(length(scenarioNums),2))
  trialDat[,1] <- output$Prev
  trialDat[,2] <- output$SevRate
  orderedResults <- trialDat[order(trialDat[,1]),]
  rm(trialDat)
  indicesToInclude <- round(seq(1,26,by=26/12))
  indicesToInclude[12] <- 26
  predictedToUse_TODOnotGreatUseSpline <- orderedResults[indicesToInclude,]
  
  
  # so no errors if exact fitting
  # maybe better is to fit a curve
  
  smoothsplineFitResults<- smooth.spline(orderedResults[,1], orderedResults[,2])
  
  xSeq <- c(0,seq(min(output$Prev)*0.9,max(output$Prev)*1.1,by=(max(output$Prev)*1.1-min(output$Prev)*0.9)/9),1)
  predictedToUse <- predict(smoothsplineFitResults, x=xSeq)
  # plot(orderedResults[,1], orderedResults[,2], xlim = c(0,1), ylim = c(0,max(predictedToUse$y)))
  # points(predictedToUse$x, predictedToUse$y,  pch=19)
  
  predictedToUse
  # , xmin = 0, xmax =0.45)
  
  
  for (scenario in 1:length(scenarioNums)) {
    fieldNew[[scenario]] <- fieldTemp[[scenario]]
    # fieldNew[[scenario]]$V5[fieldNew[[scenario]]$V4 == fieldDataID$rate] <- predictedToUse[,2]
    # fieldNew[[scenario]]$V5[fieldNew[[scenario]]$V4 == fieldDataID$prop] <- predictedToUse[,1]
    fieldNew[[scenario]]$V5[fieldNew[[scenario]]$V4 == fieldDataID$rate] <- predictedToUse$y
    fieldNew[[scenario]]$V5[fieldNew[[scenario]]$V4 == fieldDataID$prop] <- predictedToUse$x
    # check dimensions
    if (dim(fieldNew[[scenario]])[1] != dim(fieldTemp[[scenario]])[1]){
      stop(paste("   ERROR in obj 6 for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
}

#---------------------------------------------------------------------------------------------
# Age Pattern of Severe hospitalisation (0-9 years)
#---------------------------------------------------------------------------------------------
OBJ_AgePatternOfSevere<- function(modelResultID, fieldDataAll, fieldDataID,  oDirAbs, runID,  scenarioNums){
  # (scenarios 158,167,173,176)
  
  # read in predictions for particular scenario number and transform if needed for new field data and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    
    # new field data
    # field[1] = nHost[age1]/ nHost[age2]* nSevere[age2]/ nSevere[age1]
    
    # measure(20)[age2] = (nSevere[age2]*1000/(2*nHost[age2]))/(nSevere[age1]*1000/(2*nHost[age1]))
    # measure(20)[age3] = (nSevere[age3]*1000/(2*nHost[age3]))/(nSevere[age1]*1000/(2*nHost[age1]))
    
    
    # require only one survey and two ages (2,3) and measure 20
    fieldNew[[scenario]]<- outputTemp[[scenario]][outputTemp[[scenario]]$V1==1&outputTemp[[scenario]]$V3==0&(outputTemp[[scenario]]$V2==2|outputTemp[[scenario]]$V2==3),]
    # replace field new (but nHost, replace number 20)
    fieldNew[[scenario]]$V3[fieldNew[[scenario]]$V3==modelResultID$nHost] <- 20
    # transform ##TR 2020: Note change to expected severe
    fieldNew[[scenario]]$V4[fieldNew[[scenario]]$V3==20 & fieldNew[[scenario]]$V2==2] <- (outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==2 & outputTemp[[scenario]]$V3==modelResultID$expectedSevere]*1000/(2*outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==2 & outputTemp[[scenario]]$V3==modelResultID$nHost])) / (outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==1 & outputTemp[[scenario]]$V3==modelResultID$expectedSevere]*1000/(2*outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==1 & outputTemp[[scenario]]$V3==modelResultID$nHost]))
    fieldNew[[scenario]]$V4[fieldNew[[scenario]]$V3==20 & fieldNew[[scenario]]$V2==3] <- (outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==3 & outputTemp[[scenario]]$V3==modelResultID$expectedSevere]*1000/(2*outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==3 & outputTemp[[scenario]]$V3==modelResultID$nHost])) / (outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==1 & outputTemp[[scenario]]$V3==modelResultID$expectedSevere]*1000/(2*outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2==1 & outputTemp[[scenario]]$V3==modelResultID$nHost]))
    
    # check dimensions
    if (dim(fieldNew[[scenario]])[1] != dim(fieldTemp[[scenario]])[1]){
      stop(paste("   ERROR in obj 8 for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
}

#---------------------------------------------------------------------------------------------
# Direct Malaria Mortality
#---------------------------------------------------------------------------------------------
OBJ_DirectMalariaMortality <- function(modelResultID, fieldDataAll, fieldDataID,  oDirAbs, runID,  scenarioNums){
  # (scenarios 301,302,303,312,316,317,318,326,327)
  
  # read in predictions for particular scenario number and transform if needed for new field data and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    # new field data
    predTempnhost<-outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2== 1 & outputTemp[[scenario]]$V3==modelResultID$nHost]
    predTempdirDeaths <-outputTemp[[scenario]]$V4[outputTemp[[scenario]]$V2== 1 &outputTemp[[scenario]]$V3==modelResultID$expectedDirectDeaths]
    # require direct deaths episodes per 1000 person years, so multiply by 1000 and divide by 2
    predTempdirDeathsPerpersonYr<-(predTempdirDeaths/(2*predTempnhost))
    
    
    # require only one survey and 1 ages (1) and measure 16
    fieldNew[[scenario]]<- outputTemp[[scenario]][outputTemp[[scenario]]$V1==1&outputTemp[[scenario]]$V3==0 & outputTemp[[scenario]]$V2==1,]
    # replace field new (but nHost, replace number 16)
    fieldNew[[scenario]]$V3[fieldNew[[scenario]]$V3==modelResultID$nHost] <- 16
    # transform
    fieldNew[[scenario]]$V4 <- predTempdirDeathsPerpersonYr
    rm(predTempnhost, predTempdirDeaths, predTempdirDeathsPerpersonYr)	
    # check dimensions
    if (dim(fieldNew[[scenario]])[1] != dim(fieldTemp[[scenario]])[1]){
      stop(paste("   ERROR in obj 9 for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
  
  
}

#---------------------------------------------------------------------------------------------
# All cause mortality
#---------------------------------------------------------------------------------------------
OBJ_indirectMortality <- function(modelResultID, fieldDataAll, fieldDataID,  oDirAbs, runID,  scenarioNums){
  #  (scenarios 401,402,408,411,414,415,416,417,418,422,426)
  
  # read in predictions for particular scenario number and transform if needed for new field data and transform if needed for new field data
  outputTemp <- list()
  fieldTemp <- list()
  fieldNew <- list()
  for (scenario in 1:length(scenarioNums)) {
    fieldTemp[[scenario]] <- fieldDataAll[fieldDataAll$V1==scenarioNums[scenario],]
    outputTemp[[scenario]] <- read.table(paste0(oDirAbs,"/output_",runID,"/wu_5000_",scenarioNums[scenario],".txt"))
    
    # require only one survey and 1 ages (1) and measure 16
    fieldNew[[scenario]]<- outputTemp[[scenario]][outputTemp[[scenario]]$V1==1&outputTemp[[scenario]]$V3==21 & outputTemp[[scenario]]$V2==1,]
    # replace field new (but nHost, replace number 16)
    fieldNew[[scenario]]$V3[fieldNew[[scenario]]$V3==21] <- 22
    # transform
    fieldNew[[scenario]]$V4 <- outputTemp[[scenario]]$V4	
    # check dimensions
    if (dim(fieldNew[[scenario]])[1] != dim(fieldTemp[[scenario]])[1]){
      stop(paste("   ERROR in obj 10 for scenario", scenarioNums[scenario], sep=" "))
    }
    # 
  }
  return(fieldNew)
  
}


# ===========================================================================================
# MAIN FUNCTION TO PLOT ALL RESULTS - allPlots
# ===========================================================================================

createSimulatedFieldData<-function(runID,ExperimentName,GitDir,ExperimentDir,xmlDir,outputDir){
  #
  xDirAbs = paste0(ExperimentDir,ExperimentNameHetGP,xmlDir)
  oDirAbs = paste0(ExperimentDir,ExperimentNameHetGP,outputDir)
  
  #runDir<-paste(sep = "","run",runID,"/")
  
  # to delete
  # # find the value of the loss function
  # lossFunctionVal <- list.files(runDir, pattern = "^[lf]", all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = TRUE)
  
  
  # find name of parameter value and ID file
  #  param_file <- list.files(paste0(xDirAbs,runID), pattern = "[.]xmlt", all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = TRUE)
  param_file = paste0(ExperimentDir, ExperimentNameHetGP,"/PriorDraws/scenarios/scenario_buffers/small/",runID,"/params31")
  #paramID <- sub(".xmlt","", param_file) #???
  
  # to delete
  # # name the pdf for plotting output
  # pdf(file=paste(dirname,"/run",runID,"_",lossFunctionVal,"_", paramID,".pdf", sep=""),paper="a4")
  
  # -----------------------------------------------------------------
  # read in all field data
  fieldDataAll <-read.table(paste0(GitDir,"LikelihoodComp/fieldData.txt"))
  
  # get codes for field data ids
  fieldDataID <- fieldData_nameCodes()
  
  # get codes for model Results ids
  modelResultID <- modelResults_nameCodes()
  # -----------------------------------------------------------------
  # density bias for observed parasite densities in garki sites and non-garki
  readDensBiasFromFile <- 1
  # readDensBiasFromFile <- 0
  
  densBias <- getDensBiasVals(readDensBiasFromFile, runDir, paramFile)
  
  # -----------------------------------------------------------------
  # create new field data but remove the observations and make NA (replicate into simulated data)
  fieldDataAll_simulated <- fieldDataAll
  fieldDataAll_simulated$V5 <- NA
  
  # -------------------------------------------------------------------------
  # 1- INCIDENCE after intervention - 
  # (scenarios 30)
  scenarioNums <- c(30)
  tempOut <- OBJ_AgePatternIncidenceAfterIntervention(modelResultID, fieldDataAll, fieldDataID, densBias,  oDirAbs, runID, scenarioNums)
  
  for (scenario in 1:length(scenarioNums)){
    # reorder simulated output to match field data
    tempReorder <- tempOut[[scenario]][order(tempOut[[scenario]]$V3),]
    # check matches field data
    if(any(tempReorder$V1 != fieldDataAll_simulated$V2[fieldDataAll_simulated$V1==scenarioNums[scenario]])){stop(paste(' ERROR reordering (V1) wrong for scenario',scenarioNums[scenario],sep=" "))}
    if(any(tempReorder$V2 != fieldDataAll_simulated$V3[fieldDataAll_simulated$V1==scenarioNums[scenario]])){stop(paste(' ERROR reordering (V2) wrong for scenario',scenarioNums[scenario],sep=" "))}
    if(any(tempReorder$V3 != fieldDataAll_simulated$V4[fieldDataAll_simulated$V1==scenarioNums[scenario]])){stop(paste(' ERROR reordering (V3) wrong for scenario',scenarioNums[scenario],sep=" "))}
    # put simulated output into fieldDataAll_simulated
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario]]<-tempReorder$V4
    rm(tempReorder)
  }
  
  rm(tempOut)
  # -------------------------------------------------------------------------
  # 2&3- ASEXUAL - prevalence and density (2 plots produced)
  # (scenarios 24,28,29,35,34 and 31)
  scenarioNums <- c(24, 28, 29, 35, 34, 31)
  tempOut <- OBJ_AgePatternPrevalenceAndDensity(modelResultID, fieldDataAll, fieldDataID, densBias,  oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    
    if (scenario == 5){# ie 34
      # reorder simulated output to match field data
      tempReorder_1 <- tempOut[[scenario]][order(tempOut[[scenario]]$V3),]
      tempReorder <- tempReorder_1[tempReorder_1$V3!=modelResultID$totalInfs & tempReorder_1$V3!=modelResultID$totalPatentInf,]
      # check matches field data
      if(any(tempReorder$V1 != fieldDataAll_simulated$V2[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4!=fieldDataID$totalpatentInf])){stop(paste(' ERROR reordering (V1) wrong for scenario',scenarioNums[scenario],sep=" "))}
      if(any(tempReorder$V2 != fieldDataAll_simulated$V3[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4!=fieldDataID$totalpatentInf])){stop(paste(' ERROR reordering (V2) wrong for scenario',scenarioNums[scenario],sep=" "))}
      if(any(tempReorder$V3 != fieldDataAll_simulated$V4[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4!=fieldDataID$totalpatentInf])){stop(paste(' ERROR reordering (V3) wrong for scenario',scenarioNums[scenario],sep=" "))}
      # put simulated output into fieldDataAll_simulated
      # only replace nHost (0), densities(5) and nPatent(3)
      fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4!=fieldDataID$totalpatentInf] <- tempReorder$V4
      print(" ")
      print(" ")
      print(" ")
      print(scenario)
      print(fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario]])
    }else{
      # reorder simulated output to match field data
      tempReorder <- tempOut[[scenario]][order(tempOut[[scenario]]$V3),]
      # check matches field data
      if(any(tempReorder$V1 != fieldDataAll_simulated$V2[fieldDataAll_simulated$V1==scenarioNums[scenario]])){stop(paste(' ERROR reordering (V1) wrong for scenario',scenarioNums[scenario],sep=" "))}
      if(any(tempReorder$V2 != fieldDataAll_simulated$V3[fieldDataAll_simulated$V1==scenarioNums[scenario]])){stop(paste(' ERROR reordering (V2) wrong for scenario',scenarioNums[scenario],sep=" "))}
      if(any(tempReorder$V3 != fieldDataAll_simulated$V4[fieldDataAll_simulated$V1==scenarioNums[scenario]])){stop(paste(' ERROR reordering (V3) wrong for scenario',scenarioNums[scenario],sep=" "))}
      # put simulated output into fieldDataAll_simulated
      fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario]]<-tempReorder$V4
    }
    rm(tempReorder)
  }
  rm(tempOut)
  # -------------------------------------------------------------------------
  # 4- MOI - multiplicity of infection
  # (scenario 34)
  scenarioNums <- c(34)
  tempOut <- OBJ_AgePatternMOI(modelResultID, fieldDataAll, fieldDataID, oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # only replace totalpatentInf (13)
    # no need to reorder as only one outcome
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4==fieldDataID$totalpatentInf] <- tempOut[[scenario]]$V4[tempOut[[scenario]]$V3==modelResultID$totalPatentInf]
  }
  rm(tempOut)
  # -------------------------------------------------------------------------
  # 5a- ACUTE EPISODES: NDIOP AND DIELMO INCIDENCE
  # (scenarios 232 and 233)
  scenarioNums <- c(232, 233)
  tempOut <- OBJ_AgePatternIncidenceClinicalMalaria_Senegal(modelResultID, fieldDataAll, fieldDataID,  oDirAbs, runID,  scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # replace recInc (18)
    # no need to reorder as only select what is needed
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4==fieldDataID$recInc] <- tempOut[[scenario]]$V4[tempOut[[scenario]]$V3==18]
  }
  rm(tempOut)
  # -------------------------------------------------------------------------
  # 5b- ACUTE EPISODES: IDETE
  # (scenario 49)
  scenarioNums <- c(49)
  tempOut <- OBJ_AgePatternIncidenceClinicalMalaria_Idete(modelResultID, fieldDataAll, fieldDataID, oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # replace recInc (18)
    # no need to reorder as only select what is needed
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4==fieldDataID$recInc] <- tempOut[[scenario]]$V4[tempOut[[scenario]]$V3==18]
  }
  rm(tempOut)
  # -------------------------------------------------------------------------
  # 6- ACUTE EPISODES: DIELMO PYROGENIC THRESHOLD
  # (scenario 234)
  scenarioNums <- c(234)
  tempOut <- OBJ_AgePatternThresholdClinicalAttack(modelResultID, fieldDataAll, fieldDataID, densBias, oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # replace recInc (18)
    # no need to reorder 
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario] & fieldDataAll_simulated$V4==fieldDataID$pyrogenThres] <- tempOut[[scenario]]$V4[tempOut[[scenario]]$V3==19]
  }
  rm(tempOut)
  # -------------------------------------------------------------------------
  # 7- SEVERE MALARIA (prevalence vs episodes)
  # (scenarios 501,502,503,504,505,506,507,508,509,510,511,512,514,515,516,517,518,519,520,521,522,523,524,525,526,527)
  scenarioNums <- c(501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527)
  tempOut <- OBJ_SevereEpisodesVsPrevalence(modelResultID, fieldDataAll, fieldDataID, oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # replace recInc (18)
    # no need to reorder 
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario] ] <- tempOut[[scenario]]$V5
  }
  rm(tempOut)
  # ----------------------------------------------------------------------
  # 8- SEVERE MALARIA RR (relative risk)
  # (scenarios 158,167,173,176)
  scenarioNums <- c(158,167,173,176)
  tempOut <- OBJ_AgePatternOfSevere(modelResultID, fieldDataAll, fieldDataID, oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # no need to reorder 
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario]] <- tempOut[[scenario]]$V4
  }
  rm(tempOut)
  # -----------------------------------------------------------------------
  # 9- DIRECT MORTALITY
  # (scenarios 301,302,303,312,316,317,318,326,327)
  scenarioNums <- c(301,302,303,312,316,317,318,326,327)
  tempOut <- OBJ_DirectMalariaMortality(modelResultID, fieldDataAll, fieldDataID, oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # no need to reorder 
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario]] <- tempOut[[scenario]]$V4
  }
  # --------------------------------------------------------------------------------
  # 10- INDIRECT MORTALITY
  # (scenarios 401,402,408,411,414,415,416,417,418,422,426)
  scenarioNums <- c(401,402,408,411,414,415,416,417,418,422,426)
  tempOut <- OBJ_indirectMortality(modelResultID, fieldDataAll, fieldDataID, oDirAbs, runID, scenarioNums)
  for (scenario in 1:length(scenarioNums)){
    # no need to reorder 
    fieldDataAll_simulated$V5[fieldDataAll_simulated$V1==scenarioNums[scenario]] <- tempOut[[scenario]]$V4
  }
  rm(tempOut)
  # ---------------------------------------------------------------------------------------
  # print final simulated to text file
  write.csv(fieldDataAll_simulated, file = paste0(GitDir, "/simulatedFieldData/", "simulatedFieldData_", ExperimentName,runID, ".txt"), quote = FALSE)
  
  write.table(fieldDataAll_simulated, file = paste0(GitDir, "/simulatedFieldData/", "simulatedFieldData_", ExperimentName,runID, ".txt"), col.names=FALSE, row.names=FALSE, sep = "\t")
  
  
}

#for hetGP: 
runID="main_0020" 
createSimulatedFieldData(runID, ExperimentNameHetGP,GitDir,ExperimentDir,xmlDir,outputDir)
#for GPSG
runID="main_0010" 
createSimulatedFieldData(runID, ExperimentNameGPSG,GitDir,ExperimentDir,xmlDir,outputDir)
