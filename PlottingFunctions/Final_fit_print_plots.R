GitDir =  "/scicore/home/smith/reiker/GitRepos/om_fitting/"
MainExperimentDir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_hetGP_norm_sampling/"
old_param_output_dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_04_full_fit_small_hetGP/PriorDraws/Output/old_params/"
GP_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_hetGP_norm_sampling/"
GPSG_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_GPSG_norm_sampling/"

if(!require(pacman)){install.packages(pacman);require(pacman)}else{require(pacman)}
pacman::p_load(gridExtra,ggplot2,colorspace,tidyverse)

D=23;W=11;
source(paste0(GitDir,"EvalOM.parallel.R"))
source(paste0(GitDir,"PlottingFunctions/Final_fit_base_functions.R"))
# ==========================================================================================================================================================
# ==========================================================================================================================================================
# MAIN FUNCTION TO PLOT ALL RESULTS - allPlots
# ==========================================================================================================================================================
# ==========================================================================================================================================================

colnames <- list()
colnames[[length(colnames)+1]] <- paste0("elapsed")
for (w in 1:W) {colnames[[length(colnames)+1]] <- paste0("r2_",sprintf("%02i",w))}
colnames[[length(colnames)+1]] <- paste0("r2_w")
colnames[[length(colnames)+1]] <- paste0("min_lcb")
colnames[[length(colnames)+1]] <- paste0("mpe")

colnames <- unlist(colnames)

# =============================================================================
# best_hetgp_params
# =============================================================================

setwd(paste0(GP_Dir))

load(paste0(GP_Dir,'time.hetgp'))
time_hetGP = time

time_frame_hetGP <- matrix(0,nrow=length(time_hetGP),(W+4))
for (g in 1:length(time_hetGP)) {
  time_frame_hetGP[g,1] <- as.numeric(time_hetGP[[g]]$time["elapsed"])
  time_frame_hetGP[g,c(2:(2+W))] <- unlist(time_hetGP[[g]]$r2)
  time_frame_hetGP[g,(W+3)] <- time_hetGP[[g]]$obj$min_lcb
  time_frame_hetGP[g,(W+4)] <- time_hetGP[[g]]$obj$mpe
}

time_frame_hetGP <- as.data.frame(time_frame_hetGP)
colnames(time_frame_hetGP) <- colnames
time_frame_hetGP$Algorithm = "GP"
time_frame_hetGP$iteration = seq(1,dim(time_frame_hetGP)[1],1)

load("eval.state.newer")
hetgp_best = vector(length=length(evaluations$best))
for (i in 1:length(evaluations$best)){
  hetgp_best[[i]]=sum(evaluations$best[[i]]$known$Objectives)
  #print(sum(evaluations$best[[i]]$known$Objectives))
  if(hetgp_best[[i]]==0 | (i>1 && hetgp_best[[i]]>hetgp_best[[i-1]])){hetgp_best[[i]]=hetgp_best[[i-1]]}
}
time_frame_hetGP$best_known=NA
time_frame_hetGP$best_known[c(1:length(hetgp_best))]=hetgp_best

for(i in 2:dim(time_frame_hetGP)[1]){
  if(time_frame_hetGP$elapsed[i]<time_frame_hetGP$elapsed[i-1]){
    base = time_frame_hetGP$elapsed[i-1]
    for(j in (dim(time_frame_hetGP)[1]-1):i){time_frame_hetGP$best_known[j+1] = time_frame_hetGP$best_known[j]}
    time_frame_hetGP$best_known[i] = time_frame_hetGP$best_known[i-1]
    for(j in i:dim(time_frame_hetGP)[1]){time_frame_hetGP$elapsed[j] = time_frame_hetGP$elapsed[j] + base }
    break()
  }
}


best_gp_iter =time_frame_hetGP[which.min(time_frame_hetGP$best_known),"iteration"] +1 #(+1 only here)
gp_theta = evaluations$best[[best_gp_iter]]$known$parameters
gp_params = quant.to.prior(gp_theta)



# =============================================================================
#best GPSG
# =============================================================================

#gpsg dataframe
load(paste0(GPSG_Dir,'/time.gpsg'))
time_gpsg = time
time_frame_gpsg <- matrix(0,nrow=length(time_gpsg),(W+4))
for (g in 1:length(time_gpsg)) {
  time_frame_gpsg[g,1] <- as.numeric(time_gpsg[[g]]$time["elapsed"])
  time_frame_gpsg[g,c(2:(2+W))] <- unlist(time_gpsg[[g]]$r2)
  time_frame_gpsg[g,(W+3)] <- time_gpsg[[g]]$obj$min_lcb
  time_frame_gpsg[g,(W+4)] <- time_gpsg[[g]]$obj$mpe
}
time_frame_gpsg <- as.data.frame(time_frame_gpsg)
colnames(time_frame_gpsg) <- colnames
time_frame_gpsg$Algorithm = "GPSG"
time_frame_gpsg$iteration = seq(1,dim(time_frame_gpsg)[1],1)

load(paste0(GPSG_Dir,'/eval.state.newer'))
gpsg_best = vector(length=length(evaluations$best))
for (i in 1:length(evaluations$best)){
  gpsg_best[[i]]=sum(evaluations$best[[i]]$known$Objectives)
  if(gpsg_best[[i]]==0 | (i>1 && gpsg_best[[i]]>gpsg_best[[i-1]])){gpsg_best[[i]]=gpsg_best[[i-1]]}
  
}
time_frame_gpsg$best_known=NA
time_frame_gpsg$best_known[c(1:length(gpsg_best))]=gpsg_best
for(i in 2:dim(time_frame_gpsg)[1]){
  if(time_frame_gpsg$elapsed[i]<time_frame_gpsg$elapsed[i-1]){
    base = time_frame_gpsg$elapsed[i-1]
    for(j in i:dim(time_frame_gpsg)[1]){time_frame_gpsg$elapsed[j] = time_frame_gpsg$elapsed[j] + base }
    break()
  }
}


best_gpsg_iter =time_frame_gpsg[which.min(time_frame_gpsg$best_known),"iteration"] #(+1 only here)

gpsg_theta = evaluations$best[[best_gpsg_iter]]$known$parameters
gpsg_params = quant.to.prior(gpsg_theta)

OldDensityBiasGarki =4.79610772546704
OldDensityBiasnonGarki =0.177378570987455


#allPlots<-function(runID,dirname,params){
  GPrunID <- paste0("output_main_",sprintf("%04i",best_gp_iter))
  GPrunDir <- paste0(GP_Dir,"/PriorDraws/Output/small/",GPrunID,"/")

    
  GPSGrunID <- paste0("output_main_",sprintf("%04i",best_gpsg_iter))
  GPSGrunDir <- paste0(GPSG_Dir,"/PriorDraws/Output/small/",GPSGrunID,"/")

  OldrunDir <- paste0("/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_04_full_fit_small_hetGP/PriorDraws/Output/old_params/")
  
  # -----------------------------------------------------------------
  # read in all field data
  fieldDataAll <-read.table(paste0(GitDir,"/LikelihoodComp/fieldData",".txt"))
  
  # get codes for field data ids
  fieldDataID <- fieldData_nameCodes()
  
  # get codes for model Results ids
  modelResultID <- modelResults_nameCodes()
  # --------------------------------------------------------------------
  # density bias for observed parasite densities in garki sites and non-garki
  densBiasGP <- list()
  densBiasGPSG <- list()
  densBiasOld <- list()
  
  densBiasGP$nonGarki <- gp_params$DensityBiasnonGarki
  densBiasGP$garki <- gp_params$DensityBiasGarki
  
  densBiasGPSG$nonGarki <- gpsg_params$DensityBiasnonGarki
  densBiasGPSG$garki <- gpsg_params$DensityBiasGarki
  
  densBiasOld$nonGarki <- OldDensityBiasnonGarki
  densBiasOld$garki <-OldDensityBiasGarki
  # --------------------------------------------------------------------
  # set loss Vector etc
  lossVector <- array(NA,dim=c(11,3))
  lossVector_logLH <- array(NA,dim=c(11,3))
  LF <- list()
  LF$GP <- list()
  LF$GPSG <- list()
  LF$Old <- list()
  
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  #  OBJECTIVE 1
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  # (scenarios 30)
  LF$GP$obj_1 <- OBJ_AgePatternIncidenceAfterIntervention(fieldDataAll, fieldDataID, modelResultID, densBiasGP, GPrunDir)
  LF$GPSG$obj_1 <- OBJ_AgePatternIncidenceAfterIntervention(fieldDataAll, fieldDataID, modelResultID, densBiasGPSG, GPSGrunDir)
  LF$Old$obj_1 <- OBJ_AgePatternIncidenceAfterIntervention(fieldDataAll, fieldDataID, modelResultID, densBiasOld, OldrunDir)
  
  
  Obj1 = data.frame(age = LF$GP$obj_1$Data$FieldData$Matsari$age)
  Obj1$prevalence_data = LF$GP$obj_1$Data$FieldData$Matsari$prevalence
  
  Obj1$gp = LF$GP$obj_1$Data$ModelResults$Matsari$prevalence
  Obj1$gpsg = LF$GPSG$obj_1$Data$ModelResults$Matsari$prevalence
  Obj1$old = LF$Old$obj_1$Data$ModelResults$Matsari$prevalence
  
  colours = c("black","#FF7745", "#720026")
  
  #lossVector[1] <- LF$obj_1$LF$orig$total
  #lossVector_logLH[1] <- LF$obj_1$LF$logLH$total
  
  Obj1_gathered = gather(Obj1, key = "source", value = "value", -age)
  pdf(file="./OutputPlots/Model_fit_Objective_1.pdf",width=6.98,height=5.65)
  ggplot(Obj1_gathered,aes(x=age,y=value))+
    geom_point(data = Obj1_gathered[which(Obj1_gathered$source=="prevalence_data"),], aes(y=value),size=4,color="#D3D4D9")+
    geom_line(data = Obj1_gathered[which(Obj1_gathered$source!="prevalence_data"),], aes(color = source,linetype=source),size=1.2)+
    aes(group=rev(source)) + 
    scale_x_continuous(trans='log10') +
    scale_linetype_manual(name="Parameterization",values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    #facet_wrap(~site,scales = "free") +
    labs(title=paste0("Age pattern of prevalence of infection after intervention"), x="Age", y="Prevalence")+
    theme_light() + 
    theme(legend.position = c(0.01,0.86),
          legend.justification = "left",
          legend.title = element_text(face="bold",size = 12),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=15, hjust=1,colour="grey30"),
          axis.text.x=element_text(size=15, vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=15, colour="grey20",face="bold"),
          axis.title.x=element_text(size=15, colour="grey20",face="bold"),
          
          plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title=element_text(face="bold",size=15,vjust=1),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  
  lossVector[1,1] <- LF$Old$obj_1$LF$orig$total
  lossVector_logLH[1,1] <-  LF$Old$obj_1$LF$logLH$total
  
  lossVector[1,2] <- LF$GP$obj_1$LF$orig$total
  lossVector_logLH[1,2] <-  LF$GP$obj_1$LF$logLH$total
  
  lossVector[1,3] <- LF$GPSG$obj_1$LF$orig$total
  lossVector_logLH[1,3] <-  LF$GPSG$obj_1$LF$logLH$total
  
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  #  OBJECTIVES 2+3
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  # --------------------------------------------------------------------
  # ASEXUAL - prevalence and density (2 plots produced)
  # (scenarios 24,28,29,35,34 and 31)
  LF_temp <- list()
  LF_temp$GP <- list()
  LF_temp$GPSG <- list()
  LF_temp$Old <- list()
  
  
  LF_temp$GP <- OBJ_AgePatternPrevalenceAndDensity(fieldDataAll, fieldDataID, modelResultID, densBiasGP, GPrunDir)
  LF_temp$GPSG <- OBJ_AgePatternPrevalenceAndDensity(fieldDataAll, fieldDataID, modelResultID, densBiasGPSG, GPSGrunDir)
  LF_temp$Old <- OBJ_AgePatternPrevalenceAndDensity(fieldDataAll, fieldDataID, modelResultID, densBiasOld, OldrunDir)
  
  
  LF$GP$obj_2 <- LF_temp$GP$LF$obj_2
  LF$GPSG$obj_2 <- LF_temp$GPSG$LF$obj_2
  LF$Old$obj_2 <- LF_temp$Old$LF$obj_2
  
  LF$GP$obj_3 <- LF_temp$GP$LF$obj_3
  LF$GPSG$obj_3 <- LF_temp$GPSG$LF$obj_3
  LF$Old$obj_3 <- LF_temp$Old$LF$obj_3
  
  Obj2_ls=list()
  for(i in 1: length(LF_temp$GP$Data$FieldData)){
    Obj2_ls[[i]] = data.frame(age = LF_temp$GP$Data$FieldData[[i]]$age)
    Obj2_ls[[i]]$prevalence_data = LF_temp$GP$Data$FieldData[[i]]$prevalence
    Obj2_ls[[i]]$old = LF_temp$Old$Data$ModelResults[[i]]$prevalence
    Obj2_ls[[i]]$gp = LF_temp$GP$Data$ModelResults[[i]]$prevalence
    Obj2_ls[[i]]$gpsg = LF_temp$GPSG$Data$ModelResults[[i]]$prevalence
    Obj2_ls[[i]]$site = names(LF_temp$GP$Data$FieldData)[i]
  }
  Obj2 = data.frame(Obj2_ls[[1]])
  for(i in 2: length(LF_temp$GP$Data$FieldData)){
    Obj2 = rbind(Obj2,Obj2_ls[[i]])
  }
  #Obj 2 plot:prevalence
  Obj2_gathered = as.data.frame(gather(Obj2, key = "source", value = "value", -c(age,site)))
  pdf(file="./OutputPlots/Model_fit_Objective_2.pdf",height=5.13,width=7.76)
  ggplot(Obj2_gathered,aes(x=age,y=value))+
    geom_point(data = Obj2_gathered[which(Obj2_gathered$source=="prevalence_data"),], aes(y=value),size=2.1,color="#D3D4D9")+
    geom_line(data = Obj2_gathered[which(Obj2_gathered$source!="prevalence_data"),], aes(color = source,linetype=source),size=1.2)+
    aes(group=rev(source)) + 
    scale_linetype_manual(name="Parameterization",values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    facet_wrap(~site,scales = "free_x") +
    labs(title=paste0("Age patterns of prevalence"), x="Age", y="Prevalence")+
    theme_light() + 
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 14),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=14,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=14,  colour="grey20",face="bold"),
          
          plot.title=element_text(face="bold", vjust=1,size=15),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
)  
  dev.off()
  #Obj 3 plot:density

  Obj3_ls=list()
  for(i in 1: length(LF_temp$GP$Data$FieldData)){
    Obj3_ls[[i]] = data.frame(age = LF_temp$GP$Data$FieldData[[i]]$age)
    Obj3_ls[[i]]$density = LF_temp$GP$Data$FieldData[[i]]$density
    Obj3_ls[[i]]$old = LF_temp$Old$Data$ModelResults[[i]]$density
    Obj3_ls[[i]]$gp = LF_temp$GP$Data$ModelResults[[i]]$density
    Obj3_ls[[i]]$gpsg = LF_temp$GPSG$Data$ModelResults[[i]]$density
    Obj3_ls[[i]]$site = names(LF_temp$GP$Data$FieldData)[i]
  }
  Obj3 = data.frame(Obj3_ls[[1]])
  for(i in 2: length(LF_temp$GP$Data$FieldData)){
    Obj3 = rbind(Obj3,Obj3_ls[[i]])
  }
  #Obj 2 plot:prevalence
  Obj3_gathered = as.data.frame(gather(Obj3, key = "source", value = "value", -c(age,site)))
  pdf(file="./OutputPlots/Model_fit_Objective_3.pdf",height=5.82,width=9.32)
  ggplot(Obj3_gathered,aes(x=age,y=value))+
    geom_point(data = Obj3_gathered[which(Obj3_gathered$source=="density"),], aes(y=value),size=2.1,color="#D3D4D9")+
    geom_line(data = Obj3_gathered[which(Obj3_gathered$source!="density"),], aes(color = source,linetype=source),size=1.2)+
    aes(group=rev(source)) + 
    scale_linetype_manual(name="Parameterization",values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    scale_x_continuous(trans='log10') +
    facet_wrap(~site,scales = "free") +
    labs(title=paste0("Age patterns of parasite densities"), x="Age", y="Density")+
    theme_light() + 
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 14),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=14,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=14,  colour="grey20",face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"),
    
          plot.title=element_text(face="bold", vjust=1,size=15),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  
  lossVector[2,1] <- LF$Old$obj_2$orig$total
  lossVector_logLH[2,1] <-  LF$Old$obj_2$logLH$total
  
  lossVector[2,2] <- LF$GP$obj_2$orig$total
  lossVector_logLH[2,2] <-  LF$GP$obj_2$logLH$total
  
  lossVector[2,3] <- LF$GPSG$obj_2$orig$total
  lossVector_logLH[2,3] <-  LF$GPSG$obj_2$logLH$total
  
  
  
  
  lossVector[3,1] <- LF$Old$obj_3$orig$total
  lossVector_logLH[3,1] <-  LF$Old$obj_3$logLH$total
  
  lossVector[3,2] <- LF$GP$obj_3$orig$total
  lossVector_logLH[3,2] <-  LF$GP$obj_3$logLH$total
  
  lossVector[3,3] <- LF$GPSG$obj_3$orig$total
  lossVector_logLH[3,3] <-  LF$GPSG$obj_3$logLH$total
  
  
  
  
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  #  OBJECTIVES 4: Multiplicity of infection
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  # --------------------------------------------------------------------
  # MOI - multiplicity of infection
  # (scenario 34)

  LF$GP$obj_4 <- OBJ_AgePatternMOI(fieldDataAll, fieldDataID, modelResultID, GPrunDir)
  LF$GPSG$obj_4 <- OBJ_AgePatternMOI(fieldDataAll, fieldDataID, modelResultID, GPSGrunDir)
  LF$Old$obj_4 <- OBJ_AgePatternMOI(fieldDataAll, fieldDataID, modelResultID, OldrunDir)
  
  
  Obj4 = data.frame(age = LF$GP$obj_4$Data$FieldData$Navrongo$age)
  Obj4$moi = LF$GP$obj_4$Data$FieldData$Navrongo$moi
  
  Obj4$gp = LF$GP$obj_4$Data$ModelResults$moi
  Obj4$gpsg = LF$GPSG$obj_4$Data$ModelResults$moi
  Obj4$old = LF$Old$obj_4$Data$ModelResults$moi
  
  #Plot
  
  Obj4_gathered = gather(Obj4, key = "source", value = "value", -age)
  pdf(file="./OutputPlots/Model_fit_Objective_4.pdf",height=4.24,width=4.98)
  ggplot(Obj4_gathered,aes(x=age,y=value))+
    geom_point(data = Obj4_gathered[which(Obj4_gathered$source=="moi"),], aes(y=value),size=3,color="#D3D4D9")+
    geom_line(data = Obj4_gathered[which(Obj4_gathered$source!="moi"),], aes(color = source,linetype=source),size=1.2)+
    aes(group=rev(source)) + 
    #scale_x_continuous(trans='log10') +
    scale_linetype_manual(name="Parameterization",values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    #facet_wrap(~site,scales = "free") +
    labs(title=paste0("Age pattern of number of cuncurrent infections"), x="Age", y="Multiplicity of infection (MOI)")+
    theme_light() + 
    theme(legend.position = c(0.99,0.83),
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 11),
          legend.text = element_text(size = 11),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=12,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=12,colour="grey20",face="bold"),
          
          plot.title=element_text(face="bold", vjust=1, 
                                  size = 14),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  
  lossVector[4,1] <- LF$Old$obj_4$LF$orig$total
  lossVector_logLH[4,1] <-  LF$Old$obj_4$LF$logLH$total
  
  lossVector[4,2] <- LF$GP$obj_4$LF$orig$total
  lossVector_logLH[4,2] <-  LF$GP$obj_4$LF$logLH$total
  
  lossVector[4,3] <- LF$GPSG$obj_4$LF$orig$total
  lossVector_logLH[4,3] <-  LF$GPSG$obj_4$LF$logLH$total
  
  
  
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  #  OBJECTIVES 5: Severe disease
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  # --------------------------------------------------------------------
  # ACUTE EPISODES: NDIOP AND DIELMO INCIDENCE
  # (scenarios 232 and 233)

  LF$GP$obj_5a <- OBJ_AgePatternIncidenceClinicalMalaria_Senegal(fieldDataAll, fieldDataID, modelResultID, GPrunDir)
  LF$GPSG$obj_5a <- OBJ_AgePatternIncidenceClinicalMalaria_Senegal(fieldDataAll, fieldDataID, modelResultID, GPSGrunDir)
  LF$Old$obj_5a <- OBJ_AgePatternIncidenceClinicalMalaria_Senegal(fieldDataAll, fieldDataID, modelResultID, OldrunDir)
  
  Obj5a_ls=list()
  Obj5a_ls$Ndiop = data.frame(age = LF$Old$obj_5a$Data$FieldData$Ndiop$age)
  Obj5a_ls$Dielmo = data.frame(age = LF$Old$obj_5a$Data$FieldData$Dielmo$age)
  
  Obj5a_ls$Ndiop$Inc = LF$Old$obj_5a$Data$FieldData$Ndiop$incid
  Obj5a_ls$Dielmo$Inc = LF$Old$obj_5a$Data$FieldData$Dielmo$incid
  
  Obj5a_ls$Ndiop$gp = LF$GP$obj_5a$Data$ModelResults$Ndiop$incid
  Obj5a_ls$Ndiop$gpsg = LF$GPSG$obj_5a$Data$ModelResults$Ndiop$incid
  Obj5a_ls$Ndiop$old = LF$Old$obj_5a$Data$ModelResults$Ndiop$incid
  Obj5a_ls$Ndiop$site = "Ndiop"
  
  
  Obj5a_ls$Dielmo$gp = LF$GP$obj_5a$Data$ModelResults$Dielmo$incid
  Obj5a_ls$Dielmo$gpsg = LF$GPSG$obj_5a$Data$ModelResults$Dielmo$incid
  Obj5a_ls$Dielmo$old = LF$Old$obj_5a$Data$ModelResults$Dielmo$incid  
  Obj5a_ls$Dielmo$site = "Dielmo"
  
  Obj5a = data.frame(Obj5a_ls[[1]])
  for(i in 2: length(LF$GP$obj_5a$Data$FieldData)){
    Obj5a = rbind(Obj5a,Obj5a_ls[[i]])
  }
  #Obj 2 plot:prevalence
  Obj5a_gathered = as.data.frame(gather(Obj5a, key = "source", value = "value", -c(age,site)))
  pdf(file="./OutputPlots/Model_fit_Objective_5a.pdf",height=3.90,width=6.34)
  ggplot(Obj5a_gathered,aes(x=age,y=value))+
    geom_point(data = Obj5a_gathered[which(Obj5a_gathered$source=="Inc"),], aes(y=value),size=2.1,color="#D3D4D9")+
    geom_line(data = Obj5a_gathered[which(Obj5a_gathered$source!="Inc"),], aes(color = source,linetype=source),size=1.2)+
    aes(group=rev(source)) + 
    scale_linetype_manual(name="Parameterization",values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    #scale_colour_manual(values=c(gp="#00AFBB",gpsg="burlywood",old="grey10"), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    scale_x_continuous(trans='log10') +
    facet_wrap(~site,scales = "free") +
    labs(title=paste0("Age pattern of incidence of clinical malaria"), x="Age", y="Episodes / year")+
    theme_light() + 
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 14),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=14,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=14,  colour="grey20",face="bold"),
          
          plot.title=element_text(size=14,face="bold", vjust=1),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  lossVector[5,1] <- LF$Old$obj_5a$LF$orig$total
  lossVector_logLH[5,1] <-  LF$Old$obj_5a$LF$logLH$total
  
  lossVector[5,2] <- LF$GP$obj_5a$LF$orig$total
  lossVector_logLH[5,2] <-  LF$GP$obj_5a$LF$logLH$total
  
  lossVector[5,3] <- LF$GPSG$obj_5a$LF$orig$total
  lossVector_logLH[5,3] <-  LF$GPSG$obj_5a$LF$logLH$total
  
  #lossVector[5] <- LF$obj_5a$orig$total
  #lossVector_logLH[5] <- LF$obj_5a$logLH$total
  # --------------------------------------------------------------------
  # ACUTE EPISODES: IDETE
  # (scenario 49)
  
  LF$GP$obj_5b <- OBJ_AgePatternIncidenceClinicalMalaria_Idete(fieldDataAll, fieldDataID, modelResultID, GPrunDir)
  LF$GPSG$obj_5b <- OBJ_AgePatternIncidenceClinicalMalaria_Idete(fieldDataAll, fieldDataID, modelResultID, GPSGrunDir)
  LF$Old$obj_5b <- OBJ_AgePatternIncidenceClinicalMalaria_Idete(fieldDataAll, fieldDataID, modelResultID, OldrunDir)
  
  Obj5b = data.frame(age = LF$GP$obj_5b$Data$FieldData$Idete$age)
  Obj5b$inc = LF$GP$obj_5b$Data$FieldData$Idete$incid
  
  Obj5b$gp = LF$GP$obj_5b$Data$ModelResults$Idete$incid
  Obj5b$gpsg = LF$GPSG$obj_5b$Data$ModelResults$Idete$incid
  Obj5b$old = LF$Old$obj_5b$Data$ModelResults$Idete$incid

  Obj5b_gathered = gather(Obj5b, key = "source", value = "value", -age)
  accessConst<-2.7975236259999985
  
  pdf(file="./OutputPlots/Model_fit_Objective_5b.pdf",height=3.86,width=6.34)
  ggplot(Obj5b_gathered,aes(x=age,y=value))+
    geom_point(data = Obj5b_gathered[which(Obj5b_gathered$source=="inc"),], aes(y=value),size=4.5,color="#D3D4D9")+
    geom_line(data = Obj5b_gathered[which(Obj5b_gathered$source!="inc"),], aes(color = source,linetype=source),size=1.2)+
    
    geom_line(data = Obj5b_gathered[which(Obj5b_gathered$source=="gp"),], aes(y=value*accessConst), linetype="dotted",color=colours[2],size=1.2)+
    geom_line(data = Obj5b_gathered[which(Obj5b_gathered$source=="gpsg"),], aes(y=value*accessConst), linetype="dotted",color=colours[3],size=1.2)+
    geom_line(data = Obj5b_gathered[which(Obj5b_gathered$source=="old"),], aes(y=value*accessConst), linetype="dotted",color=colours[1],size=1.2)+
    
    geom_text(aes(max(age),1.6,label = "36% predicted",vjust=1.3, hjust = "inward"), color="grey40")+
    geom_text(aes(max(age),4.8,label = "total incidence",vjust=1.3, hjust = "inward"), color="grey40")+
    
    aes(group=rev(source)) + 
    #scale_x_continuous(trans='log10') +
    scale_linetype_manual(name="Parameterization",values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    #facet_wrap(~site,scales = "free") +
    labs(title=paste0("Age pattern of incidence of clinical malaria (Idete predicted)"), x="Age", y="Episodes per person-year")+
    theme_light() + 
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 14),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=14,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=14,  colour="grey20",face="bold"),
          
          plot.title=element_text(size=14,face="bold", vjust=1),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  
  
  lossVector[6,1] <- LF$Old$obj_5b$LF$orig$total
  lossVector_logLH[6,1] <-  LF$Old$obj_5b$LF$logLH$total
  
  lossVector[6,2] <- LF$GP$obj_5b$LF$orig$total
  lossVector_logLH[6,2] <-  LF$GP$obj_5b$LF$logLH$total
  
  lossVector[6,3] <- LF$GPSG$obj_5b$LF$orig$total
  lossVector_logLH[6,3] <-  LF$GPSG$obj_5b$LF$logLH$total
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  #  OBJECTIVES 6: Pyrogenic Threshold
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  # --------------------------------------------------------------------
  # ACUTE EPISODES: DIELMO PYROGENIC THRESHOLD
  # (scenario 234)

  LF$GP$obj_6 <- OBJ_AgePatternThresholdClinicalAttack(fieldDataAll, fieldDataID, modelResultID, densBiasGP, GPrunDir)
  LF$GPSG$obj_6 <- OBJ_AgePatternThresholdClinicalAttack(fieldDataAll, fieldDataID, modelResultID, densBiasGPSG, GPSGrunDir)
  LF$Old$obj_6 <- OBJ_AgePatternThresholdClinicalAttack(fieldDataAll, fieldDataID, modelResultID, densBiasOld, OldrunDir)
  
  Obj6 = data.frame(age = LF$GP$obj_6$Data$FieldData$NdiopDielmo$age)
  Obj6$Ystar = LF$GP$obj_6$Data$FieldData$NdiopDielmo$Ystar
  
  Obj6$gp = LF$GP$obj_6$Data$ModelResults$NdiopDielmo$pyrogt
  Obj6$gpsg = LF$GPSG$obj_6$Data$ModelResults$NdiopDielmo$pyrogt
  Obj6$old = LF$Old$obj_6$Data$ModelResults$NdiopDielmo$pyrogt
  
  
  Obj6_gathered = gather(Obj6, key = "source", value = "value", -age)
  
  pdf(file="./OutputPlots/Model_fit_Objective_6.pdf",height=3.41,width=5.38)
  ggplot(Obj6_gathered,aes(x=age,y=value))+
    geom_point(data = Obj6_gathered[which(Obj6_gathered$source=="Ystar"),], aes(y=value),size=3,color="#D3D4D9")+
    geom_line(data = Obj6_gathered[which(Obj6_gathered$source!="Ystar"),], aes(color = source,linetype=source),size=1.2)+
    aes(group=rev(source)) + 
    scale_x_continuous(trans='log10') +
    scale_linetype_manual(name="Parameterization", values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    #facet_wrap(~site,scales = "free") +
    labs(title=paste0("Age pattern of threshold parasite density for clinical attacks"), x="Age", y="Parasite:leucocyte ratio")+
    theme_light() + 
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 14),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=12,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=12,  colour="grey20",face="bold"),
          
          plot.title=element_text(size=12,face="bold", vjust=1),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  lossVector[7,1] <- LF$Old$obj_6$LF$orig$total
  lossVector_logLH[7,1] <-  LF$Old$obj_6$LF$logLH$total
  
  lossVector[7,2] <- LF$GP$obj_6$LF$orig$total
  lossVector_logLH[7,2] <-  LF$GP$obj_6$LF$logLH$total
  
  lossVector[7,3] <- LF$GPSG$obj_6$LF$orig$total
  lossVector_logLH[7,3] <-  LF$GPSG$obj_6$LF$logLH$total
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  #  OBJECTIVES 7: Severe disease by prevalence
  #------------------------------------------------------------------------------
  #------------------------------------------------------------------------------
  # --------------------------------------------------------------------
  # SEVERE MALARIA (prevalence vs episodes)
  # (scenarios 501,502,503,504,505,506,507,508,509,510,511,512,514,515,516,517,518,519,520,521,522,523,524,525,526,527)

  LF$GP$obj_7 <- OBJ_SevereEpisodesVsPrevalence(fieldDataAll, fieldDataID, modelResultID,  GPrunDir)
  LF$GPSG$obj_7 <- OBJ_SevereEpisodesVsPrevalence(fieldDataAll, fieldDataID, modelResultID,  GPSGrunDir)
  LF$Old$obj_7 <- OBJ_SevereEpisodesVsPrevalence(fieldDataAll, fieldDataID, modelResultID,  OldrunDir)
  
  Obj7 = list()
  Obj7$data = data.frame(prev= LF$GP$obj_7$Data$FieldData$prevalence)
  Obj7$data$episodes  = LF$GP$obj_7$Data$FieldData$MarshSnow_episodes 
  Obj7$data$source = "data"
  
  Obj7$gp = data.frame(prev= LF$GP$obj_7$Data$ModelResults$prevalence)
  Obj7$gp$episodes = LF$GP$obj_7$Data$ModelResults$sevEpPer1000personYr
  Obj7$gp$source = "gp"
  
  Obj7$gpsg = data.frame(prev= LF$GPSG$obj_7$Data$ModelResults$prevalence)
  Obj7$gpsg$episodes = LF$GPSG$obj_7$Data$ModelResults$sevEpPer1000personYr
  Obj7$gpsg$source = "gpsg"
  
  Obj7$old = data.frame(prev= LF$Old$obj_7$Data$ModelResults$prevalence)
  Obj7$old$episodes = LF$Old$obj_7$Data$ModelResults$sevEpPer1000personYr
  Obj7$old$source = "old"
  
  Obj_7_gathered =  Obj7$data
  for(i in 2:length(Obj7)){Obj_7_gathered = rbind(Obj_7_gathered,Obj7[[i]])}
  
  
  pdf(file="./OutputPlots/Model_fit_Objective_7.pdf",width=6.49,height=5.17)
  ggplot(Obj_7_gathered,aes(x=prev,y=episodes))+
    geom_line(data = Obj_7_gathered[which(Obj_7_gathered$source=="data"),], aes(y=episodes),size=1.2,color="#D3D4D9")+
    geom_point(data = Obj_7_gathered[which(Obj_7_gathered$source=="data"),], aes(y=episodes),size=2,shape=16)+
    geom_point(data = Obj_7_gathered[which(Obj_7_gathered$source!="data"),], aes(color = source,shape=source),size=4)+
    
    aes(group=rev(source)) + 
    #scale_x_continuous(trans='log10') +
    #scale_linetype_manual(values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP", gpsg="GPSG-BO", old="GA-O")) + 
    scale_shape_manual(name="Parameterization",values=c(data = 16,gp=15,gpsg=17,old=3), labels = c(data="Data",gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    scale_colour_manual(name="Parameterization",values=c(data = "#D3D4D9",gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(data="Data",gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "",shape="") + 
    #facet_wrap(~site,scales = "free") +
    labs(title=paste0(" Hospitalization rate in relation to prevalence in children (severe episodes)"), x="prevalence (0-9 year-olds)", y="Episodes / 1000 person years (0-9 year-olds)")+
    theme_light() + 
    theme(legend.position = c(0.99,0.86),
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 12),
          legend.text = element_text(size = 10),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=12,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=12,  colour="grey20",face="bold"),
          
          plot.title=element_text(face="bold", vjust=1, 
                                  size=12),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  lossVector[8,1] <- LF$Old$obj_7$LF$orig$total
  lossVector_logLH[8,1] <-  LF$Old$obj_7$LF$logLH$total
  
  lossVector[8,2] <- LF$GP$obj_7$LF$orig$total
  lossVector_logLH[8,2] <-  LF$GP$obj_7$LF$logLH$total
  
  lossVector[8,3] <- LF$GPSG$obj_7$LF$orig$total
  lossVector_logLH[8,3] <-  LF$GPSG$obj_7$LF$logLH$total
  # ---------------------------------------------------------------------
  # SEVERE MALARIA RR (relative risk)
  # (scenarios 158,167,173,176)
  LF$GP$obj_8 <- OBJ_AgePatternOfSevere(fieldDataAll, fieldDataID, modelResultID,  GPrunDir)
  LF$GPSG$obj_8 <- OBJ_AgePatternOfSevere(fieldDataAll, fieldDataID, modelResultID,  GPSGrunDir)
  LF$Old$obj_8 <- OBJ_AgePatternOfSevere(fieldDataAll, fieldDataID, modelResultID,  OldrunDir)
  
 
  Obj8_ls=list()
  Obj8_ls$Sukuta = data.frame(age = LF$Old$obj_8$Data$FieldData$Sukuta$age)
  Obj8_ls$KilifiNorth = data.frame(age = LF$Old$obj_8$Data$FieldData[["Kilifi-north"]]$age)
  Obj8_ls$KilifiSouth = data.frame(age = LF$Old$obj_8$Data$FieldData[["Kilifi-south"]]$age)
  Obj8_ls$Siaya = data.frame(age = LF$Old$obj_8$Data$FieldData$Siaya$age)
  
  Obj8_ls$Sukuta$Inc = LF$Old$obj_8$Data$FieldData$Sukuta$rate_epPer1000perYr
  Obj8_ls$KilifiNorth$Inc = LF$Old$obj_8$Data$FieldData[["Kilifi-north"]]$rate_epPer1000perYr
  Obj8_ls$KilifiSouth$Inc = LF$Old$obj_8$Data$FieldData[["Kilifi-south"]]$rate_epPer1000perYr
  Obj8_ls$Siaya$Inc = LF$Old$obj_8$Data$FieldData$Siaya$rate_epPer1000perYr
  
  #ModelPredictions
  
  Obj8_ls$Sukuta$gp = LF$GP$obj_8$Data$ModelResults$Sukuta$sevEpPer1000personYr[1:3]
  Obj8_ls$Sukuta$gpsg = LF$GPSG$obj_8$Data$ModelResults$Sukuta$sevEpPer1000personYr[1:3]
  Obj8_ls$Sukuta$old = LF$Old$obj_8$Data$ModelResults$Sukuta$sevEpPer1000personYr[1:3]
  Obj8_ls$Sukuta$site = "Sukuta"
  
  Obj8_ls$KilifiNorth$gp = LF$GP$obj_8$Data$ModelResults[["Kilifi-north"]]$sevEpPer1000personYr[1:3]
  Obj8_ls$KilifiNorth$gpsg = LF$GPSG$obj_8$Data$ModelResults[["Kilifi-north"]]$sevEpPer1000personYr[1:3]
  Obj8_ls$KilifiNorth$old = LF$Old$obj_8$Data$ModelResults[["Kilifi-north"]]$sevEpPer1000personYr[1:3]
  Obj8_ls$KilifiNorth$site = "KilifiNorth"
  
  Obj8_ls$KilifiSouth$gp = LF$GP$obj_8$Data$ModelResults[["Kilifi-north"]]$sevEpPer1000personYr[1:3]
  Obj8_ls$KilifiSouth$gpsg = LF$GPSG$obj_8$Data$ModelResults[["Kilifi-north"]]$sevEpPer1000personYr[1:3]
  Obj8_ls$KilifiSouth$old = LF$Old$obj_8$Data$ModelResults[["Kilifi-north"]]$sevEpPer1000personYr[1:3]
  Obj8_ls$KilifiSouth$site = "KilifiSouth"
  
  Obj8_ls$Siaya$gp = LF$GP$obj_8$Data$ModelResults$Siaya$sevEpPer1000personYr[1:3]
  Obj8_ls$Siaya$gpsg = LF$GPSG$obj_8$Data$ModelResults$Siaya$sevEpPer1000personYr[1:3]
  Obj8_ls$Siaya$old = LF$Old$obj_8$Data$ModelResults$Siaya$sevEpPer1000personYr[1:3]
  Obj8_ls$Siaya$site = "Siaya"

  
  Obj8 = data.frame(Obj8_ls[[1]])
  for(i in 2: length(LF$GP$obj_8$Data$FieldData)){
    Obj8 = rbind(Obj8,Obj8_ls[[i]])
  }
  #Obj 2 plot:prevalence
  Obj8_gathered = as.data.frame(gather(Obj8, key = "source", value = "value", -c(age,site)))
  pdf(file="./OutputPlots/Model_fit_Objective_8.pdf",height=5.06,width=7.50)
  ggplot(Obj8_gathered,aes(x=age,y=value))+
    geom_point(data = Obj8_gathered[which(Obj8_gathered$source=="Inc"),], aes(y=value),size=2.8,color="#D3D4D9")+
    geom_line(data = Obj8_gathered[which(Obj8_gathered$source!="Inc"),], aes(color = source,linetype=source),size=1.2)+
    aes(group=rev(source)) + 
    scale_linetype_manual(name="Parameterization",values = c(gp="solid",gpsg="solid",old="dashed"), labels = c(gp="GP-BO", gpsg="GPSG-BO", old="GA-O")) + 
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",linetype = "") + 
    #scale_x_continuous(trans='log10') +
    scale_x_continuous(breaks=c(0.5,3,7),
                      labels=c("1-11 months","1-4 years","5-9 years")) + 
    facet_wrap(~site, ncol = 2) +
    labs(title=paste0("Age pattern of hospitalisation: severe malaria"), x="Age", y="Episodes per 1000 person-years")+
    theme_light() + 
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.title = element_text(face="bold",size = 12),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  colour="grey30",angle=45, hjust = 1),
          
          axis.title.y=element_text(size=12,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=12,  colour="grey20",face="bold"),
          
          plot.title=element_text(size=14,face="bold", vjust=1),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off() #,labels=c("1-11 months","1-4 years","5-9 years")
  
  lossVector[9,1] <- LF$Old$obj_8$LF$orig$total
  lossVector_logLH[9,1] <-  LF$Old$obj_8$LF$logLH$total
  
  lossVector[9,2] <- LF$GP$obj_8$LF$orig$total
  lossVector_logLH[9,2] <-  LF$GP$obj_8$LF$logLH$total
  
  lossVector[9,3] <- LF$GPSG$obj_8$LF$orig$total
  lossVector_logLH[9,3] <-  LF$GPSG$obj_8$LF$logLH$total
  # ---------------------------------------------------------------------
  # DIRECT MORTALITY
  # (scenarios 301,302,303,312,316,317,318,326,327)

  LF$GP$obj_9 <- OBJ_DirectMalariaMortality(fieldDataAll, fieldDataID, modelResultID,  GPrunDir)
  LF$GPSG$obj_9 <- OBJ_DirectMalariaMortality(fieldDataAll, fieldDataID, modelResultID,  GPSGrunDir)
  LF$Old$obj_9 <- OBJ_DirectMalariaMortality(fieldDataAll, fieldDataID, modelResultID,  OldrunDir)
  
  Obj9 = data.frame(site = names(LF$GP$obj_9$Data$FieldData$data))
  Obj9$EIR = NA
  for(i in 1:dim(Obj9)[1]){ Obj9$EIR[i] = LF$GP$obj_9$Data$FieldData$data[[i]]$EIR}
  Obj9$dirMort = NA
  for(i in 1:dim( Obj9)[1]){ Obj9$dirMort[i] = LF$GP$obj_9$Data$FieldData$data[[i]]$dirDeathRate_epPer1000perYr}
  Obj9$uci = NA
  for(i in 1:dim( Obj9)[1]){ Obj9$uci[i] = LF$GP$obj_9$Data$FieldData$data[[i]]$dirDeathRate_epPer1000perYr_UpperCI}
  Obj9$lci = NA
  for(i in 1:dim( Obj9)[1]){ Obj9$lci[i] = LF$GP$obj_9$Data$FieldData$data[[i]]$dirDeathRate_epPer1000perYr_LowerCI}
 
   Obj9$gp = NA
  for(i in 1:dim( Obj9)[1]){ Obj9$gp[i] = LF$GP$obj_9$Data$ModelResults[[i]]$dirDeathsPer1000personYr}
   Obj9$gpsg = NA
   for(i in 1:dim( Obj9)[1]){ Obj9$gpsg[i] = LF$GPSG$obj_9$Data$ModelResults[[i]]$dirDeathsPer1000personYr}
   Obj9$old = NA
   for(i in 1:dim( Obj9)[1]){ Obj9$old[i] = LF$Old$obj_9$Data$ModelResults[[i]]$dirDeathsPer1000personYr}
   
  
  Obj9_gathered = gather(Obj9, key = "source", value = "value", -c(EIR,site,lci,uci))
  
  pdf(file="./OutputPlots/Model_fit_Objective_9.pdf",height=4.04,width=5.44)
  ggplot(Obj9_gathered,aes(x=EIR,y=value))+
    scale_x_continuous(trans='log10') +
    geom_point(data = Obj9_gathered[which(Obj9_gathered$source=="dirMort"),], aes(y=value),size=3,color="#D3D4D9")+
    geom_errorbar(aes(ymin=lci, ymax=uci), width=.02,color="grey70") +
    geom_point(data = Obj9_gathered[which(Obj9_gathered$source!="dirMort"),], aes(color = source,shape=source),size=4)+
    aes(group=rev(source)) + 
    #scale_x_continuous(trans='log10') +
    scale_shape_manual(name="Parameterization",values=c(gp=15,gpsg=17,old=3), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",shape = "") + 
    #facet_wrap(~site,scales = "free") +
    labs(title=paste0("Direct Mortality in children <5 years"), x="EIR", y="Deaths per 1000 person-years")+
    theme_light() + 
    theme(legend.position = c(0.01,0.82),
          legend.justification = "left",
          legend.title = element_text(face="bold",size = 11),
          legend.text = element_text(size = 10),
          
          axis.text.y=element_text(size=14,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=14,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=14,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=14,  colour="grey20",face="bold"),
          
          plot.title=element_text(size=14,face="bold", vjust=1),
          #text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  
  lossVector[10,1] <- LF$Old$obj_9$LF$orig$total
  lossVector_logLH[10,1] <-  LF$Old$obj_9$LF$logLH$total
  
  lossVector[10,2] <- LF$GP$obj_9$LF$orig$total
  lossVector_logLH[10,2] <-  LF$GP$obj_9$LF$logLH$total
  
  lossVector[10,3] <- LF$GPSG$obj_9$LF$orig$total
  lossVector_logLH[10,3] <-  LF$GPSG$obj_9$LF$logLH$total
  # ---------------------------------------------------------------------
  # INDIRECT MORTALITY
  # (scenarios 401,402,408,411,414,415,416,417,418,422,426)
  LF$GP$obj_10 <- OBJ_indirectMortality(fieldDataAll, fieldDataID, modelResultID,  GPrunDir)
  LF$GPSG$obj_10 <- OBJ_indirectMortality(fieldDataAll, fieldDataID, modelResultID,  GPSGrunDir)
  LF$Old$obj_10 <- OBJ_indirectMortality(fieldDataAll, fieldDataID, modelResultID,  OldrunDir)
  
  Obj10 = data.frame(site = names(LF$GP$obj_10$Data$FieldData$data))
  Obj10$EIR = NA
  for(i in 1:dim(Obj10)[1]){ Obj10$EIR[i] = LF$GP$obj_10$Data$FieldData$data[[i]]$EIR}
  Obj10$indirMort = NA
  for(i in 1:dim( Obj10)[1]){ Obj10$indirMort[i] = LF$GP$obj_10$Data$FieldData$data[[i]]$indirDeath}
  
  Obj10$gp = NA
  for(i in 1:dim( Obj10)[1]){ Obj10$gp[i] = LF$GP$obj_10$Data$ModelResults[[i]]$allCauseIMR}
  Obj10$gpsg = NA
  for(i in 1:dim( Obj10)[1]){ Obj10$gpsg[i] = LF$GPSG$obj_10$Data$ModelResults[[i]]$allCauseIMR}
  Obj10$old = NA
  for(i in 1:dim( Obj10)[1]){ Obj10$old[i] = LF$Old$obj_10$Data$ModelResults[[i]]$allCauseIMR}
  
  
  Obj10_gathered = gather(Obj10, key = "source", value = "value", -c(EIR,site))
  

  
  pdf(file="./OutputPlots/Model_fit_Objective_10.pdf",height=4.04,width=5.44)
  ggplot(Obj10_gathered,aes(x=EIR,y=value))+
    scale_x_continuous(trans='log10') +
    geom_point(data = Obj10_gathered[which(Obj10_gathered$source=="indirMort"),], aes(y=value),size=4,color="#D3D4D9")+
    geom_point(data = Obj10_gathered[which(Obj10_gathered$source!="indirMort"),], aes(color = source,shape=source),size=4)+
    aes(group=rev(source)) + 
    scale_shape_manual(name="Parameterization",values=c(gp=15,gpsg=17,old=3), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    scale_colour_manual(name="Parameterization",values=c(gp=colours[2],gpsg=colours[3],old=colours[1]), labels = c(gp="GP-BO",gpsg= "GPSG-BO",old= "GA-O")) +
    labs(colour = "",shape = "") + 
    #facet_wrap(~site,scales = "free") +
    labs(title=paste0("All-cause infant mortality rate"), x="EIR", y="Deaths per 1000 livebirths")+
    theme_light() + 
    theme(legend.position = "bottom",
          legend.justification = "right",
          legend.title = element_text(size = 14,face="bold"),
          legend.text = element_text(size = 12),
          
          axis.text.y=element_text(size=12,  hjust=1,colour="grey30"),
          axis.text.x=element_text(size=12,  vjust=1,colour="grey30"),
          
          axis.title.y=element_text(size=14,  colour="grey20",face="bold"),
          axis.title.x=element_text(size=14,  colour="grey20",face="bold"),
          
          plot.title=element_text(face="bold", vjust=1),
  #        text=element_text(family="URWHelvetica"),
          strip.text.x = element_text(size = 12, color = "white"),
          strip.background = element_rect(color = "grey70",fill="grey70"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )  
  dev.off()
  lossVector[11,1] <- LF$Old$obj_10$LF$orig$total
  lossVector_logLH[11,1] <-  LF$Old$obj_10$LF$logLH$total
  
  lossVector[11,2] <- LF$GP$obj_10$LF$orig$total
  lossVector_logLH[11,2] <-  LF$GP$obj_10$LF$logLH$total
  
  lossVector[11,3] <- LF$GPSG$obj_10$LF$orig$total
  lossVector_logLH[11,3] <-  LF$GPSG$obj_10$LF$logLH$total

  
  # ---------------------------------------------------------------------
  # LossFunction weights # TODO this is hard coded at the moment
  lossFunc_weights <- array(c(0.001, 0.001, 0.01, 0.01, 1, 1, 1, 2, 2, 1, 10),dim=c(11,1))
  # ---------------------------------------------------------------------
  # WEIGHTED LOSS FUNCTION
  weighted_LF <- list()
  weighted_LF$vector <- lossFunc_weights*as.data.frame(lossVector)
  weighted_LF$vector_new <- lossFunc_weights*as.data.frame(lossVector_logLH)
  weighted_LF$total <- colSums(weighted_LF$vector) 
  weighted_LF$total_avgOverScen <- weighted_LF$total/61
  
  ## Printing of lf into table in pdf 
  tableToPrint <- array(0,dim=c(11,7))
  tableToPrint[,c(1:3)] <- lossVector
  tableToPrint[,c(4)] <- as.numeric(lossFunc_weights)
  tableToPrint[,c(5:7)] <- as.matrix(weighted_LF$vector)

  # not the most elegant way to print to pdf but is quick fix
  plot.new()
  limUpper_y <- 50
  plot.window(xlim=c(0,1),ylim=c(0,limUpper_y)) 
  title(paste("Summary of loss functions",sep=" "),cex.main=1.3) 
  
  rowNames_descObj<-c("Age prevalence  \nafter intervention", "Age pattern \nof prevalence", "Age pattern \nof parasite density", "Multiplicity of \nInfection", "Age pattern of clinical \nincidence : Senegal", "Age pattern of clinical \nincidence : Idete", "Clinical Threshold", "Severe episodes \nwith prevalence", "Age pattern of \nsevere episodes", "Direct mortality \nwith EIR", "All-cause mortality \nwith EIR")
  rowNames_obj<-c("obj_1", "obj_2", "obj_3", "obj_4", "obj_5a", "obj_5b", "obj_6", "obj_7", "obj_8", "obj_9", "obj_10")
  rowNames_losVec<-c("lossVec_1", "lossVec_2", "lossVec_3", "lossVec_4", "lossVec_5", "lossVec_6", " lossVec_7", "lossVec_8", "lossVec_9", "lossVec_10", "lossVec_11")
  
  text(0.0,(limUpper_y-1), adj=c(0,0),lab=paste("GA total weighted loss function     ", as.character(format(weighted_LF$total_avgOverScen[1], digits = 4, nsmall = 5)),"(averaged)",sep="  "))
  text(0.0,(limUpper_y-4), adj=c(0,0),lab=paste("GP total weighted loss function     ", as.character(format(weighted_LF$total_avgOverScen[2], digits = 4, nsmall = 5)),"(averaged)",sep="  "))
  text(0.0,(limUpper_y-7), adj=c(0,0),lab=paste("GPSG total weighted loss function     ", as.character(format(weighted_LF$total_avgOverScen[3], digits = 4, nsmall = 5)),"(averaged)",sep="  "))
  
  text(0.0,(limUpper_y-14), adj=c(0,0),lab=paste("GA total weighted loss functions  ", as.character(format(weighted_LF$total[1], digits = 4, nsmall = 5)),sep=" "))
  text(0.0,(limUpper_y-17), adj=c(0,0),lab=paste("GP total weighted loss functions  ", as.character(format(weighted_LF$total[2], digits = 4, nsmall = 5)),sep=" "))
  text(0.0,(limUpper_y-20), adj=c(0,0),lab=paste("GPSG total weighted loss functions  ", as.character(format(weighted_LF$total[3], digits = 4, nsmall = 5)),sep=" "))
  
  plot.new()
  tableGrid <- data.frame(rowNames_descObj,round(lossVector, digits=6),as.character(lossFunc_weights),round(weighted_LF$vector, digits=4),round(weighted_LF$vector_new, digits=4))
  colnames(tableGrid)=c( "objective", "GA (unweighted)","GP (unweighted)","GPSG (unweighted)","weights","GA (weighted)","GP (weighted)","GP(unweighted)")
  grid.table(tableGrid)
  
  dev.off()
  write.csv(tableGrid,paste0("Fit_summary.csv"))



