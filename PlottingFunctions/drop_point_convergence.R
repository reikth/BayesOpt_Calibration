#Convergence compared to drop point

u = "reiker" #Scicore user name
GitDir = "/scicore/home/pothin/reiker/GitRepos/om_fitting/" # Local version of the Git repository
ExperimentDir = "/scicore/home/pothin/reiker/Paper_3_Model_Fitting/" #Parent folder to experiments

ExperimentNameHetGP = "2020_05_full_fit_small_hetGP_norm_sampling"
ExperimentNameGPSG = "2020_05_full_fit_small_GPSG_norm_sampling"
ExperimentNameDropPoint = "2020_10_drop_point"


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

############################################################
#
# functions
#
############################################################

D = 23 #dimensions
W = 11 #output objectives

setwd(local.directory)


############################################################
#
# Data Prep
#
############################################################
mds = read.table("master.data.small.train")
mdt = read.table("master.data.small.test")
md = rbind(mds,mdt)
no_update_fies = system("ls update.data.small* |wc -l",intern=TRUE)
for (i in 1:no_update_fies){
  ud = read.table(paste0("update.data.small",i))
  md=rbind(md,ud)
}
rm(mds,mdt,ud,no_update_fies)

weights.w = c(1,1,1,1,1,1,1,1,1,1,1) #weights are already applied in likelihood.R!
#Otherwise: weights.w = c(0.001,0.001,0.01,0.01,1,1,1,2,2,1,10)

md$OBJECTIVE12 = numeric(dim(md)[1])

for (w in 1:W) {md$OBJECTIVE12 <- md$OBJECTIVE12 + (weights.w[w]*exp(md[,((2*D)+w)]))} # for each parameter vector, compute



############################################################
#
# Parameter and objective mapping
#
############################################################


colnames <- list()
colnames[[length(colnames)+1]] <- paste0("elapsed")
for (w in 1:W) {colnames[[length(colnames)+1]] <- paste0("r2_",sprintf("%02i",w))}
colnames[[length(colnames)+1]] <- paste0("r2_w")
colnames[[length(colnames)+1]] <- paste0("min_lcb")
colnames[[length(colnames)+1]] <- paste0("mpe")

colnames <- unlist(colnames)

############################################################
#
# Data Prep
#
############################################################
#hetgp dataframe
load(paste0('time.hetgp'))
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
#manual fix because of interrupted runs: 
if(ExperimentNameHetGP == "2020_05_full_fit_small_hetGP_norm_sampling"){
  hetgp_best  = c(hetgp_best[c(1:11)], hetgp_best[11],hetgp_best[11],hetgp_best[c(12:length(hetgp_best))],hetgp_best[length(hetgp_best)]) 
}

time_frame_hetGP$best_known=NA
time_frame_hetGP$best_known=hetgp_best

for(i in 2:dim(time_frame_hetGP)[1]){
  if(time_frame_hetGP$elapsed[i]<time_frame_hetGP$elapsed[i-1]){
    base = time_frame_hetGP$elapsed[i-1]
    for(j in (dim(time_frame_hetGP)[1]-1):i){time_frame_hetGP$best_known[j+1] = time_frame_hetGP$best_known[j]}
    time_frame_hetGP$best_known[i] = time_frame_hetGP$best_known[i-1]
    for(j in i:dim(time_frame_hetGP)[1]){time_frame_hetGP$elapsed[j] = time_frame_hetGP$elapsed[j] + base }
    break()
  }
}



#gpsg dataframe
load(paste0(ExperimentDir,ExperimentNameGPSG,'/time.gpsg'))
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

load(paste0(ExperimentDir,ExperimentNameGPSG,'/eval.state.newer'))
gpsg_best = vector(length=length(evaluations$best))
for (i in 1:length(evaluations$best)){
  gpsg_best[[i]]=sum(evaluations$best[[i]]$known$Objectives)
  if(gpsg_best[[i]]==0 | (i>1 && gpsg_best[[i]]>gpsg_best[[i-1]])){gpsg_best[[i]]=gpsg_best[[i-1]]}
  
}

#manual fix because of interrupted runs: 
if(ExperimentNameGPSG == "2020_05_full_fit_small_GPSG_norm_sampling"){
  gpsg_best  = c(gpsg_best[c(1:length(gpsg_best))],rep(gpsg_best[length(gpsg_best)],3)) 
}

time_frame_gpsg$best_known=NA
time_frame_gpsg$best_known=gpsg_best[1:length(time_frame_gpsg$best_known)]
for(i in 2:dim(time_frame_gpsg)[1]){
  if(time_frame_gpsg$elapsed[i]<time_frame_gpsg$elapsed[i-1]){
    base = time_frame_gpsg$elapsed[i-1]
    for(j in i:dim(time_frame_gpsg)[1]){time_frame_gpsg$elapsed[j] = time_frame_gpsg$elapsed[j] + base }
    break()
  }
}

##Frame drop point 
load(paste0(ExperimentDir,ExperimentNameDropPoint,'/time.hetgp'))
time_hetgp_drop = time
time_frame_hetgp_drop <- matrix(0,nrow=length(time_hetgp_drop),(W+4))
for (g in 1:length(time_hetgp_drop)) {
  time_frame_hetgp_drop[g,1] <- as.numeric(time_hetgp_drop[[g]]$time["elapsed"])
  time_frame_hetgp_drop[g,c(2:(2+W))] <- unlist(time_hetgp_drop[[g]]$r2)
  time_frame_hetgp_drop[g,(W+3)] <- time_hetgp_drop[[g]]$obj$min_lcb
  time_frame_hetgp_drop[g,(W+4)] <- time_hetgp_drop[[g]]$obj$mpe
}
time_frame_hetgp_drop <- as.data.frame(time_frame_hetgp_drop)
colnames(time_frame_hetgp_drop) <- colnames
time_frame_hetgp_drop$Algorithm = "GP_drop"
time_frame_hetgp_drop$iteration = seq(1,dim(time_frame_hetgp_drop)[1],1)

load(paste0(ExperimentDir,ExperimentNameDropPoint,'/eval.state.newer'))
hetgp_drop_best = vector(length=length(evaluations$best))
for (i in 1:length(evaluations$best)){
  hetgp_drop_best[[i]]=sum(evaluations$best[[i]]$known$Objectives)
  if(hetgp_drop_best[[i]]==0 | (i>1 && hetgp_drop_best[[i]]>hetgp_drop_best[[i-1]])){hetgp_drop_best[[i]]=hetgp_drop_best[[i-1]]}
  
}



time_frame_hetgp_drop$best_known=NA
time_frame_hetgp_drop$best_known=hetgp_drop_best[1:length(time_frame_hetgp_drop$best_known)]
for(i in 2:dim(time_frame_hetgp_drop)[1]){
  if(time_frame_hetgp_drop$elapsed[i]<time_frame_hetgp_drop$elapsed[i-1]){
    base = time_frame_hetgp_drop$elapsed[i-1]
    for(j in i:dim(time_frame_hetgp_drop)[1]){time_frame_hetgp_drop$elapsed[j] = time_frame_hetgp_drop$elapsed[j] + base }
    break()
  }
}


time_frame_all = rbind(time_frame_hetGP,time_frame_gpsg,time_frame_hetgp_drop)
############################################################
#
# plotting
#
############################################################
colours = c("#FF7745", "#98D2EB", "#720026") #Sky blue: #98D2EB; #purple: #B2B1CF
pdf("./OutputPlots/Fig_2_A_covergence_plot_drop_point.pdf" , width= 5.07, height = 4.11)
ggplot(data=time_frame_all, aes(x=elapsed, y=best_known, group=Algorithm,color=Algorithm)) + 
  #geoms:
  #geom_ribbon(aes(ymin=(min_lcb), ymax=mpe,fill=Algorithm), alpha = 0.6, show.legend = T, colour = NA) + #,stat="stepribbon"
  geom_step(size=1.3) + #geom_step
  #geom_line(aes(linetype=Algorithm)) + #geom_step
  geom_point(aes(shape=Algorithm),size=3)+
  
  #previous
  geom_hline(yintercept=73, linetype="dashed", color="#D3D4D9", size=1) + #grey40
  geom_text(aes(max(time_frame_all$elapsed),73,label = 'satisfactory',vjust=1.3, hjust = "inward"), color="grey40")+
  geom_hline(yintercept=63.13401, linetype="dashed", color="#D3D4D9", size=1) + 
  geom_text(aes(max(time_frame_all$elapsed),63.13401,label = "previous fit",vjust=1.3, hjust = "inward"), color="grey40")+
  #axes:
  scale_x_continuous(labels = function(elapsed) format(elapsed, scientific = TRUE)) + 
  scale_y_continuous(limits =c(50,100)) + 
  #colors:
  scale_colour_manual(values=colours, labels=c("GP-BO", "GP-BO (drop)","GPSG-BO")) +
  scale_fill_manual(values=colours, labels=c("GP-BO", "GP-BO (drop)", "GPSG-BO")) +
  scale_shape_manual(values=c(16,18,17),labels=c("GP-BO", "GP-BO (drop)", "GPSG-BO")) +
  
  #scale_colour_discrete_qualitative(palette = "Harmonic",order = c(2:1),nmax=2) +
  #scale_fill_discrete_qualitative(palette = "Harmonic",order = c(2:1),nmax=2) +
  #Labels
  #labs(x = 'CPU time (seconds)', y = TeX("argmin $\\widehat{\\Y_{LF}}$ / LCB")) +
  labs(x = 'CPU time (seconds)', y = TeX("Current best (min $F(\\theta))$")) +
  theme_classic() + 
  
  theme(legend.position=c(0.55, 1.03),
        legend.direction="horizontal",
        legend.text=element_text(colour="grey20"),
        legend.title = element_text(),
        
        plot.margin=grid::unit(c(1.5,.5,.5,.5), "cm"),
        
        axis.text.y=element_text(size=12, family="Helvetica", hjust=1,colour="grey30"),
        axis.text.x=element_text(size=10, family="Helvetica", vjust=1,colour="grey30"),
        #axis.text.x=element_blank(),
        
        
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid=element_blank(),
        
        title=element_text(face="bold", vjust=1, 
                           family="Helvetica"),
        text=element_text(family="URWHelvetica"))

dev.off()
