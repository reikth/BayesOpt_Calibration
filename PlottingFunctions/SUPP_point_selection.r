########
#
#	Code for plotting selected points (supplement)
#	Theresa Reiker, 28.5.2020
#########

#Folder setup
require(pacman)
pacman::p_load(tidyverse,ggplot2)
u = "reiker" #Scicore user name
GitDir = "/scicore/home/smith/reiker/GitRepos/om_fitting/" # Local version of the Git repository
ExperimentDir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/" #Parent folder to experiments

ExperimentNameGP = "2020_05_full_fit_small_hetGP_norm_sampling"
ExperimentNameGPSG = "2020_05_full_fit_small_GPSG_norm_sampling"


#Plotting
local.directory = paste0(ExperimentDir,ExperimentNameGP,'/') #working directory
setwd(local.directory)
system("mkdir OutputPlots")
p = read.table("all_selected_points")

i = 1 # iteration
n = 250 #no of points per batch)

pdf("./OutputPlots/SUPP_selected_points_hetgp.pdf",width=9,height=9)
for(i in 1:(dim(p)[1]/n)){
print(i)
	x = ggplot(gather(p[c((n*(i-1)+1):(n*i)),c(1:23)]),aes(value))+
		geom_density(fill="#FF7745",aes(y=..count..))+
		facet_wrap(~key,scales ='free_x') +
		geom_rug(aes(x=value))+
		labs(title=paste0("GP-BO Iteration ",i))+
   xlim(c(-0.5,0.5)) + 
		theme_bw()
	print(x)
}
dev.off()

#GP
#iterations 1,10,20,30
selection = c(1,10,20,30)
i = selection[1]
df = gather(p[c((n*(i-1)+1):(n*i)),c(1:23)])
df$iteration = i

for(i in selection[c(1:length(selection))]){
  df_tmp = gather(p[c((n*(i-1)+1):(n*i)),c(1:23)])
  df_tmp$iteration = i
  df = rbind(df,df_tmp)
}

layered_points_gp = ggplot(df,aes(value))+
  geom_density(alpha=0.75,aes(fill=as.factor(iteration),color=as.factor(iteration),group=as.factor(iteration),y=..count..))+
  facet_wrap(~key,scales ='free_x') +
  scale_color_manual(values = c("#FAC7B7","#F58F70","#EF4E1C","#B3330C"))+
  scale_fill_manual(values = c("#FAC7B7","#F58F70","#EF4E1C","#B3330C"))+
  labs(title=paste0("GP-BO selected points"),color="iteration", fill="iteration")+
  xlim(c(-0.5,0.5)) + 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))




#GPSG
local.directory = paste0(ExperimentDir,ExperimentNameGPSG,'/') #working directory
setwd(local.directory)
system("mkdir OutputPlots")
p = read.table("all_selected_points")

i = 1 # iteration
n = 250 #no of points per batch)

pdf("./OutputPlots/SUPP_selected_points_gpsg.pdf",width=9,height=9)
for(i in 1:(dim(p)[1]/n)){
print(i)
	x = ggplot(gather(p[c((n*(i-1)+1):(n*i)),c(1:23)]),aes(value))+
		geom_density(fill="#720026",aes(y=..count..))+ #color="burlywood"
		facet_wrap(~key,scales ='free_x') +
		geom_rug(aes(x=value))+
		labs(title=paste0("GPSG-BO Iteration ",i))+
   xlim(c(-0.5,0.5)) + 
		theme_bw()
	print(x)
}
dev.off()




#GPSG
#iterations 1,10,20,23
selection = c(1,10,20,23)
i = selection[1]
df = gather(p[c((n*(i-1)+1):(n*i)),c(1:23)])
df$iteration = i

for(i in selection[c(1:length(selection))]){
  df_tmp = gather(p[c((n*(i-1)+1):(n*i)),c(1:23)])
  df_tmp$iteration = i
  df = rbind(df,df_tmp)
}

layered_points_gpsg = ggplot(df,aes(value))+
  geom_density(alpha=0.75,aes(fill=as.factor(iteration),color=as.factor(iteration),group=as.factor(iteration),y=..count..))+
  facet_wrap(~key,scales ='free_x') +
  scale_color_manual(values = c("#F6A7C1","#E92365","#9B0F3E","#610A27"))+
  scale_fill_manual(values = c("#F6A7C1","#E92365","#9B0F3E","#610A27"))+
  labs(title=paste0("GPSG-BO selected points"),color="iteration", fill="iteration")+
  xlim(c(-0.5,0.5)) + 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


pdf("./OutputPlots/SUPP_selected_points_hetgp_layers.pdf",width=11,height=9)
print(layered_points_gp)
dev.off()

pdf("./OutputPlots/SUPP_selected_points_gpsg_layers.pdf",width=11,height=9)
print(layered_points_gpsg)
dev.off()


