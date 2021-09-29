GitDir =  "/scicore/home/smith/reiker/GitRepos/om_fitting/"
MainExperimentDir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_hetGP_norm_sampling/"
old_param_output_dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_04_full_fit_small_hetGP/PriorDraws/Output/old_params/"
GP_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_hetGP_norm_sampling/"
GPSG_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_GPSG_norm_sampling/"

if(!require(pacman)){install.packages(pacman);require(pacman)}else{require(pacman)}
pacman::p_load(gridExtra,ggplot2,colorspace,tidyverse, lhs,stringr)

D=23;W=11
setwd(MainExperimentDir)
source(paste0(GitDir,"EvalOM.parallel.R"))

old_params = c(
  0.050736,
  0.03247,
  0.138161050830301,
  1514.385853233699891,
  2.03692533424484,
  10.173598698525799,
  35158523.31132510304451,
  97.334652723897705,
  2.33031045876193,
  2.53106547375805,
  0.655747311168152,
  0.916181104713054,
  6502.26335600001039,
  142601.912520000012591,
  0.177378570987455,
  0.05,
  0.736202,
  0.018777338,
  49.539046599999999,
  4.79610772546704,
  784455.599999999976717,
  1,
  0,
  0.0968,
  0.275437402,
  0.596539864,
  0,
  296.302437899999973,
  2.797523626,
  0.117383
)



sample_size = 100000

theta = randomLHS(sample_size,23)
priors_theta = apply(theta, 1,quant.to.prior) %>% unlist()


all_param_names = names(priors_theta) %>% unique() %>% as.data.frame() %>% filter(!str_detect(.,'quantile|paramNumber')) %>% pull(.)

not_fitted <- c ("sigmasq",
                 "ImmunityPenalty",
                 "ImmuneEffectorDecay",
                 "AsexualImmunityDecay",
                 "IdeteMultiplier",
                 "Estar",
                 "minuslogcompSinf")
param_names=all_param_names[!all_param_names %in% not_fitted]

param_values=data.frame(matrix(nrow=sample_size,ncol=length(param_names))) 
names(param_values) = param_names

for(i in 1:length(param_names)){
  param_values[,i] = as.numeric(priors_theta[which(names(priors_theta)==paste0(param_names[i]))])
}
param_values_gathered = gather(param_values)

param_labeller =  as.character(seq(1,23,1))
names(param_labeller)=as.character(param_names)

#create df with Col1 = "key", col2 = "GA", col3 = "GP", col2 = "GPSG"
#1. GP
load("eval.state.newer")

df3 <- param_values_gathered %>%
  group_by(key) %>%
  dplyr::summarise(Mean.SL = mean(value)) %>% as.data.frame()


for(h in 1:length(evaluations$best)){
  df3$GP = NA
    if(sum(evaluations$best[[h]]$known$Objectives) >0 ){
      best_gp = quant.to.prior(evaluations$best[[h]]$known$parameters) %>% unlist()
    }else{}


  for(i in 1:length(param_names)){
    df3[which(df3$key==paste0(param_names[i])),"GP"] = as.numeric(best_gp[which(names(best_gp)==paste0(param_names[i]))])
  }
  names(df3)[names(df3)=="GP"] =paste0("GP",h)
   
}


#2. GPSG
load(paste0(GPSG_Dir,"eval.state.newer"))
best_gpsg = quant.to.prior(evaluations$best[[length(evaluations$best)]]$known$parameters) %>% unlist()

#df3$GPSG = NA
#for(i in 1:length(param_names)){
#  df3[which(df3$key==paste0(param_names[i])),"GPSG"] = as.numeric(best_gpsg[which(names(best_gpsg)==paste0(param_names[i]))])
#}


for(h in 1:length(evaluations$best)){
  df3$GPSG = NA
  if(sum(evaluations$best[[h]]$known$Objectives) >0 ){
    best_gpsg = quant.to.prior(evaluations$best[[h]]$known$parameters) %>% unlist()
  }else{}
  
  
  for(i in 1:length(param_names)){
    df3[which(df3$key==paste0(param_names[i])),"GPSG"] = as.numeric(best_gpsg[which(names(best_gpsg)==paste0(param_names[i]))])
  }
  names(df3)[names(df3)=="GPSG"] =paste0("GPSG",h)
  
}

#3. Old params
#extract parameter order
param_order = quant.to.prior(theta[1,]) %>% unlist() %>% as.data.frame()
param_order$name = rownames(param_order)
rownames(param_order)= NULL
param_order=param_order[- grep("quantile", param_order$name),]
param_order_numbers = param_order[grep("paramNumber", param_order$name),]
param_order_names= param_order[-grep("paramNumber", param_order$name),]
param_order = cbind(param_order_numbers,param_order_names)
param_order = param_order[,c(1,4)]
names(param_order) = c("number","key")
rm(param_order_names,param_order_numbers)

old = param_order
old$GA=NA
for(i in 1:length(old_params)){
  n = old$number[i]
  old$GA[i] = old_params[n]
}
old = old[,c(2,3)]
df3 = merge(df3,old,by="key")


#Plot 

ggplot(param_values_gathered, aes(value)) + 
  #geom_histogram(aes(y = ..density.., fill = ..count..), bins = 50) +
  geom_density(color="darkblue", fill="lightblue")+ 
  facet_wrap(~factor(key,levels=unique(as.factor(param_values_gathered$key))), scales = "free", labeller = as_labeller(param_labeller)) +
  ggtitle("Prior parameter distributions") + 
  theme_bw()

logplot = ggplot(param_values_gathered, aes(value)) + 
  #geom_histogram(aes(y = ..density.., fill = ..count..), bins = 50) +
  geom_density(alpha=0.4,color=NA, fill="#D3D4D9")+ 
  scale_x_continuous(trans = 'log10') + 
  facet_wrap(~factor(key,levels=unique(as.factor(param_values_gathered$key))), scales = "free", labeller = as_labeller(param_labeller)) +
  ggtitle("Log prior parameter distributions and posterior estimates") 

GPevolution=logplot
GPSGevolution=logplot

logplot=logplot + 
  geom_vline(data=df3, aes(xintercept= GP27, colour="GP27"),size=1,show.legend=T) + 
  geom_vline(data=df3, aes(xintercept= GPSG21, colour="GPSG21"),size=1,show.legend=T)  + 
  geom_vline(data=df3, aes(xintercept= GA, colour="GA"), linetype="dashed",show.legend=T) + 
  scale_color_manual(name = "Algorithm", values = c(GP27 = "#FF7745", GPSG21 = "#720026", GA="black")) + 
  theme_classic()
all = logplot

indx <- grepl('GPSG', colnames(df3))
gpsgnames = names(df3)[indx]
indx <- grepl('GP[0-9]', colnames(df3))
gpnames = names(df3)[indx]



gpalphas = seq(from=0.1,to=1,length.out=length(gpnames))
gpsgalphas = seq(from=0.1,to=1,length.out=length(gpsgnames)) #alternative: lseq
for(g in 1:length(gpnames)){
  
  GPevolution=GPevolution + 
    geom_vline(data=df3, aes_string(xintercept= gpnames[g], color = shQuote("GP")),alpha=gpalphas[g],size=1,show.legend=T) + 
    geom_vline(data=df3, aes(xintercept= GA, colour="GA"),linetype="dashed",show.legend=T) + 

    scale_color_manual(name = "Algorithm", values = c(GP = "#FFAB8C", GA="black")) + 
    theme_classic()
  
  if(g==length(gpnames)){GPevolution = GPevolution + 
     geom_vline(data=df3, aes_string(xintercept= gpnames[g]),color="#FFCDBA", alpha=gpalphas[g],size=1,show.legend=F)
    }
}

#GPSG
for(g in 1:length(gpsgnames)){
  
  GPSGevolution=GPSGevolution + 
    geom_vline(data=df3, aes_string(xintercept= gpsgnames[g], color = shQuote("GPSG")),alpha=gpsgalphas[g],size=1,show.legend=T) + 
    geom_vline(data=df3, aes(xintercept= GA, colour="GA"),linetype="dashed",show.legend=T) + 
    
    scale_color_manual(name = "Algorithm", values = c(GPSG = "#FF9CBD", GA="black")) + 
    theme_classic()
  
  if(g==length(gpsgnames)){GPSGevolution = GPSGevolution + 
    geom_vline(data=df3, aes_string(xintercept= gpsgnames[g]),color="#720026", alpha=gpsgalphas[g],size=1,show.legend=F)
  }
}



pdf(paste0("OutputPlots/Log_prior_distributions.pdf"), width = 11.89, height = 7.61)

print(all)
print(GPevolution)
print(GPSGevolution)

dev.off()
