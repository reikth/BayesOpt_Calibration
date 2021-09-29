GitDir =  "/scicore/home/smith/reiker/GitRepos/om_fitting/"
MainExperimentDir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_hetGP_norm_sampling/"
GP_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_hetGP_norm_sampling/"
GPValidation_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_07_hetGP_synth_data_validation/"
GPSG_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_05_full_fit_small_GPSG_norm_sampling/"
GPSGValidation_Dir = "/scicore/home/smith/reiker/Paper_3_Model_Fitting/2020_07_GPSG_synth_data_validation/"

if(!require(pacman)){install.packages(pacman);require(pacman)}else{require(pacman)}
pacman::p_load(gridExtra,ggplot2,colorspace,tidyverse, lhs,stringr)

D=23;W=11;
source(paste0(GitDir,"EvalOM.parallel.R"))

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
load(paste0(GP_Dir,"eval.state.newer"))
best_gp = quant.to.prior(evaluations$best[[length(evaluations$best)]]$known$parameters) %>% unlist()

df3 <- param_values_gathered %>%
  group_by(key) %>%
  dplyr::summarise(Mean.SL = mean(value)) %>% as.data.frame()
df3$GP = NA
for(i in 1:length(param_names)){
  df3[which(df3$key==paste0(param_names[i])),"GP"] = as.numeric(best_gp[which(names(best_gp)==paste0(param_names[i]))])
}

#2. GP Validation
load(paste0(GPValidation_Dir,"eval.state.newer"))
best_gpvalid = quant.to.prior(evaluations$best[[length(evaluations$best)]]$known$parameters) %>% unlist()


for(i in 1:length(param_names)){
  df3[which(df3$key==paste0(param_names[i])),"GP_Validation"] = as.numeric(best_gpvalid[which(names(best_gpvalid)==paste0(param_names[i]))])
}


#3. GPSG
load(paste0(GPSG_Dir,"eval.state.newer"))
best_gpsg = quant.to.prior(evaluations$best[[length(evaluations$best)]]$known$parameters) %>% unlist()


df3$GPSG = NA
for(i in 1:length(param_names)){
  df3[which(df3$key==paste0(param_names[i])),"GPSG"] = as.numeric(best_gpsg[which(names(best_gpsg)==paste0(param_names[i]))])
}

#4. GPSG Validation
load(paste0(GPSGValidation_Dir,"eval.state.newer"))
best_gpsgvalid = quant.to.prior(evaluations$best[[length(evaluations$best)]]$known$parameters) %>% unlist()


for(i in 1:length(param_names)){
  df3[which(df3$key==paste0(param_names[i])),"GPSG_Validation"] = as.numeric(best_gpsgvalid[which(names(best_gpsgvalid)==paste0(param_names[i]))])
}
#Plot 

ggplot(param_values_gathered, aes(value)) + 
  #geom_histogram(aes(y = ..density.., fill = ..count..), bins = 50) +
  geom_density(color="darkblue", fill="lightblue")+ 
  facet_wrap(~key, scales = "free", labeller = as_labeller(param_labeller)) +
  ggtitle("Prior parameter distributions") + 
  theme_bw()

logplot = ggplot(param_values_gathered, aes(value)) + 
  #geom_histogram(aes(y = ..density.., fill = ..count..), bins = 50) +
  geom_density(alpha=0.4,color=NA, fill="#D3D4D9")+ 
  scale_x_continuous(trans = 'log10') + 
  facet_wrap(~factor(key,levels=unique(as.factor(param_values_gathered$key))), scales = "free", labeller = as_labeller(param_labeller)) +
  ggtitle("Data recovery validation of posterior estimates") 

logplotGP=logplot + 
  geom_vline(data=df3, aes(xintercept= GP, colour="GP"),size=1,show.legend=T) + 
  #geom_vline(data=df3, aes(xintercept= GPSG, colour="GPSG"),size=1,show.legend=T)  + 
  geom_vline(data=df3, aes(xintercept= GP_Validation, colour="GP_Validation"), linetype="dashed",show.legend=T) + 
  scale_color_manual(name = "Algorithm", values = c(GP = "#FF7745", GP_Validation="#A25EA2")) + 
  theme_classic()

logplotGPSG=logplot + 
  #geom_vline(data=df3, aes(xintercept= GP, colour="GP"),size=1,show.legend=T) + 
  geom_vline(data=df3, aes(xintercept= GPSG, colour="GPSG"),size=1,show.legend=T)  + 
  geom_vline(data=df3, aes(xintercept= GPSG_Validation, colour="GPSG_Validation"), linetype="dashed",show.legend=T) + 
  scale_color_manual(name = "Algorithm", values = c(GPSG = "#720026", GPSG_Validation="#1D005F")) + 
  theme_classic()

pdf(paste0("OutputPlots/GP_Validation_Log_Prior_distributions.pdf"), width = 11.89, height = 7.61)
print(logplotGP)
print(logplotGPSG)
dev.off()
