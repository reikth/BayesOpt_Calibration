############################################################
# EIR TO PREVALENCE
#
# Plot results requested by the Malaria Modelling Consortium.
#
# Written by A.J.Shattock - andrewjames.shattock@unibas.ch
############################################################

pacman::p_load(ggplot2,latex2exp,wrapr,tidyverse)


# ---------------------------------------------------------
# Parent function for calling plotting function
# ---------------------------------------------------------
run_results = function(o) {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Running results")
  
  # Call plotting functions according to flags
  if (o$plot$eir_error) plot_eir_error(o)
  if (o$plot$eir2prev)  plot_eir2prev(o)
  plot_prev2inc(o)
  plot_prev2nSevere(o)
  plot_prev2death(o)
  
  plot_age_incidence(o)
  plot_age_severe(o)
  plot_age_death(o)
}

# ---------------------------------------------------------
# Plot error between inputEIR and simulatedEIR
# ---------------------------------------------------------
plot_eir_error = function(o) {
  
  message(" - Plotting inputEIR-simulatedEIR errors")
  
  # ---- Figure properties ----
  
  # Text sizes
  ttl_size = 22  # Title
  lab_size = 18  # Axes labels
  stp_size = 16  # Strip sizes
  tck_size = 12  # Ticks
  lgd_size = 12  # Legend
  
  # Figure size
  fig_size = c(12, 7)  # Width & height
  
  # ---- Load data ----
  
  # Load EIR-prevalence results table
  eir_prev_pth = paste0(o$pth$results, "eir_prevalence.txt")
  eir_prev_df  = read.table(eir_prev_pth, header = TRUE, 
                            stringsAsFactors = FALSE)
  
  # Calculate error between input and simulated
  eir_df = eir_prev_df %>% 
    mutate(eir = factor(eir, levels = rev(o$sim$eir)), 
           error_abs = simulatedEIR - inputEIR, 
           error_rel = 100 * error_abs / inputEIR, 
           season    = paste("Season", season), 
           model     = paste("Model", model)) %>%
    dplyr::select(model, season, eir, error_abs, error_rel)
  
  # ---- Produce error plots ----
  
  # Two types of plot: absolute and relative
  error_types = qc(error_abs = "Absolute", 
                   error_rel = "Relative")
  
  # Also state how many burn-in years we've had
  burn_in_years = o$time$start_year - o$time$burn_in
  
  # Set legend title
  legend_title = "Input EIR"
  
  # Loop through the two types of error we're reporting
  for (error_type in names(error_types)) {
    
    # Report the error as normalised density plot
    g = ggplot(eir_df, aes_string(x = error_type, colour = "eir")) + 
      geom_density(aes(y = ..scaled.., fill = eir), alpha = 0.5) + 
      geom_vline(aes(xintercept = 0), linetype = "dashed") + 
      facet_grid(season ~ model)
    
    # Construct plot title (explaining error type and burn in period)
    title = paste0("Input vs simulated EIR (", error_types[error_type], 
                   " error, ", burn_in_years, " year burn-in)")
    
    # Set the title and axes labels
    g = g + ggtitle(title) + ylab("Normalised frequency") + 
      xlab(paste(error_types[error_type], "error"))

    # Reverse legend for both geoms (seems stupid - surely a better way)
    g = g + guides(fill   = guide_legend(title = legend_title, reverse = TRUE), 
                   colour = guide_legend(title = legend_title, reverse = TRUE))
    
    # Fix up text sizes  
    g = g + theme(plot.title   = element_text(size = ttl_size, hjust = 0.5),
                  strip.text   = element_text(size = stp_size), 
                  axis.title   = element_text(size = lab_size), 
                  axis.text.x  = element_text(size = tck_size), 
                  axis.text.y  = element_blank(),
                  legend.title = element_text(size = lgd_size), 
                  legend.text  = element_text(size = lgd_size))
    
    # Save plot to file
    ggsave(paste0(o$pth$figures, "EIR_error_", error_types[error_type], ".png"), 
           plot = g, width = fig_size[1], height = fig_size[2])
  }
}

# ---------------------------------------------------------
# Plot OpenMalaria EIR-prevalence relationship
# ---------------------------------------------------------
plot_eir2prev = function(o) {
  
  # Produce EIR-prevalence relationships for Malaria Modelling
  # Consortium. This code is designed to reproduce plots presented
  # in Melissa's four model comparison paper:
  #
  # https://www.sciencedirect.com/science/article/pii/S0140673615007254
  
  message(" - Plotting EIR-prevalence relationship")
  
  # ---- Figure properties ----
  
  # Colours (one colour per user-defined season)
  
  #colours = c("#ff7f51", "#720026") #orange and burgundy. Alternatively: 
  colours = c("black","#FF7745", "#720026")  # black for GA, orange for GP, purple for GPSG
  

  # Text sizes
  lab_size = 24  # Labels
  tck_size = 20  # Ticks
  lgd_size = 20  # Legend
  
  # Figure size
  fig_size = c(12, 7)  # Width & height
  
  # ---- Load data ----
  
  # TODO: Plot simulated EIR rather than input EIR?
  
  # Load EIR-prevalence results table
  eir_prev_pth = paste0(o$pth$results, "eir_prevalence.txt")
  eir_prev_df  = na.omit(read.table(eir_prev_pth, header = TRUE, 
                            stringsAsFactors = FALSE))
  
  # Format prevalence as percentage for pretty plotting
  eir_prev_df$prevalence = eir_prev_df$prevalence_2to10 * 100
  
  # ---- Generate plotting dataframe with lower and uppper bound ----
  
  # Remove variable we will aggregate in plotting dataframe
  plot_df = unique(eir_prev_df %>% dplyr::select(eir, access, season,model))
  
  # We want to keep grounpings of EIR, access, and seasonality
  group_by = list(eir_prev_df$eir, eir_prev_df$access, eir_prev_df$season,eir_prev_df$model)
  
  # Group by and aggreagate with mean to get best estimate
  tmp = aggregate(eir_prev_df$prevalence, by = group_by, mean) %>% 
    dplyr::rename(eir = Group.1, access = Group.2, season = Group.3, model = Group.4, prevalence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "season", "model"), all.y=T)

  # Group by and aggreagate with min and max to get lower and upper bounds
  tmp = aggregate(eir_prev_df$prevalence, by = group_by, min)%>% 
    dplyr::rename(eir = Group.1, access = Group.2, season = Group.3, model = Group.4, lower = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "season", "model"))
  tmp = aggregate(eir_prev_df$prevalence, by = group_by, max)%>% 
    dplyr::rename(eir = Group.1, access = Group.2, season = Group.3, model = Group.4, upper = x)
  plot_df =  merge(plot_df,tmp, by = c("eir", "access", "season", "model"))
  rm(tmp)
  
  # ---- Format treatment access and seasonality variables ----
  
  # Convert access to percentage strings
  plot_df$access = paste0(plot_df$access * 100, "%")
  
  # Convert seasonal indices to descriptive strings
  plot_df$season[plot_df$season == 1] = "Perennial"
  plot_df$season[plot_df$season == 2] = "Seasonal"
  
  plot_df$model[plot_df$model == "0000GA"] = "GA"
  plot_df$model[plot_df$model == "0000GP"] = "GP"
  plot_df$model[plot_df$model == "0000GPSG"] = "GPSG"
  
  plot_df$group = paste(plot_df[,"season"],plot_df[,"model"])
  
  # ---- Produce plot ----
  # Y-axis settings
  y_max   = 90 # Axis upper limit
  y_ticks = seq(0, y_max, by = 20)  # Tick marks
  
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df, aes(x = eir, y = prevalence, colour = model,group =group)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_line(size = 2,aes(linetype=season))# + 
    #facet_wrap(~model)
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours) + 
    scale_color_manual(values = colours) + 
    scale_linetype_manual(values = c("dotted","solid"))
  
  # Sort out y axes limits and ticks
  g = g + scale_y_continuous(limits = c(0, y_max), 
                             breaks = y_ticks, 
                             labels = paste0(y_ticks, "%"))
  
  # Transform x axis to log10 scale
  g = g + scale_x_continuous(limits = c( min(plot_df$eir), 
                                         max(plot_df$eir)),
                             breaks = c(1,10,100,500),
                             trans = 'log10')
  
  # Set legend title and axes labels
  g = g + labs(colour = "Parameterization", 
               linetype = "Seasonality",
               x = "Entomological Inoculation Rate (EIR)", 
               y = TeX("\\textit{Pf}PR_{2-10}"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(strip.text   = element_text(size = lab_size), 
                axis.title            = element_text(size = lab_size), 
                axis.text             = element_text(size = tck_size), 
                legend.title          = element_text(size = lgd_size,face="bold"), 
                legend.text           = element_text(size = lgd_size),
                legend.justification  = c(0, 1), 
                legend.position       = c(0.01, 1),
                legend.box.just       = "left", 
                legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
                guides(
                   linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
                   colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "EIR_prevalence_relationship.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
}



# ---------------------------------------------------------
# Plot OpenMalaria prevalence-incidence relationship
# ---------------------------------------------------------
plot_prev2inc = function(o) {
  
  # Produce prevalence-incidence relationships for Malaria Modelling
  # Consortium. This code is designed to reproduce plots presented
  # in Melissa's four model comparison paper:
  #
  # https://www.sciencedirect.com/science/article/pii/S0140673615007254
  
  message(" - Plotting prevalence- incidence relationship")
  
  # ---- Figure properties ----
  
  # Colours (one colour per user-defined season)
  #colours = c("black","#FF7745", "#720026")
  colours = c("black","#2F6089", "#508FC3","#AFCCE4")
  #colours = c("#FF7745", "#720026")  # Orange & burgundy
  
  # Y-axis settings
  y_max   = 6  # Axis upper limit
  x_max   = 75
  y_ticks = seq(0, y_max, by = 1)  # Tick marks
  x_ticks = seq(0, x_max, by = 25) 
  
  # Text sizes
  lab_size = 18  # Labels
  tck_size = 12  # Ticks
  lgd_size = 12  # Legend
  
  # Figure size
  fig_size = c(11.23, 7.13)  # Width & height
  
  # ---- Load data ----
  
  # TODO: Plot simulated EIR rather than input EIR?
  
  #load prevalence incidence results
  prev_inc_pth = paste0(o$pth$results, "prev_incidence.txt")
  prev_inc_df  = read.table(prev_inc_pth, header = TRUE, 
                            stringsAsFactors = FALSE)
  
  #merge 1st FOUR age groups to plot under 5 years of age
  #sum for all incidence measures
  tmp = prev_inc_df[(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),]
  tmp$prevalence_2to10 = NULL
  agg_by = list(tmp$sim_id,tmp$eir,tmp$season,tmp$seed,tmp$model)
  tmp = tmp %>%
    group_by(sim_id, eir, season, seed, model) %>% 
    summarise_each(sum)
  #reset age-prevalence
  tmp = tmp %>% mutate(
    age_prev=nPatent/nHost,
    age =     2.5
  )
  
  a = unique(prev_inc_df[(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),] %>% dplyr::select(sim_id,prevalence_2to10))
  tmp = merge(tmp,a,by="sim_id")

  #bind to original dataframe,  prev_inc_df 
  prev_inc_df = rbind(prev_inc_df[!(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),],tmp)
  
  rm(tmp,a)
  # Format prevalence as percentage for pretty plotting
  prev_inc_df$prevalence = prev_inc_df$prevalence_2to10 * 100
  #prev_inc_df$age_prevalence = prev_inc_df$age_prev * 100
  
  # ---- Generate plotting dataframe with lower and uppper bound ----
 
  
  # Remove variable we will aggregate in plotting dataframe
  plot_df = unique(prev_inc_df %>% dplyr::select(eir,age,access, season,model))

    # We want to keep grounpings of EIR, access, and seasonality, and model
  group_by = list(prev_inc_df$eir, prev_inc_df$access,prev_inc_df$age, prev_inc_df$season,prev_inc_df$model)
  prev_inc_df$inc_ppyo = prev_inc_df$nUncomp/prev_inc_df$nHost
  
  
  # Group by and aggreagate with mean, min and max to get lower and upper bounds
  tmp = aggregate(prev_inc_df$inc_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$inc_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$inc_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  #add prev

  tmp =   aggregate(prev_inc_df$prevalence, by = group_by, FUN = mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, prev2to10 = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  plot_df = plot_df[plot_df$age<20,]
  
  
  rm(tmp)
  plot_df$eir=NULL
  
  # ---- Format treatment access and seasonality variables ----
  
  # Convert access to percentage strings
  plot_df$access = paste0(plot_df$access * 100, "%")

  # Convert seasonal indices to descriptive strings
  plot_df$season[plot_df$season == 1] = "Perennial"
  plot_df$season[plot_df$season == 2] = "Seasonal"
  
  plot_df$model[plot_df$model == "0000GA"] = "GA"
  plot_df$model[plot_df$model == "0000GP"] = "GP-BO"
  plot_df$model[plot_df$model == "0000GPSG"] = "GPSG-BO"
  
  plot_df$group = paste(plot_df[,"season"],plot_df[,"model"])
  
  # ---- Produce plot ----
  
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df, aes(x = prev2to10,y = incidence, colour = as.factor(age),group =as.factor(age))) + 
    facet_grid(season~model) + 
    geom_line() +
    
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(age)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 


  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("0-5","5-10","10-15","15-20")) + 
    scale_color_manual(values = colours,labels=c("0-5","5-10","10-15","15-20")) 
  
  # Sort out y axes limits and ticks
  g = g + scale_y_continuous(limits = c(0, y_max), 
                             breaks = y_ticks, 
                             labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks, "%"))
  
  # Set legend title and axes labels
  g = g + labs(colour = "age group",
               fill = "age group", 
               y = TeX("Clinincal incidence (events per person per year)"), 
               x = TeX("\\textit{Pf}PR_{2-10}"))
  
  # Fix up text sizes  
    g = g + 
      theme_classic() + 
      theme(panel.spacing = unit(0, "mm"),
            panel.border = element_rect(fill=NA,color="grey60"),
            strip.text   = element_text(size = lab_size), 
            strip.background = element_rect(fill="grey90",color="grey70"),
            axis.title            = element_text(size = lab_size), 
            axis.text             = element_text(size = tck_size), 
            legend.title          = element_text(size = lgd_size,face="bold"), 
            legend.text           = element_text(size = lgd_size),
            legend.justification  = c(0, 1), 
            legend.position       = c(0.001, 0.999),
            legend.box.just       = "left", 
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
      guides(
        linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
        colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "Prevalence_incidence_relationship.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
}


# ---------------------------------------------------------
# Plot OpenMalaria prevalence-nSevere relationship
# ---------------------------------------------------------  
#Do this!!!
plot_prev2nSevere = function(o) {

  message(" - Plotting prevalence- severe incidence relationship")
  
  # ---- Figure properties ----
  
  # Colours (one colour per user-defined season)
  #colours = c("black","#FF7745", "#720026")
  colours = c("black","#2F6089", "#508FC3","#AFCCE4")
  #colours = c("#FF7745", "#720026")  # Orange & burgundy
  
  # Y-axis settings
  y_max   = 120  # Axis upper limit
  x_max   = 75
  y_ticks = seq(0, y_max, by = 20)  # Tick marks
  x_ticks = seq(0, x_max, by = 25) 
  
  # Text sizes
  lab_size = 18  # Labels
  tck_size = 12  # Ticks
  lgd_size = 12  # Legend
  
  # Figure size
  fig_size = c(11.23, 7.13)  # Width & height
  
  # ---- Load data ----
  
  # TODO: Plot simulated EIR rather than input EIR?
  
  #load prevalence incidence results
  prev_inc_pth = paste0(o$pth$results, "prev_incidence.txt")
  prev_inc_df  = read.table(prev_inc_pth, header = TRUE, 
                            stringsAsFactors = FALSE)

   
  #merge 1st FOUR age groups to plot under 5 years of age
  #sum for all incidence measures
  tmp = prev_inc_df[(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),]
  
  
  tmp$prevalence_2to10 = NULL
  agg_by = list(tmp$sim_id,tmp$eir,tmp$season,tmp$seed,tmp$model)
  tmp = tmp %>%
    group_by(sim_id, eir, season, seed, model) %>% 
    summarise_each(sum)
  #reset age-prevalence
  tmp = tmp %>% mutate(
    age_prev=nPatent/nHost,
    age =     2.5
  )
  
  a = unique(prev_inc_df[(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),] %>% dplyr::select(sim_id,prevalence_2to10))
  tmp = merge(tmp,a,by="sim_id")
  

  #bind to original dataframe,  prev_inc_df 
  prev_inc_df = rbind(prev_inc_df[!(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),],tmp)
  
  rm(tmp,a)
  # Format prevalence as percentage for pretty plotting
  prev_inc_df$prevalence = prev_inc_df$prevalence_2to10 * 100
  #prev_inc_df$age_prevalence = prev_inc_df$age_prev * 100
  
  # ---- Generate plotting dataframe with lower and uppper bound ----
  
  
  # Remove variable we will aggregate in plotting dataframe
  plot_df = unique(prev_inc_df %>% dplyr::select(eir,age,access, season,model))
  
  # We want to keep grounpings of EIR, access, and seasonality, and model
  group_by = list(prev_inc_df$eir, prev_inc_df$access,prev_inc_df$age, prev_inc_df$season,prev_inc_df$model)
  prev_inc_df$severe_ppyo = prev_inc_df$nSevere/prev_inc_df$nHost
  prev_inc_df$exp_severe_ppyo = prev_inc_df$expectedSevere/prev_inc_df$nHost
  
  
  # Group by and aggreagate with mean, min and max to get lower and upper bounds
  tmp = aggregate(prev_inc_df$severe_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, severe = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$severe_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$severe_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  #add prev
  
  tmp =   aggregate(prev_inc_df$prevalence, by = group_by, FUN = mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, prev2to10 = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  plot_df = plot_df[plot_df$age<20,]
  
  
  rm(tmp)
  plot_df$eir=NULL
  
  # ---- Format treatment access and seasonality variables ----
  
  # Convert access to percentage strings
  plot_df =plot_df %>% mutate(
    severe_pyyo = severe * 1000,
    upper = upper * 1000, 
    lower = lower * 1000
  )
  plot_df$access = paste0(plot_df$access * 100, "%")
  
  # Convert seasonal indices to descriptive strings
  plot_df$season[plot_df$season == 1] = "Perennial"
  plot_df$season[plot_df$season == 2] = "Seasonal"
  
  plot_df$model[plot_df$model == "0000GA"] = "GA"
  plot_df$model[plot_df$model == "0000GP"] = "GP-BO"
  plot_df$model[plot_df$model == "0000GPSG"] = "GPSG-BO"
  
  plot_df$group = paste(plot_df[,"season"],plot_df[,"model"])
  
  # ---- Produce plot ----
  
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df, aes(x = prev2to10,y = severe_pyyo, colour = as.factor(age),group =as.factor(age))) + 
    facet_grid(season~model) + 
    geom_line() +
    
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = as.factor(age)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("0-5","5-10","10-15","15-20")) + 
    scale_color_manual(values = colours,labels=c("0-5","5-10","10-15","15-20")) 
  
  # Sort out y axes limits and ticks
  g = g + scale_y_continuous(limits = c(0, y_max), 
                             breaks = y_ticks, 
                             labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks, "%"))
  
  # Set legend title and axes labels
  g = g + labs(colour = "age group",
               fill = "age group", 
               y = TeX("Severe incidence (events per person 1000 person-years)"), 
               x = TeX("\\textit{Pf}PR_{2-10}"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.justification  = c(0, 1), 
          legend.position       = c(0.001, 0.999),
          legend.box.just       = "left", 
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "Prevalence_severe_relationship.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
}

  
  
# ---------------------------------------------------------
# Plot OpenMalaria prevalence-deaths relationship
# ---------------------------------------------------------  

plot_prev2death = function(o) {
  # ---- Figure properties ----
  
  # Colours (one colour per user-defined season)
  #colours = c("black","#FF7745", "#720026")
  colours = c("black","#2F6089", "#508FC3","#AFCCE4")
  #colours = c("#FF7745", "#720026")  # Orange & burgundy
  
  # Y-axis settings
  y_max   = 35  # Axis upper limit
  x_max   = 80
  y_ticks = seq(0, y_max, by = 10)  # Tick marks
  x_ticks = seq(0, x_max, by = 20) 
  
  # Text sizes
  lab_size = 18  # Labels
  tck_size = 12  # Ticks
  lgd_size = 12  # Legend
  
  # Figure size
  fig_size = c(11.23, 7.13)  # Width & height
  
  # ---- Load data ----
  
  # TODO: Plot simulated EIR rather than input EIR?
  
  #load prevalence incidence results
  prev_inc_pth = paste0(o$pth$results, "prev_incidence.txt")
  prev_inc_df  = read.table(prev_inc_pth, header = TRUE, 
                            stringsAsFactors = FALSE)
  
  #merge 1st FOUR age groups to plot under 5 years of age
  #sum for all incidence measures
  tmp = prev_inc_df[(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),]
  tmp$prevalence_2to10 = NULL
  agg_by = list(tmp$sim_id,tmp$eir,tmp$season,tmp$seed,tmp$model)
  tmp = tmp %>%
    group_by(sim_id, eir, season, seed, model) %>% 
    summarise_each(sum)
  #reset age-prevalence
  tmp = tmp %>% mutate(
    age_prev=nPatent/nHost,
    age =     2.5
  )
  
  a = unique(prev_inc_df[(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5),] %>% dplyr::select(sim_id,prevalence_2to10))
  tmp = merge(tmp,a,by="sim_id")
  
  #bind to original dataframe,  prev_inc_df 
  prev_inc_df = rbind(prev_inc_df[!(prev_inc_df$age==0.25 |prev_inc_df$age==0.75 | prev_inc_df$age==1.5 | prev_inc_df$age==3.5) ,],tmp)
  
  rm(tmp,a)
  # Format prevalence as percentage for pretty plotting
  prev_inc_df$prevalence = prev_inc_df$prevalence_2to10 * 100
  #prev_inc_df$age_prevalence = prev_inc_df$age_prev * 100
  
  # ---- Generate plotting dataframe with lower and uppper bound ----
  
  
  # Remove variable we will aggregate in plotting dataframe
  plot_df = unique(prev_inc_df %>% dplyr::select(eir,age,access, season,model))
  
  # We want to keep grounpings of EIR, access, and seasonality, and model
  group_by = list(prev_inc_df$eir, prev_inc_df$access,prev_inc_df$age, prev_inc_df$season,prev_inc_df$model)
  
  prev_inc_df$sev_inc_ppyo = prev_inc_df$nSevere/prev_inc_df$nHost
  prev_inc_df$dir_death_ppyo = prev_inc_df$expectedDirectDeaths/prev_inc_df$nHost
  prev_inc_df$indir_death_ppyo = prev_inc_df$expectedIndirectDeaths/prev_inc_df$nHost
  prev_inc_df$all_death_ppyo = (prev_inc_df$expectedDirectDeaths+prev_inc_df$expectedIndirectDeaths)/prev_inc_df$nHost
  
  # Group by and aggreagate with mean, min and max to get lower and upper bounds
  
  # Group by and aggreagate with mean, min and max to get lower and upper bounds
  #indirect deaths
  tmp = aggregate(prev_inc_df$indir_death_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, indirDeath_incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$indir_death_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, indirDeath_lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$indir_death_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, indirDeath_upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))

  
  #direct deaths:
  tmp = aggregate(prev_inc_df$dir_death_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, dirDeath_incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$dir_death_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, dirDeath_lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$dir_death_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, dirDeath_upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  #all deaths:
  tmp = aggregate(prev_inc_df$all_death_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, allDeath_incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$all_death_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, allDeath_lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$all_death_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, allDeath_upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
    #add prev
  
  tmp =   aggregate(prev_inc_df$prevalence, by = group_by, FUN = mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, prev2to10 = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  plot_df = plot_df[plot_df$age<20,]
  
  
  rm(tmp)
  plot_df$eir=NULL
  
  # ---- Format treatment access and seasonality variables ----
  
  # Convert access to percentage strings
  plot_df$access = paste0(plot_df$access * 100, "%")
  #Severe incidence per 1000 per year
  plot_df = plot_df %>% 
    mutate(
      indirDeath_incidence  = indirDeath_incidence*1000,
      indirDeath_lower  =indirDeath_lower * 1000,
      indirDeath_upper = indirDeath_upper * 1000,
      dirDeath_incidence  = dirDeath_incidence*1000,
      dirDeath_lower  =dirDeath_lower * 1000,
      dirDeath_upper = dirDeath_upper * 1000,
      allDeath_incidence  = allDeath_incidence*1000,
      allDeath_lower  =allDeath_lower * 1000,
      allDeath_upper = allDeath_upper * 1000
    )
  
  # Convert seasonal indices to descriptive strings
  plot_df$season[plot_df$season == 1] = "Perennial"
  plot_df$season[plot_df$season == 2] = "Seasonal"
  
  plot_df$model[plot_df$model == "0000GA"] = "GA"
  plot_df$model[plot_df$model == "0000GP"] = "GP-BO"
  plot_df$model[plot_df$model == "0000GPSG"] = "GPSG-BO"
  
  plot_df$group = paste(plot_df[,"season"],plot_df[,"model"])
  
  # ---- Produce plot ----
  
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df, aes(x = prev2to10, y = allDeath_incidence,colour = as.factor(age),group =as.factor(age))) + 
    facet_grid(season~model) + 
    geom_line(aes(x = prev2to10, y = allDeath_incidence)) +
    geom_line(aes(x = prev2to10, y = dirDeath_incidence),linetype="dashed") +
    
    geom_ribbon(aes(ymin = allDeath_lower, ymax = allDeath_upper, fill = as.factor(age)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_ribbon(aes(ymin = dirDeath_lower, ymax = dirDeath_upper, fill = as.factor(age)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("0-5","5-10","10-15","15-20")) + 
    scale_color_manual(values = colours,labels=c("0-5","5-10","10-15","15-20")) + 
    scale_linetype_manual(values = c("solid", "dashed"), labels=c("indirect", "direct"))
  
  # Sort out y axes limits and ticks
  g = g + scale_y_continuous(limits = c(0, y_max), 
                             breaks = y_ticks, 
                             labels = paste0(y_ticks))

  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks, "%"))
  
  # Set legend title and axes labels
  g = g + labs(colour = "age group",
               fill = "age group", 
               linetype="mortality",
               y = TeX("Deaths (events per 1000 person-years)"), 
               x = TeX("\\textit{Pf}PR_{2-10}"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.justification  = c(0, 1), 
          legend.position       = c(0.001, 0.999),
          legend.box.just       = "left", 
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "Prevalence_death_relationship.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
  
  
}


# ---------------------------------------------------------
# Plot OpenMalaria clinical incidence by age
# ---------------------------------------------------------  

plot_age_incidence = function(o) {
  
  # Produce incidence by age relationships for Malaria Modelling
  # Consortium. This code is designed to reproduce plots presented
  # in Melissa's four model comparison paper:
  #
  # https://www.sciencedirect.com/science/article/pii/S0140673615007254
  
  message(" - Plotting age - incidence relationship")
  
  # ---- Figure properties ----
  
  # Colours (one colour per user-defined season)
  colours = c("black","#FF7745", "#720026") #orange and burgundy for fits
  #colours = c("black","#2F6089", "#508FC3","#AFCCE4") #blues for ages
  #colours = c("#FF7745", "#720026")  # Orange & burgundy
  
  # Y-axis settings
  y_max   = 4  # Axis upper limit
  x_max   = 20
  y_ticks = seq(0, y_max, by = 1)  # Tick marks
  x_ticks = seq(0, x_max, by = 5) 
  
  # Text sizes
  lab_size = 18  # Labels
  tck_size = 12  # Ticks
  lgd_size = 12  # Legend
  
  # Figure size
  fig_size = c(11.23, 8.5)  # Width & height
  
  # ---- Load data ----
  
  # TODO: Plot simulated EIR rather than input EIR?
  
  #load prevalence incidence results
  prev_inc_pth = paste0(o$pth$results, "prev_incidence.txt")
  prev_inc_df  = read.table(prev_inc_pth, header = TRUE, 
                            stringsAsFactors = FALSE)
  
  
  # Format prevalence as percentage for pretty plotting
  prev_inc_df$prevalence = prev_inc_df$prevalence_2to10 * 100
  #prev_inc_df$age_prevalence = prev_inc_df$age_prev * 100
  
  # ---- Generate plotting dataframe with lower and uppper bound ----
  
  
  # Remove variable we will aggregate in plotting dataframe
  plot_df = unique(prev_inc_df %>% dplyr::select(eir,age,access, season,model))
  
  # We want to keep grounpings of EIR, access, and seasonality, and model
  group_by = list(prev_inc_df$eir, prev_inc_df$access,prev_inc_df$age, prev_inc_df$season,prev_inc_df$model)
  prev_inc_df$inc_ppyo = prev_inc_df$nUncomp/prev_inc_df$nHost
  
  
  # Group by and aggreagate with mean, min and max to get lower and upper bounds
  tmp = aggregate(prev_inc_df$inc_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$inc_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$inc_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  #add prev
  
  tmp =   aggregate(prev_inc_df$prevalence, by = group_by, FUN = mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, prev2to10 = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  plot_df = plot_df[plot_df$age<20,]
  
  
  rm(tmp)
  plot_df$eir=NULL
  
  # ---- Format treatment access and seasonality variables ----
  
  # Convert access to percentage strings
  plot_df$access = paste0(plot_df$access * 100, "%")
  
  # Convert seasonal indices to descriptive strings
  plot_df$season[plot_df$season == 1] = "Perennial"
  plot_df$season[plot_df$season == 2] = "Seasonal"
  
  plot_df$model[plot_df$model == "0000GA"] = "GA"
  plot_df$model[plot_df$model == "0000GP"] = "GP-BO"
  plot_df$model[plot_df$model == "0000GPSG"] = "GPSG-BO"
  
  plot_df$group = paste(plot_df[,"season"],plot_df[,"model"])
  
  # ---- Produce plot ----
  plot_df_perennial = plot_df %>% filter(season=="Perennial")
  plot_df_seasonal = plot_df %>% filter(season=="Seasonal")
  mround <- function(x,base){
    base*round(x/base)
  }
  
  plot_df_seasonal$prev2to10round = mround(plot_df_seasonal$prev2to10,10)
  
  #select PfPR= 3,10,30,50
  plot_df_seasonal = plot_df_seasonal %>%
    mutate(
      prev_cat = case_when(
        prev2to10>=2.5 & prev2to10<=4.4   ~ "3",
        prev2to10>=9.5 & prev2to10<=10.5   ~ "10",
        prev2to10>=28 & prev2to10<=32   ~ "30",
        prev2to10>=47 & prev2to10<=53   ~ "50"
      )
    )
  plot_df_seasonal = plot_df_seasonal %>% 
    na.omit() %>%
    dplyr::select(prev_cat,model,age,incidence,lower,upper)%>%
    group_by(prev_cat,age,model)%>%
    summarise_each(funs(mean,min,max)) %>%
    dplyr::mutate(prev_cat, as.factor(prev_cat)) %>%
    dplyr::mutate(prev_cat = fct_relevel(prev_cat,sort))
  plot_df_seasonal$prev_cat2 = factor(plot_df_seasonal$prev_cat, 
                                     levels=c('3','10','30','50'),
                                     labels=c(TeX('PfPR_{2-10}=3\\%'),
                                              TeX('PfPR_{2-10}=10\\%'),
                                              TeX('PfPR_{2-10}=30\\%'),
                                              TeX('PfPR_{2-10}=50\\%')))
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df_seasonal, aes(x = age,y = incidence_mean, colour = as.factor(model),group =as.factor(model))) + 
    facet_grid(prev_cat2~model,labeller = label_parsed,scales="free_y") + 
    geom_line() +
    
    geom_ribbon(aes(ymin = lower_min, ymax = upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("GA","GP-BO","GPSG-BO")) + 
    scale_color_manual(values = colours,labels=c("GA","GP-BO","GPSG-BO")) 
  
  # Sort out y axes limits and ticks
  #g = g + scale_y_continuous(limits = c(0, y_max), 
  #                           breaks = y_ticks, 
  #                           labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks))
  
  # Set legend title and axes labels
  g = g + labs(title = "Clinical incidence by age (seasonal transmission)",
               colour = "Parameterization",
               fill = "Parameterization", 
               y = TeX("Clinincal incidence (events per person per year)"), 
               x = TeX("age"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          title                 = element_text(size = lab_size+2,face="bold"), 
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.justification  = c(0, 1), 
          legend.position       = c("bottom"),
          legend.box.just       = "right", 
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "age_incidence_relationship_seasonal.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])

  y_max=5
  #perennial
  plot_df_perennial$prev2to10round = mround(plot_df_perennial$prev2to10,10)
  
  #select PfPR= 3,10,30,50
  plot_df_perennial = plot_df_perennial %>%
    mutate(
      prev_cat = case_when(
        prev2to10>=2.5 & prev2to10<=4.4   ~ "3",
        prev2to10>=9.5 & prev2to10<=10.5   ~ "10",
        prev2to10>=28 & prev2to10<=32   ~ "30",
        prev2to10>=47 & prev2to10<=53   ~ "50"
      )
    )
  plot_df_perennial = plot_df_perennial %>% 
    na.omit() %>%
    dplyr::select(prev_cat,model,age,incidence,lower,upper)%>%
    group_by(prev_cat,age,model)%>%
    summarise_each(funs(mean,min,max)) %>%
    mutate(prev_cat, as.factor(prev_cat)) %>%
    mutate(prev_cat = fct_relevel(prev_cat,sort))
  plot_df_perennial$prev_cat2 = factor(plot_df_perennial$prev_cat, 
                                      levels=c('3','10','30','50'),
                                      labels=c(TeX('PfPR_{2-10}=3\\%'),
                                               TeX('PfPR_{2-10}=10\\%'),
                                               TeX('PfPR_{2-10}=30\\%'),
                                               TeX('PfPR_{2-10}=50\\%')))
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df_perennial, aes(x = age,y = incidence_mean, colour = as.factor(model),group =as.factor(model))) + 
    facet_grid(prev_cat2~model,labeller = label_parsed,scales="free_y") + 
    geom_line() +
    
    geom_ribbon(aes(ymin = lower_min, ymax = upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("GA","GP-BO","GPSG-BO")) + 
    scale_color_manual(values = colours,labels=c("GA","GP-BO","GPSG-BO")) 
  
  # Sort out y axes limits and ticks
  #g = g + scale_y_continuous(limits = c(0, y_max), 
  #                           breaks = y_ticks, 
  #                           labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks))
  
  # Set legend title and axes labels
  g = g + labs(title = "Clinical incidence by age (perennial transmission)",
               colour = "Parameterization",
               fill = "Parameterization", 
               y = TeX("Clinincal incidence (events per person per year)"), 
               x = TeX("age"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          title                 = element_text(size = lab_size+2,face="bold"), 
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.justification  = c(0, 1), 
          legend.position       = c("bottom"),
          legend.box.just       = "right", 
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "age_incidence_relationship_perennial.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
}


# ---------------------------------------------------------
# Plot OpenMalaria severe incidence by age
# ---------------------------------------------------------  

plot_age_severe = function(o) {
  
  # Produce severe incidence by age relationships for Malaria Modelling
  # Consortium. This code is designed to reproduce plots presented
  # in Melissa's four model comparison paper:
  #
  # https://www.sciencedirect.com/science/article/pii/S0140673615007254
  
  message(" - Plotting age - severe incidence relationship")
  
  # ---- Figure properties ----
  
  # Colours (one colour per user-defined season)
  colours = c("black","#FF7745", "#720026") #orange and burgundy for fits
  #colours = c("black","#2F6089", "#508FC3","#AFCCE4") #blues for ages
  #colours = c("#FF7745", "#720026")  # Orange & burgundy
  
  # Y-axis settings
  y_max   = 600  # Axis upper limit
  x_max   = 20
  y_ticks = seq(0, y_max, by = 100)  # Tick marks
  x_ticks = seq(0, x_max, by = 5) 
  
  # Text sizes
  lab_size = 18  # Labels
  tck_size = 12  # Ticks
  lgd_size = 12  # Legend
  
  # Figure size
  fig_size = c(11.23, 8.5)  # Width & height
  
  # ---- Load data ----
  
  # TODO: Plot simulated EIR rather than input EIR?
  
  #load prevalence incidence results
  prev_inc_pth = paste0(o$pth$results, "prev_incidence.txt")
  prev_inc_df  = read.table(prev_inc_pth, header = TRUE, 
                            stringsAsFactors = FALSE)
  
  
  # Format prevalence as percentage for pretty plotting
  prev_inc_df$prevalence = prev_inc_df$prevalence_2to10 * 100
  #prev_inc_df$age_prevalence = prev_inc_df$age_prev * 100
  
  # ---- Generate plotting dataframe with lower and uppper bound ----
  
  
  # Remove variable we will aggregate in plotting dataframe
  plot_df = unique(prev_inc_df %>% dplyr::select(eir,age,access, season,model))
  
  # We want to keep grounpings of EIR, access, and seasonality, and model
  group_by = list(prev_inc_df$eir, prev_inc_df$access,prev_inc_df$age, prev_inc_df$season,prev_inc_df$model)
  prev_inc_df$sev_inc_ppyo = prev_inc_df$nSevere/prev_inc_df$nHost
  
  
  # Group by and aggreagate with mean, min and max to get lower and upper bounds
  tmp = aggregate(prev_inc_df$sev_inc_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, severe = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$sev_inc_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$sev_inc_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  #add prev
  
  tmp =   aggregate(prev_inc_df$prevalence, by = group_by, FUN = mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, prev2to10 = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  plot_df = plot_df[plot_df$age<20,]
  
  
  rm(tmp)
  plot_df$eir=NULL
  
  # ---- Format treatment access and seasonality variables ----
  
  # Convert access to percentage strings
  plot_df$access = paste0(plot_df$access * 100, "%")
  plot_df = plot_df %>% 
    mutate(
      severe = severe*1000,
      lower =lower * 1000,
      upper = upper * 1000,
    )
  # Convert seasonal indices to descriptive strings
  plot_df$season[plot_df$season == 1] = "Perennial"
  plot_df$season[plot_df$season == 2] = "Seasonal"
  
  plot_df$model[plot_df$model == "0000GA"] = "GA"
  plot_df$model[plot_df$model == "0000GP"] = "GP-BO"
  plot_df$model[plot_df$model == "0000GPSG"] = "GPSG-BO"
  
  plot_df$group = paste(plot_df[,"season"],plot_df[,"model"])
  
  # ---- Produce plot ----
  plot_df_perennial = plot_df %>% filter(season=="Perennial")
  plot_df_seasonal = plot_df %>% filter(season=="Seasonal")
  mround <- function(x,base){
    base*round(x/base)
  }
  
  plot_df_seasonal$prev2to10round = mround(plot_df_seasonal$prev2to10,10)
  
  #select PfPR= 3,10,30,50
  plot_df_seasonal = plot_df_seasonal %>%
    mutate(
      prev_cat = case_when(
        prev2to10>=2.5 & prev2to10<=4.4   ~ "3",
        prev2to10>=9.5 & prev2to10<=10.5   ~ "10",
        prev2to10>=28 & prev2to10<=32   ~ "30",
        prev2to10>=47 & prev2to10<=53   ~ "50"
      )
    )
  plot_df_seasonal = plot_df_seasonal %>% 
    na.omit() %>%
    dplyr::select(prev_cat,model,age,severe ,lower,upper)%>%
    group_by(prev_cat,age,model)%>%
    summarise_each(funs(mean,min,max)) %>%
    mutate(prev_cat, as.factor(prev_cat)) %>%
    mutate(prev_cat = fct_relevel(prev_cat,sort))
  plot_df_seasonal$prev_cat2 = factor(plot_df_seasonal$prev_cat, 
                                      levels=c('3','10','30','50'),
                                      labels=c(TeX('PfPR_{2-10}=3\\%'),
                                               TeX('PfPR_{2-10}=10\\%'),
                                               TeX('PfPR_{2-10}=30\\%'),
                                               TeX('PfPR_{2-10}=50\\%')))
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df_seasonal, aes(x = age,y = severe_mean, colour = as.factor(model),group =as.factor(model))) + 
    facet_grid(prev_cat2~model,labeller = label_parsed,scales="free_y") + 
    geom_line() +
    
    geom_ribbon(aes(ymin = lower_min, ymax = upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("GA","GP-BO","GPSG-BO")) + 
    scale_color_manual(values = colours,labels=c("GA","GP-BO","GPSG-BO")) 
  
  # Sort out y axes limits and ticks
  #g = g + scale_y_continuous(limits = c(0, y_max), 
  #                           breaks = y_ticks, 
  #                           labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks))
  
  # Set legend title and axes labels
  g = g + labs(title = "Severe incidence by age (seasonal transmission)",
               colour = "Parameterization",
               fill = "Parameterization", 
               y = TeX("Severe incidence (events per 1000 person-years)"), 
               x = TeX("age"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          title                 = element_text(size = lab_size+2,face="bold"), 
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.justification  = c(0, 1), 
          legend.position       = c("bottom"),
          legend.box.just       = "right", 
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "age_severe_relationship_seasonal.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
  
  
  
  
  #PERENNIAL
  plot_df_perennial$prev2to10round = mround(plot_df_perennial$prev2to10,10)
  
  #select PfPR= 3,10,30,50
  plot_df_perennial = plot_df_perennial %>%
    mutate(
      prev_cat = case_when(
        prev2to10>=2.5 & prev2to10<=4.4   ~ "3",
        prev2to10>=9.5 & prev2to10<=10.5   ~ "10",
        prev2to10>=28 & prev2to10<=32   ~ "30",
        prev2to10>=47 & prev2to10<=53   ~ "50"
      )
    )
  plot_df_perennial = plot_df_perennial %>% 
    na.omit() %>%
    dplyr::select(prev_cat,model,age,severe,lower,upper)%>%
    group_by(prev_cat,age,model)%>%
    summarise_each(funs(mean,min,max)) %>%
    mutate(prev_cat, as.factor(prev_cat)) %>%
    mutate(prev_cat = fct_relevel(prev_cat,sort))
  plot_df_perennial$prev_cat2 = factor(plot_df_perennial$prev_cat, 
                                       levels=c('3','10','30','50'),
                                       labels=c(TeX('PfPR_{2-10}=3\\%'),
                                                TeX('PfPR_{2-10}=10\\%'),
                                                TeX('PfPR_{2-10}=30\\%'),
                                                TeX('PfPR_{2-10}=50\\%')))
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df_perennial, aes(x = age,y = severe_mean, colour = as.factor(model),group =as.factor(model))) + 
    facet_grid(prev_cat2~model,labeller = label_parsed,scales="free_y") + 
    geom_line() +
    
    geom_ribbon(aes(ymin = lower_min, ymax = upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("GA","GP-BO","GPSG-BO")) + 
    scale_color_manual(values = colours,labels=c("GA","GP-BO","GPSG-BO")) 
  
  # Sort out y axes limits and ticks
  #g = g + scale_y_continuous(limits = c(0, y_max), 
  #                           breaks = y_ticks, 
  #                           labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks))
  
  # Set legend title and axes labels
  g = g + labs(title = "Severe incidence by age (perennial transmission)",
               colour = "Parameterization",
               fill = "Parameterization", 
               y = TeX("Severe incidence (events per 1000 person-years)"), 
               x = TeX("age"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          title                 = element_text(size = lab_size+2,face="bold"), 
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.position       = c("bottom"),
          legend.box.just       = "right", 

          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "age_severe_relationship_perennial.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
}


  
# ---------------------------------------------------------
# Plot OpenMalaria death incidence by age
# ---------------------------------------------------------  
plot_age_death = function(o) {
  
  # Produce mortalty by age relationships for Malaria Modelling
  # Consortium. This code is designed to reproduce plots presented
  # in Melissa's four model comparison paper:
  #
  # https://www.sciencedirect.com/science/article/pii/S0140673615007254
  
  message(" - Plotting age - mortality relationship")
  
  # ---- Figure properties ----
  
  # Colours (one colour per user-defined season)
  colours = c("black","#FF7745", "#720026") #orange and burgundy for fits
  #colours = c("black","#2F6089", "#508FC3","#AFCCE4") #blues for ages
  #colours = c("#FF7745", "#720026")  # Orange & burgundy
  
  # Y-axis settings
  y_max   = 130  # Axis upper limit
  x_max   = 20
  y_ticks = seq(0, y_max, by = 25)  # Tick marks
  x_ticks = seq(0, x_max, by = 5) 
  
  # Text sizes
  lab_size = 18  # Labels
  tck_size = 12  # Ticks
  lgd_size = 12  # Legend
  
  # Figure size
  fig_size = c(11.23, 8.5)  # Width & height
  
  # ---- Load data ----
  
  # TODO: Plot simulated EIR rather than input EIR?
  
  #load prevalence incidence results
  prev_inc_pth = paste0(o$pth$results, "prev_incidence.txt")
  prev_inc_df  = read.table(prev_inc_pth, header = TRUE, 
                            stringsAsFactors = FALSE)
  
  
  # Format prevalence as percentage for pretty plotting
  prev_inc_df$prevalence = prev_inc_df$prevalence_2to10 * 100
  #prev_inc_df$age_prevalence = prev_inc_df$age_prev * 100
  
  # ---- Generate plotting dataframe with lower and uppper bound ----
  
  
  # Remove variable we will aggregate in plotting dataframe
  plot_df = unique(prev_inc_df %>% dplyr::select(eir,age,access, season,model))
  
  # We want to keep grounpings of EIR, access, and seasonality, and model
  group_by = list(prev_inc_df$eir, prev_inc_df$access,prev_inc_df$age, prev_inc_df$season,prev_inc_df$model)
  
  prev_inc_df$sev_inc_ppyo = prev_inc_df$nSevere/prev_inc_df$nHost
  prev_inc_df$dir_death_ppyo = prev_inc_df$expectedDirectDeaths/prev_inc_df$nHost
  prev_inc_df$indir_death_ppyo = prev_inc_df$expectedIndirectDeaths/prev_inc_df$nHost
  prev_inc_df$all_death_ppyo = (prev_inc_df$expectedDirectDeaths+prev_inc_df$expectedIndirectDeaths)/prev_inc_df$nHost																				 
  
  
  # Group by and aggreagate with mean, min and max to get lower and upper bounds
  
  #indirect deaths
  tmp = aggregate(prev_inc_df$indir_death_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, indirDeath_incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$indir_death_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, indirDeath_lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$indir_death_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, indirDeath_upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  
  #direct deaths:
  tmp = aggregate(prev_inc_df$dir_death_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, dirDeath_incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$dir_death_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, dirDeath_lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$dir_death_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, dirDeath_upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  #all deaths:
  tmp = aggregate(prev_inc_df$all_death_ppyo, by = group_by, mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, allDeath_incidence = x)
  plot_df = merge(plot_df,tmp, by = c("eir", "access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$all_death_ppyo, by = group_by, min)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, allDeath_lower = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  tmp = aggregate(prev_inc_df$all_death_ppyo, by = group_by, max)%>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, allDeath_upper = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  
  #add prev
  
  
  tmp =   aggregate(prev_inc_df$prevalence, by = group_by, FUN = mean) %>% 
    dplyr::rename(eir=Group.1,access = Group.2, age = Group.3, season = Group.4, model = Group.5, prev2to10 = x)
  plot_df = merge(plot_df,tmp, by = c( "eir","access", "age","season", "model"))
  plot_df = plot_df[plot_df$age<20,]
  
  
  rm(tmp)
  plot_df$eir=NULL
  
  # ---- Format treatment access and seasonality variables ----
  
  # Convert access to percentage strings
  plot_df$access = paste0(plot_df$access * 100, "%")
  
  plot_df = plot_df %>% 
    mutate(
      indirDeath_incidence  = indirDeath_incidence*1000,
      indirDeath_lower  =indirDeath_lower * 1000,
      indirDeath_upper = indirDeath_upper * 1000,
      dirDeath_incidence  = dirDeath_incidence*1000,
      dirDeath_lower  =dirDeath_lower * 1000,
      dirDeath_upper = dirDeath_upper * 1000,
      allDeath_incidence  = allDeath_incidence*1000,
      allDeath_lower  =allDeath_lower * 1000,
      allDeath_upper = allDeath_upper * 1000										
    )
  
  # Convert seasonal indices to descriptive strings
  plot_df$season[plot_df$season == 1] = "Perennial"
  plot_df$season[plot_df$season == 2] = "Seasonal"
  
  plot_df$model[plot_df$model == "0000GA"] = "GA"
  plot_df$model[plot_df$model == "0000GP"] = "GP-BO"
  plot_df$model[plot_df$model == "0000GPSG"] = "GPSG-BO"
  
  plot_df$group = paste(plot_df[,"season"],plot_df[,"model"])
  
  # ---- Produce plot ----
  plot_df_perennial = plot_df %>% filter(season=="Perennial")
  plot_df_seasonal = plot_df %>% filter(season=="Seasonal")
  mround <- function(x,base){
    base*round(x/base)
  }
  
  plot_df_seasonal$prev2to10round = mround(plot_df_seasonal$prev2to10,10)
  
  #select PfPR= 3,10,30,50
  plot_df_seasonal = plot_df_seasonal %>%
    mutate(
      prev_cat = case_when(
        prev2to10>=2.5 & prev2to10<=4.4   ~ "3",
        prev2to10>=9.5 & prev2to10<=10.5   ~ "10",
        prev2to10>=28 & prev2to10<=32   ~ "30",
        prev2to10>=47 & prev2to10<=53   ~ "50"
      )
    )
  plot_df_seasonal = plot_df_seasonal %>% 
    na.omit() %>%
    dplyr::select(prev_cat,model,age,dirDeath_incidence ,allDeath_incidence, allDeath_lower,allDeath_upper,dirDeath_lower,dirDeath_upper)%>%
    group_by(prev_cat,age,model)%>%
    summarise_each(funs(mean,min,max)) %>%
    mutate(prev_cat, as.factor(prev_cat)) %>%
    mutate(prev_cat = fct_relevel(prev_cat,sort))
  plot_df_seasonal$prev_cat2 = factor(plot_df_seasonal$prev_cat, 
                                      levels=c('3','10','30','50'),
                                      labels=c(TeX('PfPR_{2-10}=3\\%'),
                                               TeX('PfPR_{2-10}=10\\%'),
                                               TeX('PfPR_{2-10}=30\\%'),
                                               TeX('PfPR_{2-10}=50\\%')))
  
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df_seasonal, aes(x = age,y = allDeath_incidence_mean, colour = as.factor(model),group =as.factor(model))) + 
    facet_grid(prev_cat2~model,labeller = label_parsed,scales="free_y") + 
    
    geom_line(aes(x = age, y = allDeath_incidence_mean)) +
    geom_line(aes(x = age, y = dirDeath_incidence_mean), linetype="dashed") +											   
    
    
    geom_ribbon(aes(ymin = allDeath_lower_min, ymax = allDeath_upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_ribbon(aes(ymin = dirDeath_lower_min, ymax = dirDeath_upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("GA","GP-BO","GPSG-BO")) + 
    scale_color_manual(values = colours,labels=c("GA","GP-BO","GPSG-BO")) 
  
  # Sort out y axes limits and ticks
  #g = g + scale_y_continuous(limits = c(0, y_max), 
  #                           breaks = y_ticks, 
  #                           labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks))
  
  # Set legend title and axes labels
  g = g + labs(title = "Deaths by age (seasonal transmission)",
               colour = "Parameterization",
               fill = "Parameterization", 
               y = TeX("Deaths (events per 1000 person-years)"), 
               x = TeX("age"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          title                 = element_text(size = lab_size+2,face="bold"), 
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.justification  = c(0, 1), 
          legend.position       = c("bottom"),
          legend.box.just       = "right", 
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "age_death_relationship_seasonal.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
  
  
  
  
  #PERENNIAL
  plot_df_perennial$prev2to10round = mround(plot_df_perennial$prev2to10,10)
  
  
  #select PfPR= 3,10,30,50
  plot_df_perennial = plot_df_perennial %>%
    mutate(
      prev_cat = case_when(
        prev2to10>=2.5 & prev2to10<=4.4   ~ "3",
        prev2to10>=9.5 & prev2to10<=10.5   ~ "10",
        prev2to10>=28 & prev2to10<=32   ~ "30",
        prev2to10>=47 & prev2to10<=53   ~ "50"
      )
    )
  plot_df_perennial = plot_df_perennial %>% 
    na.omit() %>%
    dplyr::select(prev_cat,model,age,dirDeath_incidence ,allDeath_incidence, allDeath_lower,allDeath_upper,dirDeath_lower,dirDeath_upper)%>%
    group_by(prev_cat,age,model)%>%
    summarise_each(funs(mean,min,max)) %>%
    mutate(prev_cat, as.factor(prev_cat)) %>%
    mutate(prev_cat = fct_relevel(prev_cat,sort))
  plot_df_perennial$prev_cat2 = factor(plot_df_perennial$prev_cat, 
                                      levels=c('3','10','30','50'),
                                      labels=c(TeX('PfPR_{2-10}=3\\%'),
                                               TeX('PfPR_{2-10}=10\\%'),
                                               TeX('PfPR_{2-10}=30\\%'),
                                               TeX('PfPR_{2-10}=50\\%')))
  
  # Plot uncerainty bounds with a ribbon, best estimate with thick line
  g = ggplot(plot_df_perennial, aes(x = age,y = allDeath_incidence_mean, colour = as.factor(model),group =as.factor(model))) + 
    facet_grid(prev_cat2~model,labeller = label_parsed,scales="free_y") + 
    
    geom_line(aes(x = age, y = allDeath_incidence_mean)) +
    geom_line(aes(x = age, y = dirDeath_incidence_mean), linetype="dashed") +											   
    
    
    geom_ribbon(aes(ymin = allDeath_lower_min, ymax = allDeath_upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_ribbon(aes(ymin = dirDeath_lower_min, ymax = dirDeath_upper_max, fill = as.factor(model)), 
                linetype = 0, alpha = 0.2, show.legend = FALSE) +
    geom_point(size = 2) 
  
  
  # Manually set pretty colours and linetype
  g = g + scale_fill_manual(values = colours, labels=c("GA","GP-BO","GPSG-BO")) + 
    scale_color_manual(values = colours,labels=c("GA","GP-BO","GPSG-BO")) 
  
  # Sort out y axes limits and ticks
  #g = g + scale_y_continuous(limits = c(0, y_max), 
  #                           breaks = y_ticks, 
  #                           labels = paste0(y_ticks))
  
  #  Sort out x axes limits and ticks
  g = g + scale_x_continuous(limits = c(0, x_max),
                             breaks = x_ticks, 
                             labels = paste0(x_ticks))
  
  # Set legend title and axes labels
  g = g + labs(title = "Deaths by age (perennial transmission)",
               colour = "Parameterization",
               fill = "Parameterization", 
               y = TeX("Deaths (events per 1000 person-years)"), 
               x = TeX("age"))
  
  # Fix up text sizes  
  g = g + 
    theme_classic() + 
    theme(panel.spacing = unit(0, "mm"),
          panel.border = element_rect(fill=NA,color="grey60"),
          strip.text   = element_text(size = lab_size), 
          strip.background = element_rect(fill="grey90",color="grey70"),
          title                 = element_text(size = lab_size+2,face="bold"), 
          axis.title            = element_text(size = lab_size), 
          axis.text             = element_text(size = tck_size), 
          legend.title          = element_text(size = lgd_size,face="bold"), 
          legend.text           = element_text(size = lgd_size),
          legend.justification  = c(0, 1), 
          legend.position       = c("bottom"),
          legend.box.just       = "right", 
          legend.box.background = element_rect(fill = "transparent", colour = NA),
          legend.key            = element_rect(colour = 'grey', size = 0.2)) + 
    guides(
      linetype=guide_legend(keywidth = 3.2, keyheight = 1.4),
      colour=guide_legend(keywidth = 3.2, keyheight = 1.4))
  
  # Save plot to file
  ggsave(paste0(o$pth$figures, "age_death_relationship_perennial.png"), 
         plot = g, width = fig_size[1], height = fig_size[2])
  }


