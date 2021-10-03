#=============================================================================
# function_UncompAccess_5dayProb.R 
#
# R-code function to take fever recall data, eg DHS proportion of malaria fevers 
# treated and convert to 5 day probability of treatment for OpenMalaria.
#
# This code is an example of how MIS and DHS reported treatment of fever proportions 
# can be converted for use in OpenMalaria. Olivier Briet, September 2012. It uses
# data from Crowell et al., submitted in 2012 to PloS ONE.
#
# NOTES:
#  1) Input must be in percentages
#  2) pamd = ProportionAntiMalarialDrugs
#
#  Created by : Olivier Briet and Melissa Penny
#  Originally created : September 2012
#
#  Last modified : 21/11/2013
#			  by : Melissa Penny (melissa.penny@unibas.ch)
#=============================================================================

convert_access <- function(pamd) {
  
  # Digitized data from Crowell et al., figure 5a
  # pt = proportion treated, x = daily treatment probability, y = proportion of fever episodes treated
  pt.x <- c(0, 0.005427408, 0.01492537, 0.02170963, 0.02985075, 0.04070556, 0.05291723, 0.0651289, 0.08005427, 0.09497965, 0.1071913, 0.1180461, 0.1275441, 0.1383989, 0.1506106, 0.165536, 0.1804613, 0.1940299, 0.2075984, 0.2211669, 0.238806, 0.256445, 0.2767978, 0.2957938, 0.3107191, 0.3242877, 0.339213, 0.3622795, 0.3812754, 0.4029851, 0.431479, 0.4545455, 0.4776119, 0.4993216, 0.5237449, 0.5522388, 0.5766621, 0.5956581, 0.6200814, 0.6445047, 0.6757123, 0.7001357, 0.7177748, 0.7476255, 0.78019, 0.807327, 0.8426052, 0.8751696, 0.9104478, 0.9416554, 0.9796472, 0.9986431) 
  pt.y <- c(0.003125, 0.0359375, 0.09375, 0.15, 0.19375, 0.24375, 0.2984375, 0.353125, 0.4078125, 0.4609375, 0.5015625, 0.5359375, 0.56875, 0.6078125, 0.634375, 0.6546875, 0.6671875, 0.6921875, 0.7203125, 0.74375, 0.775, 0.7953125, 0.809375, 0.8203125, 0.8328125, 0.853125, 0.8703125, 0.8890625, 0.903125, 0.9140625, 0.91875, 0.9203125, 0.93125, 0.94375, 0.95625, 0.9640625, 0.9640625, 0.9609375, 0.9671875, 0.9796875, 0.98125, 0.98125, 0.9875, 0.9921875, 0.99375, 0.99375, 0.9953125, 0.9984375, 1, 1, 1, 1) 
  
  # Data from Crowell et al., figure 7c black line
  probtreat_daily <- c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1)
  probtreat_5daily <- c(0, 0.05880144, 0.1157784, 0.1733548, 0.2299865, 0.2818799, 0.329691, 0.3801282, 0.4195029, 0.4652633, 0.5091863, 0.5457407, 0.5781565, 0.6112154, 0.6461608, 0.668167, 0.6951173, 0.7227968, 0.74725, 0.7697838, 0.793149, 0.8094841, 0.8234381, 0.8405304, 0.8572314, 0.8700262, 0.8853289, 0.8971165, 0.908916, 0.9172805, 0.9278921, 0.9377259, 0.9457713, 0.9543433, 0.9645463, 0.9701525, 0.9760008, 0.983826, 0.9891665, 0.9944932, 1) 
  
  # Daily prob of treatment prediction based on pamd
  # Make a smoothing spline, make sure it goes through c(0,0) and c(1,1).
  sp.pamd <- smooth.spline(c(rep(0, 100), pt.y, rep(1, 100)), c(rep(0, 100), pt.x, rep(1, 100)), spar = 0.7)
  pred.pamd <- predict(sp.pamd, x = pamd / 100)
  
  # Conversion of daily prob of treatment to 5-daily prob of treatment prediction
  sp.onetofivedaily <- smooth.spline(c(rep(0, 100), probtreat_daily, rep(1, 100)), c(rep(0, 100), probtreat_5daily, rep(1, 100)), spar = 0.6)
  pred.pamd.5daily <- predict(sp.onetofivedaily, x = pred.pamd$y)
  
  # Prepare output
  pred_5daily_prob = pred.pamd.5daily$y
  
  return(pred_5daily_prob)
}

