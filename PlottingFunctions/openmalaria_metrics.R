############################################################
# OPENMALARIA METRICS
#
# Convert OpenMalaria output event indicators to human-
# understandable(ish) strings.
#
# Example usage:
# event_list = events_list()
# event_idx  = c(event_list[["nHost"]], event_list[["nPatent"]])
#
# See: https://github.com/SwissTPH/openmalaria/wiki/XmlMonitoring
############################################################

openmalaria_metrics = function() {
  
  # All event indices 
  om_event <- list()
  om_event[["nHost"]] 							            = 0
  om_event[["nInfect"]] 						            = 1
  om_event[["nExpectd"]] 						            = 2
  om_event[["nPatent"]] 						            = 3   
  om_event[["sumlogDens"]] 						          = 5
  om_event[["totalInfs"]] 					       	    = 6
  om_event[["nTransmit"]] 						          = 7
  om_event[["totalPatentInf"]] 					        = 8
  om_event[["removed"]] 						            = 9
  om_event[["sumPyrogenThresh"]] 				        = 10
  om_event[["nTreatments1"]] 					          = 11
  om_event[["nTreatments2"]] 					          = 12
  om_event[["nTreatments3"]] 					          = 13
  om_event[["nUncomp"]] 						            = 14
  om_event[["nSevere"]] 						            = 15
  om_event[["nSeq"]] 							              = 16
  om_event[["nHospitalDeaths"]] 				        = 17
  om_event[["nIndDeaths"]] 						          = 18
  om_event[["nDirDeaths"]] 						          = 19
  om_event[["nEPIVaccinations"]]  				      = 20
  om_event[["allCauseIMR"]] 					          = 21
  om_event[["nMassVaccinations"]]               = 22
  om_event[["nHospitalRecovs"]] 				        = 23
  om_event[["nHospitalSeqs"]]                   = 24
  om_event[["nIPTDoses"]] 						          = 25
  om_event[["annAvgK"]] 						            = 26
  om_event[["nNMFever"]] 						            = 27
  om_event[["removed"]] 						            = 28
  om_event[["removed"]]                         = 29
  om_event[["innoculationsPerAgeGroup"]]        = 30
  om_event[["Vector_Nv0"]] 						          = 31
  om_event[["Vector_Nv"]] 						          = 32
  om_event[["Vector_Ov"]] 						          = 33
  om_event[["Vector_Sv"]] 						          = 34
  om_event[["inputEIR"]] 						            = 35
  om_event[["simulatedEIR"]] 					          = 36
  om_event[["removed"]] 						            = 37
  om_event[["removed"]] 						            = 38
  om_event[["Clinical_RDTs"]] 					        = 39
  om_event[["Clinical_DrugUsage"]] 				      = 40
  om_event[["Clinical_FirstDayDeaths"]] 	      = 41
  om_event[["Clinical_HospitalFirstDayDeaths"]] = 42
  om_event[["nNewInfections"]] 					        = 43
  om_event[["nMassITNs"]] 						          = 44
  om_event[["nEPI_ITNs"]] 						          = 45
  om_event[["nMassIRS"]] 						            = 46
  om_event[["nMassVA"]] 						            = 47
  om_event[["Clinical_Microscopy"]] 			      = 48
  om_event[["Clinical_DrugUsageIV"]]			      = 49
  om_event[["nAddedToCohort"]] 					        = 50
  om_event[["nRemovedFromCohort"]] 				      = 51
  om_event[["nMDAs"]] 							            = 52
  om_event[["nNmfDeaths"]] 						          = 53
  om_event[["nAntibioticTreatments"]] 		      = 54
  om_event[["nMassScreenings"]] 				        = 55
  om_event[["nMassGVI"]] 						            = 56
  om_event[["nCtsIRS"]] 						            = 57
  om_event[["nCtsGVI"]] 						            = 58
  om_event[["nCtsMDA"]] 						            = 59
  om_event[["nCtsScreenings"]] 				  	      = 60
  om_event[["nSubPopRemovalTooOld"]] 			      = 61
  om_event[["nSubPopRemovalFirstEvent"]] 	      = 62
  om_event[["nLiverStageTreatments"]] 		      = 63
  om_event[["nTreatDiagnostics"]] 				      = 64
  om_event[["nMassRecruitOnly"]] 				        = 65
  om_event[["nCtsRecruitOnly"]] 				        = 66
  om_event[["nTreatDeployments"]] 				      = 67
  om_event[["sumAge"]] 							            = 68
  om_event[["nInfectByGenotype"]] 				      = 69
  om_event[["nPatentByGenotype"]] 				      = 70
  om_event[["logDensByGenotype"]] 				      = 71
  om_event[["nHostDrugConcNonZero"]] 			      = 72
  om_event[["sumLogDrugConcNonZero"]] 		      = 73
  om_event[["expectedDirectDeaths"]] 			      = 74
  om_event[["expectedHospitalDeaths"]] 		      = 75
  om_event[["expectedIndirectDeaths"]] 		      = 76
  om_event[["expectedSequelae"]] 				        = 77
  om_event[["expectedSevere"]] 					        = 78
  
  # Define some useful combinations
  om_event[["allDeaths"]] 						          = c(18, 19) 
  om_event[["sumUncompSev"]] 					          = c(14, 15) 
  om_event[["allHospitalisations"]] 			      = c(17, 23, 24) 
  
  return(om_event)
}

