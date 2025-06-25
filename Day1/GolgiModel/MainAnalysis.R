########## Main Script to run the model analysis #######

# Firstly, let isntall the "epimod" package
# for more details about the package see: https://github.com/qBioTurin/epimod

# library(devtools)
# install_github("qBioTurin/epimod",dependencies=TRUE,ref="dev-merge")

library(epimod)
#downloadContainers()

#### Preliminary steps:
# Load the R packages and function that will be exploited
library(readr)
source('./R_functions/Plot_function.R') # Script encoding the functions to generate the several plots

# Select the name of the model to run the analysis:
# 1) GolgiModel: for the last results
# 2) GolgiModelOLD: for the old results
ModelName = "GolgiModel"

### generate the list of the ordered parameters' name:
Calibration_list <- read_delim("Input/Calibration_list.csv", 
                               delim = ";", escape_double = FALSE, trim_ws = TRUE)
OptNames = c("i_Prkd1X", "i_Pi4kbX", "i_PtdIns", "i_Itga5Itgb1Fn1nMb","i_Sacm1lX", Calibration_list$init[-c(37,38)] )

#######################

# Step 1: generation of the mathematical processes stored in "ModelName".solver
model.generation(net_fname = paste0("Net/",ModelName),
                 transitions_fname = "Net/GenTransitions.cpp")
system(paste0("mv ",ModelName,".* ./Net"))

## Step 2: model analysis using the parameter configurations obtained from the calibration step.

# Movement_GM [0.1634; 0.2426]  
# Movement_MG [0.275; 0.601]
# Activation_PI4KBx [0.005,0.04] 
# FromPtdIns_toPtdIns4P [0.005,0.04] 

optim = c(2, 2, 32, 247.5, 16.875, # i_Prkd1X i_Pi4kbX i_PtdIns i_Itga5Itgb1Fn1nMb i_Sacm1lX
          .6, 0.00053*1.2, 0, #  Input_Ptdins  Input_PRKD1x Input_FN1
          0.00065*9, 0.0000165*0.8, #  ConsuptionMAPK13 ConsuptionMAPK14
          0.00005*0.8,0.001, 0.1, # Sequestration InterGeneration DissGeneration (Intergen >/= 5x DissGen)
          0.001/1.7, 0.0009,0.009,0.01, # Activation_PRKD1x Inactivation_PRDK1X Activation_PI4KBx Inactivation_PI4KBX 0.219773
          0.09, 0.00001, # FromPtdIns_toPtdIns4P FromPtdIns4P_toPtdIns
          0.01*0.5, 0.5,  # Inactivation_SACM1LX Activation_SACM1Lx
          0.05, # Docking
          0.276,0.165,  # Movement_MG Movement_GM
          0.0005*0.01 ,  # Exocytosis
          0.016, 0.007, # Endocytosis EndocytosisR [0.003333333,0.01666667]
          0.15, 0.182, #  RecyclingA RecyclingB
          0.18, # VesicleBudding
          0.0005*1.5, 0.0022/90, # CostMAPK13 CostMAPK14
          0.00025, # ConsuptionPrkd
          0.0034, # Endocytosis_A  [0.003333333,0.01666667]
          0.022,0.005, # ATAC_MN ATAC_NM
          0.022,0.005, # PPFIA1_MN PPFIA1_NM
          0.08,0.02, # PPFIA1_ATAC_MN PPFIA1_ATAC_NM
          0.08) # ExocitosiRetro

names(optim) = OptNames
optim

model.analysis(solver_fname = paste0("./Net/",ModelName,".solver"),
               parameters_fname = "./Input/Calibration_list.csv",
               functions_fname = "./R_functions/Functions.R",
               event_times = c(24,48)*60*60,
               event_function = "EventFun",
               f_time = 72*60*60,
               s_time = 20*60,
               i_time = 0,
               ini_v = optim)#,debug = T )
ModelAnalysisPlot(tracefile = paste0(ModelName,"_analysis/",ModelName,"-analysis-1.trace"))


###### What-if analysis
## 1. varying percentage of PPFIA in Nucleus and Membrane starting from 0.83 -> 0.17

model.analysis(solver_fname = "./Net/GolgiModel.solver",
               parameters_fname = "./Input/WIFanalysis_list.csv",
               functions_fname = "./R_functions/Functions.R",
               f_time =48*60*60,
               s_time = 60,
               i_time = 0,
               ini_v = optim,
               n_config = 50)

WifPlot(folder="GolgiModel_analysis/",
        model="GolgiModel",
        paramWIF = c("ATAC_MN","ATAC_NM"),
        filename = "WhatIF_ATAC_Nucleo.png",
        optim=optim)

#####################################################
save(optim,file = "Input/optim.RData")

targets = c("Mapk13","MAPK14", "PtdIns4Pb", "Sacm1lX", "Itga5Itgb1Fn1nMb", "Itga5Itgb1Fn1nMbR",
            "Ppfia1_Atacn", "Ppfia1n","VesicleFnC")
for(t in targets){
  targetName = t
  save(targetName,file = 'Input/target.RData')
  model.sensitivity(solver_fname = "./Net/GolgiModel.solver",
               parameters_fname = "./Input/Sensitivity_list.csv",
               functions_fname = "./R_functions/Functions.R",
               out_fname = t,
               f_time =72*60*60,
               s_time = 60*60,
               i_time = 0,
               event_times = c(24,48)*60*60,
               event_function = "EventFun",
               n_config = 500,
               parallel_processors = 10,
               target_value = "Target")#,debug = T )
  system(paste0("mv GolgiModel_sensitivity/prcc_",t,".RData ./Sensitivity"))
  
}

# To generate the PRCC plots, run the script PRCCplots.R in the R_functions folder
