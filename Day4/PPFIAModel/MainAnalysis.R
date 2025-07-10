########## Main Script to run the model analysis #######

# Firstly, let isntall the "epimod" package
# for more details about the package see: https://github.com/qBioTurin/epimod

# library(devtools)
# install_github("qBioTurin/epimod",dependencies=TRUE,ref="epimod-pFBA")

library(epimod)
#downloadContainers()

#### Preliminary steps:
# Load the R packages and function that will be exploited
library(readr)
source('./R_functions/Plot_function.R') # Script encoding the functions to generate the several plots

# Select the name of the model to run the analysis:
ModelName = "PPFIA1Model"

# Step 1: generation of the mathematical processes stored in "ModelName".solver
model.generation(net_fname = paste0("Net/",ModelName,".PNPRO"))

system(paste0("mv ",ModelName,".* ./Net"))

## Step 2: model analysis using the parameter configurations obtained from the calibration step.

# Movement_GM [0.1634; 0.2426]  
# Movement_MG [0.275; 0.601]
# Activation_PI4KBx [0.005,0.04] 
# FromPtdIns_toPtdIns4P [0.005,0.04] 

optim = c(2, 2, 32, 247.5, 16.875, # i_Prkd1X i_lipidkX i_lipid i_Itga5Itgb1Fn1nMb i_Sacm1lX
          .6, 0.00053*1.2, 0, #  Input_Lipid  Input_PRKD1x Input_FN1
          0.00065*9, 0.0000165*0.8, #  ConsuptionMAPKA ConsuptionMAPKB
          0.00005*0.8,0.001, 0.1, # Sequestration InterGeneration DissGeneration (Intergen >/= 5x DissGen)
          0.001/1.7, 0.0009,0.009,0.01, # Activation_PRKD1x Inactivation_PRDK1X Activation_LipidKx Inactivation_LipidKX 
          0.09, 0.00001, # Fromlipid_toP_lipid FromP_lipid_tolipid
          0.01*0.5, 0.5,  # Inactivation_SACM1LX Activation_SACM1Lx
          0.05, # Docking
          0.276,0.165,  # Movement_MG Movement_GM
          0.0005*0.01 ,  # Exocytosis
          0.016, 0.007, # Endocytosis EndocytosisR [0.003333333,0.01666667]
          0.15, 0.182, #  RecyclingA RecyclingB
          0.18, # VesicleBudding
          0.0005*1.5, 0.0022/60, # CostMAPKA CostMAPKB
          0.00025, # ConsuptionPrkd
          0.0034, # Endocytosis_A  [0.003333333,0.01666667]
          0.022,0.005, # NC_MN NC_NM
          0.022,0.005, # PPFIA1_MN PPFIA1_NM
          0.08,0.02, # PPFIA1_NC_MN PPFIA1_NC_NM
          0.08) # ExocitosiRetro


model.analysis(solver_fname = paste0("./Net/",ModelName,".solver"),
               parameters_fname = "./Input/paramsList.csv",
               functions_fname = "./R_functions/Functions.R",
               event_times = c(24,48)*60*60,
               event_function = "EventFun",
               f_time = 72*60*60,
               s_time = 20*60,
               i_time = 0,
               ini_v = optim)#,debug = T )
ModelAnalysisPlot(tracefile = paste0(ModelName,"_analysis/",ModelName,"-analysis-1.trace"))

