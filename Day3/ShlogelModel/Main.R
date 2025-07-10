########## Main Script to run the model analysis #######

# Firstly, let install the "epimod" package
# for more details about the package see: https://github.com/qBioTurin/epimod

# install.packages("devtools")
# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="epimod_pFBA")

library(epimod)
downloadContainers()


start_time <- Sys.time()
model.generation(net_fname = "./Net/Schlogl_reduced.PNPRO")
end_time <- Sys.time()-start_time

### Model Analysis
source("Rfunction/ModelAnalysisPlot.R")

##############################
## Deterministic setting
##############################

model.analysis(solver_fname = "./Schlogl_reduced.solver",
               f_time = 30, # days
               s_time = 1,
               functions_fname = "Rfunction/Functions.R",
               parameters_fname = "Input/Functions_list_Calibration.csv",
               ini_v = c(248, 0.03, 0.0001, 200, 3.5 )
)


AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./Schlogl_reduced_analysis/")

p1 <- AnalysisPlot$list.plX1$plX1.1

##### Calibration ###### 

start_time <- Sys.time()

model.calibration(parameters_fname = "Input/Functions_list_Calibration.csv",
                  functions_fname = "Rfunction/Functions.R",
                  solver_fname = "Schlogl_reduced.solver",
                  reference_data = "Input/reference_data.csv",
                  distance_measure = "msqd" ,
                  f_time = 30, # days
                  s_time = 1, # days
                  # Vectors to control the optimization
                  ini_v = c(248, 0.03, 0.0001, 200, 3 ),
                  ub_v = c(260, 0.05, 0.001, 250, 4 ),
                  lb_v = c(240, 0.02, 0.00009,150, 2 ), max.time = 1
)

end_time <- Sys.time()-start_time

source("Rfunction/CalibrationPlot.R")

plX1

