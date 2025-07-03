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
               f_time = 100, # days
               s_time = 1
)


AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./Schlogl_reduced_analysis/")

p1 <- AnalysisPlot$list.plX1$plX1.1


# changing the initial marking of the model through model.analysis

model.analysis(solver_fname = "./Schlogl_reduced.solver",
               parameters_fname = "./Input/Functions_list_ModelAnalysis.csv",
               f_time = 100, # days
               s_time = 1
)

AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./Schlogl_reduced_analysis/")
p2 <- AnalysisPlot$list.plX1$plX1.1

library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1)

## Stochastic simulations

model.analysis(solver_fname = "./Schlogl_reduced.solver",
               parameters_fname = "./Input/Functions_list_ModelAnalysis.csv",
               solver_type = "SSA",
               n_run = 100,
               parallel_processors = 2,
               f_time = 100, # days
               s_time = 1
)

AnalysisPlot = ModelAnalysisPlot(solverTraces_path = "./Schlogl_reduced_analysis",
                                 Stoch = T)
AnalysisPlot$list.plX1


### Sensitivity analysis

##########################################################
## Simple version where only the initial marking of the X1 place vary. ##
##########################################################

start_time <- Sys.time()
sensitivity<-model.sensitivity(n_config = 100,
                               parameters_fname = "Input/Functions_list_sensitivity.csv", 
                               solver_fname = "Schlogl_reduced.solver",
                               functions_fname = "Rfunction/Functions.R",
                               target_value =  "Target" ,
                               f_time = 100, # days
                               s_time = 1, # days      
                               parallel_processors = 2
)

end_time <- Sys.time()-start_time

