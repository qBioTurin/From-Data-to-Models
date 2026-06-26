InitGeneration <- function(n_file,optim_v=NULL, ini = NULL , scenario = "",optim_file = NULL)
{
  yini.names <- readRDS(n_file)
  if( !is.null(optim_file) ){
    load(optim_file)
    optim_v = optim
  }  
  y_ini <- rep(0,length(yini.names))
  names(y_ini) <- yini.names
  NvescicRab21 = 3800*0.29929029 # nCentroid 1:12 video * freq
  NvescicRab11 = 6500*0.29929029 # nCentroid 13:23 video * freq
  #un possibile range del numero di centroidi è 5000-9000
  
  PPFIA = 66 

  PPFIAN = PPFIA * 0.83
  PPFIAM = PPFIA * 0.17
  
  y_ini["PPFIA1_NCn"] =	PPFIAN /2
  y_ini["PPFIA1n"] =	PPFIAN /2
  y_ini["PPFIA1_NCm"] =	PPFIAM /2
  y_ini["PPFIA1m"] =	PPFIAM /2
  
  y_ini["Vesicle_Integrin_Rab11b"] =	 NvescicRab11
  y_ini["Vesicle_Integrin_FN"] =	 0
  
  y_ini["Vesicle_Integrin_Rab21"] =	 NvescicRab21
  y_ini["Vesicle_Integrin"] =	0
  y_ini["MAPKA"] =	5.38425526
  y_ini["Vesicle_FN"] =	269.9
  
  NC  =  94.7 # 4.4457 
  y_ini["NCn"] = 0.771915719 * NC - y_ini["PPFIA1_NCn"] 
  y_ini["NCm"] =	0.2280842805 * NC - y_ini["PPFIA1_NCm"] 
  
  y_ini["LipidKX"] = optim_v[2]
  y_ini["Lipid"] =	optim_v[3]/2
  y_ini["Integrin_FN_Mb"] =	optim_v[4]
  y_ini["SACM1LX"] =	optim_v[5]
  
  
  y_ini["MAPKA_PKD"] =		4.39252142/4
  y_ini["PKDx"] =		4.39252142/4
  y_ini["PKDX"] =   4.39252142/2
  
  y_ini["LipidKx"] =		4.39252142 - y_ini["LipidKX"]
  y_ini["P_lipid"] = y_ini["Lipid"]
  
  ##### 
  y_ini["Integrin_Ma"] =		495
  y_ini["Integrin_M_R"] =	495 - y_ini["Integrin_FN_Mb"]
  
  y_ini["MAPKB"] = 58.18
  y_ini["SACM1Lx"] = 33.75 - y_ini["SACM1LX"]
  
  if(!is.null(ini)){
    y_ini = ini
  }

  if(scenario == "MaxNC")
    y_ini = y_ini["NCn"] + y_ini["NCm"] + y_ini["PPFIA1_NCm"]+ y_ini["PPFIA1_NCn"]
  
  return( y_ini )
}


calibration.function<- function(optim_v,id){
  return(optim_v[id])
}

InterGeneration<- function(min,max){
  p = runif(1,min = min,max = max)
  return(p)
}

DissGeneration<- function(min,max){
  p = runif(1,min = min,max = max)
  return(p)
}

EventFun = function(marking,time){
  newmarking <- marking
  names(newmarking) <- readRDS("/home/docker/data/Input/NAMES.RDS")
  
  if (time == 24*60*60){
    PPFIA = 66 * 0.138#7.757895 #0.445805016432915 
  }else if(time == 48*60*60){
    PPFIA = 0.1
  }
  
  PPFIAN = PPFIA * 0.83
  PPFIAM = PPFIA * 0.17
  
  newmarking["PPFIA1_NCm"] =	PPFIAM * marking["PPFIA1_NCm"] / (marking["PPFIA1_NCn"] + marking["PPFIA1_NCm"] )
  newmarking["PPFIA1_NCn"] =	PPFIAN * marking["PPFIA1_NCn"] / (marking["PPFIA1_NCn"] + marking["PPFIA1_NCm"] )
  newmarking["PPFIA1m"] =	PPFIAM * marking["PPFIA1m"] / (marking["PPFIA1n"] + marking["PPFIA1m"] )
  newmarking["PPFIA1n"] =	PPFIAN * marking["PPFIA1n"] / (marking["PPFIA1n"] + marking["PPFIA1m"] )
  newmarking["NCn"] = marking["NCn"] + (marking["PPFIA1_NCn"] - newmarking["PPFIA1_NCn"] )#+ 20) 
  newmarking["NCm"] = marking["NCm"] + (marking["PPFIA1_NCm"] - newmarking["PPFIA1_NCm"] )#+ 20) 
  
  #return(ceiling(newmarking))
  return(newmarking) 
}

ErrorDistance <- function(reference, output) {
  ## 1. Extract reference values
  colnames(reference) = c("Place", "V")
  reference$Place <- as.character(factor(reference[, 1]))
  reference$V <- as.numeric(as.character(factor(reference[, 2])))
  
  
  ## 2. Compute total species
  output[,"LipidKTot"]  = output[,"LipidKX"] + output[,"LipidKx"]
  output[,"PKDTot"]  = output[,"PKDX"] + output[,"PKDx"] + output[,"MAPKA_PKD"]
  output[,"LipidTot"]  = output[,"P_lipid"] + output[,"Lipid"]
  output[,"SACM1LTot"]  = output[,"SACM1LX"] + output[,"SACM1Lx"]
  output[,"Integrin_MbTot"]  = output[,"Integrin_FN_Mb"] + output[,"Integrin_M_R"]
  output[,"NCnTot"]  = output[,"NCn"] + output[,"PPFIA1_NCn"]
  output[,"NCmTot"]  = output[,"NCm"] + output[,"PPFIA1_NCm"]
  
  ## 3. Check if simulation output is long enough
  if (nrow(output) < 6) return(1e19)
  
  ## 4. Extract final values for comparison
  output.value <- output[nrow(output), name.place]
  names(output.value) <- name.place
  
  ## 5. Handle special transformations for some targets
  # PtdIns4Pb (relative drop compared to initial value)
  v = output[output$Time == 0 & output$Place =="P_lipid","V"] * (reference[reference$Place == "P_lipid","V"])
  reference[reference$Place == "P_lipid", "V"] <- as.numeric(v)
  v = output[output$Time == 0 & output$Place =="PKDX","V"] * (reference[reference$Place == "PKDX","V"])
  reference[reference$Place == "PKDX", "V"] <- as.numeric(v)
  
  
  NCn = output[output$Time == 0 & output$Place == "NCnTot","V"] + output[output$Time == 0 & output$Place == "NCmTot","V"]
  reference[reference$Place == "NCnTot", "V"] = NCn * reference[reference$Place == "NCnTot", "V"]
  reference[reference$Place == "NCmTot", "V"] = NCn * reference[reference$Place == "NCmTot", "V"]
  
  ## 7. Penalize collapse of some signals at 24h
  name24h <- c(
    "PKDX", "LipidKX", "P_lipid", "PPFIA1_NCn", "PPFIA1n",
    "MAPKA", "Vesicle_FN", "Vesicle_Integrin_Rab11b", "NCn",
    "Lipid", "MAPKB", "LipidKx"
  )
  values <- output[1, name24h]
  final24h <- output[output$Time == 60*60*24, name24h]
  err24 <- sum(abs(final24h - values) / values)
  
  ## 8. Main error based on relative difference
  if (all(output.value >= -1e-5)) {
    err <- sum(abs(output.value - value) / value)
  } else {
    err <- 1e18  # Penalize negative concentrations
  }
  
  ## 9. Return total distance
  return(err + err24)
}




# function for the whatif analysis

whatif.function<- function(optim_v,id,type){
  if(type == "slow-down")
    v =optim_v[id]*(runif(n=1,min=0,max = 1) )
  else
    v=optim_v[id]*(runif(n=1,min=1,max = 3) )
  
  # return(optim_v[id]*perc)
  return( v )
  
}

sensitive.function<- function(id,n_file){
  load(n_file)
  a =  optim[id]*.4
  b = optim[id]*1.6
  v = runif(n=1,a,b) 
  # return(optim_v[id]*perc)
  return( v )
  
}


Target<-function(output)
{
  load("/home/docker/data/Input/target.RData")
  print(targetName)
  ret <- output[,targetName]
  return(as.data.frame(ret))
}
