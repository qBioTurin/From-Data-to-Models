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
  
  Ppfia = 66 #3.849674431
  
#if(is.null(PPFIAperc)){
  PpfiaN = Ppfia * 0.83
  PpfiaM = Ppfia * 0.17
# }else{
#   load(PPFIAperc)
#   PpfiaN = Ppfia * PercVariation[1]
#   PpfiaM = Ppfia * (1-PercVariation[1])
#   PercVariation = PercVariation[-1]
#   save(PercVariation,file=PPFIAperc)
# }

  y_ini["Ppfia1_Atacn"] =	PpfiaN /2
  y_ini["Ppfia1n"] =	PpfiaN /2
  y_ini["Ppfia1_Atacm"] =	PpfiaM /2
  y_ini["Ppfia1m"] =	PpfiaM /2
  
  y_ini["Vesicle_Itga5Itgb1_Rab11bC"] =	 NvescicRab11
  y_ini["Vesicle_FN_Itga5itgb1"] =	 0
  
  y_ini["Itga5Itgb1Fn1Rab21C"] =	 NvescicRab21
  y_ini["Itga5Itgb1G"] =	0
  # y_ini["Sacm1l"] =	33.75
  y_ini["Mapk13"] =	5.38425526
  y_ini["VesicleFnC"] =	269.9
  y_ini["Fn1"] = 600		#705.2216018
  
  ATAC  =  94.7 # 4.4457 
  y_ini["Atac"] = 0.771915719 * ATAC - y_ini["Ppfia1_Atacn"] 
  y_ini["Atacm"] =	0.2280842805 * ATAC - y_ini["Ppfia1_Atacm"] 
  #y_ini["Atac"] =	 ATAC
  
  if(is.null(optim_v)){
    y_ini["PrkdX"] =		runif(1,min = 1,max = 100)
    y_ini["Pi4kbX"] =		runif(1,min = 1,max = 100)
    y_ini["PtdIns"] =		runif(1,min = 1,max = 100)
    y_ini["Itga5Itgb1Fn1nMb"] =		runif(n = 1,min = 100,max = 395)
    PtdIns4P = runif(1, min = 8, max = 32)
  }else{
    #y_ini["PtdIns4P"] =	x[18]
    #y_ini["PrkdX"] =	optim_v[1]
    y_ini["Pi4kbX"] = optim_v[2]
    y_ini["PtdIns"] =	optim_v[3]/2
    y_ini["Itga5Itgb1Fn1nMb"] =	optim_v[4]
    y_ini["Sacm1lX"] =	optim_v[5]
  }
  
  y_ini["MAPK13_PKD"] =		4.39252142/4
  y_ini["PKDx"] =		4.39252142/4
  y_ini["PrkdX"] =   4.39252142/2
  
  y_ini["Pi4kbx"] =		4.39252142 - y_ini["Pi4kbX"]
  y_ini["PtdIns4Pb"] = y_ini["PtdIns"]
  
  ##### 
  y_ini["Itga5Itgb1Fn1nMa"] =		495
  y_ini["Itga5Itgb1Fn1nMbR"] =	495 - y_ini["Itga5Itgb1Fn1nMb"]
  
  y_ini["MAPK14"] = 58.18
  y_ini["Sacm1lx"] = 33.75 - y_ini["Sacm1lX"]
  
  if(!is.null(ini)){
    y_ini = ini
  }
  
  if(scenario == "MaxATAC")
    y_ini = y_ini["Atac"] + y_ini["Atacm"] + y_ini["Ppfia1_Atacm"]+ y_ini["Ppfia1_Atacn"]
  
  return( y_ini )
  #return(ceiling(y_ini) )
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
    Ppfia = 66 * 0.138#7.757895 #0.445805016432915 
  }else if(time == 48*60*60){
    Ppfia = 0.1
  }
  
  PpfiaN = Ppfia * 0.83
  PpfiaM = Ppfia * 0.17
  
  newmarking["Ppfia1_Atacm"] =	PpfiaM * marking["Ppfia1_Atacm"] / (marking["Ppfia1_Atacn"] + marking["Ppfia1_Atacm"] )
  newmarking["Ppfia1_Atacn"] =	PpfiaN * marking["Ppfia1_Atacn"] / (marking["Ppfia1_Atacn"] + marking["Ppfia1_Atacm"] )
  newmarking["Ppfia1m"] =	PpfiaM * marking["Ppfia1m"] / (marking["Ppfia1n"] + marking["Ppfia1m"] )
  newmarking["Ppfia1n"] =	PpfiaN * marking["Ppfia1n"] / (marking["Ppfia1n"] + marking["Ppfia1m"] )
  newmarking["Atac"] = marking["Atac"] + (marking["Ppfia1_Atacn"] - newmarking["Ppfia1_Atacn"] )#+ 20) 
  newmarking["Atacm"] = marking["Atacm"] + (marking["Ppfia1_Atacm"] - newmarking["Ppfia1_Atacm"] )#+ 20) 
  
  #return(ceiling(newmarking))
  return(newmarking) 
}

ErrorDistance<-function(reference, output)
{
  name.place = as.character(factor(reference[,1]))
  value = as.numeric(as.character(factor(reference[,2])))
  names(value) = name.place
  output[,"Pi4kbxTot"]  = output[,"Pi4kbX"] + output[,"Pi4kbx"]
  output[,"PKDxTot"]  = output[,"PrkdX"] + output[,"PKDx"]
  output[,"Sacm1lTot"]  = output[,"Sacm1lX"] + output[,"Sacm1lx"]
  
  # Per controllare quali posti hai nel reference che non sono contenuti nella rete:
  # name.place[! name.place %in% colnames(output)]
  
  if(length(output[,1]) < 6 ) return(1e19)
  # length(output[,1]) -> vuol dire che vai a prendere l'ultimo elemento della traccia!
  output.value = output[length(output[,1]), name.place]
  
  names(output.value) = name.place
  # Prendo il primo valore di PtdIns4P per moltiplicarlo per la percentuale data in reference per poi fare a differenza
  # con il valore finale (perche' volgiamo che sia diminuita rispetto una data percentuale)
  v = output[1,"PtdIns4Pb"] * (value["PtdIns4Pb"])
  value["PtdIns4Pb"] = as.numeric(v)
  
  # Prendo il valore di PrkdX dopo 24h per moltiplicarlo per la percentuale data in reference per poi fare a differenza
  # con il valore finale (perche' volgiamo che sia diminuita rispetto una data percentuale)
  v = output[which(output$Time==60*60*24),"PrkdX"] * (value["PrkdX"])
  value["PrkdX"] = as.numeric(v)
  
  # Prendo il primo valore di Atacm per moltiplicarlo per la percentuale data in reference per poi fare a differenza
  # con il valore finale (perche' volgiamo che sia diminuita rispetto una data percentuale)
  # ATAC = output[1,"Atacm"] + output[1,"Atacn"] 
  # v = ATAC * (value["Atacm"])
  # value["Atacm"] = as.numeric(v)
  # v = ATAC * (1-value["Atacm"])
  # value["Atacn"] = as.numeric(v)
  
  ######
  output.value["Itga5Itgb1Fn1nMb"] =	output.value["Itga5Itgb1Fn1nMb"] + output[length(output[,1]), "Itga5Itgb1Fn1nMbR"]
  
  #output.value["Atacm"] =  output.value["Atacm"]/ (output.value["Atacn"] +output.value["Atacm"])
  #output.value["Atacn"] =  output.value["Atacn"]/ (output.value["Atacn"] +output.value["Atacm"])
  
  weight = rep(1,length(name.place))
  names(weight) = name.place
  #weight[c("Atacn","Atacm")] = 2
  #weight[c("PtdIns4P")] = 4
  
  ### The last error that I want to calculate regards the fact that I dont want the trajectories of some
  #   some places to go quickly to zero!!
  #PtdIns4Pb,PrkdX,Pi4kbX  devono essere intorno al valore iniziale
  name24h = c("PrkdX","Pi4kbX","PtdIns4Pb", 
              "Ppfia1PtprfAtac","Ppfia1Ptprf", "Mapk13acH3K14",
              "VesicleFnC","Vesicle_Itga5Itgb1_Rab11bC","Atac",
              "PtdIns","Mapk14",
              "Pi4kbX", "Pi4kbx","Fn1")
  values = output[1,name24h]
  
  err24 = sum(abs(output[which(output$Time==60*60*24),name24h ] - values) /values )
  
  #errPi4kb = sum( abs( (output[,"Pi4kbX"] + output[,"Pi4kbx"]) - (output[1,"Pi4kbX"] + output[1,"Pi4kbx"]) )/ (output[1,"Pi4kbX"] + output[1,"Pi4kbx"]) )
  #errPrkd1 = sum( abs( (output[,"PrkdX"] + output[,"PKDx"]) - (output[1,"PrkdX"] + output[1,"PKDx"]) )/ (output[1,"PrkdX"] + output[1,"PKDx"]) )
  #errPtdIns =sum( abs( (output[,"PtdIns4Pb"] + output[,"PtdIns"]) - (output[1,"PtdIns4Pb"] + output[1,"PtdIns"]) )/ (output[1,"PtdIns4Pb"] + output[1,"PtdIns"]) )
  
  ###
  
  if (length(output.value[output.value < -1e-5]) == 0)
    err <- sum( ( abs(output.value[names(output.value)] - value[names(output.value)])/(value[names(output.value)]) ))
  #err <- 1/length(names(output.value))*sum( ( (output.value[names(output.value)] - value[names(output.value)]) )^2)
  else
    err <- 1e18 # NON voglio valori negativi!
  
  return(err + err24) # + 2*(errPi4kb + errPrkd1 + errPtdIns) )
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
