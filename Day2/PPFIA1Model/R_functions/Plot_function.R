library(cowplot)
library(ggplot2)
library(ggthemes)

ModelAnalysisPlot <-function(tracefile){
  # PlaceNotToPlot = c("Lysosomes",
  #                    "Sacm1lTot","MAPKA_PKD","LipidkxTot","PtdinsTrash","Fn1")
  PlaceNotToPlot = c("Fn1")
  # PlaceNotToPlot = c("Itga5Itgb1Fn1Rab21C",
  #                    "Itga5Itgb1Fn1nMb",
  #                    "Lysosomes",
  #                    "Itga5Itgb1G",
  #                    "Itga5Itgb1Fn1nMa",
  #                    "Itga5Itgb1Fn1nMbR",
  #                    "Vesicle_FN_Itga5itgb1","MAPKA_PKD")
  
  trace=read.csv(tracefile, sep = "")
  n_sim_tot<-table(trace$Time)
  n_sim <- n_sim_tot[1]
  time_delete<-as.numeric(names(n_sim_tot[n_sim_tot!=n_sim_tot[1]]))
  
  if(length(time_delete)!=0) trace = trace[which(trace$Time!=time_delete),]
  
  # Reference
  #reference <- as.data.frame(t(read.csv("./Input/reference.csv", header = FALSE, sep = "")))
  
  # New reference with multiple points!!!!
  library(readxl)
  Valori_finali <- read_excel("./Input/Valori finali.xlsx", 
                              col_names = FALSE)
  reference <- as.data.frame(Valori_finali[,2:3])
  
  name.place = as.character(factor(reference[,1]))
  value = as.numeric(as.character(factor(reference[,2])))
  names(value) = name.place
  Ref.DF = data.frame(Place = rep(name.place,n_sim), V = rep(value,n_sim) )
  ###
  
  trace[,"LipidkxTot"]  = trace[,"lipidkX"] + trace[,"Lipidkx"]
  trace[,"PKDxTot"]  = trace[,"PrkdX"] + trace[,"PKDx"] + trace[,"MAPKA_PKD"]
  trace[,"PtdInsTot"]  = trace[,"P_lipid"] + trace[,"lipid"]
  trace[,"Sacm1lTot"]  = trace[,"Sacm1lX"] + trace[,"Sacm1lx"]
  trace[,"Itga5Itgb1Fn1nMTot"]  = trace[,"Itga5Itgb1Fn1nMb"] + trace[,"Itga5Itgb1Fn1nMbR"]
  trace[,"NCTot"]  = trace[,"NCn"] + trace[,"Ppfia1_NCn"]
  trace[,"NCmTot"]  = trace[,"NCm"] + trace[,"Ppfia1_NCm"]

  if(length(unique(trace$Time)) < 72 )
    warning("Trace with number of points less than 70")
  
  trace$ID <- rep(1:n_sim[1],each = length(unique(trace$Time)) )
  trace.final <-  lapply(colnames(trace)[-which( colnames(trace)%in% c("ID","Time",PlaceNotToPlot))],function(c){
    # if(c == "Itga5Itgb1Fn1nMb")
    # {
    #   df1 = data.frame(V= trace[,c] + trace[,"Itga5Itgb1Fn1nMbR"], Time=trace$Time,ID = trace$ID,Place=c )
    #   df2 = data.frame(V= trace[,c], Time=trace$Time,ID = trace$ID,Place="Itga5Itgb1Fn1nMb_NotSum" )
    #   df = rbind(df1,df2)
    # }else{
    df = data.frame(V=trace[,c],Time=trace$Time,ID = trace$ID,Place=c )
    #}
    return( df )
  })
  
  traces <- do.call("rbind",trace.final)
  
  ##### Modifica per plottare i reference con percentuali o somme in maniera corretta:
  #Ref.DF <- rbind(Ref.DF, data.frame(Place="PKDx",V = 0))
  v = traces[traces$Time == 0 & traces$Place =="P_lipid","V"] * (Ref.DF[Ref.DF$Place == "P_lipid","V"])
  Ref.DF[Ref.DF$Place == "P_lipid", "V"] <- as.numeric(v)
  v = traces[traces$Time == 0 & traces$Place =="PrkdX","V"] * (Ref.DF[Ref.DF$Place == "PrkdX","V"])
  Ref.DF[Ref.DF$Place == "PrkdX", "V"] <- as.numeric(v)
  # v = traces[traces$Time == 0 & traces$Place =="PKDxTot","V"] - Ref.DF[Ref.DF$Place == "PrkdX", "V"]
  # Ref.DF[Ref.DF$Place == "PKDx", "V"] <- as.numeric(v)
  
  # d = data.frame(Place = c("NCn","NCm"), V = NaN )
  # Ref.DF <- rbind(Ref.DF,d)
  NCn = traces[traces$Time == 0 & traces$Place == "NCTot","V"] + traces[traces$Time == 0 & traces$Place == "NCmTot","V"]
  Ref.DF[Ref.DF$Place == "NCTot", "V"] = NCn * Ref.DF[Ref.DF$Place == "NCTot", "V"]
  Ref.DF[Ref.DF$Place == "NCmTot", "V"] = NCn * Ref.DF[Ref.DF$Place == "NCmTot", "V"]
    
  # d <- unique(traces[traces$Place %in% Ref.DF$Place,"Place"])
  # 
  # d = data.frame(Place = d, V = NaN )
  # Ref.DF <- rbind(Ref.DF,d)
  Ref.DF = Ref.DF[!Ref.DF$Place == "Fn1",]
  Ref.DF$Time = max(traces$Time)
  
  ### Green points
  BlueDF = data.frame(Time = 24,
                      Place = c("MapkA","NCn","P_lipid","PrkdX",
                                "LipidkxTot","PKDxTot","PtdInsTot","Sacm1lTot","MAPKB","VesicleFnC"),
                      V = c(5.38,4.445,
                            trace[1,"P_lipid"],trace[1,"PrkdX"],
                            4.39,4.39,32,33.75, 58.18,270)
  )
  ########
  lev =  c(
    "Ppfia1m","Ppfia1n","Ppfia1_NCm","Ppfia1_NCn","NCm","NCn","MapkA","MAPKB",
    "PKDx","MAPKA_PKD","PrkdX","Lipidkx", "lipidkX","lipid","P_lipid","LipidTrash",
    "Sacm1lx","Sacm1lX",
    "VesicleFnC",
    "Vesicle_FN_Itga5itgb1",
    "Itga5Itgb1Fn1nMb","Itga5Itgb1Fn1nMbR","Itga5Itgb1Fn1Rab21C","Lysosomes",
    "Itga5Itgb1G","Itga5Itgb1Fn1nMa","Vesicle_Itga5Itgb1_Rab11bC",
    "PKDxTot","LipidkxTot","PtdInsTot","Sacm1lTot","Itga5Itgb1Fn1nMTot",
    "NCTot","NCmTot"
  )
  lab = c("PPFIA1m","PPFIA1n","PPFIA1_NCm","PPFIA1_NCn","NCm","NCn","MapkA",
          "MAPKB","PKD","Inactive_PKD","PKD*","Lipidkx", "lipidkX","lipid","P_lipid","PI_Buffer",
          "SACM1L","SACM1L*",
          "Vesicle_FN", "Vesicle_Integrin_FN","Integrin_FN_Mb",
          "Integrin_M_R","Vesicle_Integrin_FN_Rab21","Lysosomes",
          "Endocytic_Recycling_Compartment","Integrin_Ma","Vesicle_Integrin_Rab11b",
          "Sum_PKD", "Sum_PI4KB","Sum_PI","Sum_SACM1L", "Sum_Integrins_Mb", "Sum_NCn", "Sum_NCm")

  traces$Area = NA
  traces$Area[traces$Place %in% c("MapkA","MAPKA_PKD","MAPKB","MAPKA_PKD","LipidkxTot", "PKDxTot","LipidkxTot","PtdInsTot","Sacm1lTot",
                                  "PKDx","PrkdX","Lipidkx", "lipidkX","lipid","P_lipid","PtdinsTrash","Sacm1lx","Sacm1lX","Fn1")] = "Golgi Apparatus"
  traces$Area[traces$Place %in% c("NCn","Ppfia1n","Ppfia1_NCn","NCTot")] = "Nucleus"
  traces$Area[traces$Place %in% c("VesicleFnC","Vesicle_Itga5Itgb1_Rab11bC",
                                  "Vesicle_FN_Itga5itgb1",
                                  "Itga5Itgb1Fn1Rab21C",
                                  "Itga5Itgb1G","Lysosomes")] = "Integrin recycling loop"
  traces$Area[traces$Place %in% c("Itga5Itgb1Fn1nMb","Itga5Itgb1Fn1nMbR","Itga5Itgb1Fn1nMTot","NCmTot","NCm",
                                  "Ppfia1m","Ppfia1_NCm")] = "Basal membrane"
  traces$Area[traces$Place %in% c("Itga5Itgb1Fn1nMa")] = "Apical membrane"
  
  traces$V[traces$Place == "Sacm1lTot"]  = 33.75
  
  Ref.DF = Ref.DF[!Ref.DF$Place %in% PlaceNotToPlot,]
  traces$Place =  factor(traces$Place, levels = lev, labels =lab  )
  Ref.DF$Place =  factor(Ref.DF$Place, levels = lev, labels =lab  )
    
  if(length(unique(traces$ID)) == 1)
  {
    library(dplyr)
    pl<-ggplot( )+
      facet_wrap(~ Place,ncol =  5,scales = "free")+
      geom_rect(data=traces%>%
                  select(Place,Area) %>% distinct() ,
                aes(fill = Area),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.3) +
      geom_line(data=traces,
                aes(x=Time/(60*60),y=V))+
      geom_point(data=Ref.DF,
                 aes(x=Time/(60*60),y=V,
                     col="realData"))+
      geom_vline(xintercept = c(24,48),aes(col="RemPPFIA"),linetype = "dashed")+
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=10),
            legend.title=element_text(size=14,face="bold"),
            legend.position="right",
            legend.key.size = unit(1.3, "cm"),
            legend.key.width = unit(1.3,"cm") )+
      labs(x="Hours", y="")+
      theme_bw()+
      scale_fill_manual(values = c("Golgi Apparatus" = "orange",
                                   "Apical membrane"="#FF0000",
                                   "Integrin recycling loop" = "#008080",
                                   "Basal membrane" = "#FF0000",
                                   "Nucleus"="blue"),
                        guide = guide_legend(title = "",
                                             override.aes = list(alpha = .5)))+
      scale_colour_manual(values = c( "realData"= "red" ),
                          label = c("Real Data"),name = '')
    
  }else{ # stocastico
    library(dplyr)
    meanTrace <- traces %>% group_by(Time,Place) %>%
      summarise(Vmean=mean(V)) %>%
      ungroup()
    
    pl<-ggplot( )+
      facet_wrap(~ Place,ncol =  5,scales = "free")+
      geom_line(data=traces,
                aes(x=Time/(60*60),y=V,group = ID),alpha = .2,col = "grey")+
      geom_line(data=meanTrace,
                aes(x=Time/(60*60),y=Vmean),col = "purple")+
      geom_violin(data=Ref.DF[!Ref.DF$Place %in% PlaceNotToPlot,],
                  aes(x=Time/(60*60),y=V),
                  col="red")+
      geom_vline(xintercept = c(24,48),col="blue",linetype = "dashed")+
      geom_point(data=BlueDF[!BlueDF$Place %in% PlaceNotToPlot,],
                 aes(x=Time,y=V),
                 col="green")+
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=14,face="bold"),
            legend.text=element_text(size=10),
            legend.title=element_text(size=14,face="bold"),
            legend.position="right",
            legend.key.size = unit(1.3, "cm"),
            legend.key.width = unit(1.3,"cm") )+
      labs(x="Hours", y="")
  }
  
  ggsave(pl,filename = "ModelAnalysis.png",height = 15,width = 20)
  return(pl)
}


SensitivityPlot <-function(rank=T,folder = "results_sensitivity_analysis",model,finaltime){
  load(paste0("results_sensitivity_analysis/",model,"-sensitivity.RData"))
  # Then, we read all the trajectories generated saving them in a list called
  # ListTraces. List that will be rewritten as a data frame in order to use ggplot.
  # ConfigID represents the initial condition associated to each trajectory,
  # which was generated by using the function implemented in the file Functions.R .
  
  listFile<-list.files(folder,
                       pattern = ".trace")
  
  configID<-t(sapply(1:length(listFile),
                     function(x){
                       return(c(x,config[[1]][[x]][[3]]))
                     }) )
  
  id.traces<-as.numeric(gsub("[^[:digit:].]", "",listFile) )
  
  if(rank){
    load(paste0(folder,"/ranking_",model,"-sensitivity.RData"))
    rank$idTrace = rank$id
    rank$id = as.numeric(gsub("[^[:digit:].]", "",rank$idTrace) )
  }else{
    rank <- data.frame(measure = 0 , id = id.traces)
  }
  
  # reference <- as.data.frame(read.csv("input/reference_data.csv",
  #                                     header = FALSE,
  #                                     sep = ""))
  # 
  reference <- as.data.frame(t(read.csv("Input/reference.csv", header = FALSE, sep = "")))
  name.place = as.character(factor(reference[,1]))
  value = as.numeric(as.character(factor(reference[,2])))
  names(value) = name.place
  value["P_lipid"] = NaN
  
  Ref.DF = data.frame(Place = name.place, V = value)
  
  #PlaceNotToPlot = c("Ppfia1PtprfNCm", "Ppfia1Ptprf","Ppfia1PtprfNCn","Ppfia1Ptprfn", "Sacm1l")
  PlaceNotToPlot = c("Sacm1l")
  
  ListTraces<-lapply(id.traces,
                     function(x){
                       xfile = paste0(model,"-sensitivity-",
                                      x,
                                      ".trace")
                       
                       trace.tmp=read.csv(paste0(
                         folder,"/",xfile), sep = "")
                       trace=data.frame(trace.tmp,ID = x, rank = rank[ which(rank$idTrace == xfile),1])
                       
                       if(max(trace$Time) == finaltime){
                         trace.final <-  lapply(colnames(trace)[-which( colnames(trace)%in% c("Time","ID","rank",PlaceNotToPlot))],function(c,PlaceNotToPlot){
                           
                           if(c == "Itga5Itgb1Fn1nMb")
                           {
                             df1 = data.frame(V= trace[,c] + trace[,"Itga5Itgb1Fn1nMbR"],ID = trace$ID, rank = trace$rank, Time=trace$Time,Place=c )
                             df2 = data.frame(V= trace[,c], ID = trace$ID, rank = trace$rank, Time=trace$Time,Place="Itga5Itgb1Fn1nMb_NotSum" )
                             df = rbind(df1,df2)
                           }else{
                             df = data.frame(V=trace[,c],ID = trace$ID, rank = trace$rank,Time=trace$Time,Place=c )
                           }
                           return( df )
                         },PlaceNotToPlot)
                         
                         trace.final <- do.call("rbind",trace.final)
                         return(trace.final)
                       }
                       
                     })
  
  traces <- do.call("rbind", ListTraces)
  
  ##### Modifica per plottare i reference con percentuali o somme in maniera corretta:
  v = traces[traces$Time == 0 & traces$Place =="P_lipid","V"] * (1-value["P_lipid"])
  Ref.DF[Ref.DF$Place == "P_lipid", "V"] = as.numeric(v)
  
  v = traces[traces$Time == 0 & traces$Place == "NCm","V"] * (value["NCm"])
  Ref.DF[Ref.DF$Place == "NCm", "V"] = unique(as.numeric(v))
  v = traces[traces$Time == 0 & traces$Place == "NCn","V"] * (1-value["NCm"])
  Ref.DF[Ref.DF$Place == "NCn", "V"] = unique(as.numeric(v))
  
  d <- unique(traces[traces$Place %in% Ref.DF$Place,"Place"])
  
  d = data.frame(Place = d, V = NaN )
  Ref.DF <- rbind(Ref.DF,d)
  Ref.DF$Time = max(traces$Time)
  
  pl<-ggplot( )+
    geom_line(data=traces,
              aes(x=Time,y=V,group=ID,col=rank))+
    geom_point(data=Ref.DF,
               aes(x=Time,y=V),
               col="red")+
    facet_wrap(~ Place,ncol =  5,scales = "free")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=14,face="bold"),
          legend.position="right",
          legend.key.size = unit(1.3, "cm"),
          legend.key.width = unit(1.3,"cm") )+
    labs(x="Seconds", y="V",col="Distance")
  
  ggsave(pl,filename = "Sensitivity.png",height = 10,width = 12)
  
  pl2<-ggplot( )+
    geom_histogram(data=traces[traces$Time == max(traces$Time),],
                   aes(V))+
    facet_wrap(~ Place,ncol =  5,scales = "free")+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=18),
          legend.title=element_text(size=20,face="bold"),
          legend.position="right",
          legend.key.size = unit(1.3, "cm"),
          legend.key.width = unit(1.3,"cm") )+
    labs(x="", y="V")
  
  
  ##### calcolo optim
  rank = rank[rank$id %in% unique(traces$ID),]
  rank[which.min(rank$measure),"id"] -> bestConfig
  TrNo = c("init","Endocytosis",  "Maturation_Adhesions" ,"DissociationParam","InteractionParam")
  optim =  rep(0,length(config)-length(TrNo))
  k=1
  for(x in 1:length(config) ){
    if(!config[[x]][[bestConfig]][[1]] %in% TrNo )
    {
      optim[k] = config[[x]][[bestConfig]][[3]]
      names(optim)[k] = config[[x]][[bestConfig]][[1]]
      k=k+1
    }else if(config[[x]][[bestConfig]][[1]] == "init"){
      init = config[[x]][[bestConfig]][[3]][c("PtdIns4Pa","Prkd1X","lipid")]
    }
  }
  optim = c(optim,init)
  
  return(list(pl1=pl,pl2=pl2,optim=optim))
}


WifPlot <-function(folder,model,paramWIF,filename,optim){
  load(paste0(folder,model,"-analysis.RData"))
  PlaceNotToPlot = c("Fn1")
  
  #reference <- as.data.frame(t(read.csv("./Input/reference.csv", header = FALSE, sep = "")))
  
  # New reference with multiple points!!!!
  # library(readxl)
  # Valori_finali <- read_excel("./Input/Valori finali.xlsx", 
  #                             col_names = FALSE)
  # reference <- as.data.frame(Valori_finali[,2:3])
  # 
  # name.place = as.character(factor(reference[,1]))
  # value = as.numeric(as.character(factor(reference[,2])))
  # names(value) = name.place
  # Ref.DF = data.frame(Place = name.place, V = value )
  # 
  listFile<-list.files(folder,
                       pattern = ".trace")
  
  configID<-t(sapply(1:length(listFile),
                     function(x){
                       return(c(x,config[[1]][[x]][[3]]))
                     }) )
  paramNames = sapply(1:length(config), function(x) config[[x]][[1]][[1]])
  paramNames = sub(pattern = "\t", "",paramNames)
  Indexes = which(paramNames %in% paramWIF)
  
  id.traces<-as.numeric(gsub("[^[:digit:].]", "",listFile) )
  
  ListTraces<-lapply(id.traces,
                     function(x){
                       xfile = paste0(model,"-analysis-",
                                      x,
                                      ".trace")
                       
                       trace.tmp=read.csv(paste0(
                         folder,"/",xfile), sep = "")
                       #perc = (trace.tmp$Ppfia1n[1]*2)/66
                       #perc = config[[35]][[x]][[3]]/0.08
                       
                       percSpeed = config[[ Indexes[1] ]][[ x ]][[3]]/optim[paramNames[Indexes[1]]]
                       percSlow = config[[ Indexes[2] ]][[ x ]][[3]]/optim[paramNames[Indexes[2]]]
                       trace=data.frame(trace.tmp,ID = x, perc = abs(percSpeed - percSlow) )
                       ###
                       trace[,"LipidkxTot"]  = trace[,"lipidkX"] + trace[,"Lipidkx"]
                       trace[,"PKDxTot"]  = trace[,"PrkdX"] + trace[,"PKDx"] + trace[,"MAPKA_PKD"]
                       trace[,"PtdInsTot"]  = trace[,"P_lipid"] + trace[,"lipid"]
                       trace[,"Sacm1lTot"]  = trace[,"Sacm1lX"] + trace[,"Sacm1lx"]
                       trace[,"Itga5Itgb1Fn1nMTot"]  = trace[,"Itga5Itgb1Fn1nMb"] + trace[,"Itga5Itgb1Fn1nMbR"]
                       trace[,"NCTot"]  = trace[,"NCn"] + trace[,"Ppfia1_NCn"]
                       trace[,"NCmTot"]  = trace[,"NCm"] + trace[,"Ppfia1_NCm"]
                       ###
                         trace.final <-  lapply(colnames(trace)[-which( colnames(trace)%in% c("Time","ID","percSpeed","percSlow","perc",PlaceNotToPlot))],function(c,PlaceNotToPlot){
                          
                        df = data.frame(V=trace[,c],ID = trace$ID, Perc = trace$perc,Time=trace$Time,Place=c )

                           return( df )
                         },PlaceNotToPlot)
                         
                         trace.final <- do.call("rbind",trace.final)
                         return(trace.final)
                       
                       
                     })
  
  traces <- do.call("rbind", ListTraces)
  
  lev =  c(
    "Ppfia1m","Ppfia1n","Ppfia1_NCm","Ppfia1_NCn","NCm","NCn","MapkA","MAPKB",
    "PKDx","MAPKA_PKD","PrkdX","Lipidkx", "lipidkX","lipid","P_lipid","PtdinsTrash",
    "Sacm1lx","Sacm1lX",
    "VesicleFnC",
    "Vesicle_FN_Itga5itgb1",
    "Itga5Itgb1Fn1nMb","Itga5Itgb1Fn1nMbR","Itga5Itgb1Fn1Rab21C","Lysosomes",
    "Itga5Itgb1G","Itga5Itgb1Fn1nMa","Vesicle_Itga5Itgb1_Rab11bC",
    "PKDxTot","LipidkxTot","PtdInsTot","Sacm1lTot","Itga5Itgb1Fn1nMTot",
    "NCTot","NCmTot"
  )
  lab = c("PPFIA1m","PPFIA1n","PPFIA1_NCm","PPFIA1_NCn","NCm","NCn","MapkA",
          "MAPKB","PKD","Inactive_PKD","PKD*","PI4KB","PI4KB*","PI" ,"PI4P","PI_Buffer",
          "SACM1L","SACM1L*",
          "Vesicle_FN", "Vesicle_Integrin_FN","Integrin_FN_Mb",
          "Integrin_M_R","Vesicle_Integrin_FN_Rab21","Lysosomes",
          "Endocytic_Recycling_Compartment","Integrin_Ma","Vesicle_Integrin_Rab11b",
          "Sum_PKD", "Sum_PI4KB","Sum_PI","Sum_SACM1L", "Sum_Integrins_Mb", "Sum_NCn", "Sum_NCm")
  
  traces$Area = NA
  traces$Area[traces$Place %in% c("MapkA","MAPKA_PKD","MAPKB","MAPKA_PKD","LipidkxTot", "PKDxTot","LipidkxTot","PtdInsTot","Sacm1lTot",
                                  "PKDx","PrkdX","Lipidkx", "lipidkX","lipid","P_lipid","PtdinsTrash","Sacm1lx","Sacm1lX","Fn1")] = "Golgi Apparatus"
  traces$Area[traces$Place %in% c("NCn","Ppfia1n","Ppfia1_NCn","NCTot")] = "Nucleus"
  traces$Area[traces$Place %in% c("VesicleFnC","Vesicle_Itga5Itgb1_Rab11bC",
                                  "Vesicle_FN_Itga5itgb1",
                                  "Itga5Itgb1Fn1Rab21C",
                                  "Itga5Itgb1G","Lysosomes")] = "Integrin recycling loop"
  traces$Area[traces$Place %in% c("Itga5Itgb1Fn1nMb","Itga5Itgb1Fn1nMbR","Itga5Itgb1Fn1nMTot","NCmTot","NCm",
                                  "Ppfia1m","Ppfia1_NCm")] = "Basal membrane"
  traces$Area[traces$Place %in% c("Itga5Itgb1Fn1nMa")] = "Apical membrane"
  
  traces$V[traces$Place == "Sacm1lTot"]  = 33.75
  
  # Ref.DF = Ref.DF[!Ref.DF$Place %in% PlaceNotToPlot,]
  # v = traces[traces$Time == 0 & traces$Place =="P_lipid","V"] * (Ref.DF[Ref.DF$Place == "P_lipid","V"])
  # Ref.DF[Ref.DF$Place == "P_lipid", "V"] <- as.numeric(v)
  # v = traces[traces$Time == 0 & traces$Place =="PrkdX","V"] * (Ref.DF[Ref.DF$Place == "PrkdX","V"])
  # Ref.DF[Ref.DF$Place == "PrkdX", "V"] <- as.numeric(v)
  # 
  traces$Place =  factor(traces$Place, levels = lev, labels =lab  )
  #Ref.DF$Place =  factor(Ref.DF$Place, levels = lev, labels =lab  )
  #Ref.DF$Time = max(traces$Time)
  
  library(dplyr)
  pl<-ggplot( )+
    facet_wrap(~ Place,ncol =  5,scales = "free")+
    geom_rect(data=traces%>%
                select(Place,Area) %>% distinct() ,
              aes(fill = Area),
              xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.3) +
    geom_line(data=traces,
               aes(x=Time/(60*60),y=V, col = Perc, group = Perc))+
    # geom_point(data=Ref.DF,
    #            aes(x=Time/(60*60),y=V),
    #            col="red")+
    geom_vline(xintercept = c(24),aes(col="RemPPFIA"),linetype = "dashed")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=14,face="bold"),
          legend.position="right",
          legend.key.size = unit(1.3, "cm"),
          legend.key.width = unit(1.3,"cm") )+
    labs(x="Hours", y="", col = "Costants\ difference")+
    theme_bw()+
    scale_fill_manual(values = c("Golgi Apparatus" = "orange",
                                 "Apical membrane"="#FF0000",
                                 "Integrin recycling loop" = "#008080",
                                 "Basal membrane" = "#FF0000",
                                 "Nucleus"="blue"),
                      guide = guide_legend(title = "",
                                           override.aes = list(alpha = .5)))
  
  ggsave(pl,filename = filename,height = 15,width = 20)
  
  return(pl)
}
