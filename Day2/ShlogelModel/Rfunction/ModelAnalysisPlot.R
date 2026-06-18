library(dplyr)
library(ggplot2)

ModelAnalysisPlot=function(solverTraces_path, Stoch = F, print=T){
  
  traces = list.files(path = solverTraces_path, pattern = ".trace")
  i = 0
  list.plX1 <- list()
  list.plAll <- list()
  list.plAllMean <- list()
  list.plX1dens <- list()
  
  for(tr in traces){
    
    i = i + 1
    solverName_path = paste(solverTraces_path, "/", tr, sep = "")
    
    trace <-read.csv(solverName_path,sep = "")
    n_sim_tot<-table(trace$Time)
    n_sim <- n_sim_tot[1]
    time_delete<-as.numeric(names(n_sim_tot[n_sim_tot!=n_sim_tot[1]]))
    
    if(length(time_delete)!=0) trace = trace[which(trace$Time!=time_delete),]
    
    trace$ID <- rep(1:n_sim[1],each = length(unique(trace$Time)) )
    
    trace.final <-  lapply(colnames(trace)[-which( colnames(trace)%in% c("ID","Time"))],function(c){
      return(data.frame(V=trace[,c], ID = trace$ID,Time=trace$Time,Compartment=c ) )
    })
    trace.final <- do.call("rbind",trace.final)
    
    name <- paste("plX1", i, sep=".")
    
    plX1<-ggplot( )+
      geom_line(data=trace,
                aes(x=Time,y=X1,group=ID), 
                color = "red", size = 1)+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=20,face="bold"),
            legend.text=element_text(size=18),
            legend.title=element_text(size=20,face="bold"),
            legend.position="right",
            legend.key.size = unit(1.3, "cm"),
            legend.key.width = unit(1.3,"cm") )+
      labs(title = paste("p_ini: ", trace.final$V[1]), x="Time", y="X1")
    
    list.plX1[[name]] = plX1
    
    if(Stoch){
      
      meanTrace <- trace %>% group_by(Time) %>%
        summarise(X1=mean(X1))
      
      meanTrace.final <- lapply(colnames(meanTrace)[-which(colnames(meanTrace)=="Time")],function(c){
        return(data.frame(V=unlist(meanTrace[,c]), Time=meanTrace$Time,Compartment=c ) )
      })
      meanTrace.final <- do.call("rbind",meanTrace.final)
      
      name <- paste("plAll",i,sep=".")
      
      plAll <-ggplot( )+
        geom_line(data=trace.final,
                  aes(x=Time,y=V,group=ID))+
        geom_line(data=meanTrace.final,
                  aes(x=Time,y=V,col="Mean"),
                  linetype="dashed")+
        facet_grid(~Compartment)+
        theme(axis.text=element_text(size=15),
              axis.title=element_text(size=15,face="bold"),
              legend.text=element_text(size=18),
              legend.title=element_text(size=18,face="bold"),
              legend.position="bottom",
              legend.key.size = unit(1, "cm"),
              legend.key.width = unit(1,"cm") )+
        labs(x="Time", y="Population")
      
      list.plAll[[name]] = plAll

      name <- paste("plAllMean",i,sep=".")
      
      plAllMean <-ggplot( )+
        geom_line(data=meanTrace.final,
                  aes(x=Time,y=V,col=Compartment),
                  linetype="dashed")+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20,face="bold"),
              legend.text=element_text(size=18),
              legend.title=element_text(size=20,face="bold"),
              legend.position="right",
              legend.key.size = unit(1.3, "cm"),
              legend.key.width = unit(1.3,"cm") )+
        labs(x="Time", y="Mean Population")
      
      list.plAllMean[[name]] = plAllMean

      name <- paste("plX1dens",i,sep=".")
      
      plX1dens<-ggplot(trace[trace$Time==max(trace$Time),])+
        geom_histogram(aes(X1))+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20,face="bold"),
              legend.text=element_text(size=18),
              legend.title=element_text(size=20,face="bold"),
              legend.position="right",
              legend.key.size = unit(1.3, "cm"),
              legend.key.width = unit(1.3,"cm") )
      
      list.plX1dens[[name]] = plX1dens
      
      list.plX1[[i]] <- list.plX1[[i]]+
        geom_line(data=meanTrace,
                  aes(x=Time,y=X1,col="Mean"),
                  linetype="dashed")+
        labs(x="Time", y="X1",col="")
      
      ListReturn<-list(list.plX1 = list.plX1, list.plX1dens = list.plX1dens, 
                       list.plAll = list.plAll,list.plAllMean = list.plAllMean)
      
    }
    else{
      
      list.plAll <- list()
      name <- paste("plAll",i,sep=".")
      
      plAll <-ggplot( )+
        geom_line(data=trace.final,
                  aes(x=Time,y=V,col=Compartment))+
        theme(axis.text=element_text(size=18),
              axis.title=element_text(size=20,face="bold"),
              legend.text=element_text(size=18),
              legend.title=element_text(size=20,face="bold"),
              legend.position="right",
              legend.key.size = unit(1.3, "cm"),
              legend.key.width = unit(1.3,"cm") )+
        labs(title = paste("p_ini: ", trace.final$V[1]), x="Time", y="X1")
      
      list.plAll[[name]] = plAll
      
      ListReturn<-list(list.plX1 = list.plX1,  list.plAll = list.plAll)
    }
  }
  
  if(print){
    for(j in 1:length(ListReturn))
      print(ListReturn[j])
  }
  
  return(ListReturn)
}
