library(dplyr)
library(ggplot2)
library(tidyr)

trLev = c("Interaction_PPFIA1ATAC","Dissociation_PPFIA1ATAC",
       "PPFIA1_MN",
       "PPFIA1_NM","ATAC_MN","ATAC_NM","PPFIA1_ATAC_MN",
       "PPFIA1_ATAC_NM","TranscriptionalReg_MAPK13", "TranscriptionalReg_MAPK14",
       "Consumption_MAPK13","Consumption_MAPK14","Sequestration",
       "Input_PKD","Activation_PKD","Inactivation_PKD","Consumption_PKD",
       "Activation_PI4KB","Inactivation_PI4KB*","Input_PI",
       "Phosphorylation","Dephosphorylation","Activation_SACM1L","Inactivation_SACM1L*",
       "FN_Vesicle_Budding","Exocytosis","Docking","Endocytosis","Endocytosis_R",
       "Endocytosis_A","Movement_MG","Recycling_A","Recycling_B","Movement_GM","Retrieval")

trLab = c("InterGeneration","DissGeneration",
          "PPFIA1_MN","PPFIA1_NM","ATAC_MN","ATAC_NM","PPFIA1_ATAC_MN", "PPFIA1_ATAC_NM",
          "CostMAPK13","CostMAPK14" , "ConsuptionMAPK13","ConsuptionMAPK14","Sequestration",
          "Input_PKDx","Activation_PRKD1x","Inactivation_PRDK1X", "Consumption_PKD",
          "Activation_PI4KBx", "Inactivation_PI4KBX",
          "Input_Ptdins","FromPtdIns_toPtdIns4P" ,"FromPtdIns4P_toPtdIns",
          "Activation_SACM1Lx","Inactivation_SACM1LX", 
          "VesicleBudding","Exocytosis","Docking","Endocytosis","EndocytosisR",
          "Endocytosis_A","Movement_MG","RecyclingA",           
          "RecyclingB","Movement_GM","ExocitosiRetro")    
          

names(trLev) = trLab

listRData<-list.files("Sensitivity",
                     pattern = ".RData")


for(i in listRData){
  load(paste0("Sensitivity/",i))
  place = sub(".RData",replacement = "",x = sub("prcc_",replacement = "",x = i))
  prcc$PRCC$Time = seq(0,72*60*60,60*60)
  dataPRCC = prcc$PRCC %>% 
    gather(key = "Transitions", value = "V",-Time)
  dataPRCC$Transitions = sub("\t","",dataPRCC$Transitions)
  dataPRCC$Transitions = trLev[dataPRCC$Transitions]
  
  transName = dataPRCC %>%
    filter(V > 0.21 | V < -0.21) %>%
    select(Transitions) %>%
    distinct()

  dataPRCC = dataPRCC[dataPRCC$Transitions %in% transName$Transitions,]
  pl = ggplot(dataPRCC) +
    geom_rect(data= data.frame(Time=c(0,60*60),V=c(0.2,-0.2)),
              fill = "yellow", alpha = .5,
              xmin = -Inf, xmax = Inf,
              ymin = -0.2, ymax = 0.2)+
    geom_line(
              aes(x=Time/(60*60),y=V,col=Transitions,group=Transitions))+
    theme_bw()+
    labs(x = "Hours",y = "PRCC",title = paste0("Place: ",place))+
    # geom_text(data = dataPRCC[which(dataPRCC$Time == max(dataPRCC$Time) ),],
    #           aes( x = Time/(60*60), y =V,label = Transitions)
    #           )+
    geom_vline(xintercept = c(24,48),col="blue",linetype = "dashed")+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=14,face="bold"))
  pl
  ggsave(pl,filename = paste0("Sensitivity/PRCC",place,".png"),height = 8,width = 18) 
  
  
}
