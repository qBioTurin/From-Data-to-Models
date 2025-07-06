library(readr)
name.generation = function(path,out)
{
  Names <- read_table2(path)
  
  transStart<- which(Names[,1] == "#TRANSITION")
  
  NAMES = unlist(Names[1:(transStart-1),1])
  NAMES = unname(NAMES)
  
  saveRDS(NAMES,file=paste0(out,"NAMES.RDS"))
  ##################################################
  # 
  # library(readxl)
  # Data <- read_excel("input/Data.xlsx", na = "NA")
  # saveRDS(Data,file="./input/Data.RDS")
  # ############################
  # 
  # reference <- Data[,c("Place name","final perc")]
  # NA.pos <- which(is.na(reference[,2]))
  # reference<-as.data.frame(reference[-NA.pos,])
  # 
  # reference$`Place name`<- which(NAMES%in%reference$`Place name` )
  # write.table(t(reference),"./input/reference.csv",row.names = F, col.names= F, sep =" ")
  
}

name.generation(path = "./PPFIA1Model.PlaceTransition",out = "Input/")
# reference.generation = function(path="input/Data.csv")
# {
#   Data <- read_excel(path, na = "NA")
#   saveRDS(Data,file="./input/Data.RDS")
#   ############################
# 
#   reference <- Data[,c("Place name","final perc")]
#   NA.pos <- which(is.na(reference[,2]))
#   reference<-as.data.frame(reference[-NA.pos,])
# 
#   reference$`Place name`<- which(NAMES%in%reference$`Place name` )
#   write.table(t(reference),"./input/reference.csv",row.names = F, col.names= F, sep =" ")
# 
# }