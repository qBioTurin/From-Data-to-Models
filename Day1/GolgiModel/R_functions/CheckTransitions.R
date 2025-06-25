#### Check the transitions!!
library(readr)
transition <- read_table2("./Net/GolgiModel.PlaceTransition",
                             skip = 25)
transList <- read_delim("./Input/Sensitivity_list_Prova0312.csv",
                               delim = ";", escape_double = FALSE, trim_ws = TRUE)

transition$`#TRANSITION` [! transition$`#TRANSITION` %in% transList$init]
transList$init[! transList$init %in% transition$`#TRANSITION`]
## Mettere stessa velocità di movimento Nucleo -> membrana
##########