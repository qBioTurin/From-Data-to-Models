#### Check the transitions!!
library(readr)
transition <- read_table2("./Net/PPFIA1Model.PlaceTransition",
                             skip = 28)
transList <- read_delim("./Input/paramsList.csv",
                               delim = ";", escape_double = FALSE, trim_ws = TRUE)

transition$`#TRANSITION` [! transition$`#TRANSITION` %in% transList$init]
transList$init[! transList$init %in% transition$`#TRANSITION`]
## Mettere stessa velocità di movimento Nucleo -> membrana
##########