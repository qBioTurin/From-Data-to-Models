install.packages("devtools")
devtools::install_github("qBioTurin/ORCA", ref="main", dependencies=TRUE)

ORCA::ORCA.run()
