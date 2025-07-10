
#remotes::install_github("https://github.com/qBioTurin/epimod_FBAfunctions", ref="main",dependencies = T,force = T)

model = epimodFBAfunctions::FBA4Greatmod.generation(fba_mat = "MetabolicModel/Recon3D_301.mat")

epimodFBAfunctions::showFBAmethods()

epimodFBAfunctions::writeFBAfile(theObject = model,fba_fname = "Recon3D",dest_dir = "MetabolicModel/"  )

ExR <- getExchangesR(model)


library(readxl)
Genes4FBA <- read_excel("Genes4FBA.xlsx")


grep(pattern = "5587.1" ,x = genesFromGeneAss)


model@react_id[grep(pattern = "PIK4" ,x = model@react_id)]
