# From Data to Models:  Analyzing and Integrating Biological Data into Mechanistic Models





## Golgi_model
   ![PN_Retrieval](https://user-images.githubusercontent.com/81301099/161936773-a740437c-599a-4607-86ec-6d5b03c38d5a.png)
    
   Fig. 1) PN representation of PPFIA1-dependent apico-basolateral polarity pathway.
  
  <p> The model consist in 27 places and 37 transitions.
  
  Our simplified representation of a single cell encompasses five cellular compartments, as described in Fig. 1:
  1. basal membrane; 
  2. cytosol;
  3. nucleus;
  4. Golgi Apparatus; 
  5. apical membrane.
  
These compartments partake in three functional and interconnected phenomena, highlighted in the PN with colored boxes (Fig. 1):
  1. the integrin α5β1 recycling loop, coupled with FN secretion, which spans the entirety of the cell, from the basal membrane through the cytoplasm up to the apical membrane (see Fig. 1, teal colored box);
  2. ATAC-dependent modulation of MAPK13 and MAPK14 gene expression (see Fig. 1, blue box);
  3. regulation of Golgi secretion through the control of PI4P levels on TGN membranes (see Fig. 1, yellow box). 
                                                                                                                     
  
<p> The model simulates the PPFIA1-dependent pathways occurring in a single representative EC as PPFIA1 is artificially removed from the cell to transition from a physiological condition to one where polarity is absent. <p\>
 
  
<p> Model construction and analysis was performed using the GreatMod framework, which provides an intuitive graphical interface for model construction based on PNs and extensions, and automatic derivation of the low-level mathematical processes (deterministic or stochastic) characterizing the system dynamics. <p\>

The GolgiMod_model folder contains the main R script "MainAnalysis.R", where it is possible to find the commands exploited during the analysis, and several folders.
  
Folders:

1.  **Net** contains the PNPRO file corresponding to the PN exploited to model the PPFIA1-dependent apico-basolateral polarity pathway, the C++ code regarding the general transitions, and the solver (which should be generated before starting the analysis);
2.  **R\_functions** stores all the R scripts to generate the plots and to run the analysis;
3.  **Input** contains all the csv files (e.g., reference, parameter lists) necessary to run the analysis;
4.  **Sensitivity** stores the RDatas storing the PRCC values for a subset of places and their respective plot (saved as png).
  
<p> A detailed description of the model and results is reported in Dr. Dora Tortarolo's PhD Thesis. <p\>
    
