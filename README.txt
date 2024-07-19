kpop: A kernel balancing approach for reducing specification assumption in survey weighting

Erin Hartman: ekhartmant@berkeley.edu
Chad Hazlett: chazlett@ucla.edu
Ciara Sterbenz: cster@ucla.edu

Replication for results in 2016 Application (Section 6) and Simulation (Section 5.2)

#######
data_prep_rep.R

	This file loads all necessary data from a folder named "data" (requires downloading Pew survey and election returns manually through dataverse, instructions commented within), cleans data (dropping units that done report vote preference etc), predicting and filling in missing age and education for units in pew, recodes/bins all auxiliary variables, and creates relevant interactions. It then builds and runs a 3-way lasso multinomial logit in CCES (using 80\% training) and projects this model back onto the Pew and CCES data to create our modeled outcome. All of this is described in section 6.1. It produces a .Rout file "data_prep_out.Rout" which describes the data cleaning process step by step. It produces two data files, "pew_prepared.rds" and "cces_prepared.rds" in a folder named" generated_data." The exact data cleaned and used by the authors may not be directly reproducible due to seed selection, but are additionally supplied as "pew_prepared_orig.rds" and "cces_prepared_orig.rds" in the "generated_data" folder.
	
######
generated_data (folder) 

	folder where cleaned data are saved and loaded from in the remainder of the files
	
######
data (folder)

	folder where the raw data should be saved. Note that while Pew and election returns are publicly available, they require a simple free login bypass to download. Users are thus required to download Pew and election returns manually to this folder.

######	
app_kpop_run.R

	This file runs kpop for all specifications and saves the results to a weights file "app_update_[CURRENTDATE].Rdata" to the "weights" folder that is then loaded below in the results .Rmd that produces all figures and plots in the paper and appendix related to the application in Section 6. Note that given the large size of the CCES data (40k row matrix), this file is time consuming to run as taking the SVD and balancing on columns of a large matrix is computationally intensive. This file does not analyze or produce results, but simply saves the kpop weights to a file which is then loaded in the results file "app_res_rep.Rmd" 

######	
weights (folder)

	folder where the kpop weights produced by "app_kpop_run.R" and "sims_rep.R" are saved. Original kpop weights produced by authors for the application in Section 6 are in the file "app_update_june2023-06-28.Rdata." Weights across a range of choices of b are in the file "app_allBPOPW_2023-07-01.Rdata." To save computational time, there is an additional file "Kscree_svdsA_reduced.Rdata" which saves the results of the SVD of several different kernel matrices for the scree plot in Appendix B.1 which is loaded in "app_kpop_run.R." All three of these files are loaded by "app_res_rep.Rmd". Finally, this folder contains the simulation results file in "sims_kpop_2023-04-18_nsims1000.RData" loaded by "sims_res_rep.Rmd" 


######
app_res_rep.Rmd

	This file loads both the prepared data from the "generated_data" folder and the necessary weights from the "weights" folder to produce all the results in Section 6 and the Appendix, including 1) the kpop weights file, 2) an additional kpop weights file for various choices of b, and 3) an additional file containing various SVD results for different kernel matrices to save some computation time in producing the scree plot in Appendix B.1. This Markdown file produces all tables and figures in the paper and appendix. 
	It also can produces more in-depth results/adjustments for interested replicators. Users may also adjust the function call to request: a different set of weighting estimators (additional available methods include kpop with forced ebalance convergence and post stratification on all variables), a different set of auxiliary variable margins, the mean absolute error on margins not weighted by the CCES population targets (which themselves are weighted by CCES's internal weights), or change the population target from the 3-way multinomial logit modeled vote difference in the weighted CCES to the actual self-reported vote difference in the weighted CCES.
	Note that this file is written for a user to run inline in RStudio including walkthrough comments.  As such, the knitted pdf may have aesthetic deficiencies (e.g. tables not quite fitting). The authors recommend an interested replicator explore the results of the application inline.

######
sims_rep.R

	This script will run the simulations presented in Section 5.2. NB: This is computationally intensive and take very long to run. It loads CCES data from "cces_prepared_orig.rds" in the "generated_data" folder. The file will printout updates as the simulations progress. Given the long runtime, we recommend running the file from the command-line with a "nohup R CMD BATCH sims_rep.R & " command and checking the .Rout file with cat to check progress from the printouts. Results of the simulations are saved to a file "sims_kpop_[CURRENTDATE]_nsims1000.RData"
	
	
######
sims_res_rep.Rmd

	All simulation results in section 5.2 and all of Appendix D are produced in this markdown file. It loads the prepared CCES data from the "generated_data" folder and the simulation weights file output from the script "sims_rep.R" from the "weights" folder.  The original simulation results used by the authors are in file "sims_kpop_2023-04-18_nsims1000.RData" with the original prepared CCES data in "cces_prepared_orig.rds."
