
## Mammographic density / BMI / breast cancer MR project


### Background 


 to add


### This repository

Main analysis scripts and metadata (see details below):

```
├── set_paths.R
├── 00v1_MD_data_processing.Rmd
├── 00v2_MD_data_processing.Rmd
├── 01_process_gwas_summary.Rmd
├── 02_mr_BMI_to_MD.Rmd
├── 03_mr_mediators_to_BC.Rmd
├── 04_mvmr_run_analysis.Rmd
├── 05_mvmr_create_plots.Rmd * not added yet
├── 06_mediation_analysis.Rmd  
├── functions.R
├── functions_mvmr.R
├── metadata
│   └── data_lookup.csv
│   └── data_lookup_BCAC.csv

```


### How to reproduce this analysis, i.e. Main Workflow


1.  Script `set_paths.R` is imported by all other scripts, and is used to set the environment of where the project is run.

2. MD data processing 

`00v1_MD_data_processing.Rmd` - original MD GWAS files processing (adjusted for BMI)

`00v1_MD_data_processing.Rmd` - re-run MD GWAS files processing (unadjusted)


2. Processing all datasets into 'tidy' format

`01_process_gwas_summary.Rmd` is required for processing data that comes as text files (i.e. GWAS summary stats) from the IEU GWAS pipeline or other sources. This script has to be run to convert raw data into `outcome` data frames and to extract instruments from each GWAS (in `exposure` format) and save them to be used directly in MR analysis in subsequent scripts. The names of raw files, tidy outcome data frames, and exposure instruments are all get saved in the metadata file `data_lookup.csv` upon generation. (NB the metadata file has to contain raw file names and the desired output prefixes before running this Rmd).

3. Rmd `02_mr_BMI_to_MD.Rmd` runs univariable MR of Childhood and Adult body size on MD (a version of MD data to use is specified via metadata file `data_lookup.csv`). The code has to be run interactively per trait category. The results merged by trait category will be stored in `Results` directory outside the codebase. 


4. Rmd `03_mr_mediators-to-BC.Rmd` is used to run univariable MR of all mediators (and BMI) on Breast cancer (`ieu-a-1126`). The code has to be run interactively per trait category, the results are stored in `Results` outside the codebase. After the analysis, forest plots can be created for each trait category. To recreate the plots, don't need to rerun the full analysis, can just read in the merged files. The plots will be saved in the codebase in `figures/`. 

5. Rmd `04_mvmr_run_analysis.Rmd` runs four types of MVMR with each mediator, but first we run MVMR with BMIs only (Analysis 0).  

	*	 (Analysis 0) multivariable MR: childhood BMI + adult BMI -> Breast Cancer 
	*	 (Analysis 1) multivariable MR: childhood BMI + adult BMI -> mediator
	*	 (Analysis 2) multivariable MR: childhood BMI + mediator -> Breast Cancer 
	*	 (Analysis 3) multivariable MR: adult BMI + mediator -> Breast Cancer 
	*	 (Analysis 4) multivariable MR: childhood BMI + adult BMI + mediator -> Breast Cancer 



05_create_plot.Rmd

- S1 figure


	The code has to be run interactively per trait category, the results are stored in `Results/trair_category/` outside the codebase. The analysis is structured as a large for-loop that will perform all four MVMR for each mediator in the selected trait category when individually specified T/F outside the loop. After all analyses have been performed, the mediators within each trait category are collated into a single dataframe and saved in `Results/trait_category/merged/` directory.

  MVMR analysis is performed using modified code from `TwoSampleMR` package and also `MVMR` package for comparison. MVMR sensitivity tests are done using `MVMR` package and phenotypic correlations estimated using `metaCCA` package. 

6. Rmd `05_mvmr_create_plots.Rmd` creates forest plots from all MVMR analyses for each trait category separately and by trait category. The plots are saved to `figures` within the codebase. The code also creates summary plot of all direct Childhood BMI estimates from Analyses 2 & 3.


7. Rmd `06_mediation_analysis.Rmd` contains a workflow for performing MR mediation analysis using Difference and Product method, with CIs calculation using Delta and Propagation of Errors methods. The script also contains methods for plotting the calculated indirect estimates (+CIs) as forest plots. 














