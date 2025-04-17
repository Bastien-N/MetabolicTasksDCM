# MetabolicTasksDCM
Scripts and data related to the manuscript titled "Metabolic task analysis reveals distinct metabotypes in end-stage dilated cardiomyopathy"

The manuscript is available as a preprint on Biorxiv (https://doi.org/10.1101/2025.03.28.645901)

## Run order
### Subscripts
R scripts and function in the R/Subscripts folder will be sourced as required.
Matlab scripts in the Matlab/Subscripts folder *must be added to the Matlab Path* before running any of the other Matlab scripts.
### Main scripts
The first scripts which should be run are:
- 1st the R scripts numbered "0", to preprocess the expression data.
- 2nd the Matlab script "N0_task_analysis_essential_rxns.m", to generate essential reactions required for CellFie

Then, the other Metlab scripts can be run.
Finally, the remaining R scripts can be run in their numbered order.