Metabolic tasks scripts:
The scripts function as their equivalents from the CobraToolbox. 
The only added functionnalities are:
- fixing an important issue with getting essential tasks.
- the 'saveUsed' parameter, which allows to return all the reactions carying 
flux to perform a task during the essential reaction check. 
- the ability to use the EQU field from the RAVEN toolbox when designing tasks.
- the ability to allow any input/output by using "ALLMETS"

Metabolic task lists:
- Latest version is MetabolicTasks_2022_3.xlsx
- MetabolicTasks_2022_RAVEN_3.xlsx can be used with RAVEN's checkTasks function (may require closing off exchange reactions prior). It is also usefull to visually inspect tasks as it uses metNames instead of metID.


Cellfie_2022:
Computes the Cellfie scores for metabolic tasks from gene expression data, based on the essential reactions needed for the task in a reference model.

cobra_to_Raven_tasks.R: R script which can be used to convert our COBRA tasklists to a RAVEN tasklist. Could be modified to allow conversion other COBRA tasklists