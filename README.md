# Scripts for CAGI 5 ID challenge assessment
 

## REQUIREMENTS
```
R:
>install.packages('ROCR')
>install.packages('plotrix')
>install.packages('gridExtra)
>install.packages('ggplot2')
```
## USAGE

to run all analyses, you have just to run the following command: Rscript ./src/CAGI_assessment_main.R
      This script will compute all statistics and make plots
      This script expects a folder structure like this to run:
      ./src <- contains all scripts
      ./results <- will contain all performance tables and plots
      ./data <- has to contain 3 folders: experimental_value, submissions, template 
      
      These are Input needed:
            * experimental values file  in ./data/experimental_value
            * submission template in ./data/template
            * submission files in ./data/submissions

      Running the script these Output files will be generated in ./results
           * AUC of each submission
           * barplot with amount of patients correctly predicted by n groups for each phenotype
           * heatmap of AUC of all submissions on each phenotype
           * barplot with number of correctly predicted variants by each submission
           * barplot with amount of correctly predicted variants by n groups
           * pie chart with number of patiets with and without variants
           * barplot with number of Causative, Putative causative and Contributing factor variants
           * barplot with number of patients with known features for each phenotype 
           * table with assessment performed on each patient
           * tables with n of correct predictions, n of correct variants, n of correct prediction and variant (all patients or only patients with variant)
           * tables with number of correct predictions for each patient on each class
           * tables with submission performance on each phenotype (AUC, MCC, ACC, F1, TPR, PPV, TNR, NPV, FNR, TP, TN, FN, FP)
           * barplots with submission performance(TP, TN, FP, FN) computed on submitted p-values (P = "*" are not considered) and experimental value, computed with all or only patients with variants
           * table with submission performance on all phenotypes (AUC, MCC, ACC, F1, TPR, PPV, TNR, NPV, FNR, TP, TN, FN, FP)
           * table with (for each submission): total number of correctly predicted variants, total predicted (different) variants, Correctly predicted variants/Experimental variants, Correctly predicted variants/Predicted Variants 
           * table with rank of AUC on all submissions and phenotypes
           * table with patient, variant, number of groups with correct answer
           
