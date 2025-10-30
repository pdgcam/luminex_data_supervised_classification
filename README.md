Data structure: 
Univariate: Singe variable used for predict 
Multivariate: 2 or more variables used for prediction 

Classification targets: 
Flavivirus versus not (binomial)
Any DENV vs no DENV (binomial)
ZIKV vs no ZIKV (binomial) 
Any DENV vs ZIKV given flavivirus infection (binomial)
CHIKV vs DENV
Any DENV vs ZIKV vs no infection (multinomial) (kept as placeholder for now) 
DENV serotype given infection (multinomial)
DENV serotype vs neg 

Analysis: 
Longitudinal: Ratio of pre- and post-infection timepoints per sample)
Cross-sectional:  Last infection per sample 

Models: 
Binomial Multivariate : GLM, Random Forest, Support Vector Machine 
Binomial Univariate : GLM, Decision Tree, Support Vector Machine  
Multinomial Multivariate: Random Forest, Naive Bayes

Metrics: 
AUROC
AUPRC 
Brier
Stratified Brier (for multinomial) 



Further detaials about scripts: 


1) PreProcessing.R : 
Takes dataset with logged MFI values, Target infection and metadata (patient ID, days since infection) 
Use HI titres to calculate log mean HI ratio between each blood draw, if ratio < 1.6, classify sample as 'no infection'
Use only samples classified as 'no infection' for model fitting, else remove samples
Returns all_processed_dfs, which contains ratio_data (used for longitudinal analysis) and cross_sectional_data (used for cross-sectional analysis) 
Ratio data = Post infection / Pre infection ratio for each individual patient 
Cross_sectional_data = MFI titre value at last blood draw for each individual patient 

Note: If data already in correct format (for longitudinal / cross-sectional analysis, do not need to run PreProcessing.R, can fit models directly) 



2) Longitudinal_ModelFitting.R
Takes all_processed_dfs and extract ratio data
Runs select_targets function to get targets / classification questions to be evaluated
Runs bionomial / multinomial model fitting and evaulation (using train_binary_models / train_multinomial_models functions respectively)
Optionally can use train_multiple_targets function to train multiple targets simulataneously 
Runs train_multiple_targets_univariate function to repeat  bionomial analysis using only a single variable as input



3) CrossSectional_ModelFitting.R
Runs same function as Longitudinal_ModelFitting but on cross_sectional_data (can be extracted from all_processed_dfs)

4) PostProcessing.R
Imports modeling results and plots AUCROC, AUPRC, Brier and Strat Brier 


5) PreProcessing_withCHIKV.R and CHIKV_dengue_Longitudinal_ModelFitting.R 
Longitudinal binomial analysis for CHIKV vs dengue 



