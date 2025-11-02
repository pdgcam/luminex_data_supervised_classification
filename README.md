# Luminex Classification Analysis

A machine learning pipeline for predicting flavivirus infections using serological data.

---

## Data Structure

### Analysis Types
- **Univariate**: Single variable used for prediction
- **Multivariate**: Two or more variables used for prediction

### Analysis Approaches
- **Longitudinal**: Ratio of pre- and post-infection timepoints per sample
- **Cross-sectional**: Last draw per sample

---

## Classification Targets

### Binomial Classifications
- Flavivirus vs. not flavivirus
- Any DENV vs. no DENV
- ZIKV vs. no ZIKV
- Any DENV vs. ZIKV (given flavivirus infection)
- CHIKV vs. DENV
  
### Multinomial Classifications
- Any DENV vs. ZIKV vs. no infection 
- DENV serotype classification (given infection)
- DENV serotype vs. negative


---

## Models

### Binomial Models
- **Multivariate**: GLM, Random Forest, Support Vector Machine
- **Univariate**: GLM, Decision Tree, Support Vector Machine

### Multinomial Models
- **Multivariate**: Random Forest, Naive Bayes

---

## Evaluation Metrics

- **AUROC** - Area Under the Receiver Operating Characteristic curve
- **AUPRC** - Area Under the Precision-Recall Curve
- **Brier Score** - Probabilistic prediction accuracy
- **Stratified Brier Score** - For multinomial classifications

---

## Pipeline Scripts

### 1. `PreProcessing.R`

**Purpose**: Prepares raw serological data for model fitting

**Input**:
- Dataset with logged MFI values
- Target infection labels
- Metadata (patient ID, days since infection)

**Process**:
- Calculates log mean HI ratio between blood draws using HI titres
- Classifies samples as "no infection" if ratio < 1.6
- Retains only "no infection" samples for model fitting
- Generates two datasets:
  - `ratio_data` - For longitudinal analysis (post/pre infection ratios per sample)
  - `cross_sectional_data` - For cross-sectional analysis (MFI at last blood draw per sample)

**Output**: `all_processed_dfs` containing both analysis-ready datasets

**Note**: If data is already formatted correctly, skip directly to model fitting scripts.

---

### 2. `Longitudinal_ModelFitting.R`

**Purpose**: Trains and evaluates models on longitudinal (ratio) data

**Process**:
1. Extracts `ratio_data` from `all_processed_dfs`
2. Selects classification targets using [`select_targets()`](https://github.com/pdgcam/luminex_data_supervised_classification/blob/157dc37d8f1428e468f60d98c171c1801f78b9be/Functions.R#L80) function
3. Trains and evaluates models:
   - Binomial models: [`train_binary_models()`](https://github.com/pdgcam/luminex_data_supervised_classification/blob/1ad414c5ab848c14bfaac784e8273f3822aef94c/Functions.R#L170)
   - Multinomial models: [`train_multinomial_models()`](https://github.com/pdgcam/luminex_data_supervised_classification/blob/1ad414c5ab848c14bfaac784e8273f3822aef94c/Functions.R#L421)
   - Multiple targets: [`train_multiple_targets()`](https://github.com/pdgcam/luminex_data_supervised_classification/blob/1ad414c5ab848c14bfaac784e8273f3822aef94c/Functions.R#L732)
 (optional, if want to run on multiple targets simulataneously)
4. Performs univariate analysis: [`train_multiple_targets_univariate()`](https://github.com/pdgcam/luminex_data_supervised_classification/blob/1ad414c5ab848c14bfaac784e8273f3822aef94c/Functions.R#L845)


---

### 3. `CrossSectional_ModelFitting.R`

**Purpose**: Trains and evaluates models on cross-sectional data

**Process**: Identical to `Longitudinal_ModelFitting.R`, but operates on `cross_sectional_data` extracted from `all_processed_dfs`

---

### 4. `PostProcessing.R`

**Purpose**: Visualises model performance metrics

**Process**:
- Imports modeling results
- Generates plots for:
  - AUROC
  - AUPRC
  - Brier Score
  - Stratified Brier Score

---

### 5. `PreProcessing_withCHIKV.R` & `CHIKV_dengue_Longitudinal_ModelFitting.R`

**Purpose**: Pipeline for CHIKV vs. dengue longitudinal binomial analysis

---

## 📦 Requirements

- R with required packages for GLM, Random Forest, SVM, Naive Bayes, and Decision Trees
- Serological data with MFI values and infection labels
- Metadata including patient IDs and infection timepoints
