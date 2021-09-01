# SARSCoV2modelingNHP
The code provided is this folder is an original mlxtran code describing mechanistic models developed to fit individual viral load dynamics in NHP vaccine studies.

## Mlxtran code
Each text file correspond to the mlxtran code describing a specific model and need to be uploaded within the Monolix software (developed by Lixoft) as model file to be used.
In this folder, 9 model files are given. They can be clustered into 3 groups 
1. The basic model used in the project to fit viral load dynamics without any adjustment of parameters for covariates. However, adjustment for categorical covariates defined in the datasets could be added directly in the Monolix software [*MlxtranCode_MechanisticModel_SARSCoV2_VL_NoExternalCov.txt*]
2. The basic model with two model parameters adjusted for 1 or 2 time-varying covariates (beta and delta) as described below:
   - The parameter beta adjusted for 1 covariate [*MlxtranCode_MechanisticModel_SARSCoV2_VL_TimeVaryingCovBeta.txt*]
   - The parameter delta adjusted for 1 covariate [*MlxtranCode_MechanisticModel_SARSCoV2_VL_TimeVaryingCovDelta.txt*]
   - The parameter beta adjuste for 2 covariates [*MlxtranCode_MechanisticModel_SARSCoV2_VL_2TimeVaryingCovBeta.txt*]
   - The parmeter beta adjusted for 1 covariate and delta for the other [*MlxtranCode_MechanisticModel_SARSCoV2_VL_1TimeVaryingCovBeta_1TimeVaryingCovDelta.txt*]. 
3. The basic model with two model parameters adjusted for 1 or 2 baseline covariates (covariates measured at the initial time of the study) (beta and delta) as described below:
   - The parameter beta adjusted for 1 covariate [*MlxtranCode_MechanisticModel_SARSCoV2_VL_BaselineCovBeta.txt*]
   - The parameter delta adjusted for 1 covariate [*MlxtranCode_MechanisticModel_SARSCoV2_VL_BaselineCovDelta.txt*]
   - The parameter beta adjuste for 2 covariates [*MlxtranCode_MechanisticModel_SARSCoV2_VL_2BaselineCovBeta.txt*]
   - The parmeter beta adjusted for 1 covariate and delta for the other [*MlxtranCode_MechanisticModel_SARSCoV2_VL_1BaselineCovBeta_1BaselineCovDelta.txt*]. 

## [Monolix] Data formatting
These files have be written considering datasets, uploaded on Monolix for the model parameter estimation, structured as follows: 
1) one column of Ids
2) one column for time of observations
3) one column for observations gathering four types observations: both gRNA viral load and sgRNA viral loads measured in two distinct compartments (Trachea and Nasopharynx). 
4) one column for observation id (four in this specific case).
5) one column for censorship
6) one column for individual weights, defined as regressor
Additional columns can be considered for categorical covariates.

For the models adjusted for time-varying covariates, the following elements must be added:

7) one column **for each** covariate providing values of slopes used to perform linear interpolation and obtain continuous dynamics of the covariate(s). These columns are then used as regressor by Monolix
8) one column **for each** covariate providing values of intercepts used to perform linear interpolation and obtain continuous dynamics of the covariate(s). These columns are then used as regressor by Monolix


For the models adjusted for baseline covariates, the following elements must be added:

7) one column **for each** covariate providing individual baseline values. These columns are then used as regressor by Monolix




