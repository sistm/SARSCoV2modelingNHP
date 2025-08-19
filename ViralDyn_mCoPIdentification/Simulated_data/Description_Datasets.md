
# **Description of the dataset **

Author: Marie Alexandre  
Date: 2025-08-19

These files have been written considering datasets uploaded on Monolix
for the model parameter estimation and are structured as follows:

| Column.Name | Type | Description | Model.Type |
|:---|:---|:---|:---|
| AnimalID | character | Unique ID of NHPs included in the preclinical trial | Column used for all models |
| Time | integer | Time in days since the exposure to Wuhan SARS-CoV-2 | Column used for all models |
| log10VL | numeric | Viral load measurement collected in log10 scale (log10 cp/mL) | Column used for all models |
| obs | integer | Categorical variable added for Monolix and defining at each row the type of data observed: `1` for gRNA viral load in Trachea, `2` for sgRNA viral load in Trachea, `3` for gRNA viral load in Nasopharynx, `4` for sgRNA viral load in Nasopharynx. | Column used for all models |
| Group | character | Categorical variable defining the immunological group to which NHP belongs at the time of exposure (Wuhan infection). Three group categories are possible here: (1) `Naive` - NHPs free of any previous immunization, (2) `Convalescent` - NHPs have been previoulsy infected by Wuhan SARS-CoV-2, (3) `Conv-CD40` - NHPs have been previously infected by Wuhan SARS-CoV-2 and then received CD40.RBD vaccine. | Column used for all models |
| Type | character | Categorical variable defining the type of data is observed: (1) `Total` for genomic viral load, (2) `sg` for subgenomic viral load | Column used for all models |
| organ | character | Categorical variable defining in which URT compartment measurements have been performed: (1) `Trachea` for viral loads collected in the trachea, (2) `Naso` for viral loads collected in the nasopharynx. | Column used for all models |
| Weight | numeric | Weight of NHPs, in Kgs, used as regressor in Monolix | Column used for all models |
| censored | binary | Binary variable defining whether the observation is left-censored (`1`) or not (`0`). Limits of detection are fixed at 476 and 749 cp/mL for genomic and subgenomic viral loads, respectively. | Column used for all models |
| Marker | numeric | Baseline value of the marker (ECL here) in natural scale, used as regressor in Monolix | Column used only for model adjusted for baseline covariates. |
| MARKER.Slope | numeric | One column for each potential marker used as covariate providing the `value of Slopes` of the immune marker (e.g. ECL) calculated and used to perform linear interpolation in Monolix to obtain continous dynamics of the covariate(s). Column(s) used as regressors in Monolix. | Column(s) used only for model adjusted for for time-varying covariates. |
| MARKER.Intercept | numeric | One column for each potential marker used as covariate providing the `value of intercepts` of the immune marker (e.g. ECL) calculated and used to perform linear interpolation in Monolix to obtain continous dynamics of the covariate(s). Column(s) used as regressors in Monolix. | Column(s) used only for model adjusted for for time-varying covariates. |
