
# **Description of the dataset **

Author: Marie Alexandre  
Date: 2025-08-12

The file *Simulated_dataset_VL_BAb_NAb_Delta_PostExpo.txt* contains data
(viral and antibody dynamics) that have been simulated from our original
dataset.

| Column.Name | Type | Description |
|:---|:---|:---|
| SubjectID | character | Unique ID of NHPs included in the preclinical trial |
| Group | character | Categorical variable defining the immunological group to which NHP belongs at the time of exposure (Delta infection). Six group categories are possible here: (1) `Naive` - NHPs free of any previous immunization, (2) `CD40.RBDv - Naive` - NHPs are only vaccinated with CD40.RBDv vaccine, (3) `Conv` - NHPs have been previously infected by Wuhan SARS-CoV-2, (4) `CD40.RBDv - Conv` - NHPs have been previously infected by Wuhan SARS-CoV-2 and then received CD40.RBDv vaccine, (5) `CD40.PanCoV - Conv` - NHPS have been previously infected by Wuhan SARS-CoV-2 and then received CD40.PanCoV vaccine, (6) `mRNA - Conv` - NHPs have been previously infected by Wuhan SARS-CoV-2 and then received original pfizer mRNA vaccine |
| StatutImm | character | Categorical variable defining whether NHPs have already been immunized (`Convalescent`) or not (`Naive`). |
| DayPostExpo | Integer | Time in days since the exposure to Delta SARS-CoV-2 |
| Value | numeric | value of measurements collected in natural scale (cp/mL for viral loads, AU/mL for antibodies) |
| Censored | Binary | Binary variable defining whether the observation is left-censored (`1`) or not (`0`). Limits of detection are fixed at 476 and 749 cp/mL for genomic and subgenomic viral loads, respectively, and at 0.1 AU/mL for neutralizing antibodies (no LOD for binding antibodies). |
| CensoredValue | numeric | Observations accounting for censoring, i.e. either directly measured for uncensored data, or imputed to LOD value for left-censored data. |
| MeasureType | Character | Categorical variable defining the type of data is observed: (1) `gRNA` for genomic viral load, (2) `sgRNA` for subgenomic viral load, (3) `IgG CoV-2 RBD L452R T478K` for binding IgG against Delta RBD and (4) `ACE2 CoV-2 RBD L452R T478K` for antibodies inhibiting binding of the RBD domain to the ACE2 receptor. |
| SampleType | character | Categorical variable defining in which fluid measurements have been performed: (1) `Tracheal fluid` for viral loads collected in the trachea, (2) `Nasal fluid` for viral loads collected in the nasopharynx, and (3) `Serum` for antibodies (whether binding or neutralizing) collected in blood samples |
| Observation | numeric | Value of measurements collected in log10 scale (log10 cp/mL for viral loads and log10 AU/mL for antibodies) |
| CensObservations | numeric | Observations in log10 scale accounting for censoring, i.e. either directly measured for uncensored data, or imputed to LOD value (in log10 scale) for left-censored data. Column used as observation in Monolix for model estimation. |
| obsid | integrer | categorical variable added for Monolix and defining at each row the type of data observed: `1` for gRNA viral load in Trachea, `2` for sgRNA viral load in Trachea, `3` for gRNA viral load in Nasopharynx, `4` for sgRNA viral load in Nasopharynx, `5` for binding antiodies and `6` for neutralizing antibodies. |
| Inoc | numeric | Concentration of virions inoculated at exposure, in virions/mL. Used by Monolix as initial condition for inoculum compartment. |
| Inoc_VoC | character | Strain of SARS-CoV-2 virus used for infection. |
| Weight | numeric | Weight of NHPs, in Kgs |
| BAb0 | numeric | Value in AU/mL of binding antibodies collected at the time of exposure. |
| NAb0 | numeric | Value in AU/mL of neutralizing antibodies collected at the time of exposure. |
