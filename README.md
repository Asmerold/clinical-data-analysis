# Analysis of Biochemical Relapse After Radiotherapy in Localized Prostate Cancer

## Overview
This R project processes longitudinal clinical data to define PSA nadir, identify Phoenix-defined biochemical relapse, and construct time-to-event endpoints (BCR, MFS, OS) for patients with localized prostate cancer treated with radiotherapy.

## Key Features
- Data cleaning, outlier filtering, date parsing.
- PSA nadir & Phoenix relapse (+2 ng/mL above nadir).
- Time-to-event endpoints: BCR, MFS, OS, CSS.
- Competing risks: Fine-Gray model (tidycmprsk).
- Landmark analysis (12 months).
- Cohort creation for localized prostate cancer.

## Visualizations
- Kaplan-Meier curves
- Cox regression forest plots
- Fine-Gray cumulative incidence curves

## Files in this repository
- `df4 (5).R` – main analysis script (data preparation, survival endpoints, landmark)
- `forest_plots R.R` – script for forest plots (Cox models)
- `*.png` – example figures (Kaplan-Meier, forest plots, etc.)

## Dependencies
R packages: dplyr, tidyr, lubridate, survival, survminer, tidycmprsk, ggsurvfit, tidyverse, ggplot2.

## Author
Vadim Kulikov, MD  
Clinical Data Analyst portfolio project
