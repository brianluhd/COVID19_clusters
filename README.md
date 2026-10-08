# Validation of a Diabetes Subtype Classification Model Using Data from U.S. Adults Before and After the COVID-19 Pandemic 

Analysis code for [Lu et al., 2026.](https://doi.org/10.3390/metabo16030204)

The fixed effect components of the accompanying linear mixed-effect models for scaling insulin-derived HOMA values to the c-peptide derived HOMA1 values used in the study can be found in exports/HOMA_regression.

# Model usage

**C-peptide and insulin values have different clinical implications and are not interchangeable.**
**As such, the derived HOMA values also have different clinical interpretations and are not clinically interchangeable.** 
**The scaled HOMA values are intended to enable diabetes subtype assignments and should not be used directly for clinical interpretations or clinical decisions.**

The HOMA1-B scaling model requires a tibble input with columns “insulin_homa1b”, “race”, and “gender”. 
Allowed input values for “race” are “Non-Hispanic White”, “Non-Hispanic Black”, and “Other Race - Including Multi-Racial”. 
Allowed input values for “gender” are “Female” and “Male”.

The HOMA1-IR scaling model requires a tibble input with columns “insulin_homa1ir”, “race”, and “gender”. 
Allowed input values for “race” are “Non-Hispanic White”, “Non-Hispanic Black”, and “Other Race - Including Multi-Racial”. 
Allowed input values for “gender” are “Female” and “Male”.

Insulin-derived HOMA values were calculated as original described by [Matthews et al. 1985](https://doi.org/10.1007/bf00280883).
C-peptide-derived HOMA values were calculated as described in [Lu et al., 2026.](https://doi.org/10.3390/metabo16030204).

The scaled c-peptide HOMA1 values can be predicted with

```
library(dplyr)
library(lmerTest)

homa1b_fe_model <- readRDS("exports/HOMA_regression/homa1b_fe_model.RDS")
homa1ir_fe_model <- readRDS("exports/HOMA_regression/homa1ir_fe_model.RDS")

data <- data %>%
  mutate(homa1b_pred = exp(homa1b_fe_model$predict(homa1b_fe_model, newdata = data)), 
         homa1ir_pred = exp(homa1ir_fe_model$predict(homa1ir_fe_model, newdata = data)))

```