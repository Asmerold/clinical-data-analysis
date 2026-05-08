library(survival)
library(dplyr)
library(lubridate)

options(scipen = 999)

# 1. First add the necessary variables (year and age) if they are not already in the database
df_final <- df_final %>%
  mutate(
    rt_year = year(rt_start_date),
    age_at_diag = as.numeric(difftime(diag_date, birth_date, units = "days")) / 365.25
  )

cat("\n============================================\n")
cat("1. BIOCHEMICAL RELAPSE-FREE SURVIVAL (BCR)\n")
cat("============================================\n")
# Убрали ttn_months!
cox_bcr <- coxph(Surv(bcr_months_15y, bcr_event_15y) ~ phenotype +
                   gleason_group + age_at_diag + psa_baseline_static +
                   rt_year + adt_status + max_dose_prostate,
                 data = df_loc)
print(summary(cox_bcr))

cox_bcr_separate <- coxph(Surv(bcr_months_15y, bcr_event_15y) ~ depth_cat + speed_cat +
                   gleason_group + age_at_diag + psa_baseline_static +
                   rt_year + adt_status + max_dose_prostate,
                 data = df_loc)
print(summary(cox_bcr_separate))

cox_bcr_mult <- coxph(Surv(bcr_months_15y, bcr_event_15y) ~ depth_cat * speed_cat +
                            gleason_group + age_at_diag + psa_baseline_static +
                            rt_year + adt_status + max_dose_prostate,
                          data = df_loc)
print(summary(cox_bcr_mult))


cat("\n============================================\n")
cat("2. OVERALL SURVIVAL (OS)\n")
cat("============================================\n")
cox_os <- coxph(Surv(os_months_15y, os_event_15y) ~ phenotype +
                  gleason_group + age_at_diag + psa_baseline_static +
                  rt_year + adt_status + max_dose_prostate,
                data = df_loc)
print(summary(cox_os))


cat("\n============================================\n")
cat("3. CANCER-SPECIFIC SURVIVAL (CSS)\n")
cat("============================================\n")
# Removed adt_status (to prevent HR from going into millions)!
cox_css <- coxph(Surv(os_months_15y, css_event_15y) ~ phenotype +
                   gleason_group + age_at_diag + psa_baseline_static +
                   rt_year + max_dose_prostate,
                 data = df_loc)
print(summary(cox_css))
# ---------------------------------------------------------------------------------------------------
# MULTIVARIABLE FINE-GRAY MODEL (Cancer-specific mortality)
# ---------------------------------------------------------------------------------------------------
fg_multi <- tidycmprsk::crr(Surv(os_months_15y, fg_status_15y) ~ phenotype +
                              gleason_group + age_at_diag + psa_baseline_static +
                              rt_year + max_dose_prostate,
                            data = df_loc)

# Print nice results
fg_multi %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  print()
