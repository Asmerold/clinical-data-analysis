library(dplyr)
library(lubridate)
library(survival)
library(survminer)   # For Kaplan-Mayer and Forest Plots
library(tidycmprsk)  # For Fine-Gray
library(ggsurvfit)  # For figures Fine-Gray (ggcuminc)
library(tidyverse)

# ==================================================================================================
# 1-2. Upload and cleaning files
# ==================================================================================================
raw_df <- read.csv("G:/Мой диск/статьи 2026/статья по надиру/df3.csv",
                   sep=";", dec=",", fileEncoding = "windows-1251")

clean_num <- function(x) { as.numeric(gsub(",", ".", as.character(x))) }

df <- raw_df %>%
  mutate(
    across(contains("date"), dmy),
    psa_value = clean_num(psa_value),
    psa_baseline_static = clean_num(psa_baseline_static),
    dose_prostate = clean_num(dose_prostate),
    dose_nodes = clean_num(dose_nodes),
    prostate_volume = clean_num(prostate_volume),
    orch_flag = as.logical(orch_flag)
  ) %>%
  arrange(patient_id, event_date)

# ==================================================================================================
# 3. CHRONOLOGICAL CALCULATION OF NADIR AND PHOENIX
# ==================================================================================================
psa_logic <- df %>%
  group_by(patient_id) %>%
  mutate(rt_start = min(event_date[dose_prostate > 0], na.rm = TRUE)) %>%
  filter(event_date >= rt_start & !is.na(psa_value)) %>%
  mutate(
    running_nadir = cummin(psa_value),
    is_relapse = psa_value >= (running_nadir + 2)
  ) %>%
  summarise(
    rt_start_date = first(rt_start),
    relapse_date_phoenix = min(event_date[is_relapse], na.rm = TRUE),
    psa_nadir_true = min(psa_value[event_date <= coalesce(relapse_date_phoenix, max(event_date))], na.rm = TRUE),
    nadir_date = event_date[psa_value == psa_nadir_true & event_date <= coalesce(relapse_date_phoenix, max(event_date))][1],
    .groups = "drop"
  )

# ==================================================================================================
# 4-5. BUILD MASTER TABLE AND REMOVE INFINITES
# ==================================================================================================
master_table <- df %>%
  group_by(patient_id) %>%
  summarise(
    birth_date = first(birth_date),
    diag_date = first(diag_date),
    t_stage = first(t_stage),
    n_stage = first(n_stage),
    m_stage = first(m_stage),
    gleason_score = first(gleason_score),
    psa_baseline_static = first(psa_baseline_static),
    prostate_volume = max(prostate_volume, na.rm = TRUE), # Keep volume
    risk_group = first(risk_group),
    adt_type_static = first(adt_type_static),
    metastasis_date = first(metastasis_date),
    death_date = first(death_date),
    max_dose_prostate = max(dose_prostate, na.rm = TRUE),
    max_dose_nodes = max(dose_nodes, na.rm = TRUE),
    max_dose_segment = max(dose_segment, na.rm = TRUE),

    # SUPER‑HT logic
    adt_start = min(event_date[adt_ag == "s" | adt_aa == "s" | adt_gt2 == "s" | orch_flag == TRUE], na.rm = TRUE),
    adt_end = max(event_date[adt_ag == "f" | adt_aa == "f" | adt_gt2 == "f"], na.rm = TRUE),
    is_orch = any(orch_flag == TRUE),

    last_status = last(status_text[!is.na(status_text)]),
    last_contact = max(event_date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  inner_join(psa_logic, by = "patient_id") %>%
  mutate(across(where(is.Date), ~ { .x[is.infinite(.x)] <- NA; .x })) %>%
  mutate(across(where(is.numeric), ~ { .x[is.infinite(.x)] <- NA; .x }))

master_table <- master_table %>%
  mutate(last_status = str_replace_all(last_status, ",", ""))

# ==================================================================================================
# 6. BASIC CALCULATIONS, CLINICAL GROUPS AND FLAGS
# ==================================================================================================
df_base <- master_table %>%
  mutate(
    # HT status and PSA density
    adt_status = ifelse(!is.na(adt_start) | (!is.na(adt_type_static) & adt_type_static != ""), 1, 0),
    psa_density = psa_baseline_static / prostate_volume,

    # Salvage HT flag (started AFTER relapse)
    adt_after_relapse = ifelse(!is.na(adt_start) & !is.na(relapse_date_phoenix) & adt_start > relapse_date_phoenix, 1, 0),

    ttn_months = as.numeric(difftime(nadir_date, rt_start_date, units = "days")) / 30.44,

    # OS (Overall survival)
    os_event = ifelse(last_status %in% c("cancer_death", "other_death"), 1, 0),
    os_months = as.numeric(difftime(coalesce(death_date, last_contact), rt_start_date, units = "days")) / 30.44,

    # BCR (Biochemical relapse)
    bcr_event = ifelse(!is.na(relapse_date_phoenix), 1, 0),
    bcr_months = case_when(
      bcr_event == 1 ~ as.numeric(difftime(relapse_date_phoenix, rt_start_date, units = "days")) / 30.44,
      TRUE ~ os_months
    ),

    # MFS (Metastasis‑free survival)
    mfs_event = ifelse(!is.na(metastasis_date), 1, 0),
    mfs_months = case_when(
      mfs_event == 1 ~ as.numeric(difftime(metastasis_date, rt_start_date, units = "days")) / 30.44,
      TRUE ~ os_months
    ),

    # CSS and FG Status
    css_event_km = ifelse(last_status == "cancer_death", 1, 0),
    fg_status = case_when(
      last_status == "cancer_death" ~ "1_cancer",
      last_status == "other_death" ~ "2_other",
      TRUE ~ "0_censored"
    ),

    # Phenotypes and groups
    depth_cat = ifelse(psa_nadir_true < 0.5, "Deep", "Shallow"),
    speed_cat = ifelse(ttn_months < 12, "Fast", "Slow"),
    phenotype = paste(depth_cat, speed_cat, sep = " & "),

    gleason_cat = case_when(
      gleason_score <= 6 ~ "<= 6",
      gleason_score == 7 ~ "7",
      gleason_score == 8 ~ "8",
      gleason_score >= 9 ~ "9-10",
      TRUE ~ NA_character_
    ),

    clinical_group = case_when(
      t_stage <= 2 & n_stage == 0 & m_stage == 0 ~ "1_Localized",
      (t_stage <= 2 & n_stage == 1 & m_stage == 0) | (t_stage >= 3 & m_stage == 0) ~ "2_Locally_Advanced",
      m_stage == 1 ~ "3_Metastatic",
      TRUE ~ "1_Localized"
    )
  ) %>%
  filter(os_months > 0 & ttn_months >= 0) %>%
  distinct(birth_date, diag_date, .keep_all = TRUE)

# ==================================================================================================
# 7. 15‑YEAR AND 12‑MONTH LANDMARK VARIABLES
# ==================================================================================================
max_followup <- 180
lm_time <- 12

df_final <- df_base %>%
  mutate(
    # --- 15 YEAR CENSORING ---
    os_event_15y = ifelse(os_months > max_followup, 0, os_event),
    os_months_15y = ifelse(os_months > max_followup, max_followup, os_months),

    bcr_event_15y = ifelse(bcr_months > max_followup, 0, bcr_event),
    bcr_months_15y = ifelse(bcr_months > max_followup, max_followup, bcr_months),

    mfs_event_15y = ifelse(mfs_months > max_followup, 0, mfs_event),
    mfs_months_15y = ifelse(mfs_months > max_followup, max_followup, mfs_months),

    css_event_15y = ifelse(os_months > max_followup, 0, css_event_km),
    fg_status_15y = ifelse(os_months > max_followup, "0_censored", fg_status),

    # --- 12-MONTH LANDMARK ---
    # OS/CSS/FG Landmark
    os_months_lm12 = ifelse(os_months >= lm_time, os_months_15y - lm_time, NA_real_),
    os_event_lm12 = ifelse(os_months >= lm_time, os_event_15y, NA_real_),
    css_event_lm12 = ifelse(os_months >= lm_time, css_event_15y, NA_real_),
    fg_status_lm12 = ifelse(os_months >= lm_time, fg_status_15y, NA_character_),

    # BCR Landmark
    bcr_months_lm12 = ifelse(bcr_months >= lm_time, bcr_months_15y - lm_time, NA_real_),
    bcr_event_lm12 = ifelse(bcr_months >= lm_time, bcr_event_15y, NA_real_),

    # MFS Landmark
    mfs_months_lm12 = ifelse(mfs_months >= lm_time, mfs_months_15y - lm_time, NA_real_),
    mfs_event_lm12 = ifelse(mfs_months >= lm_time, mfs_event_15y, NA_real_)
  ) %>%

  mutate(
    # 1. Those where order is critical (reference)
    phenotype = factor(phenotype, levels = c("Deep & Slow", "Deep & Fast", "Shallow & Fast", "Shallow & Slow")),
    gleason_cat = factor(gleason_cat, levels = c("<= 6", "7", "8", "9-10")),
    adt_status = factor(adt_status, levels = c(0, 1)),
    fg_status_lm12 = factor(fg_status_lm12, levels = c("0_censored", "1_cancer", "2_other")),

    # 2. All other factors (convert to factor without specific order)
    across(c(adt_after_relapse, is_orch, depth_cat, speed_cat,
             clinical_group, fg_status, fg_status_15y), as.factor)
  )

df_final <- df_final %>%
  mutate(
    # Year of radiotherapy start
    rt_year = year(rt_start_date),

    # Age at diagnosis (days / 365.25)
    age_at_diag = as.numeric(difftime(diag_date, birth_date, units = "days")) / 365.25
  )

df_final <- df_final %>%
  filter(ttn_months <= 36,           # Remove those who took years to reach nadir (noise)
         prostate_volume < 400,
         max_dose_prostate < 100,
         max_dose_prostate >= 56,
         !is.na(t_stage),
         !is.na(m_stage),
         !is.na(n_stage))


save(df_final, file = "G:/Мой диск/статьи 2026/статья по надиру/df_final.RData")

# Final check
cat("\n=== ULTIMATE DATABASE READY ===")
cat("\nTotal columns:", ncol(df_final))
cat("\nTotal patients:", nrow(df_final))
cat("\nPatients with salvage HT (adt_after_relapse=1):", sum(df_final$adt_after_relapse == "1", na.rm=T))



# ====================================================================================================
# LOCALIZED PROSTATE CANCER
# ====================================================================================================
# ==============================================================================
# FINAL FORMATION OF LOCALIZED CANCER GROUP (df_loc)
# ==============================================================================

df_loc <- df_final %>%
  filter(
    # 1. Only localized cancer
    clinical_group == "1_Localized",

    # 2. Clinical sample purity
    ttn_months <= 36,               # Nadir within 3 years
    psa_baseline_static <= 40,      # Your new threshold (exclude extremely high PSA)
    max_dose_prostate >= 60,        # Only radical doses
    max_dose_prostate < 100,        # Remove dose entry errors
    prostate_volume < 400,          # Remove giant glands (noise)
    
    # 3. Remove missing values in critical variables for Cox
    !is.na(t_stage),
    !is.na(n_stage),
    !is.na(m_stage),
    !is.na(gleason_score)
  ) %>%
  mutate(
    # 4. Gleason grouping (merge 8‑10 for model stability)
    gleason_group = case_when(
      as.numeric(as.character(gleason_score)) <= 6 ~ "\u2264 6",
      as.numeric(as.character(gleason_score)) == 7 ~ "7",
      as.numeric(as.character(gleason_score)) >= 8 ~ "8-10"
    ),
    gleason_group = factor(gleason_group, levels = c("\u2264 6", "7", "8-10")),

    # 5. Set reference level for phenotypes
    phenotype = factor(phenotype, levels = c("Deep & Slow", "Deep & Fast", "Shallow & Fast", "Shallow & Slow"))
  )

# Remove duplicates
df_loc <- df_loc %>% distinct(patient_id, .keep_all = TRUE)

# Check final size
cat("\n=== df_loc FORMATION COMPLETED ===")
cat("\nNumber of patients in df_loc:", nrow(df_loc))
cat("\nNumber of events (relapses):", sum(df_loc$bcr_event_15y == 1))
cat("\n======================================\n")

# Save the ready group into our .RData file (overwrites it)
save(df_final, df_loc, file = "G:/Мой диск/статьи 2026/статья по надиру/final_data.RData")
