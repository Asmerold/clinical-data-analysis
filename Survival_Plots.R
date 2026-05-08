library(dplyr)
library(lubridate)
library(survival)
library(survminer)
library(tidycmprsk)
library(survRM2)
library(ggsurvfit)
library(cmprsk)
library(broom)

# ----------------------------------------
# LOCALIZED OVERALL SURVIVAL
# ----------------------------------------

fit_os <- survfit(Surv(os_months_15y, os_event_15y) ~ phenotype, data = df_loc)
sd_os <- survdiff(Surv(os_months_15y, os_event_15y) ~ phenotype, data = df_loc)

# Print the test result (including p-value) to the console
print(sd_os)

# 1. Pairwise comparisons (WHO DIFFERS FROM WHOM)
# This function returns a matrix of p-values
pairwise_os_loc <- pairwise_survdiff(Surv(os_months_15y, os_event_15y) ~ phenotype, data = df_loc)
print("Pairwise differences OS (Localized):")
print(pairwise_os_loc)


# 1. "Bulletproof" function to calculate RMST
calc_rmst_pair <- function(data, group1, group2) {
  # Prepare data without missing values
  sub_data <- data %>%
    filter(phenotype %in% c(group1, group2)) %>%
    filter(!is.na(os_months_15y) & !is.na(os_event_15y)) %>%
    mutate(arm = ifelse(phenotype == group1, 1, 0))

  # AUTOMATIC HORIZON (to avoid tau errors)
  max_t1 <- max(sub_data$os_months_15y[sub_data$arm == 1], na.rm = TRUE)
  max_t0 <- max(sub_data$os_months_15y[sub_data$arm == 0], na.rm = TRUE)
  safe_tau <- min(max_t1, max_t0, 180) # Take 15 years (180) or the maximum available

  # Protection for cases where there are no events
  if(safe_tau <= 0 || nrow(sub_data) == 0) {
    return(data.frame(Comparison = paste(group1, "vs", group2),
                      Horizon_Months = 0,
                      Diff_Months = NA,
                      CI_95 = "Unable to calculate",
                      p_value = "NA"))
  }
  # Calculate RMST
  res <- rmst2(time = sub_data$os_months_15y,
               status = sub_data$os_event_15y,
               arm = sub_data$arm,
               tau = safe_tau)

  # Correctly extract data (unadjusted.result)
  diff <- res$unadjusted.result[1, "Est."]
  low  <- res$unadjusted.result[1, "lower .95"]
  high <- res$unadjusted.result[1, "upper .95"]
  p    <- res$unadjusted.result[1, "p"]

  # Return a ready row
  return(data.frame(Comparison = paste(group1, "vs", group2),
                    Horizon_Months = round(safe_tau, 1), # Show the horizon up to which we calculate
                    Diff_Months = round(diff, 1),
                    CI_95 = paste0("(", round(low, 1), " to ", round(high, 1), ")"),
                    p_value = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))))
}


# 3. Compare our reference (Deep & Slow) with all others
pair1 <- calc_rmst_pair(df_loc, "Deep & Slow", "Deep & Fast")
pair2 <- calc_rmst_pair(df_loc, "Deep & Slow", "Shallow & Fast")
pair3 <- calc_rmst_pair(df_loc, "Deep & Slow", "Shallow & Slow")

# 4. Combine everything into one table
rmst_summary_table <- bind_rows(pair1, pair2, pair3)

# Print the result
print("=== Summary table of clinical benefit (OS Gain) ===")
print(rmst_summary_table)

# ====================================================================================================
# KM PLOTS
# ====================================================================================================
ggsurvplot(fit_os,
           data = df_loc,
           pval = TRUE,
           censor = FALSE,
           break.time.by = 36,
           xlim = c(0, 180),
           palette = c("orangered2", "royalblue2", "limegreen", "gold2"),
           title = "Overall Survival by Phenotype Localized Cancer",
           risk.table = FALSE,
           ggtheme = theme_bw())

print(fit_os)

summary(fit_os, times = c(12,24,36,48,60,72,84,96,108,120))
summary(fit_os, times = c(12,36,60,120, 180))

# ----------------------------------------
# LOCALIZED BIOCHEMICAL RELAPSE-FREE SURVIVAL
# ----------------------------------------

fit_bcr <- survfit(Surv(bcr_months_15y, bcr_event_15y) ~ phenotype, data = df_loc)

# 1. Pairwise comparisons (WHO DIFFERS FROM WHOM)
# This function returns a matrix of p-values
pairwise_bcr_loc <- pairwise_survdiff(Surv(bcr_months_15y, bcr_event_15y) ~ phenotype, data = df_loc)
print("Pairwise differences BCR (Localized):")
print(pairwise_bcr_loc)


# 1. "Bulletproof" function to calculate RMST
calc_rmst_pair <- function(data, group1, group2) {
  # Prepare data without missing values
  sub_data <- data %>%
    filter(phenotype %in% c(group1, group2)) %>%
    filter(!is.na(bcr_months_15y) & !is.na(bcr_event_15y)) %>%
    mutate(arm = ifelse(phenotype == group1, 1, 0))

  # AUTOMATIC HORIZON (to avoid tau errors)
  max_t1 <- max(sub_data$bcr_months_15y[sub_data$arm == 1], na.rm = TRUE)
  max_t0 <- max(sub_data$bcr_months_15y[sub_data$arm == 0], na.rm = TRUE)
  safe_tau <- min(max_t1, max_t0, 180) # Take 15 years (180) or the maximum available

   # Protection for cases where there are no events
  if(safe_tau <= 0 || nrow(sub_data) == 0) {
    return(data.frame(Comparison = paste(group1, "vs", group2),
                      Horizon_Months = 0,
                      Diff_Months = NA,
                      CI_95 = "Unable to calculate",
                      p_value = "NA"))
  }

  # Calculate RMST
  res <- rmst2(time = sub_data$bcr_months_15y,
               status = sub_data$bcr_event_15y,
               arm = sub_data$arm,
               tau = safe_tau)

  # Correctly extract data (unadjusted.result)
  diff <- res$unadjusted.result[1, "Est."]
  low  <- res$unadjusted.result[1, "lower .95"]
  high <- res$unadjusted.result[1, "upper .95"]
  p    <- res$unadjusted.result[1, "p"]

  # Return a ready row
  return(data.frame(Comparison = paste(group1, "vs", group2),
                    Horizon_Months = round(safe_tau, 1), # Show the horizon up to which we calculate
                    Diff_Months = round(diff, 1),
                    CI_95 = paste0("(", round(low, 1), " to ", round(high, 1), ")"),
                    p_value = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))))
}

# 3. Compare our reference (Deep & Slow) with all others
pair1 <- calc_rmst_pair(df_loc, "Deep & Slow", "Deep & Fast")
pair2 <- calc_rmst_pair(df_loc, "Deep & Slow", "Shallow & Fast")
pair3 <- calc_rmst_pair(df_loc, "Deep & Slow", "Shallow & Slow")

# 4. Combine everything into one table
rmst_summary_table <- bind_rows(pair1, pair2, pair3)

# Print the result
print("=== Summary table of clinical benefit (BCR Gain) ===")
print(rmst_summary_table)

# ====================================================================================================
# KM PLOTS
# ====================================================================================================
ggsurvplot(fit_bcr,
           data = df_loc,
           pval = TRUE,
           censor = FALSE,
           break.time.by = 36,
           xlim = c(0, 180),
           palette = c("orangered2", "royalblue2", "limegreen", "gold2"),
           title = "Biochemical Recurrence Free Survival by Phenotype Localized Cancer",
           risk.table = FALSE,
           ggtheme = theme_bw())

print(fit_bcr)

summary(fit_bcr, times = c(12,24,36,48,60,72,84,96,108,120))
summary(fit_bcr, times = c(12,36,60,120, 180))


# ----------------------------------------
# LOCALIZED CANCER-SPECIFIC SURVIVAL
# ----------------------------------------

fit_css <- survfit(Surv(os_months_15y,css_event_15y) ~ phenotype, data = df_loc)

# Pairwise comparisons (WHO DIFFERS FROM WHOM)
# This function returns a matrix of p-values
pairwise_css_localized <- pairwise_survdiff(Surv(os_months_15y, css_event_15y) ~ phenotype, data = df_loc)
print("Pairwise differences CSS (Localized):")
print(pairwise_css_localized)

# ====================================================================================================
# RMST + PAIRWISE COMPARISONS
# ====================================================================================================

# 1. "Bulletproof" function to calculate RMST
calc_rmst_pair <- function(data, group1, group2) {
  # Prepare data without missing values
  sub_data <- data %>%
    filter(phenotype %in% c(group1, group2)) %>%
    filter(!is.na(os_months_15y) & !is.na(css_event_15y)) %>%
    mutate(arm = ifelse(phenotype == group1, 1, 0))

  # AUTOMATIC HORIZON (to avoid tau errors)
  max_t1 <- max(sub_data$os_months_15y[sub_data$arm == 1], na.rm = TRUE)
  max_t0 <- max(sub_data$os_months_15y[sub_data$arm == 0], na.rm = TRUE)
  safe_tau <- min(max_t1, max_t0, 180) # Take 15 years (180) or the maximum available

  # Calculate RMST
  res <- rmst2(time = sub_data$os_months_15y,
               status = sub_data$css_event_15y,
               arm = sub_data$arm,
               tau = safe_tau)

  # Correctly extract data (unadjusted.result)
  diff <- res$unadjusted.result[1, "Est."]
  low  <- res$unadjusted.result[1, "lower .95"]
  high <- res$unadjusted.result[1, "upper .95"]
  p    <- res$unadjusted.result[1, "p"]

  # Return a ready row
  return(data.frame(Comparison = paste(group1, "vs", group2),
                    Horizon_Months = round(safe_tau, 1), # Show the horizon up to which we calculate
                    Diff_Months = round(diff, 1),
                    CI_95 = paste0("(", round(low, 1), " to ", round(high, 1), ")"),
                    p_value = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))))
}


# 3. Compare our reference (Deep & Slow) with all others
pair1 <- calc_rmst_pair(df_loc, "Deep & Slow", "Deep & Fast")
pair2 <- calc_rmst_pair(df_loc, "Deep & Slow", "Shallow & Fast")
pair3 <- calc_rmst_pair(df_loc, "Deep & Slow", "Shallow & Slow")

# 4. Combine everything into one table
rmst_summary_table <- bind_rows(pair1, pair2, pair3)

# Print the result
print("=== Summary table of clinical benefit (CSS Gain) ===")
print(rmst_summary_table)

# ====================================================================================================
# KM PLOTS
# ====================================================================================================
ggsurvplot(fit_css,
           data = df_loc,
           pval = TRUE,
           censor = FALSE,
           break.time.by = 36,
           xlim = c(0, 180),
           palette = c("orangered2", "royalblue2", "limegreen", "gold2"),
           title = "Cancer-Specific Survival by Phenotype Localized Cancer",
           risk.table = TRUE,
           ggtheme = theme_bw())

print(fit_css)

summary(fit_css, times = c(12,24,36,48,60,72,84,96,108,120))
summary(fit_css, times = c(12,36,60,120,180))

# ============================================================================
# FINE - GRAY TEST
# ============================================================================

# === 0. FIX THE REFERENCE (MUST ASSIGN WITH <- ) ===
df_loc <- df_loc %>%
  mutate(phenotype = factor(phenotype, levels = c("Deep & Slow", "Deep & Fast", "Shallow & Fast", "Shallow & Slow")))

cat("\n====================================================================\n")
cat("1. PAIRWISE DIFFERENCES (UNIVARIABLE FINE-GRAY)\n")
cat("====================================================================\n")
# Instead of a matrix of p-values we obtain sHR (how many times higher the risk of dying from cancer)
fg_uni <- tidycmprsk::crr(Surv(os_months_15y, fg_status_15y) ~ phenotype, data = df_loc)

fg_pairwise_table <- fg_uni %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    estimate = round(estimate, 2),
    CI = paste0("(", round(conf.low, 2), " - ", round(conf.high, 2), ")"),
    p.value = round(p.value, 4)
  ) %>% select(term, estimate, CI, p.value)

print(as.data.frame(fg_pairwise_table))


cat("\n====================================================================\n")
cat("2. SUMMARY TABLE OF CLINICAL BENEFIT (RMST ДЛЯ CSS)\n")
cat("====================================================================\n")
# Function to calculate lost life months (using css_event_15y)
calc_rmst_pair_css <- function(data, group1, group2) {
  sub_data <- data %>%
    filter(phenotype %in% c(group1, group2)) %>%
    filter(!is.na(os_months_15y) & !is.na(css_event_15y)) %>%
    mutate(arm = ifelse(phenotype == group1, 1, 0))

  max_t1 <- max(sub_data$os_months_15y[sub_data$arm == 1], na.rm = TRUE)
  max_t0 <- max(sub_data$os_months_15y[sub_data$arm == 0], na.rm = TRUE)
  safe_tau <- min(max_t1, max_t0, 180)

  res <- rmst2(time = sub_data$os_months_15y, status = sub_data$css_event_15y, arm = sub_data$arm, tau = safe_tau)

  diff <- res$unadjusted.result[1, "Est."]
  low  <- res$unadjusted.result[1, "lower .95"]
  high <- res$unadjusted.result[1, "upper .95"]
  p    <- res$unadjusted.result[1, "p"]

  return(data.frame(Comparison = paste(group1, "vs", group2),
                    Horizon = round(safe_tau, 1),
                    Loss_Months = round(diff, 1),
                    CI_95 = paste0("(", round(low, 1), " to ", round(high, 1), ")"),
                    p_val = round(p, 4)))
}

pair1_css <- calc_rmst_pair_css(df_loc, "Deep & Slow", "Deep & Fast")
pair2_css <- calc_rmst_pair_css(df_loc, "Deep & Slow", "Shallow & Fast")
pair3_css <- calc_rmst_pair_css(df_loc, "Deep & Slow", "Shallow & Slow")

rmst_css_table <- bind_rows(pair1_css, pair2_css, pair3_css)
print(rmst_css_table)


cat("\n====================================================================\n")
cat("3. CUMULATIVE CANCER MORTALITY BY YEAR (Proportions)\n")
cat("====================================================================\n")
# Calculate exact percentages (Cumulative Incidence)
cuminc_css <- tidycmprsk::cuminc(Surv(os_months_15y, fg_status_15y) ~ phenotype, data = df_loc)

cif_results_tidy <- cuminc_css %>%
  tidy(times = c(12, 36, 60, 120, 180)) %>%
  filter(outcome == "1_cancer") %>%
  mutate(
    Mortality_Percent = paste0(round(estimate * 100, 1), "%"),
    CI_Percent = paste0("(", round(conf.low * 100, 1), "% - ", round(conf.high * 100, 1), "%)")
  ) %>%
  select(strata, time, Mortality_Percent, CI_Percent)

print(as.data.frame(cif_results_tidy))

# ============================================================================
# CIF PLOT (Journal quality)
# ============================================================================
my_cols_fg <- c("Deep & Slow" = "#4DAF4A",
                "Deep & Fast" = "#377EB8",
                "Shallow & Fast" = "#E41A1C",
                "Shallow & Slow" = "#FF7F00")

p_css <- ggcuminc(cuminc_css, outcome = "1_cancer", linewidth = 1.2) +
  scale_color_manual(values = my_cols_fg) +
  scale_x_continuous(breaks = seq(0, 180, by = 36), limits = c(0, 180)) +

  # Add risk table (highly appreciated in Q1)
  add_risktable(risktable_stats = "n.risk", size = 3) +

  # Add P-value to the plot
  annotate("text", x = 5, y = 0.25, label = p_label, size = 4.5, fontface = "italic", hjust = 0) +

  labs(
    title = "Cancer-Specific Mortality by Phenotype",
    subtitle = "Accounting for competing risks of death from other causes",
    x = "Months after Radiotherapy",
    y = "Cumulative Incidence of Cancer Death"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_css)



# ====================================================================================================
# THREE PANEL FIGURE
# ====================================================================================================
library(ggplot2)
library(ggpubr)
library(survminer)
library(ggsurvfit)
library(tidycmprsk)

# 1. COLORS (Fixes the grey line problem)
# Named vector for Fine-Gray
my_cols_fg <- c("Deep & Slow" = "#4DAF4A",
                "Deep & Fast" = "#377EB8",
                "Shallow & Fast" = "#E41A1C",
                "Shallow & Slow" = "#FF7F00")
# Unnamed vector for KM (so survminer doesn't get confused by names)
my_cols_km <- unname(my_cols_fg)
line_size <- 1.0

# --- PLOT A: OS (Overall survival) ---
fit_os <- survfit(Surv(os_months_15y, os_event_15y) ~ phenotype, data = df_loc)

p_os <- ggsurvplot(fit_os, data = df_loc,
                   censor = FALSE, pval = TRUE,
                   palette = my_cols_km, size = line_size,
                   xlim = c(0, 180), break.time.by = 36,
                   legend.title = "Nadir Phenotype",
                   legend.labs = c("Deep & Slow", "Deep & Fast", "Shallow & Fast", "Shallow & Slow"),

                   ggtheme = theme_bw())$plot +
  labs(title = "A. Overall Survival", y = "Survival Probability")


# --- PLOT B: BCR (Recurrence‑free survival) ---
fit_bcr <- survfit(Surv(bcr_months_15y, bcr_event_15y) ~ phenotype, data = df_loc)

p_bcr <- ggsurvplot(fit_bcr, data = df_loc,
                    censor = FALSE, pval = TRUE,
                    pval.coord = c(5, 0.55),
                    palette = my_cols_km, size = line_size,
                    xlim = c(0, 180), break.time.by = 36,
                    legend.title = "Nadir Phenotype",
                    legend.labs = c("Deep & Slow", "Deep & Fast", "Shallow & Fast", "Shallow & Slow"),

                    ggtheme = theme_bw())$plot +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  labs(title = "B. Recurrence-free Survival", y = "Survival Probability")


# --- PLOT C: CSS (Cancer‑specific mortality / Fine-Gray) ---
cuminc_css <- tidycmprsk::cuminc(Surv(os_months_15y, fg_status_15y) ~ phenotype, data = df_loc)
p_val_raw <- cuminc_css$cmprsk$Tests[1, "pv"]
p_label <- ifelse(p_val_raw < 0.0001, "p < 0.0001", paste("p =", round(p_val_raw, 4)))

p_css <- ggcuminc(cuminc_css, outcome = "1_cancer", linewidth = line_size) +
  scale_color_manual(name = "Nadir Phenotype", values = my_cols_fg) +
  scale_x_continuous(breaks = seq(0, 180, by = 36), limits = c(0, 180)) +
  annotate("text", x = 5, y = 0.35, label = p_label, size = 5, hjust = 0) +
  labs(title = "C. Cancer-Specific Mortality", y = "Cumulative Incidence of Death") +
  theme_bw()


# --- COMBINE USING ggarrange ---
final_panel <- ggarrange(p_os, p_bcr, p_css,
                         ncol = 3, nrow = 1,
                         common.legend = TRUE,
                         legend = "bottom")

# Print
print(final_panel)

# Save (you can change width and height to make the image wider or taller)
ggsave("Figure2_Survival_Final.pdf", final_panel, width = 16, height = 5)
