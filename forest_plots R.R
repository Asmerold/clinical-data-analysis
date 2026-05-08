# ===================================================================================================
# FOREST PLOTS
# ===================================================================================================
library(dplyr)
library(broom)
library(forestplot)

# --------------------------------------------------------------------------------------------------
# BCR
# --------------------------------------------------------------------------------------------------

# Extract data
res_df <- tidy(cox_bcr, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    var_name = case_when(
      term == "phenotypeDeep & Fast" ~ "Deep & Fast",
      term == "phenotypeShallow & Fast" ~ "Shallow & Fast",
      term == "phenotypeShallow & Slow" ~ "Shallow & Slow",
      term == "gleason_group7" ~ "Gleason 7",
      term == "gleason_group8-10" ~ "Gleason 8-10",
      term == "age_at_diag" ~ "Age at Diagnosis",
      term == "psa_baseline_static" ~ "Baseline PSA",
      term == "rt_year" ~ "Year of Radiotherapy",
      term == "adt_status1" ~ "ADT Status (Yes vs No)",
      term == "max_dose_prostate" ~ "Total Dose (Prostate)",
      TRUE ~ term
    ),
    hr_ci = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    p_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

# Create header and reference rows
header <- data.frame(var_name="Variable", hr_ci="Hazard Ratio (95% CI)", p_val="p-value", estimate=NA, conf.low=NA, conf.high=NA)

pheno_header <- data.frame(var_name="Nadir Phenotype", hr_ci="", p_val="", estimate=NA, conf.low=NA, conf.high=NA)
pheno_ref <- data.frame(var_name="Deep & Slow (Reference)", hr_ci="1.00", p_val="", estimate=1, conf.low=1, conf.high=1)

gleason_header <- data.frame(var_name="Gleason Score", hr_ci="", p_val="", estimate=NA, conf.low=NA, conf.high=NA)
gleason_ref <- data.frame(var_name="Gleason <= 6 (Reference)", hr_ci="1.00", p_val="", estimate=1, conf.low=1, conf.high=1)

# Assemble final table for plotting
plot_df <- bind_rows(
  header,
  pheno_header, pheno_ref, filter(res_df, grepl("phenotype", term)),
  gleason_header, gleason_ref, filter(res_df, grepl("gleason_cat", term)),
  filter(res_df, !grepl("phenotype|gleason_cat", term))
)

# Prepare text
tabletext <- cbind(plot_df$var_name, plot_df$hr_ci, plot_df$p_val)
# Indicate bold rows (headers)
is_summary_vec <- plot_df$var_name %in% c("Variable", "Nadir Phenotype", "Gleason Score")

# Draw Forest Plot (Cochrane Style)
forestplot(
  labeltext = tabletext,
  mean = plot_df$estimate,
  lower = plot_df$conf.low,
  upper = plot_df$conf.high,
  is.summary = is_summary_vec,
  xlog = TRUE,
  zero = 1,
  boxsize = 0.2,
  lineheight = unit(0.8, "cm"),
  xticks = c(0.25, 0.5, 1, 2, 4, 8),
  clip = c(0.2, 8.0),
  ci.vertices = TRUE,
  ci.vertices.height = 0.05,
  col = fpColors(box = "royalblue", line = "darkblue", zero = "gray50"),
  hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4)),
  txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", cex = 0.9),
                   ticks = gpar(fontfamily = "sans", cex = 0.8),
                   xlab  = gpar(fontfamily = "sans", cex = 0.9)),
  xlab = "Hazard Ratio (log scale)"
)

# --------------------------------------------------------------------------------------------------
# CSS
# --------------------------------------------------------------------------------------------------

# Extract data
res_df <- tidy(fg_multi, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    var_name = case_when(
      term == "phenotypeDeep & Fast" ~ "Deep & Fast",
      term == "phenotypeShallow & Fast" ~ "Shallow & Fast",
      term == "phenotypeShallow & Slow" ~ "Shallow & Slow",
      term == "gleason_group7" ~ "Gleason 7",
      term == "gleason_group8-10" ~ "Gleason 8-10",
      term == "age_at_diag" ~ "Age at Diagnosis",
      term == "psa_baseline_static" ~ "Baseline PSA",
      term == "rt_year" ~ "Year of Radiotherapy",
      term == "adt_status1" ~ "ADT Status (Yes vs No)",
      term == "max_dose_prostate" ~ "Total Dose (Prostate)",
      TRUE ~ term
    ),
    hr_ci = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    p_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

# Create header and reference rows
header <- data.frame(var_name="Variable", hr_ci="Hazard Ratio (95% CI)", p_val="p-value", estimate=NA, conf.low=NA, conf.high=NA)

pheno_header <- data.frame(var_name="Nadir Phenotype", hr_ci="", p_val="", estimate=NA, conf.low=NA, conf.high=NA)
pheno_ref <- data.frame(var_name="Deep & Slow (Reference)", hr_ci="1.00", p_val="", estimate=1, conf.low=1, conf.high=1)

gleason_header <- data.frame(var_name="Gleason Score", hr_ci="", p_val="", estimate=NA, conf.low=NA, conf.high=NA)
gleason_ref <- data.frame(var_name="Gleason <= 6 (Reference)", hr_ci="1.00", p_val="", estimate=1, conf.low=1, conf.high=1)

# Assemble final table for plotting
plot_df <- bind_rows(
  header,
  pheno_header, pheno_ref, filter(res_df, grepl("phenotype", term)),
  gleason_header, gleason_ref, filter(res_df, grepl("gleason_cat", term)),
  filter(res_df, !grepl("phenotype|gleason_cat", term))
)

# Prepare text
tabletext <- cbind(plot_df$var_name, plot_df$hr_ci, plot_df$p_val)
# Indicate bold rows (headers)
is_summary_vec <- plot_df$var_name %in% c("Variable", "Nadir Phenotype", "Gleason Score")

# Draw Forest Plot (Cochrane Style)
forestplot(
  labeltext = tabletext,
  mean = plot_df$estimate,
  lower = plot_df$conf.low,
  upper = plot_df$conf.high,
  is.summary = is_summary_vec,
  xlog = TRUE,
  zero = 1,
  boxsize = 0.2,
  lineheight = unit(0.8, "cm"),
  xticks = c(0.25, 0.5, 1, 2, 4, 8, 16),
  clip = c(0.2, 8.0),
  ci.vertices = TRUE,
  ci.vertices.height = 0.05,
  col = fpColors(box = "royalblue", line = "darkblue", zero = "gray50"),
  hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4)),
  txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", cex = 0.9),
                   ticks = gpar(fontfamily = "sans", cex = 0.8),
                   xlab  = gpar(fontfamily = "sans", cex = 0.9)),
  xlab = "Hazard Ratio (log scale)"
)

# --------------------------------------------------------------------------------------------------
# OS
# --------------------------------------------------------------------------------------------------

# Extract data
res_df <- tidy(cox_os, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    var_name = case_when(
      term == "phenotypeDeep & Fast" ~ "Deep & Fast",
      term == "phenotypeShallow & Fast" ~ "Shallow & Fast",
      term == "phenotypeShallow & Slow" ~ "Shallow & Slow",
      term == "gleason_group7" ~ "Gleason 7",
      term == "gleason_group8-10" ~ "Gleason 8-10",
      term == "age_at_diag" ~ "Age at Diagnosis",
      term == "psa_baseline_static" ~ "Baseline PSA",
      term == "rt_year" ~ "Year of Radiotherapy",
      term == "adt_status1" ~ "ADT Status (Yes vs No)",
      term == "max_dose_prostate" ~ "Total Dose (Prostate)",
      TRUE ~ term
    ),
    hr_ci = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    p_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )

# Create header and reference rows
header <- data.frame(var_name="Variable", hr_ci="Hazard Ratio (95% CI)", p_val="p-value", estimate=NA, conf.low=NA, conf.high=NA)

pheno_header <- data.frame(var_name="Nadir Phenotype", hr_ci="", p_val="", estimate=NA, conf.low=NA, conf.high=NA)
pheno_ref <- data.frame(var_name="Deep & Slow (Reference)", hr_ci="1.00", p_val="", estimate=1, conf.low=1, conf.high=1)

gleason_header <- data.frame(var_name="Gleason Score", hr_ci="", p_val="", estimate=NA, conf.low=NA, conf.high=NA)
gleason_ref <- data.frame(var_name="Gleason <= 6 (Reference)", hr_ci="1.00", p_val="", estimate=1, conf.low=1, conf.high=1)

# Assemble final table for plotting
plot_df <- bind_rows(
  header,
  pheno_header, pheno_ref, filter(res_df, grepl("phenotype", term)),
  gleason_header, gleason_ref, filter(res_df, grepl("gleason_cat", term)),
  filter(res_df, !grepl("phenotype|gleason_cat", term))
)

# Prepare text
tabletext <- cbind(plot_df$var_name, plot_df$hr_ci, plot_df$p_val)
# Indicate bold rows (headers)
is_summary_vec <- plot_df$var_name %in% c("Variable", "Nadir Phenotype", "Gleason Score")

# Draw Forest Plot (Cochrane Style)
forestplot(
  labeltext = tabletext,
  mean = plot_df$estimate,
  lower = plot_df$conf.low,
  upper = plot_df$conf.high,
  is.summary = is_summary_vec,
  xlog = TRUE,
  zero = 1,
  boxsize = 0.2,
  lineheight = unit(0.8, "cm"),
  xticks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32),
  clip = c(0.2, 8.0),
  ci.vertices = TRUE,
  ci.vertices.height = 0.05,
  col = fpColors(box = "royalblue", line = "darkblue", zero = "gray50"),
  hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:4)),
  txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", cex = 0.9),
                   ticks = gpar(fontfamily = "sans", cex = 0.8),
                   xlab  = gpar(fontfamily = "sans", cex = 0.9)),
  xlab = "Hazard Ratio (log scale)"
)

# --------------------------------------------------------------------------------------------------
# COMBINED PLOTS
# --------------------------------------------------------------------------------------------------

library(dplyr)
library(broom)
library(forestplot)

# 1. Collect data
get_model_data <- function(model, label) {
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(outcome = label)
}

res_bcr <- get_model_data(cox_bcr, "BCR")
res_css <- get_model_data(cox_css, "CSS")
res_os  <- get_model_data(cox_os, "OS")
all_res <- bind_rows(res_bcr, res_css, res_os)

# 2. Row names (Variable list)
row_names <- c(
  "Variable",
  "Nadir Phenotype",
  "   Deep & Slow (Reference)",
  "   Deep & Fast",
  "   Shallow & Fast",
  "   Shallow & Slow",
  "Gleason Score",
  "   Gleason \u2264 6 (Reference)",
  "   Gleason 7",
  "   Gleason 8-10",
  "Age at Diagnosis",
  "Baseline PSA",
  "Year of Radiotherapy",
  "Total Dose (Prostate)"
)

# 3. Initialize matrices
n_rows <- length(row_names)
means <- matrix(NA, nrow = n_rows, ncol = 3)
lows  <- matrix(NA, nrow = n_rows, ncol = 3)
highs <- matrix(NA, nrow = n_rows, ncol = 3)

# 4. Mapping
map_row <- c(
  "phenotypeDeep & Fast" = 4,
  "phenotypeShallow & Fast" = 5,
  "phenotypeShallow & Slow" = 6,
  "gleason_group7" = 9,
  "gleason_group8-10" = 10,
  "age_at_diag" = 11,
  "psa_baseline_static" = 12,
  "rt_year" = 13,
  "max_dose_prostate" = 14
)

# 5. Fill matrices
fill_matrices <- function(outcome_idx, outcome_label) {
  subset_data <- all_res %>% filter(outcome == outcome_label)
  means[3, outcome_idx] <<- 1.0; lows[3, outcome_idx] <<- 1.0; highs[3, outcome_idx] <<- 1.0
  means[8, outcome_idx] <<- 1.0; lows[8, outcome_idx] <<- 1.0; highs[8, outcome_idx] <<- 1.0
  for(term_name in names(map_row)) {
    val <- subset_data %>% filter(term == term_name)
    if(nrow(val) > 0) {
      row_idx <- map_row[term_name]
      means[row_idx, outcome_idx] <<- val$estimate
      lows[row_idx, outcome_idx]  <<- val$conf.low
      highs[row_idx, outcome_idx] <<- val$conf.high
    }
  }
}

fill_matrices(1, "BCR")
fill_matrices(2, "CSS")
fill_matrices(3, "OS")

# 6.  Set bold rows and lines (FIXED)
# The table has only 2 visual columns: Name and Plot. Therefore columns = 1:2.
is_summary_vec <- row_names %in% c("Variable", "Nadir Phenotype", "Gleason Score")

# 7. DRAW
forestplot(
  labeltext = row_names,
  mean = means,
  lower = lows,
  upper = highs,

  is.summary = is_summary_vec,
  xlog = TRUE,

  # Colors
  col = fpColors(box = c("#E41A1C", "#377EB8", "#4DAF4A"),
                 line = c("#E41A1C", "#377EB8", "#4DAF4A")),

  legend = c("BCR (Phoenix)", "CSS (Cancer-Specific)", "OS (Overall Survival)"),

  # Display parameters
  zero = 1,
  boxsize = 0.15,
  ci.vertices = TRUE,
  ci.vertices.height = 0.05,
  xticks = c(0.1, 0.5, 1, 2, 5, 10),
  clip = c(0.1, 10),

  # FIXED LINES (now columns = 1:2)
  hrzl_lines = list("2" = gpar(lwd = 1, columns = 1:2),
                    "7" = gpar(lwd = 0.5, lty = 2, columns = 1:2)),

  txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", cex = 0.8),
                   legend = gpar(fontfamily = "sans", cex = 0.8),
                   ticks = gpar(fontfamily = "sans", cex = 0.7)),

  xlab = "Hazard Ratio (95% CI) - Log Scale",
  title = "Impact of PSA Nadir Phenotypes on Outcomes"
)
