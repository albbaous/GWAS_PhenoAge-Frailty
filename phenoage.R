# Install packages if not already installed
install.packages("tidyverse")
install.packages("tibble")
install.packages("readr")

# Load required libraries
library(tidyverse)
library(tibble)
library(readr)
library(mice)
library(jsonlite)
library(tidyr)

setwd('/Users/user/Desktop/Biostatistics2')

# Read in dataset
df <- read_csv("cohort_data_phenoage.csv")

# Define mapping of UKB fields -> human-readable labels
biomarkers <- tribble(
  ~column,                    ~label,
  "participant.eid",          "eid",
  "participant.p21003_i0",    "Age",
  "participant.p31",          "Sex",
  "participant.p22009_a1",    "PC1",
  "participant.p22009_a2",    "PC2",
  "participant.p22009_a3",    "PC3",
  "participant.p22009_a4",    "PC4",
  "participant.p22009_a5",    "PC5",
  "participant.p22009_a6",    "PC6",
  "participant.p22009_a7",    "PC7",
  "participant.p22009_a8",    "PC8",
  "participant.p22009_a9",    "PC9",
  "participant.p22009_a10",   "PC10",
  # PhenoAge-specific biomarkers
  "participant.p30600_i0",    "Albumin",
  "participant.p30700_i0",    "Creatinine",
  "participant.p30740_i0",    "Glucose",
  "participant.p30710_i0",    "CRP",
  "participant.p30180_i0",    "Lymphocyte_percent",
  "participant.p30044_i0",    "Mean_cell_volume",
  "participant.p30070_i0",    "Red_cell_distribution_width",
  "participant.p30610_i0",    "Alkaline_phosphatase",
  "participant.p30000_i0",    "White_blood_cell_count"
)

# Rename columns using mapping
df <- df %>%
  rename_with(~ biomarkers$label[match(.x, biomarkers$column)],
              .cols = intersect(names(df), biomarkers$column))

# Inspect renamed dataset
glimpse(df)

# -------------------------------
# Initial sample size
# -------------------------------
cat("Number of individuals before all filtering:", nrow(df), "\n")

# -------------------------------
# Filter: remove rows with NA in any of the 10 PC columns
# -------------------------------
df <- df %>%
  filter(
    !is.na(PC1) & !is.na(PC2) & !is.na(PC3) & !is.na(PC4) &
      !is.na(PC5) & !is.na(PC6) & !is.na(PC7) & !is.na(PC8) &
      !is.na(PC9) & !is.na(PC10)
  )

# -------------------------------
# Define PhenoAge biomarker columns
# -------------------------------
phenoage_markers <- c(
  "Albumin",
  "Creatinine",
  "Glucose",
  "CRP",
  "Lymphocyte_percent",
  "Mean_cell_volume",
  "Red_cell_distribution_width",
  "Alkaline_phosphatase",
  "White_blood_cell_count",
  "Age" 
)

# -------------------------------
# Filter: keep rows with at least 1 biomarker present and not all but one missing
# -------------------------------
df <- df %>%
  filter(
    rowSums(is.na(df[phenoage_markers])) < length(phenoage_markers) &
      rowSums(is.na(df[setdiff(phenoage_markers, "Age")])) < (length(phenoage_markers) - 1)
  )

# -------------------------------
# Final sample size
# -------------------------------
cat("Number of individuals after all filtering:", nrow(df), "\n")

# -------------------------------
# Missing summary for PhenoAge biomarkers
# -------------------------------
missing_summary_post_filter <- df %>%
  summarise(across(
    all_of(phenoage_markers),
    list(
      missing_pct = ~ mean(is.na(.)) * 100,
      missing_count = ~ sum(is.na(.))
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Biomarker", "Metric"),
    names_pattern = "^(.*)_(missing_pct|missing_count)$"
  ) %>%
  pivot_wider(names_from = Metric, values_from = value) %>%
  arrange(desc(missing_pct))

# Print updated missing summary
print(missing_summary_post_filter)

# -------------------------------
# Imputation including PCs
# -------------------------------

# Define columns for imputation: Age, Sex, PhenoAge biomarkers
imputation_cols <- c("Age", "Sex", phenoage_markers)

# Include the PCs for retention (no imputation needed, but keep them)
pcs <- paste0("PC", 1:10)

# Prepare data for imputation: only impute biomarkers + Age/Sex
df_imputation <- df %>%
  select(all_of(imputation_cols))

# Create predictor matrix (1 = predictor, 0 = not predictor)
predictor_matrix <- matrix(1, nrow = ncol(df_imputation), ncol = ncol(df_imputation))
colnames(predictor_matrix) <- colnames(df_imputation)
rownames(predictor_matrix) <- colnames(df_imputation)

# Exclude biomarkers as predictors for each other (set within-block to 0)
predictor_matrix[which(colnames(predictor_matrix) %in% phenoage_markers),
                 which(colnames(predictor_matrix) %in% phenoage_markers)] <- 0

# Perform imputation using PMM (fast, robust)
imputed_data <- mice(
  df_imputation,
  method = "pmm",
  m = 5,
  predictorMatrix = predictor_matrix,
  seed = 123
)

# Extract completed dataset
df_completed_imputed <- complete(imputed_data, 1)

# Keep eid and PCs at the front, followed by imputed biomarkers + Age/Sex
df_completed <- df_completed_imputed %>%
  select(everything()) %>%            # this keeps imputed columns first
  bind_cols(df %>% select(eid, all_of(pcs))) %>%   # temporarily bind eid+PCs
  relocate(eid, all_of(pcs), .before = 1)          # move eid+PCs to front

# Inspect completed dataset
glimpse(df_completed)

# -------------------------------
# PhenoAge calculation
# -------------------------------

# Gompertz hazard slope
gamma <- 0.0076927

# Linear predictor weights for the 9 biomarkers + chronological age
# Replace x_i with the column names in your dataset
weights <- c(
  "Albumin" = -0.0336,
  "Creatinine" = 0.0095,
  "Glucose" = 0.1953,
  "CRP" = 0.0954,
  "Lymphocyte_percent" = -0.0120,
  "Mean_cell_volume" = 0.0268,
  "Red_cell_distribution_width" = 0.3306,
  "Alkaline_phosphatase" = 0.0019,
  "White_blood_cell_count" = 0.0554,
  "Age" = 0.0804  # chronological age
)

# Linear predictor constant
xb_constant <- -19.9067

# Function to compute mortality risk and PhenoAge for a single row
compute_phenoage <- function(row, weights, xb_constant, gamma) {
  
  # Compute linear predictor xb
  xb <- xb_constant + sum(row[names(weights)] * weights)
  
  # Mortality risk formula
  mortality_risk <- 1 - exp(-exp(xb) * ((exp(120 * gamma) - 1) / gamma))
  
  # PhenoAge formula
  phenoage <- 141.50225 + log(-0.00553 * log(1 - mortality_risk)) / 0.090165
  
  return(phenoage)
}

# Apply the function row-wise to your dataframe
df_completed$PhenoAge <- apply(df_completed[, names(weights)], 1, compute_phenoage,
                               weights = weights,
                               xb_constant = xb_constant,
                               gamma = gamma)

# Inspect PhenoAge
head(df_completed)

# -------------------------------
# Create PHENOTYPE file '.pheno'
# -------------------------------

df_phen <- df_completed %>%
  # Rename 'eid' to 'IID'
  rename(IID = eid) %>%
  # Add 'FID' column (same as IID)
  mutate(FID = IID) %>%
  # Select columns in order: FID, IID, PhenoAge
  select(FID, IID, PhenoAge)

# Inspect
head(df_phen)

# Save as tab-delimited .pheno file
write.table(
  df_phen,
  file = "ukb_phenotype_data.pheno",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# -------------------------------
# Create COVARIATE file '.cov'
# -------------------------------

df_cov <- df_completed %>%
  # Keep columns for covariates: eid, Age, Sex, PCs
  select(eid, Age, Sex, starts_with("PC")) %>%
  # Rename 'eid' to 'IID'
  rename(IID = eid) %>%
  # Add FID column (same as IID)
  mutate(FID = IID) %>%
  # Reorder columns: FID, IID, then all others
  select(FID, IID, everything())

# Inspect
head(df_cov)

# Save as tab-delimited .cov file
write.table(
  df_cov,
  file = "ukb_covariates.cov",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
