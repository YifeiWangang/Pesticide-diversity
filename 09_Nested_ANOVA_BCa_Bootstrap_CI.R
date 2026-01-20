# ==============================================================================
# Title: Advanced ARG Profiling (Nested ANOVA & BCa Bootstrap CI)
# Description: Statistical analysis of active ARG abundance (Box-Cox transformed) 
#              evaluating multi-pesticide interactions and nesting effects.
# ==============================================================================

# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(tidyverse)
library(boot)

# 2. Data Loading and Matrix Preprocessing
# ------------------------------------------------------------------------------
# - Columns 1-3: Metadata (Column 1: sample ID, Column 2: Diversity, Column 3: Mixture (specific pesticide combination))
# - Column 4: active ARG abundance (Box-Cox transformed)
# - Rows: Each row represents a unique Sample coupled with its Grouping info.
df <- read_delim(file.choose(), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# 3. Data Cleaning and Binary Feature Engineering
# ------------------------------------------------------------------------------
# Extract dummy variables (A, B, C, D) from the 'Mixture' column 
# to allow the linear model to identify interaction terms.
stat <- df %>%
  mutate(
    A = ifelse(str_detect(Mixture, "A"), 1, 0),
    B = ifelse(str_detect(Mixture, "B"), 1, 0),
    C = ifelse(str_detect(Mixture, "C"), 1, 0),
    D = ifelse(str_detect(Mixture, "D"), 1, 0),
    Diversity = factor(Diversity)
  )

# 4. Nested ANOVA (Diversity & Mixture nesting)
# ------------------------------------------------------------------------------
# Evaluate the significance of diversity gradients and specific 
# pesticide compositions nested within those levels.
mod_aov <- aov(ARGs_bc ~ Diversity + Mixture %in% Diversity, data = stat)

cat("--- Nested ANOVA Summary ---\n")
print(summary(mod_aov))

# 5. Bootstrap Coefficient Analysis (10,000 Iterations)
# ------------------------------------------------------------------------------

# A. Define the fitting function: Full interaction model to capture synergy and antagonism
boot_fn <- function(data, indices) {
  d <- data[indices, ]
  fit <- lm(ARGs_bc ~ A * B * C * D, data = d)
  return(coef(fit))
}

# B. Run Bootstrap (BCa method is more accurate for skewed data distributions)
set.seed(123)
cat("\nRunning BCa Bootstrap Analysis (R=10000)...\n")
boot_results <- boot(stat, statistic = boot_fn, R = 10000)

# C. Extract Confidence Intervals (Preserving logic with 'try' to prevent interruption)
# 
boot_ci_output <- t(sapply(1:length(boot_results$t0), function(i) {
  
  # Prioritize the Bias-Corrected and Accelerated (BCa) method
  ci <- try(boot.ci(boot_results, type = "bca", index = i), silent = TRUE)
  
  # If BCa fails (e.g., due to zero variance in resamples), fallback to Percentile method
  if (inherits(ci, "try-error") || is.null(ci)) {
    ci <- try(boot.ci(boot_results, type = "perc", index = i), silent = TRUE)
  }
  
  # Structure the output for Lower and Upper bounds
  if (inherits(ci, "try-error")) {
    return(c(NA, NA))
  } else if (!is.null(ci$bca)) {
    return(ci$bca[4:5]) # Return BCa values
  } else {
    return(ci$percent[4:5]) # Return Percentile values
  }
}))

# 6. Final Result Compilation
# ------------------------------------------------------------------------------
# Create the dataframe using standard names first to avoid reference errors
final_results_table <- data.frame(
  Term      = names(boot_results$t0),
  Estimate  = as.numeric(boot_results$t0),
  LCI       = boot_ci_output[, 1],
  UCI       = boot_ci_output[, 2]
)

# Apply significance classification
final_results_table <- final_results_table %>%
  mutate(
    Significance = case_when(
      is.na(LCI) ~ "Incalculable",
      LCI > 0    ~ "Significant Positive",
      UCI < 0    ~ "Significant Negative",
      TRUE       ~ "Not Significant"
    )
  )

# Rename columns for formal publication formatting
colnames(final_results_table) <- c("Term", "Beta (Estimate)", "95% CI Lower", "95% CI Upper", "Significance")

cat("\n--- Compiled Statistical Results (A*B*C*D Interaction) ---\n")
print(final_results_table, digits = 4)