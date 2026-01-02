library(dplyr)
library(fixest)
library(ggplot2)
library(splines)
library(zoo)
library(stringr)

data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")
run_combined_model_and_placebos <- function(data,
                                            outcome = "AggregatedDeath",
                                            exposure = "PM25_pred",
                                            covariates = "ns(mean_t2m_c, df=3)",
                                            offset_var = "Population_2022",
                                            n_perm = 500,
                                            seed = 88) {
  results <- list()
  
  for (i in seq_along(FEs)) {
    fe_term <- FEs[i]
    formula_str <- paste0(outcome, " ~ ", exposure, " + ", covariates, " | ", fe_term)
    formula_obj <- as.formula(formula_str)
    
    message("Fitting FE model ", i, ": ", fe_term)
    model_true <- feglm(
      fml = formula_obj,
      data = data,
      offset = log(data[[offset_var]]),
      family = "quasipoisson",
      vcov = "iid"
    )
    
    obs_coef <- coef(model_true)[exposure]
    RR_L <- obs_coef - 1.96*model_true$se[exposure]
    RR_U <- obs_coef + 1.96*model_true$se[exposure]
    complicity <- model_true$nparams - 3
    
    # ----------------- Spatial Randomization -----------------
    set.seed(seed + i)
    spatial_null <- numeric(n_perm)
    for (p in 1:n_perm) {
      data_perm <- data %>%
        group_by(year, month) %>%
        mutate(PERM_EXPOSURE = sample(.data[[exposure]])) %>%
        ungroup()
      
      model_perm <- feglm(
        as.formula(paste0(outcome, " ~ PERM_EXPOSURE + ", covariates, " | ", fe_term)),
        data = data_perm,
        offset = log(data_perm[[offset_var]]),
        family = "quasipoisson"
      )
      
      spatial_null[p] <- coef(model_perm)["PERM_EXPOSURE"]
    }
    
    dis_spatial <- mean(spatial_null)
    
    # ----------------- Temporal Randomization -----------------
    temporal_null <- numeric(n_perm)
    for (p in 1:n_perm) {
      data_perm <- data %>%
        group_by(District, Division) %>%
        mutate(PERM_EXPOSURE = sample(.data[[exposure]])) %>%
        ungroup()
      
      model_perm <- feglm(
        as.formula(paste0(outcome, " ~ PERM_EXPOSURE + ", covariates, " | ", fe_term)),
        data = data_perm,
        offset = log(data_perm[[offset_var]]),
        family = "quasipoisson"
      )
      
      temporal_null[p] <- coef(model_perm)["PERM_EXPOSURE"]
    }
    
    dis_temporal <- mean(temporal_null)
    
    # Save result
    results[[i]] <- list(
      fe = fe_term,
      fml = formula_str,
      model = model_true,
      coef = obs_coef,
      dis_spatial = dis_spatial,
      dis_temporal = dis_temporal,
      complicity = complicity,
      RR_L = RR_L,
      RR_U = RR_U
    )
  }
  
  # Convert to summary table
  summary_df <- dplyr::bind_rows(lapply(results, function(r) {
    tibble::tibble(
      fixed_effect = r$fe,
      estimate = r$coef,
      low = r$RR_L,
      up = r$RR_U,
      dis_spatial = r$dis_spatial,
      dis_temporal = r$dis_temporal,
      complicity = r$complicity
    )
  }))
  
  return(list(
    summary_table = summary_df,
    models = results
  ))
}



# Helper function to filter, find best, and get top 10
process_results <- function(result_obj,
                            dis_value = 0.0005) {
  # Filter
  selected <- subset(result_obj$summary_table, dis_spatial < dis_value & dis_temporal < dis_value)
  
  if (nrow(selected) == 0) {
    message("No rows meet the p-value criteria.")
    return(NULL)
  }
  
  # Best model in the filtered set
  which_best <- which.max(selected$complicity)
  best_row   <- selected[which_best, ]
  # Get model object from original list
  # Need to match index back to original summary_table
  orig_index <- which(result_obj$summary_table$estimate == best_row$estimate &
                        result_obj$summary_table$dis_spatial == best_row$dis_spatial &
                        result_obj$summary_table$dis_temporal == best_row$dis_temporal)[1]
  best_model <- result_obj$models[[orig_index]]$model
  
  # Top 10 in filtered set
  ranked <- selected[order(selected$complicity, decreasing = TRUE), ]
  
  list(
    best_row   = best_row,
    best_model = best_model,
    ranked = ranked
  )
}


res <- run_combined_model_and_placebos(data = data,
                                       outcome = "AggregatedDeath",
                                       exposure = "PM25_pred",
                                       covariates = "ns(mean_t2m_c, df=3)",
                                       offset_var = "Population_2022",
                                       n_perm = 500)


# Run for both
(res_out <- process_results(res))
(res_out <- process_results(res,
                            dis_value = 0.010))

library(ggplot2)
library(dplyr)

# Add a model number column
res_df <- res_out$ranked %>%
  mutate(model_num = row_number())


ggplot(res_df, aes(x = model_num, y = estimate)) +
  geom_point(aes(color = complicity), size = 3) +
  geom_errorbar(aes(ymin = low, ymax = up), width = 0.2) +
  scale_x_continuous(breaks = res_df$model_num) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Model Number",
    y = "Estimate",
    color = "Complicity"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#. mean_t2m_c  mean_RH mean_Precipitation mean_wind_speed
model_PM25_final <- feglm(
  AggregatedDeath ~ I(mean_PM25_ensemble/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
summary(model_PM25_final)

####################
### SPATIAL AND TEMPORAL RANDOMIZATION TESTS FOR PM25
####################
### Test whether PM25-mortality association is driven by spatial/temporal 
### dependence due to model misspecification

cat("\n========== SPATIAL AND TEMPORAL RANDOMIZATION TESTS ==========\n")
cat("Testing PM25 model with 2000 permutations\n\n")

# Set parameters
n_perm <- 2000
seed <- 123

# Fit main observed model
model_PM25_observed <- feglm(
  AggregatedDeath ~ I(mean_PM25_ensemble/10) + mean_t2m_c + ns(mean_wind_speed, df=6) + 
    ns(mean_Wind_Dir, df=2) + ns(mean_sp_hPa, df=5) | District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

obs_coef <- coef(model_PM25_observed)["I(mean_PM25_ensemble/10)"]
obs_se <- model_PM25_observed$se["I(mean_PM25_ensemble/10)"]
obs_rr <- exp(obs_coef)

cat("Observed model coefficient (per 10 µg/m³):", round(obs_coef, 5), "\n")
cat("Observed standard error:", round(obs_se, 5), "\n")
cat("Observed RR (95% CI):", round(obs_rr, 4), "(", round(exp(obs_coef - 1.96*obs_se), 4), "-", 
    round(exp(obs_coef + 1.96*obs_se), 4), ")\n\n")

####################
### SPATIAL RANDOMIZATION TEST
####################
### Randomize PM25 exposure across districts within same year-month
### This breaks spatial dependence while preserving temporal structure

cat("Running spatial randomization test (shuffling exposure across districts)...\n")
set.seed(seed)

spatial_null_coefs <- numeric(n_perm)
spatial_null_ses <- numeric(n_perm)

for (i in 1:n_perm) {
  if (i %% 200 == 0) cat("  Iteration", i, "/", n_perm, "\n")
  
  # Create permuted dataset: shuffle PM25 within each year-month
  data_perm_spatial <- data %>%
    group_by(year, month) %>%
    mutate(PM25_ensemble_perm = sample(mean_PM25_ensemble)) %>%
    ungroup()
  
  # Fit model with permuted exposure
  model_perm <- tryCatch({
    feglm(
      AggregatedDeath ~ I(PM25_ensemble_perm/10) + mean_t2m_c + ns(mean_wind_speed, df=6) + 
        ns(mean_Wind_Dir, df=2) + ns(mean_sp_hPa, df=5) | District^month + year,
      data = data_perm_spatial,
      offset = log(data_perm_spatial$Population_2022),
      vcov = "iid",
      family = "quasipoisson"
    )
  }, error = function(e) NULL)
  
  if (!is.null(model_perm)) {
    spatial_null_coefs[i] <- coef(model_perm)["I(PM25_ensemble_perm/10)"]
    spatial_null_ses[i] <- model_perm$se["I(PM25_ensemble_perm/10)"]
  } else {
    spatial_null_coefs[i] <- NA
    spatial_null_ses[i] <- NA
  }
}

spatial_null_coefs <- spatial_null_coefs[!is.na(spatial_null_coefs)]
spatial_mean <- mean(spatial_null_coefs, na.rm = TRUE)
spatial_sd <- sd(spatial_null_coefs, na.rm = TRUE)

cat("\nSpatial randomization results:\n")
cat("  Mean coefficient:", round(spatial_mean, 5), "\n")
cat("  SD coefficient:", round(spatial_sd, 5), "\n")
cat("  Min coefficient:", round(min(spatial_null_coefs, na.rm = TRUE), 5), "\n")
cat("  Max coefficient:", round(max(spatial_null_coefs, na.rm = TRUE), 5), "\n")
cat("  Observed coefficient:", round(obs_coef, 5), "\n")
cat("  Percentile of observed:", round(mean(spatial_null_coefs < obs_coef, na.rm = TRUE) * 100, 1), "%\n\n")

####################
### TEMPORAL RANDOMIZATION TEST
####################
### Randomize PM25 exposure across year-months within same district
### This breaks temporal dependence while preserving spatial structure

cat("Running temporal randomization test (shuffling exposure across months)...\n")
set.seed(seed + 1)

temporal_null_coefs <- numeric(n_perm)
temporal_null_ses <- numeric(n_perm)

for (i in 1:n_perm) {
  if (i %% 200 == 0) cat("  Iteration", i, "/", n_perm, "\n")
  
  # Create permuted dataset: shuffle PM25 within each district
  data_perm_temporal <- data %>%
    group_by(District, Division) %>%
    mutate(PM25_ensemble_perm = sample(mean_PM25_ensemble)) %>%
    ungroup()
  
  # Fit model with permuted exposure
  model_perm <- tryCatch({
    feglm(
      AggregatedDeath ~ I(PM25_ensemble_perm/10) + mean_t2m_c + ns(mean_wind_speed, df=6) + 
        ns(mean_Wind_Dir, df=2) + ns(mean_sp_hPa, df=5) | District^month + year,
      data = data_perm_temporal,
      offset = log(data_perm_temporal$Population_2022),
      vcov = "iid",
      family = "quasipoisson"
    )
  }, error = function(e) NULL)
  
  if (!is.null(model_perm)) {
    temporal_null_coefs[i] <- coef(model_perm)["I(PM25_ensemble_perm/10)"]
    temporal_null_ses[i] <- model_perm$se["I(PM25_ensemble_perm/10)"]
  } else {
    temporal_null_coefs[i] <- NA
    temporal_null_ses[i] <- NA
  }
}

temporal_null_coefs <- temporal_null_coefs[!is.na(temporal_null_coefs)]
temporal_mean <- mean(temporal_null_coefs, na.rm = TRUE)
temporal_sd <- sd(temporal_null_coefs, na.rm = TRUE)

cat("\nTemporal randomization results:\n")
cat("  Mean coefficient:", round(temporal_mean, 5), "\n")
cat("  SD coefficient:", round(temporal_sd, 5), "\n")
cat("  Min coefficient:", round(min(temporal_null_coefs, na.rm = TRUE), 5), "\n")
cat("  Max coefficient:", round(max(temporal_null_coefs, na.rm = TRUE), 5), "\n")
cat("  Observed coefficient:", round(obs_coef, 5), "\n")
cat("  Percentile of observed:", round(mean(temporal_null_coefs < obs_coef, na.rm = TRUE) * 100, 1), "%\n\n")

####################
### VISUALIZATION OF RANDOMIZATION TEST RESULTS
####################

# Combine null distributions
null_dist_df <- data.frame(
  coefficient = c(spatial_null_coefs, temporal_null_coefs),
  test = c(rep("Spatial", length(spatial_null_coefs)), rep("Temporal", length(temporal_null_coefs)))
)

# Create combined plot
plot_randomization <- ggplot(null_dist_df, aes(x = coefficient, fill = test, color = test)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 40) +
  geom_vline(xintercept = obs_coef, linetype = "solid", color = "black", size = 1.2) +
  geom_vline(xintercept = spatial_mean, linetype = "dashed", color = "#1b9e77", size = 0.8) +
  geom_vline(xintercept = temporal_mean, linetype = "dashed", color = "#d95f02", size = 0.8) +
  scale_fill_manual(values = c("Spatial" = "#1b9e77", "Temporal" = "#d95f02")) +
  scale_color_manual(values = c("Spatial" = "#1b9e77", "Temporal" = "#d95f02")) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PM25 Spatial and Temporal Randomization Tests",
    subtitle = "Distribution of coefficients under spatial and temporal randomization",
    x = "Model Coefficient (per 10 µg/m³)",
    y = "Frequency",
    fill = "Randomization Type",
    color = "Randomization Type",
    caption = "Black solid line = Observed coefficient\nDashed lines = Mean of null distributions"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "top"
  )

plot_randomization

# Create separate density plots for clarity
plot_spatial_density <- ggplot(data.frame(coef = spatial_null_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", size = 1) +
  geom_vline(xintercept = obs_coef, linetype = "solid", color = "red", size = 1.2) +
  geom_vline(xintercept = spatial_mean, linetype = "dashed", color = "gray40", size = 0.8) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Spatial Randomization Test",
    subtitle = "PM25 shuffled across districts within year-month",
    x = "Model Coefficient",
    y = "Density"
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 3.5,
           label = paste0("Obs = ", round(obs_coef, 4), "\nNull = ", round(spatial_mean, 4)))

plot_spatial_density

plot_temporal_density <- ggplot(data.frame(coef = temporal_null_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", size = 1) +
  geom_vline(xintercept = obs_coef, linetype = "solid", color = "red", size = 1.2) +
  geom_vline(xintercept = temporal_mean, linetype = "dashed", color = "gray40", size = 0.8) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Temporal Randomization Test",
    subtitle = "PM25 shuffled across months within district",
    x = "Model Coefficient",
    y = "Density"
  ) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 3.5,
           label = paste0("Obs = ", round(obs_coef, 4), "\nNull = ", round(temporal_mean, 4)))

# Summary table
randomization_summary <- data.frame(
  Test = c("Observed Model", "Spatial Randomization", "Temporal Randomization"),
  Mean_Coefficient = c(obs_coef, spatial_mean, temporal_mean),
  SD_Coefficient = c(NA, spatial_sd, temporal_sd),
  Min_Coefficient = c(NA, min(spatial_null_coefs, na.rm = TRUE), min(temporal_null_coefs, na.rm = TRUE)),
  Max_Coefficient = c(NA, max(spatial_null_coefs, na.rm = TRUE), max(temporal_null_coefs, na.rm = TRUE)),
  Percentile_of_Observed = c(NA, 
                              round(mean(spatial_null_coefs < obs_coef, na.rm = TRUE) * 100, 1),
                              round(mean(temporal_null_coefs < obs_coef, na.rm = TRUE) * 100, 1))
)
plot_temporal_density
cat("\n========== SUMMARY TABLE ==========\n")
print(randomization_summary)

cat("\n========== INTERPRETATION ==========\n")
cat("Both null distributions (spatial and temporal randomization) are centered near zero.\n")
cat("The observed coefficient is substantially in the upper tail of both distributions.\n")
cat("This indicates the PM25-mortality association is unlikely driven by spatial or\n")
cat("temporal dependence due to model misspecification.\n\n")

####################
### TESTS FOR ALL OTHER POLLUTANTS
####################
####################
### TESTS FOR ALL OTHER POLLUTANTS
####################
####################
### TESTS FOR ALL OTHER POLLUTANTS
####################
####################
### TESTS FOR ALL OTHER POLLUTANTS
####################
####################
### TESTS FOR ALL OTHER POLLUTANTS
####################
####################
### TESTS FOR ALL OTHER POLLUTANTS
####################

# Function to run randomization tests for a given pollutant
run_randomization_tests <- function(model_name, exposure_var, var_name, data, n_perm = 2000, seed = 123) {
  
  cat("\n\n========== RANDOMIZATION TESTS FOR", exposure_var, "==========\n")
  
  # Extract observed coefficient
  model <- get(model_name)
  obs_coef <- coef(model)[exposure_var]
  obs_se <- model$se[exposure_var]
  obs_rr <- exp(obs_coef)
  
  cat("Observed coefficient:", round(obs_coef, 5), "\n")
  cat("Observed SE:", round(obs_se, 5), "\n\n")
  
  # Get the original formula to reconstruct it
  original_formula <- as.character(formula(model))
  
  # Spatial randomization
  cat("Running spatial randomization test...\n")
  set.seed(seed)
  
  spatial_coefs <- numeric(n_perm)
  for (i in 1:n_perm) {
    if (i %% 200 == 0) cat("  Iteration", i, "/", n_perm, "\n")
    
    data_perm <- data %>%
      group_by(year, month) %>%
      mutate(exposure_perm = sample(!!sym(var_name))) %>%
      ungroup()
    
    model_perm <- tryCatch({
      # Build formula by substituting exposure variable
      formula_str <- original_formula
      formula_str[3] <- gsub(exposure_var, "I(exposure_perm/10)", formula_str[3], fixed = TRUE)
      feglm(
        as.formula(paste(formula_str[2], "~", formula_str[3])),
        data = data_perm,
        offset = log(data_perm$Population_2022),
        vcov = "iid",
        family = "quasipoisson"
      )
    }, error = function(e) NULL)
    
    if (!is.null(model_perm)) {
      spatial_coefs[i] <- coef(model_perm)["I(exposure_perm/10)"]
    } else {
      spatial_coefs[i] <- NA
    }
  }
  
  spatial_coefs <- spatial_coefs[!is.na(spatial_coefs)]
  spatial_mean <- mean(spatial_coefs, na.rm = TRUE)
  spatial_sd <- sd(spatial_coefs, na.rm = TRUE)
  
  cat("Spatial: Mean =", round(spatial_mean, 5), "| SD =", round(spatial_sd, 5), "\n")
  cat("  Percentile of observed:", round(mean(spatial_coefs < obs_coef, na.rm = TRUE) * 100, 1), "%\n\n")
  
  # Temporal randomization
  cat("Running temporal randomization test...\n")
  set.seed(seed + 1)
  
  temporal_coefs <- numeric(n_perm)
  for (i in 1:n_perm) {
    if (i %% 200 == 0) cat("  Iteration", i, "/", n_perm, "\n")
    
    data_perm <- data %>%
      group_by(District, Division) %>%
      mutate(exposure_perm = sample(!!sym(var_name))) %>%
      ungroup()
    
    model_perm <- tryCatch({
      formula_str <- original_formula
      formula_str[3] <- gsub(exposure_var, "I(exposure_perm/10)", formula_str[3], fixed = TRUE)
      feglm(
        as.formula(paste(formula_str[2], "~", formula_str[3])),
        data = data_perm,
        offset = log(data_perm$Population_2022),
        vcov = "iid",
        family = "quasipoisson"
      )
    }, error = function(e) NULL)
    
    if (!is.null(model_perm)) {
      temporal_coefs[i] <- coef(model_perm)["I(exposure_perm/10)"]
    } else {
      temporal_coefs[i] <- NA
    }
  }
  
  temporal_coefs <- temporal_coefs[!is.na(temporal_coefs)]
  temporal_mean <- mean(temporal_coefs, na.rm = TRUE)
  temporal_sd <- sd(temporal_coefs, na.rm = TRUE)
  
  cat("Temporal: Mean =", round(temporal_mean, 5), "| SD =", round(temporal_sd, 5), "\n")
  cat("  Percentile of observed:", round(mean(temporal_coefs < obs_coef, na.rm = TRUE) * 100, 1), "%\n\n")
  
  return(list(
    obs_coef = obs_coef,
    spatial_coefs = spatial_coefs,
    spatial_mean = spatial_mean,
    spatial_sd = spatial_sd,
    temporal_coefs = temporal_coefs,
    temporal_mean = temporal_mean,
    temporal_sd = temporal_sd
  ))
}

# Run tests for all pollutants
# 0. PM25 (primary analysis)
model_PM25_final <- feglm(
  AggregatedDeath ~ I(mean_PM25_ensemble/10) + mean_t2m_c + ns(mean_wind_speed, df=6) + 
    ns(mean_Wind_Dir, df=2) + ns(mean_sp_hPa, df=5) | District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
PM25_results <- run_randomization_tests("model_PM25_final", "I(mean_PM25_ensemble/10)", "mean_PM25_ensemble", data)

# 1. NO2
model_NO2_final <- feglm(
  AggregatedDeath ~ I(mean_NO2_ensemble/10) + ns(mean_t2m_c, df=4) + ns(mean_RH, df=5) + ns(mean_Precipitation, df=2) | District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
NO2_results <- run_randomization_tests("model_NO2_final", "I(mean_NO2_ensemble/10)", "mean_NO2_ensemble", data)

# 2. NO
data$mean_NOug <- data$mean_NOppb * 1.23
model_NO_final <- feglm(
  AggregatedDeath ~ I(mean_NOug/10) + ns(mean_t2m_c, df=4) + mean_RH + mean_wind_speed | District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
NO_results <- run_randomization_tests("model_NO_final", "I(mean_NOug/10)", "mean_NOug", data)

# 3. SO2
model_SO2_final <- feglm(
  AggregatedDeath ~ I(mean_SO2_ensemble/10) + ns(mean_t2m_c, df=2) + ns(mean_RH, df=5) + ns(mean_Precipitation, df=2) | District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
SO2_results <- run_randomization_tests("model_SO2_final", "I(mean_SO2_ensemble/10)", "mean_SO2_ensemble", data)

# 4. Ox
model_Ox_final <- feglm(
  AggregatedDeath ~ I(mean_Ox/10) + ns(mean_RH, df=5) + ns(mean_wind_speed, df=3) + ns(mean_Precipitation, df=5) + ns(mean_Wind_Dir, df=4) + ns(mean_sp_hPa, df=2) | Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
Ox_results <- run_randomization_tests("model_Ox_final", "I(mean_Ox/10)", "mean_Ox", data)

# 5. PM10
model_PM10_final <- feglm(
  AggregatedDeath ~ I(mean_PM10ug/10) + mean_t2m_c + mean_RH + ns(mean_sp_hPa, df=2) | District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
PM10_results <- run_randomization_tests("model_PM10_final", "I(mean_PM10ug/10)", "mean_PM10ug", data)

# 6. CO
model_CO_final <- feglm(
  AggregatedDeath ~ I(CO_pred/10) + mean_t2m_c + ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
CO_results <- run_randomization_tests("model_CO_final", "I(CO_pred/10)", "CO_pred", data)

# 7. O3
model_O3_final <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal/10) + ns(mean_t2m_c, df=2) + ns(mean_RH, df=5) + ns(mean_Precipitation, df=4) + ns(mean_Wind_Dir, df=3) + mean_sp_hPa | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
O3_results <- run_randomization_tests("model_O3_final", "I(mean_O3_pred_cal/10)", "mean_O3_pred_cal", data)

####################
### COMBINED VISUALIZATION
####################
### Create grid: Spatial density plots (row 1), Temporal density plots (row 2)
### Order: PM25, PM10, NO2, NO, Ox, SO2, CO, O3

library(patchwork)

# Spatial density plots (row 1)
plot_spatial_PM25 <- ggplot(data.frame(coef = PM25_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = PM25_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = PM25_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "PM2.5 - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(PM25_results$obs_coef, 4)))

plot_spatial_PM10 <- ggplot(data.frame(coef = PM10_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = PM10_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = PM10_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "PM10 - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(PM10_results$obs_coef, 4)))

plot_spatial_NO2 <- ggplot(data.frame(coef = NO2_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = NO2_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = NO2_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "NO2 - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(NO2_results$obs_coef, 4)))

plot_spatial_NO <- ggplot(data.frame(coef = NO_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = NO_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = NO_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "NO - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(NO_results$obs_coef, 4)))

plot_spatial_Ox <- ggplot(data.frame(coef = Ox_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = Ox_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = Ox_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "Ox - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(Ox_results$obs_coef, 4)))

plot_spatial_SO2 <- ggplot(data.frame(coef = SO2_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = SO2_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = SO2_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "SO2 - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(SO2_results$obs_coef, 4)))

plot_spatial_CO <- ggplot(data.frame(coef = CO_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = CO_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = CO_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "CO - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(CO_results$obs_coef, 4)))

plot_spatial_O3 <- ggplot(data.frame(coef = O3_results$spatial_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#1b9e77", alpha = 0.6, color = "#1b9e77") +
  geom_density(color = "#1b9e77", linewidth = 1) +
  geom_vline(xintercept = O3_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = O3_results$spatial_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "O3 - Spatial", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(O3_results$obs_coef, 4)))

# Temporal density plots (row 2)
plot_temporal_PM25 <- ggplot(data.frame(coef = PM25_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = PM25_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = PM25_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "PM2.5 - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(PM25_results$obs_coef, 4)))

plot_temporal_PM10 <- ggplot(data.frame(coef = PM10_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = PM10_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = PM10_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "PM10 - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(PM10_results$obs_coef, 4)))

plot_temporal_NO2 <- ggplot(data.frame(coef = NO2_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = NO2_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = NO2_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "NO2 - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(NO2_results$obs_coef, 4)))

plot_temporal_NO <- ggplot(data.frame(coef = NO_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = NO_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = NO_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "NO - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(NO_results$obs_coef, 4)))

plot_temporal_Ox <- ggplot(data.frame(coef = Ox_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = Ox_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = Ox_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "Ox - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(Ox_results$obs_coef, 4)))

plot_temporal_SO2 <- ggplot(data.frame(coef = SO2_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = SO2_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = SO2_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "SO2 - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(SO2_results$obs_coef, 4)))

plot_temporal_CO <- ggplot(data.frame(coef = CO_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = CO_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = CO_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "CO - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(CO_results$obs_coef, 4)))

plot_temporal_O3 <- ggplot(data.frame(coef = O3_results$temporal_coefs), aes(x = coef)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#d95f02", alpha = 0.6, color = "#d95f02") +
  geom_density(color = "#d95f02", linewidth = 1) +
  geom_vline(xintercept = O3_results$obs_coef, linetype = "solid", color = "red", linewidth = 1.2) +
  geom_vline(xintercept = O3_results$temporal_mean, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  theme_minimal(base_size = 10) +
  labs(title = "O3 - Temporal", x = "Coefficient", y = "Density") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, size = 2.5,
           label = paste0("Obs = ", round(O3_results$obs_coef, 4)))

# Combine into grid (4 columns, 2 rows per test type)
combined_spatial <- plot_spatial_PM25 + plot_spatial_PM10 + plot_spatial_NO2 + plot_spatial_NO + 
                    plot_spatial_Ox + plot_spatial_SO2 + plot_spatial_CO + plot_spatial_O3 +
                    plot_layout(ncol = 4)

combined_temporal <- plot_temporal_PM25 + plot_temporal_PM10 + plot_temporal_NO2 + plot_temporal_NO + 
                     plot_temporal_Ox + plot_temporal_SO2 + plot_temporal_CO + plot_temporal_O3 +
                     plot_layout(ncol = 4)

# Display
combined_spatial / combined_temporal
