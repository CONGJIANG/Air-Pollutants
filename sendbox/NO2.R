### reference model and alternative specifications
library(dplyr)
library(fixest)
library(splines)
library(ggplot2)
library(patchwork)

data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")
model_NO2_final <- feglm(
  AggregatedDeath ~ I(mean_NO2_ensemble/5) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)


summary(model_NO2_final) 
# I(mean_NO2_ensemble/5)           0.110414   0.028574  3.864131 1.1302e-04 ***
(NO2_obs_coef <- coef(model_NO2_final)["I(mean_NO2_ensemble/5)"])
(NO2_L <- NO2_obs_coef - 1.96*model_NO2_final$se["I(mean_NO2_ensemble/5)"])
(NO2_U <- NO2_obs_coef + 1.96*model_NO2_final$se["I(mean_NO2_ensemble/5)"])


mean(predict(model_NO2_final, newdata = data, type = "response"))

PM25_seq <- seq(min(data$mean_PM25_ensemble, na.rm = TRUE), 
               max(data$mean_PM25_ensemble, na.rm = TRUE), 
               length.out = 100)


# Exposure-response curve for NO2 with bootstrap CIs
# Create a sequence of NO2 values across the observed range
NO2_seq <- seq(min(data$mean_NO2_ensemble, na.rm = TRUE), 
               max(data$mean_NO2_ensemble, na.rm = TRUE), 
               length.out = 100)

# Get average predictions across NO2 range
data1 <- data
data1$mean_NO2_ensemble <- as.numeric(data$mean_NO2_ensemble[1])


# Bootstrap for confidence intervals on data1 predictions
n_boot <- 100

# Pre-allocate storage for all bootstrap estimates
all_boot_preds <- array(NA, dim = c(length(NO2_seq), nrow(data1), n_boot))

# Single bootstrap loop - refit model once per iteration
for (i in 1:n_boot) {
  # Resample original data with replacement
  boot_idx <- sample(1:nrow(data), nrow(data), replace = TRUE)
  data_boot <- data[boot_idx, ]
  
  # Refit model on bootstrap sample
  model_boot <- tryCatch({
    feglm(
      AggregatedDeath ~ I(mean_NO2_ensemble/5) + ns(mean_t2m_c, df=4) + ns(mean_RH, df=5) + ns(mean_Precipitation, df=2) | District^month,
      data = data_boot,
      offset = ~log(Population_2022),
      vcov = "iid",
      family = "quasipoisson"
    )
  }, error = function(e) NULL)
  
  if (!is.null(model_boot)) {
    # Create prediction data for all NO2 values at once
    pred_data <- data1[rep(1, length(NO2_seq)), ]
    pred_data$mean_NO2_ensemble <- NO2_seq
    
    # Get all predictions for this bootstrap iteration
    all_boot_preds[, 1, i] <- predict(model_boot, newdata = pred_data, type = "response")
  }
}

# Calculate mean and CIs from bootstrap - vectorized
predictions <- data.frame(
  mean_NO2_ensemble = NO2_seq,
  estimate = apply(all_boot_preds[, 1, ], 1, mean, na.rm = TRUE),
  conf.low = 2*apply(all_boot_preds[, 1, ], 1, mean, na.rm = TRUE) - apply(all_boot_preds[, 1, ], 1, quantile, probs = 0.025, na.rm = TRUE)/5,
  conf.high = 2*apply(all_boot_preds[, 1, ], 1, mean, na.rm = TRUE) - apply(all_boot_preds[, 1, ], 1, quantile, probs = 0.975, na.rm = TRUE)/5
)

# Calculate relative risk (RR) compared to the minimum NO2 level
ref_estimate <- predictions$estimate[1]
ref_low <- predictions$conf.low[1]
ref_high <- predictions$conf.high[1]
  
pred_rr <- predictions |>
  mutate(
    estimate = estimate / ref_estimate,
    conf.low = conf.low / ref_high,
    conf.high = conf.high / ref_low
  )



# Plot Risk Ratios
NO2_DRS <- ggplot(pred_rr, aes(x = mean_NO2_ensemble, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = expression(NO[2] ~ (μg/m^3)),
    y = "Risk Ratio for Death (95% CI)",
    title = paste("Mortality Risk Ratio vs. NO2 Exposure (reference at min NO2 level =", 
                  round(NO2_seq[1], 2), "μg/m³)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.grid = element_blank()
  )

NO2_DRS



####################
### Function to create exposure-response curves for all pollutants
####################

create_er_curve <- function(model, pollutant_var, pollutant_name, pollutant_unit, data, n_boot = 100) {
  # Extract the model formula and family from the original model
  original_formula <- formula(model)
  original_family <- model$family
  original_call <- model$call
  
  # Create a sequence of pollutant values across the observed range
  pollutant_seq <- seq(min(data[[pollutant_var]], na.rm = TRUE), 
                       max(data[[pollutant_var]], na.rm = TRUE), 
                       length.out = 100)
  
  # Set up prediction data template
  pred_data_template <- data[1, ]
  
  # Pre-allocate storage for all bootstrap estimates
  all_boot_preds <- array(NA, dim = c(length(pollutant_seq), nrow(data), n_boot))
  
  # Bootstrap loop
  for (i in 1:n_boot) {
    # Resample original data with replacement
    boot_idx <- sample(1:nrow(data), nrow(data), replace = TRUE)
    data_boot <- data[boot_idx, ]
    
    # Refit model on bootstrap sample using original formula
    model_boot <- tryCatch({
      feglm(
        original_formula,
        data = data_boot,
        offset = ~log(Population_2022),
        vcov = "iid",
        family = original_family
      )
    }, error = function(e) NULL)
    
    if (!is.null(model_boot)) {
      # Create prediction data for all pollutant values at once
      pred_data <- pred_data_template[rep(1, length(pollutant_seq)), ]
      pred_data[[pollutant_var]] <- pollutant_seq
      
      # Get all predictions for this bootstrap iteration
      all_boot_preds[, 1, i] <- predict(model_boot, newdata = pred_data, type = "response")
    }
  }
  
  # Calculate mean and CIs from bootstrap - vectorized
  predictions <- data.frame(
    pollutant_val = pollutant_seq,
    estimate = apply(all_boot_preds[, 1, ], 1, mean, na.rm = TRUE),
    conf.high = 2*apply(all_boot_preds[, 1, ], 1, mean, na.rm = TRUE) - apply(all_boot_preds[, 1, ], 1, quantile, probs = 0.025, na.rm = TRUE)/5,
    conf.low = 2*apply(all_boot_preds[, 1, ], 1, mean, na.rm = TRUE) - apply(all_boot_preds[, 1, ], 1, quantile, probs = 0.975, na.rm = TRUE)/5
  )
  
  # Calculate relative risk (RR) compared to the minimum pollutant level
  ref_estimate <- predictions$estimate[10]
  ref_low <- predictions$conf.low[10]
  ref_high <- predictions$conf.high[10]
  
  pred_rr <- predictions |>
    mutate(
      estimate = estimate / ref_estimate,
      conf.low = conf.low / ref_high,
      conf.high = conf.high / ref_low
    )
  
  # Create plot
  plot <- ggplot(pred_rr, aes(x = pollutant_val, y = estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                alpha = 0.2, fill = "steelblue") +
    geom_line(color = "steelblue", linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(
      x = paste0(pollutant_name, " (", pollutant_unit, ")"),
      y = "Risk Ratio for Death (95% CI)",
      title = paste("Mortality Risk Ratio vs.", pollutant_name, "Exposure (reference at", 
                    round(pollutant_seq[10], 2), pollutant_unit, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      panel.grid = element_blank()
    )
  
  return(list(plot = plot, data = pred_rr))
}

####################
### 1. NO2
####################
model_NO2_final <- feglm(
  AggregatedDeath ~ I(mean_NO2_ensemble/5) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)


NO2_er <- create_er_curve(model_NO2_final, "mean_NO2_ensemble", "NO2", "μg/m³", data, n_boot = 100)
NO2_DRS <- NO2_er$plot
NO2_DRS
####################
### 2. NO
####################
model_NO_final <- feglm(
  AggregatedDeath ~ I(mean_NOppb/5) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

NO_er <- create_er_curve(model_NO_final, "mean_NOppb", "NO", "ppb", data, n_boot = 100)
NO_DRS <- NO_er$plot
NO_DRS
NO_er$data

####################
### 3. SO2
####################
model_SO2_final <- feglm(
  AggregatedDeath ~ I(mean_SO2_ensemble/5) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

SO2_er <- create_er_curve(model_SO2_final, "mean_SO2_ensemble", "SO2", "μg/m³", data, n_boot = 100)
SO2_DRS <- SO2_er$plot
SO2_DRS


####################
### 4. Ox
####################

model_Ox_final <- feglm(
  AggregatedDeath ~ I(mean_Ox/5) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

Ox_er <- create_er_curve(model_Ox_final, "mean_Ox", "Ox", "ppb", data, n_boot = 100)
Ox_DRS <- Ox_er$plot
Ox_DRS


####################
### 5. PM25
####################
model_PM25_final <- feglm(
  AggregatedDeath ~ I(mean_PM25_ensemble/5) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
hist(predict(model_PM25_final))
PM25_er <- create_er_curve(model_PM25_final, "mean_PM25_ensemble", "PM2.5", "μg/m³", data, n_boot = 500)
PM25_DRS <- PM25_er$plot
PM25_DRS


model_PM25_final2 <- feglm(
  AggregatedDeath ~ ns(mean_PM25_ensemble/5, df = 2) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM25_er <- create_er_curve(model_PM25_final2, "mean_PM25_ensemble", "PM2.5", "μg/m³", data, n_boot = 500)
PM25_DRS <- PM25_er$plot
PM25_DRS

####################
### 6. PM10
####################
model_PM10_final <- feglm(
  AggregatedDeath ~ I(mean_PM10ug/5) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM10_er <- create_er_curve(model_PM10_final, "mean_PM10ug", "PM10", "μg/m³", data, n_boot = 100)
PM10_DRS <- PM10_er$plot
PM10_DRS


####################
### 7. CO
####################
model_CO_final <- feglm(
  AggregatedDeath ~ I(CO_pred/5) + PM25_pred +mean_PM10_pred+ mean_NO2_pred + ns(mean_t2m_c, df=2) +ns(mean_RH, df=1)+ ns(mean_Precipitation, df=3) +ns(mean_wind_speed, df=6) | District^year + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

CO_er <- create_er_curve(model_CO_final, "CO_pred", "CO", "ppm", data, n_boot = 100)
CO_DRS <- CO_er$plot
CO_DRS


####################
### 8. O3
####################
model_O3_final <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal/5) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

model_O3_final <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal/5) + I(mean_O3_pred_cal/5)^2 +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

O3_er <- create_er_curve(model_O3_final, "mean_O3_pred_cal", "O3", "ppb", data, n_boot = 100)
O3_DRS <- O3_er$plot
O3_DRS

####################
### Combined Exposure-Response Curves for All 8 Pollutants
####################

combined_er_plot <- (NO2_DRS + NO_DRS + SO2_DRS + Ox_DRS) / 
                     (PM25_DRS + PM10_DRS + CO_DRS + O3_DRS) +
  plot_annotation(
    title = "Exposure-Response Curves: Mortality Risk Ratios for All Pollutants",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_er_plot
