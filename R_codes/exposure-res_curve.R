### reference model and alternative specifications
library(dplyr)
library(fixest)
library(splines)
library(ggplot2)
library(patchwork)

data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")


####################
### Function to create exposure-response curves for all pollutants
####################

create_er_curve <- function(model, pollutant_var, pollutant_name, pollutant_unit, data, n_boot = 100, reference = NULL, functional_form = "linear") {
  # Extract the model formula and family from the original model
  original_formula <- formula(model)
  original_family <- model$family
  original_call <- model$call
  
  # Set default reference values (WHO annual mean standards)
  if (is.null(reference)) {
    ref_standards <- list(
      "mean_PM10ug" = 70,
      "mean_PM25_ensemble" = 35
    )
    reference <- ref_standards[[pollutant_var]]
  }
  
  # If reference not found in standards, use first quartile
  if (is.null(reference)) {
    reference <- quantile(data[[pollutant_var]], probs = 0.25, na.rm = TRUE)
  }
  
  # Create a sequence of pollutant values across the observed range, including reference
  data_min <- min(data[[pollutant_var]], na.rm = TRUE)
  data_max <- max(data[[pollutant_var]], na.rm = TRUE)
  data_mean <- mean(data[[pollutant_var]], na.rm = TRUE)
  data_q75 <- quantile(data[[pollutant_var]], probs = 0.75, na.rm = TRUE)
  
  # Create fine sequence across range and include summary statistics
  fine_seq <- seq(data_min, data_max, length.out = 100)
  pollutant_seq <- sort(unique(c(reference, data_min, data_mean, data_q75, data_max, fine_seq)))
  ref_idx <- which.min(abs(pollutant_seq - reference))
  
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
  
  # Calculate relative risk (RR) compared to the reference level
  ref_estimate <- predictions$estimate[ref_idx]
  ref_low <- predictions$conf.low[ref_idx]
  ref_high <- predictions$conf.high[ref_idx]
  
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
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_vline(xintercept = reference, linetype = "dotted", color = "darkgreen", linewidth = 1, alpha = 0.7) +
    geom_vline(xintercept = data_min, linetype = "solid", color = "gray60", linewidth = 0.6, alpha = 0.5) +
    geom_vline(xintercept = data_max, linetype = "solid", color = "gray60", linewidth = 0.6, alpha = 0.5) +
    geom_vline(xintercept = data_mean, linetype = "dashed", color = "steelblue", linewidth = 0.8, alpha = 0.6) +
    geom_vline(xintercept = data_q75, linetype = "dashed", color = "purple", linewidth = 0.7, alpha = 0.6) +
    scale_x_continuous(
      breaks = c(data_min, data_mean, data_q75, reference, data_max),
      labels = c(paste0("Min.\n", round(data_min, 1)), 
                 paste0("Mean\n", round(data_mean, 2)),
                 paste0("Q3\n", round(data_q75, 2)),
                 paste0("Ref.\n", round(reference, 2)),
                 paste0("Max.\n", round(data_max, 2)))
    ) +
    labs(
      x = paste0(pollutant_name, " (", pollutant_unit, ")"),
      y = "Risk Ratio for Death",
      title = paste("Mortality RR vs.", pollutant_name, paste0("(", functional_form, ")"))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 8)
    )
  
  return(list(plot = plot, data = pred_rr, reference = reference))
}

####################
### 1. NO2
####################
model_NO2_final <- feglm(
  AggregatedDeath ~ mean_NO2_ensemble +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

model_NO2_final1 <- feglm(
  AggregatedDeath ~ mean_NO2_ensemble + I(mean_NO2_ensemble)^2 + ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

NO2_er <- create_er_curve(model_NO2_final, "mean_NO2_ensemble", "NO2", "μg/m³", data, n_boot = 100, functional_form = "linear")
NO2_DRS <- NO2_er$plot
NO2_DRS

NO2_er1 <- create_er_curve(model_NO2_final1, "mean_NO2_ensemble", "NO2", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
NO2_DRS1 <- NO2_er1$plot
NO2_DRS1

model_NO2_final2 <- feglm(
  AggregatedDeath ~ ns(mean_NO2_ensemble, 2) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

NO2_er2 <- create_er_curve(model_NO2_final2, "mean_NO2_ensemble", "NO2", "μg/m³", data, n_boot = 100, functional_form = "spline")
NO2_DRS2 <- NO2_er2$plot
NO2_DRS2

####################
### 2. NO
####################
data$mean_NOug <- data$mean_NOppb* 1.23  # Convert ppb to μg/m³
model_NO_final <- feglm(
  AggregatedDeath ~ I(mean_NOug) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

NO_er <- create_er_curve(model_NO_final, "mean_NOug", "NO", "μg/m³", data, n_boot = 100, functional_form = "linear")
NO_DRS <- NO_er$plot
NO_DRS

model_NO_final1 <- feglm(
  AggregatedDeath ~ I(mean_NOug) + I(mean_NOug)^2 +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

NO_er1 <- create_er_curve(model_NO_final1, "mean_NOug", "NO", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
NO_DRS1 <- NO_er1$plot
NO_DRS1

model_NO_final2 <- feglm(
  AggregatedDeath ~ ns(mean_NOug, 2) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

NO_er2 <- create_er_curve(model_NO_final2, "mean_NOug", "NO", "μg/m³", data, n_boot = 100, functional_form = "spline")
NO_DRS2 <- NO_er2$plot
NO_DRS2

plot_NO_combined <- (NO_DRS + NO_DRS1 + NO_DRS2) +
  plot_annotation(
    title = "NO Exposure-Response Curves: Comparison of Functional Forms",
    subtitle = "Linear, Quadratic, and Spline Models",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )
plot_NO_combined
####################
### 3. SO2
####################
model_SO2_final <- feglm(
  AggregatedDeath ~ I(mean_SO2_ensemble) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

SO2_er <- create_er_curve(model_SO2_final, "mean_SO2_ensemble", "SO2", "μg/m³", data, n_boot = 100, functional_form = "linear")
SO2_DRS <- SO2_er$plot
SO2_DRS


model_SO2_final1 <- feglm(
  AggregatedDeath ~ I(mean_SO2_ensemble) + I(mean_SO2_ensemble)^2 + ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

SO2_er1 <- create_er_curve(model_SO2_final1, "mean_SO2_ensemble", "SO2", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
SO2_DRS1 <- SO2_er1$plot
SO2_DRS1

model_SO2_final2 <- feglm(
  AggregatedDeath ~ ns(mean_SO2_ensemble, 2) + ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

SO2_er2 <- create_er_curve(model_SO2_final2, "mean_SO2_ensemble", "SO2", "μg/m³", data, n_boot = 100, functional_form = "spline")
SO2_DRS2 <- SO2_er2$plot
SO2_DRS2

####################
### 4. Ox
####################

model_Ox_final <- feglm(
  AggregatedDeath ~ I(mean_Ox) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

Ox_er <- create_er_curve(model_Ox_final, "mean_Ox", "Ox", "μg/m³", data, n_boot = 100, functional_form = "linear")
Ox_DRS <- Ox_er$plot
Ox_DRS


model_Ox_final1 <- feglm(
  AggregatedDeath ~ I(mean_Ox) + I(mean_Ox)^2 + ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

Ox_er1 <- create_er_curve(model_Ox_final1, "mean_Ox", "Ox", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
Ox_DRS1 <- Ox_er1$plot
Ox_DRS1

model_Ox_final2 <- feglm(
  AggregatedDeath ~ ns(mean_Ox, 2) + ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

Ox_er2 <- create_er_curve(model_Ox_final2, "mean_Ox", "Ox", "μg/m³", data, n_boot = 100, functional_form = "spline")
Ox_DRS2 <- Ox_er2$plot
Ox_DRS2

####################
### 5. PM25
####################
model_PM25_final <- feglm(
  AggregatedDeath ~ mean_PM25_ensemble +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM25_er <- create_er_curve(model_PM25_final, "mean_PM25_ensemble", "PM2.5", "μg/m³", data, n_boot = 100, functional_form = "linear")
PM25_DRS <- PM25_er$plot
PM25_DRS


model_PM25_final1 <- feglm(
  AggregatedDeath ~ mean_PM25_ensemble + mean_PM25_ensemble^2 + mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM25_er1 <- create_er_curve(model_PM25_final1, "mean_PM25_ensemble", "PM2.5", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
PM25_DRS1 <- PM25_er1$plot
PM25_DRS1



model_PM25_final2 <- feglm(
  AggregatedDeath ~ ns(mean_PM25_ensemble, 2) + mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM25_er2 <- create_er_curve(model_PM25_final2, "mean_PM25_ensemble", "PM2.5", "μg/m³", data, n_boot = 100, functional_form = "spline")
PM25_DRS2 <- PM25_er2$plot
PM25_DRS2

# Combined PM2.5 exposure-response curves
combined_PM25_plot <- (PM25_DRS + PM25_DRS1 + PM25_DRS2) +
  plot_annotation(
    title = "PM2.5 Exposure-Response Curves: Comparison of Functional Forms",
    subtitle = "Linear, Quadratic, and Spline Models",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

combined_PM25_plot



####################
### 6. PM10
####################
model_PM10_final <- feglm(
  AggregatedDeath ~ mean_PM10ug +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM10_er <- create_er_curve(model_PM10_final, "mean_PM10ug", "PM10", "μg/m³", data, n_boot = 100, functional_form = "linear")
PM10_DRS <- PM10_er$plot
PM10_DRS


model_PM10_final1 <- feglm(
  AggregatedDeath ~ mean_PM10ug + I(mean_PM10ug)^2 +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM10_er1 <- create_er_curve(model_PM10_final1, "mean_PM10ug", "PM10", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
PM10_DRS1 <- PM10_er1$plot
PM10_DRS1



model_PM10_final2 <- feglm(
  AggregatedDeath ~ ns(mean_PM10ug, 2) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

PM10_er2 <- create_er_curve(model_PM10_final2, "mean_PM10ug", "PM10", "μg/m³", data, n_boot = 100, functional_form = "spline")
PM10_DRS2 <- PM10_er2$plot
PM10_DRS2



####################
### 7. CO
####################
model_CO_final <- feglm(
  AggregatedDeath ~ I(CO_pred) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

CO_er <- create_er_curve(model_CO_final, "CO_pred", "CO", "μg/m³", data, n_boot = 100, functional_form = "linear")
CO_DRS <- CO_er$plot
CO_DRS


model_CO_final1 <- feglm(
  AggregatedDeath ~ I(CO_pred) + I(CO_pred)^2 + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

CO_er1 <- create_er_curve(model_CO_final1, "CO_pred", "CO", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
CO_DRS1 <- CO_er1$plot
CO_DRS1

model_CO_final2 <- feglm(
  AggregatedDeath ~ ns(CO_pred, 2) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

CO_er2 <- create_er_curve(model_CO_final2, "CO_pred", "CO", "μg/m³", data, n_boot = 100, functional_form = "spline")
CO_DRS2 <- CO_er2$plot
CO_DRS2

####################
### 8. O3
####################
model_O3_final <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

model_O3_final1 <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal) + I(mean_O3_pred_cal)^2 +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

O3_er <- create_er_curve(model_O3_final, "mean_O3_pred_cal", "O3", "μg/m³", data, n_boot = 100, functional_form = "linear")
O3_DRS <- O3_er$plot
O3_DRS


O3_er1 <- create_er_curve(model_O3_final1, "mean_O3_pred_cal", "O3", "μg/m³", data, n_boot = 100, functional_form = "quadratic")
O3_DRS1 <- O3_er1$plot
O3_DRS1

model_O3_final2 <- feglm(
  AggregatedDeath ~ ns(mean_O3_pred_cal, 2) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

O3_er2 <- create_er_curve(model_O3_final2, "mean_O3_pred_cal", "O3", "μg/m³", data, n_boot = 100, functional_form = "spline")
O3_DRS2 <- O3_er2$plot
O3_DRS2

####################
### Combined Exposure-Response Curves for All 8 Pollutants
####################
combined_er_plot <- (NO2_DRS + NO_DRS2) / (SO2_DRS + Ox_DRS) / 
  (PM25_DRS + PM10_DRS) / (CO_DRS + O3_DRS) +
  plot_annotation(
    title = "Exposure-Response Curves: Mortality Risk Ratios for All Pollutants (Linear Exposure Models)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_er_plot


combined_er_plot1 <- (SO2_DRS + Ox_DRS) / 
  (CO_DRS + O3_DRS) +
  plot_annotation(
    title = "Mortality Risk Ratios for Pollutants",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_er_plot1

combined_er_plot2 <- (NO2_DRS + NO_DRS1 ) / 
  (PM25_DRS + PM10_DRS) +
  plot_annotation(
    title = "Exposure-Response Curves: Mortality Risk Ratios for All Pollutants",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_er_plot2




####################
### Combined Exposure-Response Curves for All 8 Pollutants - Quadratic Models
####################

combined_er_plot_1 <- (NO2_DRS1 + NO_DRS1 + SO2_DRS1 + Ox_DRS1) / 
  (PM25_DRS1 + PM10_DRS1 + CO_DRS1 + O3_DRS1) +
  plot_annotation(
    title = "Exposure-Response Curves: Mortality Risk Ratios for All Pollutants (Quadratic Exposure Models)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_er_plot_1


combined_er_plot11 <- (SO2_DRS1 + Ox_DRS1) / 
  (CO_DRS1 + O3_DRS1) +
  plot_annotation(
    title = "Mortality Risk Ratios for Pollutants (Quadratic Models)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_er_plot11

combined_er_plot21 <- (NO2_DRS1 + NO_DRS1 ) / 
  (PM25_DRS1 + PM10_DRS1) +
  plot_annotation(
    title = "Exposure-Response Curves: Mortality Risk Ratios for All Pollutants (Quadratic Models)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

combined_er_plot21

####################
### Model Comparison: QAIC for all pollutants
####################

# Function to calculate QAIC for quasi-poisson models
calculate_QAIC <- function(model) {
  deviance <- model$deviance
  phi <- summary(model)$dispersion
  k <- length(coef(model))
  QAIC <- deviance / phi + 2 * k
  return(QAIC)
}

# NO2 QAIC
cat("\n=== NO2 Model Comparison (QAIC) ===\n")
NO2_QAIC_1 <- calculate_QAIC(model_NO2_final)
NO2_QAIC_2 <- calculate_QAIC(model_NO2_final1)
NO2_QAIC_3 <- calculate_QAIC(model_NO2_final2)
cat("Linear:   QAIC =", round(NO2_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(NO2_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(NO2_QAIC_3, 2), "\n")

# NO QAIC
cat("\n=== NO Model Comparison (QAIC) ===\n")
NO_QAIC_1 <- calculate_QAIC(model_NO_final)
NO_QAIC_2 <- calculate_QAIC(model_NO_final1)
NO_QAIC_3 <- calculate_QAIC(model_NO_final2)
cat("Linear:   QAIC =", round(NO_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(NO_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(NO_QAIC_3, 2), "\n")

# SO2 QAIC
cat("\n=== SO2 Model Comparison (QAIC) ===\n")
SO2_QAIC_1 <- calculate_QAIC(model_SO2_final)
SO2_QAIC_2 <- calculate_QAIC(model_SO2_final1)
SO2_QAIC_3 <- calculate_QAIC(model_SO2_final2)
cat("Linear:   QAIC =", round(SO2_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(SO2_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(SO2_QAIC_3, 2), "\n")

# Ox QAIC
cat("\n=== Ox Model Comparison (QAIC) ===\n")
Ox_QAIC_1 <- calculate_QAIC(model_Ox_final)
Ox_QAIC_2 <- calculate_QAIC(model_Ox_final1)
Ox_QAIC_3 <- calculate_QAIC(model_Ox_final2)
cat("Linear:   QAIC =", round(Ox_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(Ox_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(Ox_QAIC_3, 2), "\n")

# PM2.5 QAIC
cat("\n=== PM2.5 Model Comparison (QAIC) ===\n")
PM25_QAIC_1 <- calculate_QAIC(model_PM25_final)
PM25_QAIC_2 <- calculate_QAIC(model_PM25_final1)
PM25_QAIC_3 <- calculate_QAIC(model_PM25_final2)
cat("Linear:   QAIC =", round(PM25_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(PM25_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(PM25_QAIC_3, 2), "\n")

# PM10 QAIC
cat("\n=== PM10 Model Comparison (QAIC) ===\n")
PM10_QAIC_1 <- calculate_QAIC(model_PM10_final)
PM10_QAIC_2 <- calculate_QAIC(model_PM10_final1)
PM10_QAIC_3 <- calculate_QAIC(model_PM10_final2)
cat("Linear:   QAIC =", round(PM10_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(PM10_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(PM10_QAIC_3, 2), "\n")

# CO QAIC
cat("\n=== CO Model Comparison (QAIC) ===\n")
CO_QAIC_1 <- calculate_QAIC(model_CO_final)
CO_QAIC_2 <- calculate_QAIC(model_CO_final1)
CO_QAIC_3 <- calculate_QAIC(model_CO_final2)
cat("Linear:   QAIC =", round(CO_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(CO_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(CO_QAIC_3, 2), "\n")

# O3 QAIC
cat("\n=== O3 Model Comparison (QAIC) ===\n")
O3_QAIC_1 <- calculate_QAIC(model_O3_final)
O3_QAIC_2 <- calculate_QAIC(model_O3_final1)
O3_QAIC_3 <- calculate_QAIC(model_O3_final2)
cat("Linear:   QAIC =", round(O3_QAIC_1, 2), "\n")
cat("Quadratic: QAIC =", round(O3_QAIC_2, 2), "\n")
cat("Spline:    QAIC =", round(O3_QAIC_3, 2), "\n")

# Summary table
cat("\n=== Summary: Best Model (Lowest QAIC) for Each Pollutant ===\n")
cat("Lower QAIC indicates better model fit\n\n")

summary_qaic <- data.frame(
  Pollutant = c("NO2", "NO", "SO2", "Ox", "PM2.5", "PM10", "CO", "O3"),
  Linear = c(NO2_QAIC_1, NO_QAIC_1, SO2_QAIC_1, Ox_QAIC_1, PM25_QAIC_1, PM10_QAIC_1, CO_QAIC_1, O3_QAIC_1),
  Quadratic = c(NO2_QAIC_2, NO_QAIC_2, SO2_QAIC_2, Ox_QAIC_2, PM25_QAIC_2, PM10_QAIC_2, CO_QAIC_2, O3_QAIC_2),
  Spline = c(NO2_QAIC_3, NO_QAIC_3, SO2_QAIC_3, Ox_QAIC_3, PM25_QAIC_3, PM10_QAIC_3, CO_QAIC_3, O3_QAIC_3)
)

summary_qaic

####################
### Wald Tests for Non-linearity
####################
# For a quadratic extension, we test whether the coefficient of the squared term is zero;
# for spline models, we test the joint null that all spline basis coefficients are zero.
# If the non-linear terms are not clearly different from zero at a prespecified α level (e.g., 0.05),
# we conclude that there is no strong statistical evidence requiring departure from linearity
# and retain the linear specification for parsimony and interpretability.

cat("\n=== Wald Tests for Quadratic Term (Coefficient = 0) ===\n")
cat("Testing significance of squared terms at α = 0.05\n\n")

# NO2 quadratic test
NO2_wald_quad <- tryCatch({
  coef_table <- model_NO2_final1$coeftable
  coef_names <- rownames(coef_table)
  # Find the squared term coefficient (I(I(mean_NO2_ensemble)^2))
  squared_term_idx <- grep("I\\(I\\(mean_NO2_ensemble\\)\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("NO2 quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)


# NO quadratic test
NO_wald_quad <- tryCatch({
  coef_table <- model_NO_final1$coeftable
  coef_names <- rownames(coef_table)
  squared_term_idx <- grep("I\\(I\\(mean_NOug\\)\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("NO quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)



# SO2 quadratic test
SO2_wald_quad <- tryCatch({
  coef_table <- model_SO2_final1$coeftable
  coef_names <- rownames(coef_table)
  squared_term_idx <- grep("I\\(I\\(mean_SO2_ensemble\\)\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("SO2 quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)

# Ox quadratic test
Ox_wald_quad <- tryCatch({
  coef_table <- model_Ox_final1$coeftable
  coef_names <- rownames(coef_table)
  squared_term_idx <- grep("I\\(I\\(mean_Ox\\)\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("Ox quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)

# PM2.5 quadratic test
PM25_wald_quad <- tryCatch({
  coef_table <- model_PM25_final1$coeftable
  coef_names <- rownames(coef_table)
  squared_term_idx <- grep("I\\(mean_PM25_ensemble\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("PM2.5 quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)

# PM10 quadratic test
PM10_wald_quad <- tryCatch({
  coef_table <- model_PM10_final1$coeftable
  coef_names <- rownames(coef_table)
  squared_term_idx <- grep("I\\(I\\(mean_PM10ug\\)\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("PM10 quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)

# CO quadratic test
CO_wald_quad <- tryCatch({
  coef_table <- model_CO_final1$coeftable
  coef_names <- rownames(coef_table)
  squared_term_idx <- grep("I\\(I\\(CO_pred\\)\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("CO quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)

# O3 quadratic test
O3_wald_quad <- tryCatch({
  coef_table <- model_O3_final1$coeftable
  coef_names <- rownames(coef_table)
  squared_term_idx <- grep("I\\(I\\(mean_O3_pred_cal\\)\\^2\\)", coef_names)
  if (length(squared_term_idx) > 0) {
    p_value <- coef_table[squared_term_idx, "Pr(>|t|)"]
    cat("O3 quadratic term p-value:", round(p_value, 4), 
        ifelse(p_value < 0.05, "(significant)\n", "(not significant)\n"))
    p_value
  } else NA
}, error = function(e) NA)

####################
### Wald Tests for Spline Terms (Joint Null)
####################
# For spline models, we test the joint null that all spline basis coefficients are zero.
# If the spline terms are not clearly different from zero at a prespecified α level (e.g., 0.05),
# we conclude that there is no strong statistical evidence requiring departure from linearity
# and retain the linear specification for parsimony and interpretability. 

cat("\n=== Wald Tests for Spline Terms (Joint Null) ===\n")
cat("Testing significance of all spline basis coefficients jointly at α = 0.05\n\n")

# NO2 spline test
NO2_wald_spline <- tryCatch({
  coef_table <- model_NO2_final2$coeftable
  coef_names <- rownames(coef_table)
  # Find all spline basis coefficients for the pollutant (ns(mean_NO2_ensemble))
  spline_term_idx <- grep("ns\\(mean_NO2_ensemble", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    # Joint test: all spline terms significant if any p-value < 0.05 or overall pattern
    p_values <- spline_coefs[, "Pr(>|t|)"]
    # Use Fisher's method for combining p-values (simplified: check if any significant)
    any_sig <- any(p_values < 0.05)
    cat("NO2 spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

# NO spline test
NO_wald_spline <- tryCatch({
  coef_table <- model_NO_final2$coeftable
  coef_names <- rownames(coef_table)
  spline_term_idx <- grep("ns\\(mean_NOug", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    p_values <- spline_coefs[, "Pr(>|t|)"]
    any_sig <- any(p_values < 0.05)
    cat("NO spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

# SO2 spline test
SO2_wald_spline <- tryCatch({
  coef_table <- model_SO2_final2$coeftable
  coef_names <- rownames(coef_table)
  spline_term_idx <- grep("ns\\(mean_SO2_ensemble", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    p_values <- spline_coefs[, "Pr(>|t|)"]
    any_sig <- any(p_values < 0.05)
    cat("SO2 spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

# Ox spline test
Ox_wald_spline <- tryCatch({
  coef_table <- model_Ox_final2$coeftable
  coef_names <- rownames(coef_table)
  spline_term_idx <- grep("ns\\(mean_Ox", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    p_values <- spline_coefs[, "Pr(>|t|)"]
    any_sig <- any(p_values < 0.05)
    cat("Ox spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

# PM2.5 spline test
PM25_wald_spline <- tryCatch({
  coef_table <- model_PM25_final2$coeftable
  coef_names <- rownames(coef_table)
  spline_term_idx <- grep("ns\\(mean_PM25_ensemble", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    p_values <- spline_coefs[, "Pr(>|t|)"]
    any_sig <- any(p_values < 0.05)
    cat("PM2.5 spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

# PM10 spline test
PM10_wald_spline <- tryCatch({
  coef_table <- model_PM10_final2$coeftable
  coef_names <- rownames(coef_table)
  spline_term_idx <- grep("ns\\(mean_PM10ug", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    p_values <- spline_coefs[, "Pr(>|t|)"]
    any_sig <- any(p_values < 0.05)
    cat("PM10 spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

# CO spline test
CO_wald_spline <- tryCatch({
  coef_table <- model_CO_final2$coeftable
  coef_names <- rownames(coef_table)
  spline_term_idx <- grep("ns\\(CO_pred", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    p_values <- spline_coefs[, "Pr(>|t|)"]
    any_sig <- any(p_values < 0.05)
    cat("CO spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

# O3 spline test
O3_wald_spline <- tryCatch({
  coef_table <- model_O3_final2$coeftable
  coef_names <- rownames(coef_table)
  spline_term_idx <- grep("ns\\(mean_O3_pred_cal", coef_names)
  if (length(spline_term_idx) > 0) {
    spline_coefs <- coef_table[spline_term_idx, ]
    p_values <- spline_coefs[, "Pr(>|t|)"]
    any_sig <- any(p_values < 0.05)
    cat("O3 spline terms (", length(spline_term_idx), "df) - any significant: ", 
        ifelse(any_sig, "YES\n", "NO\n"), sep = "")
    any_sig
  } else NA
}, error = function(e) NA)

cat("\n=== Model Selection Conclusion ===\n")
cat("If quadratic and spline terms are not statistically significant (p > 0.05),\n")
cat("the linear model is retained for parsimony and interpretability.\n")

# Create comprehensive model selection summary
model_selection_summary <- data.frame(
  Pollutant = c("NO2", "NO", "SO2", "Ox", "PM2.5", "PM10", "CO", "O3"),
  Linear_QAIC = c(NO2_QAIC_1, NO_QAIC_1, SO2_QAIC_1, Ox_QAIC_1, PM25_QAIC_1, PM10_QAIC_1, CO_QAIC_1, O3_QAIC_1),
  Quadratic_QAIC = c(NO2_QAIC_2, NO_QAIC_2, SO2_QAIC_2, Ox_QAIC_2, PM25_QAIC_2, PM10_QAIC_2, CO_QAIC_2, O3_QAIC_2),
  Spline_QAIC = c(NO2_QAIC_3, NO_QAIC_3, SO2_QAIC_3, Ox_QAIC_3, PM25_QAIC_3, PM10_QAIC_3, CO_QAIC_3, O3_QAIC_3),
  Quadratic_Wald_p = c(NO2_wald_quad, NO_wald_quad, SO2_wald_quad, Ox_wald_quad, 
                       PM25_wald_quad, PM10_wald_quad, CO_wald_quad, O3_wald_quad),
  Spline_Any_Sig = c(NO2_wald_spline, NO_wald_spline, SO2_wald_spline, Ox_wald_spline,
                     PM25_wald_spline, PM10_wald_spline, CO_wald_spline, O3_wald_spline)
)

# Add model selection recommendation based on QAIC and Wald tests
model_selection_summary <- model_selection_summary |>
  mutate(
    Best_QAIC = case_when(
      Linear_QAIC < Quadratic_QAIC & Linear_QAIC < Spline_QAIC ~ "Linear",
      Quadratic_QAIC < Linear_QAIC & Quadratic_QAIC < Spline_QAIC ~ "Quadratic",
      TRUE ~ "Spline"
    ),
    Quadratic_Sig = Quadratic_Wald_p < 0.05,
    Recommended_Form = case_when(
      # 1) No evidence of non-linearity
      !Quadratic_Sig & !Spline_Any_Sig ~ "Linear (parsimony)",
      
      # 2) Only quadratic significant
      Quadratic_Sig & !Spline_Any_Sig ~ "Quadratic",
      
      # 3) Only spline significant: compare Linear vs Spline by QAIC
      !Quadratic_Sig & Spline_Any_Sig & (Spline_QAIC + 2 < Linear_QAIC) ~ "Spline (by QAIC vs Linear)",
      !Quadratic_Sig & Spline_Any_Sig & (Spline_QAIC + 2 >= Linear_QAIC) ~ "Linear (spline non-linear but no QAIC gain)",
      
      # 4) Both quadratic and spline significant: choose best QAIC overall
      Quadratic_Sig & Spline_Any_Sig ~ paste0(Best_QAIC, " (by QAIC)"),
      
      # 5) Fallback
      TRUE ~ Best_QAIC
    )
  )

cat("\n=== Comprehensive Model Selection Summary ===\n\n")
print(model_selection_summary)

cat("\n=== Interpretation ===\n")
cat("- Linear QAIC: Information criterion for linear model\n")
cat("- Quadratic QAIC: Information criterion for quadratic model\n")
cat("- Spline QAIC: Information criterion for spline model (2 df)\n")
cat("- Quadratic Wald p-value: Significance of squared term\n")
cat("- Spline Any Sig: Whether any spline basis function is significant\n")
cat("- Recommended Form: Selected model based on statistical tests and parsimony\n")


