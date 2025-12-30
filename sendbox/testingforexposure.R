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
model_PM25_final

model_PM25_final1 <- feglm(
  AggregatedDeath ~ mean_PM25_ensemble + mean_PM25_ensemble^2 + mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)


model_PM25_final1


model_PM25_final2 <- feglm(
  AggregatedDeath ~ ns(mean_PM25_ensemble, 2) + mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

model <- model_PM25_final2


# Function to calculate QAIC for quasi-poisson models
calculate_QAIC <- function(model) {
  deviance <- model$deviance
  phi <- summary(model)$dispersion
  k <- length(coef(model))
  QAIC <- deviance / phi + 2 * k
  return(QAIC)
}

# Calculate QAIC for each model
QAIC_1 <- calculate_QAIC(model_PM25_final)
QAIC_2 <- calculate_QAIC(model_PM25_final1)
QAIC_3 <- calculate_QAIC(model_PM25_final2)

# Summary of model comparisons
cat("\n=== Model Comparison Summary (QAIC) ===\n")
cat("Model 1 (Linear): QAIC =", round(QAIC_1, 2), "\n")
cat("Model 2 (Quadratic): QAIC =", round(QAIC_2, 2), "\n")
cat("Model 3 (Spline ns=2): QAIC =", round(QAIC_3, 2), "\n")
cat("\nModel with lowest QAIC is preferred (lower is better)\n")



