library(dplyr)
library(fixest)
library(ggplot2)
library(splines)
library(zoo)
library(stringr)

data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")

# Standardize district names
district_mapping <- c(
  "Barishal" = "Barisal",
  "Bogura" = "Bogra",
  "Brahmanbaria" = "Brahamanbaria",
  "Chattogram" = "Chittagong",
  "Cumilla" = "Comilla",
  "Jashore" = "Jessore"
)
data$District <- str_replace_all(data$District, district_mapping)

# Create lag variables for lag analysis
# Sort by district and date to ensure proper lagging
data <- data |>
  arrange(District, year, month) |>
  group_by(District) |>
  mutate(
    # Lag 2 (2 months prior)
    NO2_lag2 = lag(mean_NO2_ensemble, 2),
    NOug_lag2 = lag(mean_NOppb, 2) * 1.23,
    SO2_lag2 = lag(mean_SO2_ensemble, 2),
    Ox_lag2 = lag(mean_Ox, 2),
    PM25_lag2 = lag(mean_PM25_ensemble, 2),
    PM10_lag2 = lag(mean_PM10ug, 2),
    CO_lag2 = lag(CO_pred, 2),
    O3_lag2 = lag(mean_O3_pred_cal, 2),
    
    # Lag 4 (4 months prior)
    NO2_lag4 = lag(mean_NO2_ensemble, 4),
    NOug_lag4 = lag(mean_NOppb, 4) * 1.23,
    SO2_lag4 = lag(mean_SO2_ensemble, 4),
    Ox_lag4 = lag(mean_Ox, 4),
    PM25_lag4 = lag(mean_PM25_ensemble, 4),
    PM10_lag4 = lag(mean_PM10ug, 4),
    CO_lag4 = lag(CO_pred, 4),
    O3_lag4 = lag(mean_O3_pred_cal, 4),
    
    # Lag 6 (6 months prior)
    NO2_lag6 = lag(mean_NO2_ensemble, 6),
    NOug_lag6 = lag(mean_NOppb, 6) * 1.23,
    SO2_lag6 = lag(mean_SO2_ensemble, 6),
    Ox_lag6 = lag(mean_Ox, 6),
    PM25_lag6 = lag(mean_PM25_ensemble, 6),
    PM10_lag6 = lag(mean_PM10ug, 6),
    CO_lag6 = lag(CO_pred, 6),
    O3_lag6 = lag(mean_O3_pred_cal, 6),
    
    # Lag 8 (8 months prior)
    NO2_lag8 = lag(mean_NO2_ensemble, 8),
    NOug_lag8 = lag(mean_NOppb, 8) * 1.23,
    SO2_lag8 = lag(mean_SO2_ensemble, 8),
    Ox_lag8 = lag(mean_Ox, 8),
    PM25_lag8 = lag(mean_PM25_ensemble, 8),
    PM10_lag8 = lag(mean_PM10ug, 8),
    CO_lag8 = lag(CO_pred, 8),
    O3_lag8 = lag(mean_O3_pred_cal, 8),
    
    # Cumulative 0-3 months (rolling average of current + prior 2 months)
    NO2_cum03 = rollmean(mean_NO2_ensemble, k = 3, fill = NA, align = "right"),
    NOug_cum03 = rollmean(mean_NOppb * 1.23, k = 3, fill = NA, align = "right"),
    SO2_cum03 = rollmean(mean_SO2_ensemble, k = 3, fill = NA, align = "right"),
    Ox_cum03 = rollmean(mean_Ox, k = 3, fill = NA, align = "right"),
    PM25_cum03 = rollmean(mean_PM25_ensemble, k = 3, fill = NA, align = "right"),
    PM10_cum03 = rollmean(mean_PM10ug, k = 3, fill = NA, align = "right"),
    CO_cum03 = rollmean(CO_pred, k = 3, fill = NA, align = "right"),
    O3_cum03 = rollmean(mean_O3_pred_cal, k = 3, fill = NA, align = "right"),
    
    # Cumulative 0-6 months (rolling average of current + prior 5 months)
    NO2_cum06 = rollmean(mean_NO2_ensemble, k = 6, fill = NA, align = "right"),
    NOug_cum06 = rollmean(mean_NOppb * 1.23, k = 6, fill = NA, align = "right"),
    SO2_cum06 = rollmean(mean_SO2_ensemble, k = 6, fill = NA, align = "right"),
    Ox_cum06 = rollmean(mean_Ox, k = 6, fill = NA, align = "right"),
    PM25_cum06 = rollmean(mean_PM25_ensemble, k = 6, fill = NA, align = "right"),
    PM10_cum06 = rollmean(mean_PM10ug, k = 6, fill = NA, align = "right"),
    CO_cum06 = rollmean(CO_pred, k = 6, fill = NA, align = "right"),
    O3_cum06 = rollmean(mean_O3_pred_cal, k = 6, fill = NA, align = "right"),
    
    # Cumulative 0-12 months (rolling average of current + prior 11 months)
    NO2_cum012 = rollmean(mean_NO2_ensemble, k = 12, fill = NA, align = "right"),
    NOug_cum012 = rollmean(mean_NOppb * 1.23, k = 12, fill = NA, align = "right"),
    SO2_cum012 = rollmean(mean_SO2_ensemble, k = 12, fill = NA, align = "right"),
    Ox_cum012 = rollmean(mean_Ox, k = 12, fill = NA, align = "right"),
    PM25_cum012 = rollmean(mean_PM25_ensemble, k = 12, fill = NA, align = "right"),
    PM10_cum012 = rollmean(mean_PM10ug, k = 12, fill = NA, align = "right"),
    CO_cum012 = rollmean(CO_pred, k = 12, fill = NA, align = "right"),
    O3_cum012 = rollmean(mean_O3_pred_cal, k = 12, fill = NA, align = "right")
  ) |>
  ungroup()

####################
### LAG ANALYSIS FOR ALL POLLUTANTS
### Lag 0, 2, 4, 6, 8, 0-3, 0-6, 0-12
####################

####################
### 1. NO2 LAG ANALYSIS
####################

# Lag 0 (current month) - same as main model
model_NO2_lag0 <- feglm(
  AggregatedDeath ~ I(mean_NO2_ensemble/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_NO2_lag2 <- feglm(
  AggregatedDeath ~ I(NO2_lag2/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_NO2_lag4 <- feglm(
  AggregatedDeath ~ I(NO2_lag4/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_NO2_lag6 <- feglm(
  AggregatedDeath ~ I(NO2_lag6/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_NO2_lag8 <- feglm(
  AggregatedDeath ~ I(NO2_lag8/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_NO2_cum03 <- feglm(
  AggregatedDeath ~ I(NO2_cum03/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_NO2_cum06 <- feglm(
  AggregatedDeath ~ I(NO2_cum06/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_NO2_cum012 <- feglm(
  AggregatedDeath ~ I(NO2_cum012/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract NO2 lag results
NO2_lag_results <- data.frame(
  Pollutant = "NO2",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_NO2_lag0)["I(mean_NO2_ensemble/10)"],
    coef(model_NO2_lag2)["I(NO2_lag2/10)"],
    coef(model_NO2_lag4)["I(NO2_lag4/10)"],
    coef(model_NO2_lag6)["I(NO2_lag6/10)"],
    coef(model_NO2_lag8)["I(NO2_lag8/10)"],
    coef(model_NO2_cum03)["I(NO2_cum03/10)"],
    coef(model_NO2_cum06)["I(NO2_cum06/10)"],
    coef(model_NO2_cum012)["I(NO2_cum012/10)"]
  ),
  SE = c(
    model_NO2_lag0$se["I(mean_NO2_ensemble/10)"],
    model_NO2_lag2$se["I(NO2_lag2/10)"],
    model_NO2_lag4$se["I(NO2_lag4/10)"],
    model_NO2_lag6$se["I(NO2_lag6/10)"],
    model_NO2_lag8$se["I(NO2_lag8/10)"],
    model_NO2_cum03$se["I(NO2_cum03/10)"],
    model_NO2_cum06$se["I(NO2_cum06/10)"],
    model_NO2_cum012$se["I(NO2_cum012/10)"]
  )
)

NO2_lag_results$RR <- exp(NO2_lag_results$Estimate)
NO2_lag_results$RR_L <- exp(NO2_lag_results$Estimate - 1.96 * NO2_lag_results$SE)
NO2_lag_results$RR_U <- exp(NO2_lag_results$Estimate + 1.96 * NO2_lag_results$SE)

NO2_lag_results

####################
### 2. NO LAG ANALYSIS
####################
data$mean_NOug <- data$mean_NOppb* 1.23  # Convert ppb to μg/m³

# Lag 0 (current month)
model_NO_lag0 <- feglm(
  AggregatedDeath ~ I(mean_NOug/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_NO_lag2 <- feglm(
  AggregatedDeath ~ I(NOug_lag2/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_NO_lag4 <- feglm(
  AggregatedDeath ~ I(NOug_lag4/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_NO_lag6 <- feglm(
  AggregatedDeath ~ I(NOug_lag6/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_NO_lag8 <- feglm(
  AggregatedDeath ~ I(NOug_lag8/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_NO_cum03 <- feglm(
  AggregatedDeath ~ I(NOug_cum03/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_NO_cum06 <- feglm(
  AggregatedDeath ~ I(NOug_cum06/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_NO_cum012 <- feglm(
  AggregatedDeath ~ I(NOug_cum012/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract NO lag results
NO_lag_results <- data.frame(
  Pollutant = "NO",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_NO_lag0)["I(mean_NOug/10)"],
    coef(model_NO_lag2)["I(NOug_lag2/10)"],
    coef(model_NO_lag4)["I(NOug_lag4/10)"],
    coef(model_NO_lag6)["I(NOug_lag6/10)"],
    coef(model_NO_lag8)["I(NOug_lag8/10)"],
    coef(model_NO_cum03)["I(NOug_cum03/10)"],
    coef(model_NO_cum06)["I(NOug_cum06/10)"],
    coef(model_NO_cum012)["I(NOug_cum012/10)"]
  ),
  SE = c(
    model_NO_lag0$se["I(mean_NOug/10)"],
    model_NO_lag2$se["I(NOug_lag2/10)"],
    model_NO_lag4$se["I(NOug_lag4/10)"],
    model_NO_lag6$se["I(NOug_lag6/10)"],
    model_NO_lag8$se["I(NOug_lag8/10)"],
    model_NO_cum03$se["I(NOug_cum03/10)"],
    model_NO_cum06$se["I(NOug_cum06/10)"],
    model_NO_cum012$se["I(NOug_cum012/10)"]
  )
)

NO_lag_results$RR <- exp(NO_lag_results$Estimate)
NO_lag_results$RR_L <- exp(NO_lag_results$Estimate - 1.96 * NO_lag_results$SE)
NO_lag_results$RR_U <- exp(NO_lag_results$Estimate + 1.96 * NO_lag_results$SE)

NO_lag_results

####################
### 3. SO2 LAG ANALYSIS
####################

# Lag 0 (current month)
model_SO2_lag0 <- feglm(
  AggregatedDeath ~ I(mean_SO2_ensemble/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_SO2_lag2 <- feglm(
  AggregatedDeath ~ I(SO2_lag2/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_SO2_lag4 <- feglm(
  AggregatedDeath ~ I(SO2_lag4/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_SO2_lag6 <- feglm(
  AggregatedDeath ~ I(SO2_lag6/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_SO2_lag8 <- feglm(
  AggregatedDeath ~ I(SO2_lag8/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_SO2_cum03 <- feglm(
  AggregatedDeath ~ I(SO2_cum03/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_SO2_cum06 <- feglm(
  AggregatedDeath ~ I(SO2_cum06/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_SO2_cum012 <- feglm(
  AggregatedDeath ~ I(SO2_cum012/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract SO2 lag results
SO2_lag_results <- data.frame(
  Pollutant = "SO2",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_SO2_lag0)["I(mean_SO2_ensemble/10)"],
    coef(model_SO2_lag2)["I(SO2_lag2/10)"],
    coef(model_SO2_lag4)["I(SO2_lag4/10)"],
    coef(model_SO2_lag6)["I(SO2_lag6/10)"],
    coef(model_SO2_lag8)["I(SO2_lag8/10)"],
    coef(model_SO2_cum03)["I(SO2_cum03/10)"],
    coef(model_SO2_cum06)["I(SO2_cum06/10)"],
    coef(model_SO2_cum012)["I(SO2_cum012/10)"]
  ),
  SE = c(
    model_SO2_lag0$se["I(mean_SO2_ensemble/10)"],
    model_SO2_lag2$se["I(SO2_lag2/10)"],
    model_SO2_lag4$se["I(SO2_lag4/10)"],
    model_SO2_lag6$se["I(SO2_lag6/10)"],
    model_SO2_lag8$se["I(SO2_lag8/10)"],
    model_SO2_cum03$se["I(SO2_cum03/10)"],
    model_SO2_cum06$se["I(SO2_cum06/10)"],
    model_SO2_cum012$se["I(SO2_cum012/10)"]
  )
)

SO2_lag_results$RR <- exp(SO2_lag_results$Estimate)
SO2_lag_results$RR_L <- exp(SO2_lag_results$Estimate - 1.96 * SO2_lag_results$SE)
SO2_lag_results$RR_U <- exp(SO2_lag_results$Estimate + 1.96 * SO2_lag_results$SE)

SO2_lag_results

####################
### 4. OX LAG ANALYSIS
####################

# Lag 0 (current month)
model_Ox_lag0 <- feglm(
  AggregatedDeath ~ I(mean_Ox/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_Ox_lag2 <- feglm(
  AggregatedDeath ~ I(Ox_lag2/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_Ox_lag4 <- feglm(
  AggregatedDeath ~ I(Ox_lag4/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_Ox_lag6 <- feglm(
  AggregatedDeath ~ I(Ox_lag6/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_Ox_lag8 <- feglm(
  AggregatedDeath ~ I(Ox_lag8/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_Ox_cum03 <- feglm(
  AggregatedDeath ~ I(Ox_cum03/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_Ox_cum06 <- feglm(
  AggregatedDeath ~ I(Ox_cum06/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_Ox_cum012 <- feglm(
  AggregatedDeath ~ I(Ox_cum012/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract Ox lag results
Ox_lag_results <- data.frame(
  Pollutant = "Ox",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_Ox_lag0)["I(mean_Ox/10)"],
    coef(model_Ox_lag2)["I(Ox_lag2/10)"],
    coef(model_Ox_lag4)["I(Ox_lag4/10)"],
    coef(model_Ox_lag6)["I(Ox_lag6/10)"],
    coef(model_Ox_lag8)["I(Ox_lag8/10)"],
    coef(model_Ox_cum03)["I(Ox_cum03/10)"],
    coef(model_Ox_cum06)["I(Ox_cum06/10)"],
    coef(model_Ox_cum012)["I(Ox_cum012/10)"]
  ),
  SE = c(
    model_Ox_lag0$se["I(mean_Ox/10)"],
    model_Ox_lag2$se["I(Ox_lag2/10)"],
    model_Ox_lag4$se["I(Ox_lag4/10)"],
    model_Ox_lag6$se["I(Ox_lag6/10)"],
    model_Ox_lag8$se["I(Ox_lag8/10)"],
    model_Ox_cum03$se["I(Ox_cum03/10)"],
    model_Ox_cum06$se["I(Ox_cum06/10)"],
    model_Ox_cum012$se["I(Ox_cum012/10)"]
  )
)

Ox_lag_results$RR <- exp(Ox_lag_results$Estimate)
Ox_lag_results$RR_L <- exp(Ox_lag_results$Estimate - 1.96 * Ox_lag_results$SE)
Ox_lag_results$RR_U <- exp(Ox_lag_results$Estimate + 1.96 * Ox_lag_results$SE)

Ox_lag_results

####################
### 5. PM25 LAG ANALYSIS
####################

# Lag 0 (current month)
model_PM25_lag0 <- feglm(
  AggregatedDeath ~ I(mean_PM25_ensemble/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_PM25_lag2 <- feglm(
  AggregatedDeath ~ I(PM25_lag2/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_PM25_lag4 <- feglm(
  AggregatedDeath ~ I(PM25_lag4/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_PM25_lag6 <- feglm(
  AggregatedDeath ~ I(PM25_lag6/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_PM25_lag8 <- feglm(
  AggregatedDeath ~ I(PM25_lag8/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_PM25_cum03 <- feglm(
  AggregatedDeath ~ I(PM25_cum03/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_PM25_cum06 <- feglm(
  AggregatedDeath ~ I(PM25_cum06/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_PM25_cum012 <- feglm(
  AggregatedDeath ~ I(PM25_cum012/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract PM25 lag results
PM25_lag_results <- data.frame(
  Pollutant = "PM25",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_PM25_lag0)["I(mean_PM25_ensemble/10)"],
    coef(model_PM25_lag2)["I(PM25_lag2/10)"],
    coef(model_PM25_lag4)["I(PM25_lag4/10)"],
    coef(model_PM25_lag6)["I(PM25_lag6/10)"],
    coef(model_PM25_lag8)["I(PM25_lag8/10)"],
    coef(model_PM25_cum03)["I(PM25_cum03/10)"],
    coef(model_PM25_cum06)["I(PM25_cum06/10)"],
    coef(model_PM25_cum012)["I(PM25_cum012/10)"]
  ),
  SE = c(
    model_PM25_lag0$se["I(mean_PM25_ensemble/10)"],
    model_PM25_lag2$se["I(PM25_lag2/10)"],
    model_PM25_lag4$se["I(PM25_lag4/10)"],
    model_PM25_lag6$se["I(PM25_lag6/10)"],
    model_PM25_lag8$se["I(PM25_lag8/10)"],
    model_PM25_cum03$se["I(PM25_cum03/10)"],
    model_PM25_cum06$se["I(PM25_cum06/10)"],
    model_PM25_cum012$se["I(PM25_cum012/10)"]
  )
)

PM25_lag_results$RR <- exp(PM25_lag_results$Estimate)
PM25_lag_results$RR_L <- exp(PM25_lag_results$Estimate - 1.96 * PM25_lag_results$SE)
PM25_lag_results$RR_U <- exp(PM25_lag_results$Estimate + 1.96 * PM25_lag_results$SE)

PM25_lag_results

####################
### 6. PM10 LAG ANALYSIS
####################

# Lag 0 (current month)
model_PM10_lag0 <- feglm(
  AggregatedDeath ~ I(mean_PM10ug/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_PM10_lag2 <- feglm(
  AggregatedDeath ~ I(PM10_lag2/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_PM10_lag4 <- feglm(
  AggregatedDeath ~ I(PM10_lag4/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_PM10_lag6 <- feglm(
  AggregatedDeath ~ I(PM10_lag6/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_PM10_lag8 <- feglm(
  AggregatedDeath ~ I(PM10_lag8/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_PM10_cum03 <- feglm(
  AggregatedDeath ~ I(PM10_cum03/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_PM10_cum06 <- feglm(
  AggregatedDeath ~ I(PM10_cum06/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_PM10_cum012 <- feglm(
  AggregatedDeath ~ I(PM10_cum012/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract PM10 lag results
PM10_lag_results <- data.frame(
  Pollutant = "PM10",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_PM10_lag0)["I(mean_PM10ug/10)"],
    coef(model_PM10_lag2)["I(PM10_lag2/10)"],
    coef(model_PM10_lag4)["I(PM10_lag4/10)"],
    coef(model_PM10_lag6)["I(PM10_lag6/10)"],
    coef(model_PM10_lag8)["I(PM10_lag8/10)"],
    coef(model_PM10_cum03)["I(PM10_cum03/10)"],
    coef(model_PM10_cum06)["I(PM10_cum06/10)"],
    coef(model_PM10_cum012)["I(PM10_cum012/10)"]
  ),
  SE = c(
    model_PM10_lag0$se["I(mean_PM10ug/10)"],
    model_PM10_lag2$se["I(PM10_lag2/10)"],
    model_PM10_lag4$se["I(PM10_lag4/10)"],
    model_PM10_lag6$se["I(PM10_lag6/10)"],
    model_PM10_lag8$se["I(PM10_lag8/10)"],
    model_PM10_cum03$se["I(PM10_cum03/10)"],
    model_PM10_cum06$se["I(PM10_cum06/10)"],
    model_PM10_cum012$se["I(PM10_cum012/10)"]
  )
)

PM10_lag_results$RR <- exp(PM10_lag_results$Estimate)
PM10_lag_results$RR_L <- exp(PM10_lag_results$Estimate - 1.96 * PM10_lag_results$SE)
PM10_lag_results$RR_U <- exp(PM10_lag_results$Estimate + 1.96 * PM10_lag_results$SE)

PM10_lag_results

####################
### 7. CO LAG ANALYSIS
####################

# Lag 0 (current month)
model_CO_lag0 <- feglm(
  AggregatedDeath ~ I(CO_pred/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_CO_lag2 <- feglm(
  AggregatedDeath ~ I(CO_lag2/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_CO_lag4 <- feglm(
  AggregatedDeath ~ I(CO_lag4/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_CO_lag6 <- feglm(
  AggregatedDeath ~ I(CO_lag6/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_CO_lag8 <- feglm(
  AggregatedDeath ~ I(CO_lag8/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_CO_cum03 <- feglm(
  AggregatedDeath ~ I(CO_cum03/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_CO_cum06 <- feglm(
  AggregatedDeath ~ I(CO_cum06/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_CO_cum012 <- feglm(
  AggregatedDeath ~ I(CO_cum012/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract CO lag results
CO_lag_results <- data.frame(
  Pollutant = "CO",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_CO_lag0)["I(CO_pred/10)"],
    coef(model_CO_lag2)["I(CO_lag2/10)"],
    coef(model_CO_lag4)["I(CO_lag4/10)"],
    coef(model_CO_lag6)["I(CO_lag6/10)"],
    coef(model_CO_lag8)["I(CO_lag8/10)"],
    coef(model_CO_cum03)["I(CO_cum03/10)"],
    coef(model_CO_cum06)["I(CO_cum06/10)"],
    coef(model_CO_cum012)["I(CO_cum012/10)"]
  ),
  SE = c(
    model_CO_lag0$se["I(CO_pred/10)"],
    model_CO_lag2$se["I(CO_lag2/10)"],
    model_CO_lag4$se["I(CO_lag4/10)"],
    model_CO_lag6$se["I(CO_lag6/10)"],
    model_CO_lag8$se["I(CO_lag8/10)"],
    model_CO_cum03$se["I(CO_cum03/10)"],
    model_CO_cum06$se["I(CO_cum06/10)"],
    model_CO_cum012$se["I(CO_cum012/10)"]
  )
)

CO_lag_results$RR <- exp(CO_lag_results$Estimate)
CO_lag_results$RR_L <- exp(CO_lag_results$Estimate - 1.96 * CO_lag_results$SE)
CO_lag_results$RR_U <- exp(CO_lag_results$Estimate + 1.96 * CO_lag_results$SE)

CO_lag_results

####################
### 8. O3 LAG ANALYSIS
####################

# Lag 0 (current month)
model_O3_lag0 <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 2 (2 months prior)
model_O3_lag2 <- feglm(
  AggregatedDeath ~ I(O3_lag2/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 4 (4 months prior)
model_O3_lag4 <- feglm(
  AggregatedDeath ~ I(O3_lag4/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 6 (6 months prior)
model_O3_lag6 <- feglm(
  AggregatedDeath ~ I(O3_lag6/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Lag 8 (8 months prior)
model_O3_lag8 <- feglm(
  AggregatedDeath ~ I(O3_lag8/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-3 months
model_O3_cum03 <- feglm(
  AggregatedDeath ~ I(O3_cum03/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-6 months
model_O3_cum06 <- feglm(
  AggregatedDeath ~ I(O3_cum06/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-12 months
model_O3_cum012 <- feglm(
  AggregatedDeath ~ I(O3_cum012/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract O3 lag results
O3_lag_results <- data.frame(
  Pollutant = "O3",
  Lag = c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12"),
  Estimate = c(
    coef(model_O3_lag0)["I(mean_O3_pred_cal/10)"],
    coef(model_O3_lag2)["I(O3_lag2/10)"],
    coef(model_O3_lag4)["I(O3_lag4/10)"],
    coef(model_O3_lag6)["I(O3_lag6/10)"],
    coef(model_O3_lag8)["I(O3_lag8/10)"],
    coef(model_O3_cum03)["I(O3_cum03/10)"],
    coef(model_O3_cum06)["I(O3_cum06/10)"],
    coef(model_O3_cum012)["I(O3_cum012/10)"]
  ),
  SE = c(
    model_O3_lag0$se["I(mean_O3_pred_cal/10)"],
    model_O3_lag2$se["I(O3_lag2/10)"],
    model_O3_lag4$se["I(O3_lag4/10)"],
    model_O3_lag6$se["I(O3_lag6/10)"],
    model_O3_lag8$se["I(O3_lag8/10)"],
    model_O3_cum03$se["I(O3_cum03/10)"],
    model_O3_cum06$se["I(O3_cum06/10)"],
    model_O3_cum012$se["I(O3_cum012/10)"]
  )
)

O3_lag_results$RR <- exp(O3_lag_results$Estimate)
O3_lag_results$RR_L <- exp(O3_lag_results$Estimate - 1.96 * O3_lag_results$SE)
O3_lag_results$RR_U <- exp(O3_lag_results$Estimate + 1.96 * O3_lag_results$SE)

O3_lag_results

# Combine all lag results
all_lag_results <- rbind(NO2_lag_results, NO_lag_results, SO2_lag_results, Ox_lag_results, PM25_lag_results, PM10_lag_results, CO_lag_results, O3_lag_results)

all_lag_results

# Calculate RR and 95% CIs for plotting
all_lag_results$RR <- exp(all_lag_results$Estimate)
all_lag_results$RR_L <- exp(all_lag_results$Estimate - 1.96 * all_lag_results$SE)
all_lag_results$RR_U <- exp(all_lag_results$Estimate + 1.96 * all_lag_results$SE)

####################
### LAG EFFECTS PLOTS
####################

# Define lag order for plotting
lag_order <- c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12")

# Create individual pollutant plots
plot_lag_NO2 <- all_lag_results |>
  filter(Pollutant == "NO2") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "NO2", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_NO <- all_lag_results |>
  filter(Pollutant == "NO") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "NO", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_SO2 <- all_lag_results |>
  filter(Pollutant == "SO2") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "SO2", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_Ox <- all_lag_results |>
  filter(Pollutant == "Ox") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "Ox", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_PM25 <- all_lag_results |>
  filter(Pollutant == "PM25") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "PM2.5", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_PM10 <- all_lag_results |>
  filter(Pollutant == "PM10") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "PM10", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_CO <- all_lag_results |>
  filter(Pollutant == "CO") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "CO", x = "Lag", y = "Risk Ratio (95% CI) per 10 ppm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_O3 <- all_lag_results |>
  filter(Pollutant == "O3") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "O3", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine all lag plots (2x4 grid)
library(patchwork)

combined_lag_plot <- (plot_lag_PM25 + plot_lag_PM10 + plot_lag_NO2 + plot_lag_NO) / 
                     (plot_lag_Ox + plot_lag_SO2 + plot_lag_CO + plot_lag_O3) +
  plot_annotation(
    title = "Lag Effects of Air Pollutants on Mortality",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 12),
                  plot.caption = element_text(size = 10, hjust = 0))
  )

combined_lag_plot


#################################################################################
#################################################################################
# ABSOLUTE MORTALITY INCREASE PER 100,000 POPULATION - LAG ANALYSIS
#################################################################################
#################################################################################

cat("\n========== ABSOLUTE MORTALITY INCREASE FOR LAG EFFECTS ==========\n")

# Calculate baseline monthly mortality rate per 100,000 people
baseline_rate_100k <- with(data,
  mean(AggregatedDeath / Population_2022, na.rm = TRUE) * 1e5
)

cat("Baseline monthly mortality rate:", sprintf("%.2f", baseline_rate_100k), "per 100,000 population\n\n")

# Convert lag results to absolute mortality increase
absolute_lag_results <- all_lag_results |>
  mutate(
    Delta_Deaths_100k = baseline_rate_100k * (RR - 1),
    Delta_Deaths_100k_L = baseline_rate_100k * (RR_L - 1),
    Delta_Deaths_100k_U = baseline_rate_100k * (RR_U - 1)
  ) |>
  select(Pollutant, Lag, RR, RR_L, RR_U, Delta_Deaths_100k, Delta_Deaths_100k_L, Delta_Deaths_100k_U)

cat("Absolute mortality increase per 100,000 population per month (all lags):\n")
print(absolute_lag_results)

# Create visualization by pollutant showing lag effects with absolute mortality scale
lag_order <- c("Lag 0", "Lag 2", "Lag 4", "Lag 6", "Lag 8", "Lag 0-3", "Lag 0-6", "Lag 0-12")

# Plot for each pollutant
plot_abs_lag_NO2 <- absolute_lag_results |>
  filter(Pollutant == "NO2") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "NO2", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_abs_lag_NO <- absolute_lag_results |>
  filter(Pollutant == "NO") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "NO", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_abs_lag_Ox <- absolute_lag_results |>
  filter(Pollutant == "Ox") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "Ox", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_abs_lag_SO2 <- absolute_lag_results |>
  filter(Pollutant == "SO2") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "SO2", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_abs_lag_PM25 <- absolute_lag_results |>
  filter(Pollutant == "PM25") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "PM2.5", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_abs_lag_PM10 <- absolute_lag_results |>
  filter(Pollutant == "PM10") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "PM10", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_abs_lag_CO <- absolute_lag_results |>
  filter(Pollutant == "CO") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "CO", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_abs_lag_O3 <- absolute_lag_results |>
  filter(Pollutant == "O3") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), width = 0.2, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal(base_size = 11) +
  labs(title = "O3", x = "Lag", y = "Deaths per 100,000 per Month (95% CI)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine all absolute lag plots (2x4 grid)
combined_abs_lag_plot <- (plot_abs_lag_PM25 + plot_abs_lag_PM10 + plot_abs_lag_NO2 + plot_abs_lag_NO) / 
                         (plot_abs_lag_Ox + plot_abs_lag_SO2 + plot_abs_lag_CO + plot_abs_lag_O3) +
  plot_annotation(
    title = "Absolute Mortality Increase (per 100,000 population per month) - Lag Effects",
    subtitle = "Associated with 10-unit increase in pollutant concentration",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

cat("\n========== ABSOLUTE LAG EFFECTS PLOTS ==========\n")
combined_abs_lag_plot

#################################################################################