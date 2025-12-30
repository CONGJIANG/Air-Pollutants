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
    # Lag 1 (1 month prior)
    NO2_lag1 = lag(mean_NO2_ensemble, 1),
    NOug_lag1 = lag(mean_NOppb, 1) * 1.23,
    SO2_lag1 = lag(mean_SO2_ensemble, 1),
    Ox_lag1 = lag(mean_Ox, 1),
    PM25_lag1 = lag(mean_PM25_ensemble, 1),
    PM10_lag1 = lag(mean_PM10ug, 1),
    CO_lag1 = lag(CO_pred, 1),
    O3_lag1 = lag(mean_O3_pred_cal, 1),
    
    # Lag 2 (2 months prior)
    NO2_lag2 = lag(mean_NO2_ensemble, 2),
    NOug_lag2 = lag(mean_NOppb, 2) * 1.23,
    SO2_lag2 = lag(mean_SO2_ensemble, 2),
    Ox_lag2 = lag(mean_Ox, 2),
    PM25_lag2 = lag(mean_PM25_ensemble, 2),
    PM10_lag2 = lag(mean_PM10ug, 2),
    CO_lag2 = lag(CO_pred, 2),
    O3_lag2 = lag(mean_O3_pred_cal, 2),
    
    # Lag 3 (3 months prior)
    NO2_lag3 = lag(mean_NO2_ensemble, 3),
    NOug_lag3 = lag(mean_NOppb, 3) * 1.23,
    SO2_lag3 = lag(mean_SO2_ensemble, 3),
    Ox_lag3 = lag(mean_Ox, 3),
    PM25_lag3 = lag(mean_PM25_ensemble, 3),
    PM10_lag3 = lag(mean_PM10ug, 3),
    CO_lag3 = lag(CO_pred, 3),
    O3_lag3 = lag(mean_O3_pred_cal, 3),
    
    # Lag 4 (4 months prior)
    NO2_lag4 = lag(mean_NO2_ensemble, 4),
    NOug_lag4 = lag(mean_NOppb, 4) * 1.23,
    SO2_lag4 = lag(mean_SO2_ensemble, 4),
    Ox_lag4 = lag(mean_Ox, 4),
    PM25_lag4 = lag(mean_PM25_ensemble, 4),
    PM10_lag4 = lag(mean_PM10ug, 4),
    CO_lag4 = lag(CO_pred, 4),
    O3_lag4 = lag(mean_O3_pred_cal, 4),
    
    # Cumulative 0-12 months (rolling average of current + prior 11 months)
    NO2_cum012 = rollmean(mean_NO2_ensemble, k = 12, fill = NA, align = "right"),
    NOug_cum012 = rollmean(mean_NOppb * 1.23, k = 12, fill = NA, align = "right"),
    SO2_cum012 = rollmean(mean_SO2_ensemble, k = 12, fill = NA, align = "right"),
    Ox_cum012 = rollmean(mean_Ox, k = 12, fill = NA, align = "right"),
    PM25_cum012 = rollmean(mean_PM25_ensemble, k = 12, fill = NA, align = "right"),
    PM10_cum012 = rollmean(mean_PM10ug, k = 12, fill = NA, align = "right"),
    CO_cum012 = rollmean(CO_pred, k = 12, fill = NA, align = "right"),
    O3_cum012 = rollmean(mean_O3_pred_cal, k = 12, fill = NA, align = "right"),
    
    # Cumulative 0-24 months (rolling average of current + prior 23 months)
    NO2_cum024 = rollmean(mean_NO2_ensemble, k = 24, fill = NA, align = "right"),
    NOug_cum024 = rollmean(mean_NOppb * 1.23, k = 24, fill = NA, align = "right"),
    SO2_cum024 = rollmean(mean_SO2_ensemble, k = 24, fill = NA, align = "right"),
    Ox_cum024 = rollmean(mean_Ox, k = 24, fill = NA, align = "right"),
    PM25_cum024 = rollmean(mean_PM25_ensemble, k = 24, fill = NA, align = "right"),
    PM10_cum024 = rollmean(mean_PM10ug, k = 24, fill = NA, align = "right"),
    CO_cum024 = rollmean(CO_pred, k = 24, fill = NA, align = "right"),
    O3_cum024 = rollmean(mean_O3_pred_cal, k = 24, fill = NA, align = "right")
  ) |>
  ungroup()

####################
### LAG ANALYSIS FOR ALL POLLUTANTS
### Lag 0 (concurrent), Lag 1, Lag 0-12, and Lag 0-24
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

# Lag 1 (1 month prior)
model_NO2_lag1 <- feglm(
  AggregatedDeath ~ I(NO2_lag1/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
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

# Lag 3 (3 months prior)
model_NO2_lag3 <- feglm(
  AggregatedDeath ~ I(NO2_lag3/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
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

# Cumulative 0-12 months
model_NO2_cum012 <- feglm(
  AggregatedDeath ~ I(NO2_cum012/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_NO2_cum024 <- feglm(
  AggregatedDeath ~ I(NO2_cum024/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract NO2 lag results
NO2_lag_results <- data.frame(
  Pollutant = "NO2",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_NO2_lag0)["I(mean_NO2_ensemble/10)"],
    coef(model_NO2_lag1)["I(NO2_lag1/10)"],
    coef(model_NO2_lag2)["I(NO2_lag2/10)"],
    coef(model_NO2_lag3)["I(NO2_lag3/10)"],
    coef(model_NO2_lag4)["I(NO2_lag4/10)"],
    coef(model_NO2_cum012)["I(NO2_cum012/10)"],
    coef(model_NO2_cum024)["I(NO2_cum024/10)"]
  ),
  SE = c(
    model_NO2_lag0$se["I(mean_NO2_ensemble/10)"],
    model_NO2_lag1$se["I(NO2_lag1/10)"],
    model_NO2_lag2$se["I(NO2_lag2/10)"],
    model_NO2_lag3$se["I(NO2_lag3/10)"],
    model_NO2_lag4$se["I(NO2_lag4/10)"],
    model_NO2_cum012$se["I(NO2_cum012/10)"],
    model_NO2_cum024$se["I(NO2_cum024/10)"]
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

# Lag 1 (1 month prior)
model_NO_lag1 <- feglm(
  AggregatedDeath ~ I(NOug_lag1/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
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

# Lag 3 (3 months prior)
model_NO_lag3 <- feglm(
  AggregatedDeath ~ I(NOug_lag3/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
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

# Cumulative 0-12 months
model_NO_cum012 <- feglm(
  AggregatedDeath ~ I(NOug_cum012/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_NO_cum024 <- feglm(
  AggregatedDeath ~ I(NOug_cum024/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract NO lag results
NO_lag_results <- data.frame(
  Pollutant = "NO",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_NO_lag0)["I(mean_NOug/10)"],
    coef(model_NO_lag1)["I(NOug_lag1/10)"],
    coef(model_NO_lag2)["I(NOug_lag2/10)"],
    coef(model_NO_lag3)["I(NOug_lag3/10)"],
    coef(model_NO_lag4)["I(NOug_lag4/10)"],
    coef(model_NO_cum012)["I(NOug_cum012/10)"],
    coef(model_NO_cum024)["I(NOug_cum024/10)"]
  ),
  SE = c(
    model_NO_lag0$se["I(mean_NOug/10)"],
    model_NO_lag1$se["I(NOug_lag1/10)"],
    model_NO_lag2$se["I(NOug_lag2/10)"],
    model_NO_lag3$se["I(NOug_lag3/10)"],
    model_NO_lag4$se["I(NOug_lag4/10)"],
    model_NO_cum012$se["I(NOug_cum012/10)"],
    model_NO_cum024$se["I(NOug_cum024/10)"]
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

# Lag 1 (1 month prior)
model_SO2_lag1 <- feglm(
  AggregatedDeath ~ I(SO2_lag1/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
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

# Lag 3 (3 months prior)
model_SO2_lag3 <- feglm(
  AggregatedDeath ~ I(SO2_lag3/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
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

# Cumulative 0-12 months
model_SO2_cum012 <- feglm(
  AggregatedDeath ~ I(SO2_cum012/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_SO2_cum024 <- feglm(
  AggregatedDeath ~ I(SO2_cum024/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract SO2 lag results
SO2_lag_results <- data.frame(
  Pollutant = "SO2",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_SO2_lag0)["I(mean_SO2_ensemble/10)"],
    coef(model_SO2_lag1)["I(SO2_lag1/10)"],
    coef(model_SO2_lag2)["I(SO2_lag2/10)"],
    coef(model_SO2_lag3)["I(SO2_lag3/10)"],
    coef(model_SO2_lag4)["I(SO2_lag4/10)"],
    coef(model_SO2_cum012)["I(SO2_cum012/10)"],
    coef(model_SO2_cum024)["I(SO2_cum024/10)"]
  ),
  SE = c(
    model_SO2_lag0$se["I(mean_SO2_ensemble/10)"],
    model_SO2_lag1$se["I(SO2_lag1/10)"],
    model_SO2_lag2$se["I(SO2_lag2/10)"],
    model_SO2_lag3$se["I(SO2_lag3/10)"],
    model_SO2_lag4$se["I(SO2_lag4/10)"],
    model_SO2_cum012$se["I(SO2_cum012/10)"],
    model_SO2_cum024$se["I(SO2_cum024/10)"]
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

# Lag 1 (1 month prior)
model_Ox_lag1 <- feglm(
  AggregatedDeath ~ I(Ox_lag1/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
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

# Lag 3 (3 months prior)
model_Ox_lag3 <- feglm(
  AggregatedDeath ~ I(Ox_lag3/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
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

# Cumulative 0-12 months
model_Ox_cum012 <- feglm(
  AggregatedDeath ~ I(Ox_cum012/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_Ox_cum024 <- feglm(
  AggregatedDeath ~ I(Ox_cum024/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract Ox lag results
Ox_lag_results <- data.frame(
  Pollutant = "Ox",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_Ox_lag0)["I(mean_Ox/10)"],
    coef(model_Ox_lag1)["I(Ox_lag1/10)"],
    coef(model_Ox_lag2)["I(Ox_lag2/10)"],
    coef(model_Ox_lag3)["I(Ox_lag3/10)"],
    coef(model_Ox_lag4)["I(Ox_lag4/10)"],
    coef(model_Ox_cum012)["I(Ox_cum012/10)"],
    coef(model_Ox_cum024)["I(Ox_cum024/10)"]
  ),
  SE = c(
    model_Ox_lag0$se["I(mean_Ox/10)"],
    model_Ox_lag1$se["I(Ox_lag1/10)"],
    model_Ox_lag2$se["I(Ox_lag2/10)"],
    model_Ox_lag3$se["I(Ox_lag3/10)"],
    model_Ox_lag4$se["I(Ox_lag4/10)"],
    model_Ox_cum012$se["I(Ox_cum012/10)"],
    model_Ox_cum024$se["I(Ox_cum024/10)"]
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

# Lag 1 (1 month prior)
model_PM25_lag1 <- feglm(
  AggregatedDeath ~ I(PM25_lag1/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
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

# Lag 3 (3 months prior)
model_PM25_lag3 <- feglm(
  AggregatedDeath ~ I(PM25_lag3/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
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

# Cumulative 0-12 months
model_PM25_cum012 <- feglm(
  AggregatedDeath ~ I(PM25_cum012/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_PM25_cum024 <- feglm(
  AggregatedDeath ~ I(PM25_cum024/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract PM25 lag results
PM25_lag_results <- data.frame(
  Pollutant = "PM25",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_PM25_lag0)["I(mean_PM25_ensemble/10)"],
    coef(model_PM25_lag1)["I(PM25_lag1/10)"],
    coef(model_PM25_lag2)["I(PM25_lag2/10)"],
    coef(model_PM25_lag3)["I(PM25_lag3/10)"],
    coef(model_PM25_lag4)["I(PM25_lag4/10)"],
    coef(model_PM25_cum012)["I(PM25_cum012/10)"],
    coef(model_PM25_cum024)["I(PM25_cum024/10)"]
  ),
  SE = c(
    model_PM25_lag0$se["I(mean_PM25_ensemble/10)"],
    model_PM25_lag1$se["I(PM25_lag1/10)"],
    model_PM25_lag2$se["I(PM25_lag2/10)"],
    model_PM25_lag3$se["I(PM25_lag3/10)"],
    model_PM25_lag4$se["I(PM25_lag4/10)"],
    model_PM25_cum012$se["I(PM25_cum012/10)"],
    model_PM25_cum024$se["I(PM25_cum024/10)"]
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

# Lag 1 (1 month prior)
model_PM10_lag1 <- feglm(
  AggregatedDeath ~ I(PM10_lag1/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
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

# Lag 3 (3 months prior)
model_PM10_lag3 <- feglm(
  AggregatedDeath ~ I(PM10_lag3/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
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

# Cumulative 0-12 months
model_PM10_cum012 <- feglm(
  AggregatedDeath ~ I(PM10_cum012/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_PM10_cum024 <- feglm(
  AggregatedDeath ~ I(PM10_cum024/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract PM10 lag results
PM10_lag_results <- data.frame(
  Pollutant = "PM10",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_PM10_lag0)["I(mean_PM10ug/10)"],
    coef(model_PM10_lag1)["I(PM10_lag1/10)"],
    coef(model_PM10_lag2)["I(PM10_lag2/10)"],
    coef(model_PM10_lag3)["I(PM10_lag3/10)"],
    coef(model_PM10_lag4)["I(PM10_lag4/10)"],
    coef(model_PM10_cum012)["I(PM10_cum012/10)"],
    coef(model_PM10_cum024)["I(PM10_cum024/10)"]
  ),
  SE = c(
    model_PM10_lag0$se["I(mean_PM10ug/10)"],
    model_PM10_lag1$se["I(PM10_lag1/10)"],
    model_PM10_lag2$se["I(PM10_lag2/10)"],
    model_PM10_lag3$se["I(PM10_lag3/10)"],
    model_PM10_lag4$se["I(PM10_lag4/10)"],
    model_PM10_cum012$se["I(PM10_cum012/10)"],
    model_PM10_cum024$se["I(PM10_cum024/10)"]
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

# Lag 1 (1 month prior)
model_CO_lag1 <- feglm(
  AggregatedDeath ~ I(CO_lag1/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
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

# Lag 3 (3 months prior)
model_CO_lag3 <- feglm(
  AggregatedDeath ~ I(CO_lag3/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
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

# Cumulative 0-12 months
model_CO_cum012 <- feglm(
  AggregatedDeath ~ I(CO_cum012/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_CO_cum024 <- feglm(
  AggregatedDeath ~ I(CO_cum024/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract CO lag results
CO_lag_results <- data.frame(
  Pollutant = "CO",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_CO_lag0)["I(CO_pred/10)"],
    coef(model_CO_lag1)["I(CO_lag1/10)"],
    coef(model_CO_lag2)["I(CO_lag2/10)"],
    coef(model_CO_lag3)["I(CO_lag3/10)"],
    coef(model_CO_lag4)["I(CO_lag4/10)"],
    coef(model_CO_cum012)["I(CO_cum012/10)"],
    coef(model_CO_cum024)["I(CO_cum024/10)"]
  ),
  SE = c(
    model_CO_lag0$se["I(CO_pred/10)"],
    model_CO_lag1$se["I(CO_lag1/10)"],
    model_CO_lag2$se["I(CO_lag2/10)"],
    model_CO_lag3$se["I(CO_lag3/10)"],
    model_CO_lag4$se["I(CO_lag4/10)"],
    model_CO_cum012$se["I(CO_cum012/10)"],
    model_CO_cum024$se["I(CO_cum024/10)"]
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

# Lag 1 (1 month prior)
model_O3_lag1 <- feglm(
  AggregatedDeath ~ I(O3_lag1/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
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

# Lag 3 (3 months prior)
model_O3_lag3 <- feglm(
  AggregatedDeath ~ I(O3_lag3/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
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

# Cumulative 0-12 months
model_O3_cum012 <- feglm(
  AggregatedDeath ~ I(O3_cum012/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Cumulative 0-24 months
model_O3_cum024 <- feglm(
  AggregatedDeath ~ I(O3_cum024/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

# Extract O3 lag results
O3_lag_results <- data.frame(
  Pollutant = "O3",
  Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12", "Lag 0-24"),
  Estimate = c(
    coef(model_O3_lag0)["I(mean_O3_pred_cal/10)"],
    coef(model_O3_lag1)["I(O3_lag1/10)"],
    coef(model_O3_lag2)["I(O3_lag2/10)"],
    coef(model_O3_lag3)["I(O3_lag3/10)"],
    coef(model_O3_lag4)["I(O3_lag4/10)"],
    coef(model_O3_cum012)["I(O3_cum012/10)"],
    coef(model_O3_cum024)["I(O3_cum024/10)"]
  ),
  SE = c(
    model_O3_lag0$se["I(mean_O3_pred_cal/10)"],
    model_O3_lag1$se["I(O3_lag1/10)"],
    model_O3_lag2$se["I(O3_lag2/10)"],
    model_O3_lag3$se["I(O3_lag3/10)"],
    model_O3_lag4$se["I(O3_lag4/10)"],
    model_O3_cum012$se["I(O3_cum012/10)"],
    model_O3_cum024$se["I(O3_cum024/10)"]
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
lag_order <- c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 0-12")

# Create individual pollutant plots
plot_lag_NO2 <- all_lag_results |>
  filter(Pollutant == "NO2", Lag != "Lag 0-24") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "NO2", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_NO <- all_lag_results |>
  filter(Pollutant == "NO", Lag != "Lag 0-24") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "NO", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_SO2 <- all_lag_results |>
  filter(Pollutant == "SO2", Lag != "Lag 0-24") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "SO2", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_Ox <- all_lag_results |>
  filter(Pollutant == "Ox", Lag != "Lag 0-24") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "Ox", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_PM25 <- all_lag_results |>
  filter(Pollutant == "PM25", Lag != "Lag 0-24") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "PM2.5", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_PM10 <- all_lag_results |>
  filter(Pollutant == "PM10", Lag != "Lag 0-24") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "PM10", x = "Lag", y = "Risk Ratio (95% CI) per 10 µg/m³") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_CO <- all_lag_results |>
  filter(Pollutant == "CO", Lag != "Lag 0-24") |>
  mutate(Lag = factor(Lag, levels = lag_order)) |>
  ggplot(aes(x = Lag, y = RR)) +
  geom_point(size = 3, color = "darkblue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, color = "darkblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  labs(title = "CO", x = "Lag", y = "Risk Ratio (95% CI) per 10 ppm") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_lag_O3 <- all_lag_results |>
  filter(Pollutant == "O3", Lag != "Lag 0-24") |>
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

combined_lag_plot <- (plot_lag_NO2 + plot_lag_NO + plot_lag_Ox + plot_lag_SO2) / 
                     (plot_lag_PM25 + plot_lag_PM10 + plot_lag_CO + plot_lag_O3) +
  plot_annotation(
    title = "Lag Effects of Air Pollutants on Mortality",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 12),
                  plot.caption = element_text(size = 10, hjust = 0))
  )

combined_lag_plot

