### reference model and alternative specifications
library(dplyr)
library(fixest)
library(ggplot2)
library(splines)
library(zoo)

data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")


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


summary(model_NO2_final) 
# I(mean_NO2_ensemble/5)           0.110414   0.028574  3.864131 1.1302e-04 ***
(NO2_obs_coef <- coef(model_NO2_final)["I(mean_NO2_ensemble/5)"])
(NO2_L <- NO2_obs_coef - 1.96*model_NO2_final$se["I(mean_NO2_ensemble/5)"])
(NO2_U <- NO2_obs_coef + 1.96*model_NO2_final$se["I(mean_NO2_ensemble/5)"])

# Exposure-response curve for NO2
# Create a sequence of NO2 values across the observed range
NO2_seq <- seq(min(data$mean_NO2_ensemble, na.rm = TRUE), 
               max(data$mean_NO2_ensemble, na.rm = TRUE), 
               length.out = 100)

# Create prediction dataset with NO2 varying and other variables at their mean
pred_data_NO2 <- data.frame(
  mean_NO2_ensemble = NO2_seq,
  mean_t2m_c = mean(data$mean_t2m_c, na.rm = TRUE),
  mean_RH = mean(data$mean_RH, na.rm = TRUE),
  mean_Precipitation = mean(data$mean_Precipitation, na.rm = TRUE),
  District = data$District[1],
  month = data$month[1],
  Population_2022 = 1
)

# Get predictions on log scale
pred_NO2 <- predict(model_NO2_final, newdata = pred_data_NO2)

# Get coefficient and SE for confidence interval
coef_se <- model_NO2_final$se["I(mean_NO2_ensemble/5)"]

# Create data frame with RR and 95% CI
NO2_er_df <- data.frame(
  NO2 = NO2_seq,
  RR = exp(pred_NO2),
  RR_L = exp(pred_NO2 - 1.96 * coef_se),
  RR_U = exp(pred_NO2 + 1.96 * coef_se)
)

ggplot(NO2_er_df, aes(x = NO2, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: NO2 and Death Risk",
    x = "NO2 (µg/m³)",
    y = "Risk Ratio (95% CI)"
  )

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

summary(model_NO_final)

(NO_obs_coef <- coef(model_NO_final)["I(mean_NOppb/5)"])
(NO_L <- NO_obs_coef - 1.96*model_NO_final$se["I(mean_NOppb/5)"])
(NO_U <- NO_obs_coef + 1.96*model_NO_final$se["I(mean_NOppb/5)"])

# Exposure-response curve for NO
NO_seq <- seq(min(data$mean_NOppb, na.rm = TRUE), 
              max(data$mean_NOppb, na.rm = TRUE), 
              length.out = 100)

pred_data_NO <- data.frame(
  mean_NOppb = NO_seq,
  mean_t2m_c = mean(data$mean_t2m_c, na.rm = TRUE),
  mean_RH = mean(data$mean_RH, na.rm = TRUE),
  mean_wind_speed = mean(data$mean_wind_speed, na.rm = TRUE),
  District = data$District[1],
  month = data$month[1],
  Population_2022 = 1
)

pred_NO <- predict(model_NO_final, newdata = pred_data_NO)
coef_se_NO <- model_NO_final$se["I(mean_NOppb/5)"]

NO_er_df <- data.frame(
  NO = NO_seq,
  RR = exp(pred_NO),
  RR_L = exp(pred_NO - 1.96 * coef_se_NO),
  RR_U = exp(pred_NO + 1.96 * coef_se_NO)
)

ggplot(NO_er_df, aes(x = NO, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: NO and Death Risk",
    x = "NO (ppb)",
    y = "Risk Ratio (95% CI)"
  )


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

summary(model_SO2_final)

(SO2_obs_coef <- coef(model_SO2_final)["I(mean_SO2_ensemble/5)"])
(SO2_L <- SO2_obs_coef - 1.96*model_SO2_final$se["I(mean_SO2_ensemble/5)"])
(SO2_U <- SO2_obs_coef + 1.96*model_SO2_final$se["I(mean_SO2_ensemble/5)"])

# Exposure-response curve for SO2
SO2_seq <- seq(min(data$mean_SO2_ensemble, na.rm = TRUE), 
               max(data$mean_SO2_ensemble, na.rm = TRUE), 
               length.out = 100)

pred_data_SO2 <- data.frame(
  mean_SO2_ensemble = SO2_seq,
  mean_t2m_c = mean(data$mean_t2m_c, na.rm = TRUE),
  mean_RH = mean(data$mean_RH, na.rm = TRUE),
  mean_Precipitation = mean(data$mean_Precipitation, na.rm = TRUE),
  District = data$District[1],
  month = data$month[1],
  year = data$year[1],
  Population_2022 = 1
)

pred_SO2 <- predict(model_SO2_final, newdata = pred_data_SO2)
coef_se_SO2 <- model_SO2_final$se["I(mean_SO2_ensemble/5)"]

SO2_er_df <- data.frame(
  SO2 = SO2_seq,
  RR = exp(pred_SO2),
  RR_L = exp(pred_SO2 - 1.96 * coef_se_SO2),
  RR_U = exp(pred_SO2 + 1.96 * coef_se_SO2)
)

ggplot(SO2_er_df, aes(x = SO2, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: SO2 and Death Risk",
    x = "SO2 (µg/m³)",
    y = "Risk Ratio (95% CI)"
  )


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

summary(model_Ox_final) # I(mean_O3_pred/5)                0.010250
(Ox_obs_coef <- coef(model_Ox_final)["I(mean_Ox/5)"])
(Ox_L <- Ox_obs_coef - 1.96*model_Ox_final$se["I(mean_Ox/5)"])
(Ox_U <- Ox_obs_coef + 1.96*model_Ox_final$se["I(mean_Ox/5)"])

# Exposure-response curve for Ox
Ox_seq <- seq(min(data$mean_Ox, na.rm = TRUE), 
              max(data$mean_Ox, na.rm = TRUE), 
              length.out = 100)

pred_data_Ox <- data.frame(
  mean_Ox = Ox_seq,
  mean_RH = mean(data$mean_RH, na.rm = TRUE),
  mean_wind_speed = mean(data$mean_wind_speed, na.rm = TRUE),
  mean_Precipitation = mean(data$mean_Precipitation, na.rm = TRUE),
  mean_Wind_Dir = mean(data$mean_Wind_Dir, na.rm = TRUE),
  mean_sp_hPa = mean(data$mean_sp_hPa, na.rm = TRUE),
  Division = data$Division[1],
  year = data$year[1],
  month = data$month[1],
  Population_2022 = 1
)

pred_Ox <- predict(model_Ox_final, newdata = pred_data_Ox)
coef_se_Ox <- model_Ox_final$se["I(mean_Ox/5)"]

Ox_er_df <- data.frame(
  Ox = Ox_seq,
  RR = exp(pred_Ox),
  RR_L = exp(pred_Ox - 1.96 * coef_se_Ox),
  RR_U = exp(pred_Ox + 1.96 * coef_se_Ox)
)

ggplot(Ox_er_df, aes(x = Ox, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: Ox and Death Risk",
    x = "Ox (ppb)",
    y = "Risk Ratio (95% CI)"
  )


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
summary(model_PM25_final) # 0.001566
(PM25_obs_coef <- coef(model_PM25_final)["I(mean_PM25_ensemble/5)"])
(PM25_L <- PM25_obs_coef - 1.96*model_PM25_final$se["I(mean_PM25_ensemble/5)"])
(PM25_U <- PM25_obs_coef + 1.96*model_PM25_final$se["I(mean_PM25_ensemble/5)"])

# Exposure-response curve for PM2.5
PM25_seq <- seq(min(data$mean_PM25_ensemble, na.rm = TRUE), 
                max(data$mean_PM25_ensemble, na.rm = TRUE), 
                length.out = 100)

pred_data_PM25 <- data.frame(
  mean_PM25_ensemble = PM25_seq,
  mean_t2m_c = mean(data$mean_t2m_c, na.rm = TRUE),
  mean_wind_speed = mean(data$mean_wind_speed, na.rm = TRUE),
  mean_Wind_Dir = mean(data$mean_Wind_Dir, na.rm = TRUE),
  mean_sp_hPa = mean(data$mean_sp_hPa, na.rm = TRUE),
  District = data$District[1],
  month = data$month[1],
  year = data$year[1],
  Population_2022 = 1
)

pred_PM25 <- predict(model_PM25_final, newdata = pred_data_PM25)
coef_se_PM25 <- model_PM25_final$se["I(mean_PM25_ensemble/5)"]

PM25_er_df <- data.frame(
  PM25 = PM25_seq,
  RR = exp(pred_PM25),
  RR_L = exp(pred_PM25 - 1.96 * coef_se_PM25),
  RR_U = exp(pred_PM25 + 1.96 * coef_se_PM25)
)

ggplot(PM25_er_df, aes(x = PM25, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: PM2.5 and Death Risk",
    x = "PM2.5 (µg/m³)",
    y = "Risk Ratio (95% CI)"
  )



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

summary(model_PM10_final) # 0.001566
(PM10_obs_coef <- coef(model_PM10_final)["I(mean_PM10ug/5)"])
(PM10_L <- PM10_obs_coef - 2*model_PM10_final$se["I(mean_PM10ug/5)"])
(PM10_U <- PM10_obs_coef + 2*model_PM10_final$se["I(mean_PM10ug/5)"])

# Exposure-response curve for PM10
PM10_seq <- seq(min(data$mean_PM10ug, na.rm = TRUE), 
                max(data$mean_PM10ug, na.rm = TRUE), 
                length.out = 100)

pred_data_PM10 <- data.frame(
  mean_PM10ug = PM10_seq,
  mean_t2m_c = mean(data$mean_t2m_c, na.rm = TRUE),
  mean_RH = mean(data$mean_RH, na.rm = TRUE),
  mean_sp_hPa = mean(data$mean_sp_hPa, na.rm = TRUE),
  District = data$District[1],
  month = data$month[1],
  Population_2022 = 1
)

pred_PM10 <- predict(model_PM10_final, newdata = pred_data_PM10)
coef_se_PM10 <- model_PM10_final$se["I(mean_PM10ug/5)"]

PM10_er_df <- data.frame(
  PM10 = PM10_seq,
  RR = exp(pred_PM10),
  RR_L = exp(pred_PM10 - 1.96 * coef_se_PM10),
  RR_U = exp(pred_PM10 + 1.96 * coef_se_PM10)
)

ggplot(PM10_er_df, aes(x = PM10, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: PM10 and Death Risk",
    x = "PM10 (µg/m³)",
    y = "Risk Ratio (95% CI)"
  )


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

summary(model_CO_final)
(CO_obs_coef <- coef(model_CO_final)["I(CO_pred/5)"])
(CO_L <- CO_obs_coef - 1.96*model_CO_final$se["I(CO_pred/5)"])
(CO_U <- CO_obs_coef + 1.96*model_CO_final$se["I(CO_pred/5)"])

# Exposure-response curve for CO
CO_seq <- seq(min(data$CO_pred, na.rm = TRUE), 
              max(data$CO_pred, na.rm = TRUE), 
              length.out = 100)

pred_data_CO <- data.frame(
  CO_pred = CO_seq,
  PM25_pred = mean(data$PM25_pred, na.rm = TRUE),
  mean_PM10_pred = mean(data$mean_PM10_pred, na.rm = TRUE),
  mean_NO2_pred = mean(data$mean_NO2_pred, na.rm = TRUE),
  mean_t2m_c = mean(data$mean_t2m_c, na.rm = TRUE),
  mean_RH = mean(data$mean_RH, na.rm = TRUE),
  mean_Precipitation = mean(data$mean_Precipitation, na.rm = TRUE),
  mean_wind_speed = mean(data$mean_wind_speed, na.rm = TRUE),
  District = data$District[1],
  year = data$year[1],
  month = data$month[1],
  Population_2022 = 1
)

pred_CO <- predict(model_CO_final, newdata = pred_data_CO)
coef_se_CO <- model_CO_final$se["I(CO_pred/5)"]

CO_er_df <- data.frame(
  CO = CO_seq,
  RR = exp(pred_CO),
  RR_L = exp(pred_CO - 1.96 * coef_se_CO),
  RR_U = exp(pred_CO + 1.96 * coef_se_CO)
)

ggplot(CO_er_df, aes(x = CO, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: CO and Death Risk",
    x = "CO (ppm)",
    y = "Risk Ratio (95% CI)"
  )


####################
### 8.O3
####################
model_O3_final <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal/5) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = ~log(Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

summary(model_O3_final) # I(mean_O3_pred/5)                0.010250
(O3_obs_coef <- coef(model_O3_final)["I(mean_O3_pred_cal/5)"])
(O3_L <- O3_obs_coef - 1.96*model_O3_final$se["I(mean_O3_pred_cal/5)"])
(O3_U <- O3_obs_coef + 1.96*model_O3_final$se["I(mean_O3_pred_cal/5)"])

# Exposure-response curve for O3
O3_seq <- seq(min(data$mean_O3_pred_cal, na.rm = TRUE), 
              max(data$mean_O3_pred_cal, na.rm = TRUE), 
              length.out = 100)

pred_data_O3 <- data.frame(
  mean_O3_pred_cal = O3_seq,
  mean_t2m_c = mean(data$mean_t2m_c, na.rm = TRUE),
  mean_RH = mean(data$mean_RH, na.rm = TRUE),
  mean_Precipitation = mean(data$mean_Precipitation, na.rm = TRUE),
  mean_Wind_Dir = mean(data$mean_Wind_Dir, na.rm = TRUE),
  mean_sp_hPa = mean(data$mean_sp_hPa, na.rm = TRUE),
  District = data$District[1],
  year = data$year[1],
  month = data$month[1],
  Population_2022 = 1
)

pred_O3 <- predict(model_O3_final, newdata = pred_data_O3)
coef_se_O3 <- model_O3_final$se["I(mean_O3_pred_cal/5)"]

O3_er_df <- data.frame(
  O3 = O3_seq,
  RR = exp(pred_O3),
  RR_L = exp(pred_O3 - 1.96 * coef_se_O3),
  RR_U = exp(pred_O3 + 1.96 * coef_se_O3)
)

ggplot(O3_er_df, aes(x = O3, y = RR)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = RR_L, ymax = RR_U), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Exposure-Response Curve: O3 and Death Risk",
    x = "O3 (ppb)",
    y = "Risk Ratio (95% CI)"
  )




library(ggplot2)





library(ggplot2)

coef_df <- data.frame(
  Pollutant = c("NO2", "CO", "O3", "PM2.5", "SO2", "NO", "PM10", "Ox"),
  Estimate = c(
    coef(model_NO2_final)["I(mean_NO2_ensemble/5)"],
    coef(model_CO_final)["I(CO_pred/5)"],
    coef(model_O3_final)["I(mean_O3_pred_cal/5)"],
    coef(model_PM25_final)["I(mean_PM25_ensemble/5)"],
    coef(model_SO2_final)["I(mean_SO2_ensemble/5)"],
    coef(model_NO_final)["I(mean_NOppb/5)"],
    coef(model_PM10_final)["I(mean_PM10ug/5)"],
    coef(model_Ox_final)["I(mean_Ox/5)"]
  ),
  SE = c(
    model_NO2_final$se["I(mean_NO2_ensemble/5)"],
    model_CO_final$se["I(CO_pred/5)"],
    model_O3_final$se["I(mean_O3_pred_cal/5)"],
    model_PM25_final$se["I(mean_PM25_ensemble/5)"],
    model_SO2_final$se["I(mean_SO2_ensemble/5)"],
    model_NO_final$se["I(mean_NOppb/5)"],
    model_PM10_final$se["I(mean_PM10ug/5)"],
    model_Ox_final$se["I(mean_Ox/5)"]
  )
)


# Compute 95% CI on log-scale
coef_df$Lower <- coef_df$Estimate - 1.96 * coef_df$SE
coef_df$Upper <- coef_df$Estimate + 1.96 * coef_df$SE

# Exponentiate for Risk Ratios
coef_df$RR    <- exp(coef_df$Estimate)
coef_df$RR_L  <- exp(coef_df$Lower)
coef_df$RR_U  <- exp(coef_df$Upper)


# Plot
ggplot(coef_df, aes(x = Pollutant, y = Estimate)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 14) +
  labs(x = "Pollutant", y = "Coefficient Estimate (95% CI)")

# Plot (RR with Pollutant on y-axis)
ggplot(coef_df, aes(x = Pollutant, y = RR)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.15, color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(y = "Risk Ratio (95% CI) for Death", x = "Pollutant")

