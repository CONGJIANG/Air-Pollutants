library(dplyr)
library(fixest)
library(ggplot2)
library(splines)
library(zoo)
library(tidyr)
library(readxl)
library(patchwork)

data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")

####################
### 1. NO2
####################
model_NO2_final <- feglm(
  AggregatedDeath ~ I(mean_NO2_ensemble/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)


summary(model_NO2_final) #I(mean_NO2_ensemble/5)           0.051140   0.026509  1.929193 0.05376838 . 
# I(mean_NO2_ensemble/5)           0.110414   0.028574  3.864131 1.1302e-04 ***
(NO2_obs_coef <- coef(model_NO2_final)["I(mean_NO2_ensemble/10)"])
(NO2_L <- NO2_obs_coef - 1.96*model_NO2_final$se["I(mean_NO2_ensemble/10)"])
(NO2_U <- NO2_obs_coef + 1.96*model_NO2_final$se["I(mean_NO2_ensemble/10)"])

####################
### 2. NO
####################
data$mean_NOug <- data$mean_NOppb* 1.23  # Convert ppb to μg/m³
model_NO_final <- feglm(
  AggregatedDeath ~ I(mean_NOug/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

summary(model_NO_final)

(NO_obs_coef <- coef(model_NO_final)["I(mean_NOug/10)"])
(NO_L <- NO_obs_coef - 1.96*model_NO_final$se["I(mean_NOug/10)"])
(NO_U <- NO_obs_coef + 1.96*model_NO_final$se["I(mean_NOug/10)"])


####################
### 3. SO2
####################
model_SO2_final <- feglm(
  AggregatedDeath ~ I(mean_SO2_ensemble/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

summary(model_SO2_final)

(SO2_obs_coef <- coef(model_SO2_final)["I(mean_SO2_ensemble/10)"])
(SO2_L <- SO2_obs_coef - 1.96*model_SO2_final$se["I(mean_SO2_ensemble/10)"])
(SO2_U <- SO2_obs_coef + 1.96*model_SO2_final$se["I(mean_SO2_ensemble/10)"])


####################
### 4. Ox
####################

model_Ox_final <- feglm(
  AggregatedDeath ~ I(mean_Ox/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

summary(model_Ox_final) # I(mean_O3_pred/5)                0.010250
(Ox_obs_coef <- coef(model_Ox_final)["I(mean_Ox/10)"])
(Ox_L <- Ox_obs_coef - 1.96*model_Ox_final$se["I(mean_Ox/10)"])
(Ox_U <- Ox_obs_coef + 1.96*model_Ox_final$se["I(mean_Ox/10)"])


####################
### 5. PM25
####################
model_PM25_final <- feglm(
  AggregatedDeath ~ I(mean_PM25_ensemble/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| District^month + year,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)
summary(model_PM25_final) # 0.001566
(PM25_obs_coef <- coef(model_PM25_final)["I(mean_PM25_ensemble/10)"])
(PM25_L <- PM25_obs_coef - 1.96*model_PM25_final$se["I(mean_PM25_ensemble/10)"])
(PM25_U <- PM25_obs_coef + 1.96*model_PM25_final$se["I(mean_PM25_ensemble/10)"])

(exp(PM25_obs_coef)*100 -100)
summary(data$mean_PM25_ensemble)

exp(0.001696)*88
exp(0.001696)*75.19
exp(0.001696)*12.55
exp(0.001696)*48.50
exp(0.001696)*125.66
exp(0.001696)*265.04
####################
### 6. PM10
####################
model_PM10_final <- feglm(
  AggregatedDeath ~ I(mean_PM10ug/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

summary(model_PM10_final) # 0.001566
(PM10_obs_coef <- coef(model_PM10_final)["I(mean_PM10ug/10)"])
(PM10_L <- PM10_obs_coef - 1.96*model_PM10_final$se["I(mean_PM10ug/10)"])
(PM10_U <- PM10_obs_coef + 1.96*model_PM10_final$se["I(mean_PM10ug/10)"])


####################
### 7. CO
####################
model_CO_final <- feglm(
  AggregatedDeath ~ I(CO_pred/10) + mean_t2m_c +ns(mean_wind_speed, df=5) | District^year + month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

summary(model_CO_final)
(CO_obs_coef <- coef(model_CO_final)["I(CO_pred/10)"])
(CO_L <- CO_obs_coef - 1.96*model_CO_final$se["I(CO_pred/10)"])
(CO_U <- CO_obs_coef + 1.96*model_CO_final$se["I(CO_pred/10)"])


####################
### 8.O3
####################
model_O3_final <- feglm(
  AggregatedDeath ~ I(mean_O3_pred_cal/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District^year+month,
  data = data,
  offset = log(data$Population_2022),
  vcov = "iid",
  family = "quasipoisson"
)

summary(model_O3_final) # I(mean_O3_pred/5)                0.010250
(O3_obs_coef <- coef(model_O3_final)["I(mean_O3_pred_cal/10)"])
(O3_L <- O3_obs_coef - 1.96*model_O3_final$se["I(mean_O3_pred_cal/10)"])
(O3_U <- O3_obs_coef + 1.96*model_O3_final$se["I(mean_O3_pred_cal/10)"])




library(ggplot2)

coef_df <- data.frame(
  Pollutant = factor(c("O3", "CO", "Ox", "SO2", "NO", "NO2", "PM10", "PM2.5"), 
                     levels = c("O3", "CO", "Ox", "SO2", "NO", "NO2", "PM10", "PM2.5")),
  Estimate = c(
    coef(model_O3_final)["I(mean_O3_pred_cal/10)"],
    coef(model_CO_final)["I(CO_pred/10)"],
    coef(model_Ox_final)["I(mean_Ox/10)"],
    coef(model_SO2_final)["I(mean_SO2_ensemble/10)"],
    coef(model_NO_final)["I(mean_NOug/10)"],
    coef(model_NO2_final)["I(mean_NO2_ensemble/10)"],
    coef(model_PM10_final)["I(mean_PM10ug/10)"],
    coef(model_PM25_final)["I(mean_PM25_ensemble/10)"]
  ),
  SE = c(
    model_O3_final$se["I(mean_O3_pred_cal/10)"],
    model_CO_final$se["I(CO_pred/10)"],
    model_Ox_final$se["I(mean_Ox/10)"],
    model_SO2_final$se["I(mean_SO2_ensemble/10)"],
    model_NO_final$se["I(mean_NOug/10)"],
    model_NO2_final$se["I(mean_NO2_ensemble/10)"],
    model_PM10_final$se["I(mean_PM10ug/10)"],
    model_PM25_final$se["I(mean_PM25_ensemble/10)"]
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
  geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.32, color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_text(aes(label = paste0(sprintf("%.2f", RR), " (", sprintf("%.2f", RR_L), "-", sprintf("%.2f", RR_U), ")")), 
            vjust = -1.5, hjust = 0.0, size = 3.2) +
  coord_flip(ylim = c(0.9999999999, 1.652)) +
  theme_minimal(base_size = 14) +
  labs(y = "Risk Ratio (95% CI) for Death", x = "Pollutant")

#################################################################################
#################################################################################
# ABSOLUTE MORTALITY INCREASE PER 100,000 POPULATION PER MONTH
#################################################################################
#################################################################################

cat("\n========== ABSOLUTE MORTALITY INCREASE CALCULATIONS ==========\n")

# Calculate baseline monthly mortality rate per 100,000 people
baseline_rate_100k <- with(data,
  mean(AggregatedDeath / Population_2022, na.rm = TRUE) * 1e5
)

cat("Baseline monthly mortality rate:", sprintf("%.2f", baseline_rate_100k), "per 100,000 population\n\n")

# Create dataframe for absolute mortality increase
absolute_mortality_df <- data.frame(
  Pollutant = coef_df$Pollutant,
  Coef = coef_df$Estimate,
  SE = coef_df$SE,
  RR = coef_df$RR,
  RR_L = coef_df$RR_L,
  RR_U = coef_df$RR_U
)

# Calculate absolute mortality increase per 100,000 for 10-unit increase in pollutant
absolute_mortality_df$Delta_Deaths_100k <- 
  baseline_rate_100k * (absolute_mortality_df$RR - 1)

absolute_mortality_df$Delta_Deaths_100k_L <- 
  baseline_rate_100k * (absolute_mortality_df$RR_L - 1)

absolute_mortality_df$Delta_Deaths_100k_U <- 
  baseline_rate_100k * (absolute_mortality_df$RR_U - 1)

cat("Absolute mortality increase per 100,000 population per month:\n")
print(absolute_mortality_df[, c("Pollutant", "RR", "Delta_Deaths_100k", "Delta_Deaths_100k_L", "Delta_Deaths_100k_U")])

# Create visualization for absolute mortality increase
plot_absolute_mortality <- ggplot(absolute_mortality_df, aes(x = Pollutant, y = Delta_Deaths_100k)) +
  geom_point(size = 3, color = "darkred") +
  geom_errorbar(aes(ymin = Delta_Deaths_100k_L, ymax = Delta_Deaths_100k_U), 
                width = 0.32, color = "darkred") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_text(aes(label = paste0(sprintf("%.2f", Delta_Deaths_100k), 
                               " (", sprintf("%.2f", Delta_Deaths_100k_L), 
                               "-", sprintf("%.2f", Delta_Deaths_100k_U), ")")), 
            vjust = -1.2, hjust = 0.0, size = 3.0) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Absolute Mortality Increase per 100,000 Population per Month",
    subtitle = "Associated with 10-unit increase in pollutant concentration",
    y = "Additional Deaths per 100,000 per Month (95% CI)",
    x = "Pollutant"
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor.x = element_blank()
  )

cat("\n========== FOREST PLOT: ABSOLUTE MORTALITY INCREASE ==========\n")
plot_absolute_mortality

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Load original data with sex and age information
data_full <- read_excel("/Users/cjiang/Downloads/Aggregated data for mortality.xlsx")

# Prepare exposure data: aggregate by month/year (same exposure for all sex/age groups in same month)
exposure_data <- data |>
  select(year, month, District,
         mean_NO2_ensemble, mean_NOppb, mean_SO2_ensemble, mean_Ox, mean_PM25_ensemble, 
         mean_PM10ug, mean_CO_pred, CO_pred, mean_O3_pred_cal, 
         mean_t2m_c, mean_RH, mean_Precipitation, mean_wind_speed, 
         mean_Wind_Dir, mean_sp_hPa, Population_2022) |>
  distinct()
exposure_data$mean_NOug <- exposure_data$mean_NOppb * 1.23  # Convert ppb to μg/m³
# Prepare mortality data: pivot to long format with age groups
data_subgroup <- data_full |>
  mutate(
    year = Year,
    month = Month,
    District = Distirct  # Fix typo in Excel column name
  ) |>
  select(year, month, District, Division, Sex, Age0to5, Age6to14, Age15to49, `Age50+`) |>
  pivot_longer(
    cols = c(Age0to5, Age6to14, Age15to49, `Age50+`),
    names_to = "Age_Category",
    values_to = "Deaths"
  ) |>
  mutate(
    Age_Group = case_when(
      Age_Category %in% c("Age0to5", "Age6to14", "Age15to49") ~ "Age <50",
      Age_Category == "Age50+" ~ "Age 50+",
      TRUE ~ NA_character_
    )
  ) |>
  select(year, month, District, Division, Sex, Age_Group, Deaths) |>
  group_by(year, month, District, Division, Sex, Age_Group) |>
  summarise(Outcome = sum(Deaths, na.rm = TRUE), .groups = "drop") |>
  # Merge with exposure data by month/year/district
  left_join(
    exposure_data,
    by = c("year", "month", "District")
  ) |>
  filter(!is.na(Sex), !is.na(Age_Group), !is.na(District), !is.na(Outcome), 
         !is.na(Population_2022), !is.na(mean_NO2_ensemble))

cat("Data subgroup created with", nrow(data_subgroup), "rows (after removing districts without RDS data)\n")
cat("Unique sex values:", unique(data_subgroup$Sex), "\n")
cat("Unique age groups:", unique(data_subgroup$Age_Group), "\n")
cat("Sample rows:\n")
print(head(data_subgroup[, c("year", "month", "District", "Sex", "Age_Group", "Outcome", "mean_NO2_ensemble")], 10))
# Function to fit subgroup models
fit_subgroup_models <- function(formula_str, data, pollutant_name, divisor = 10) {
  results <- list()
  
  # By Sex
  for (sex in unique(data$Sex)) {
    if (!is.na(sex)) {
      data_subset <- data |> filter(Sex == sex)
      cat("  Fitting", pollutant_name, "for Sex =", sex, "with", nrow(data_subset), "rows\n")
      
      model <- tryCatch({
        feglm(as.formula(formula_str), data = data_subset, 
              offset = log(data_subset$Population_2022), 
              vcov = "iid", family = "quasipoisson")
      }, error = function(e) {
        cat("    Error:", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(model)) {
        # Extract coefficient - try multiple patterns
        coef_name <- NA
        for (var in names(coef(model))) {
          if (grepl(paste0("/", divisor), var)) {
            coef_name <- var
            break
          }
        }
        
        if (!is.na(coef_name)) {
          results[[paste0(pollutant_name, "_Sex_", sex)]] <- list(
            RR = exp(coef(model)[coef_name]),
            L = exp(coef(model)[coef_name] - 1.96 * model$se[coef_name]),
            U = exp(coef(model)[coef_name] + 1.96 * model$se[coef_name]),
            Group = paste0(pollutant_name, " - ", sex),
            Category = "Sex"
          )
          cat("    Success! RR =", sprintf("%.3f", exp(coef(model)[coef_name])), "\n")
        }
      }
    }
  }
  
  # By Age Group
  age_groups <- data |> pull(Age_Group) |> unique()
  
  for (age_grp in age_groups) {
    if (!is.na(age_grp)) {
      data_subset <- data |> filter(Age_Group == age_grp)
      cat("  Fitting", pollutant_name, "for Age_Group =", age_grp, "with", nrow(data_subset), "rows\n")
      
      model <- tryCatch({
        feglm(as.formula(formula_str), data = data_subset, 
              offset = log(data_subset$Population_2022), 
              vcov = "iid", family = "quasipoisson")
      }, error = function(e) {
        cat("    Error:", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(model)) {
        coef_name <- NA
        for (var in names(coef(model))) {
          if (grepl(paste0("/", divisor), var)) {
            coef_name <- var
            break
          }
        }
        
        if (!is.na(coef_name)) {
          results[[paste0(pollutant_name, "_Age_", age_grp)]] <- list(
            RR = exp(coef(model)[coef_name]),
            L = exp(coef(model)[coef_name] - 1.96 * model$se[coef_name]),
            U = exp(coef(model)[coef_name] + 1.96 * model$se[coef_name]),
            Group = paste0(pollutant_name, " - ", age_grp),
            Category = "Age Group"
          )
          cat("    Success! RR =", sprintf("%.3f", exp(coef(model)[coef_name])), "\n")
        }
      }
    }
  }
  
  return(results)
}

formulas <- list(
  NO2 = "Outcome ~ I(mean_NO2_ensemble/10) +ns(mean_t2m_c, df=4) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| Division + month",
  NO = "Outcome ~ I(mean_NOug/10) +ns(mean_t2m_c, df=4) +mean_RH +mean_wind_speed| Division + month",
  SO2 = "Outcome ~ I(mean_SO2_ensemble/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5)+ns(mean_Precipitation, df=2)| District^month",
  Ox = "Outcome ~ I(mean_Ox/10) +ns(mean_RH, df=5) +ns(mean_wind_speed, df=3)+ns(mean_Precipitation, df=5)+ns(mean_Wind_Dir, df=4)+ns(mean_sp_hPa, df=2)| Division + year + month",
  PM25 = "Outcome ~ I(mean_PM25_ensemble/10) +mean_t2m_c  +ns(mean_wind_speed, df=6)+ns(mean_Wind_Dir, df=2)+ns(mean_sp_hPa, df=5)| Division + month",
  PM10 = "Outcome ~ I(mean_PM10ug/10) +mean_t2m_c +mean_RH +ns(mean_sp_hPa, df=2)| District + month",
  CO = "Outcome ~ I(CO_pred/10) + mean_t2m_c +ns(mean_wind_speed, df=2) | Division + year",
  O3 = "Outcome ~ I(mean_O3_pred_cal/10) +ns(mean_t2m_c, df=2) +ns(mean_RH, df=5) +ns(mean_Precipitation, df=4)+ns(mean_Wind_Dir, df=3)+mean_sp_hPa| District + year + month + District^year"
)

# Fit subgroup models for each pollutant using merged data
subgroup_results <- list()
for (pol in names(formulas)) {
  cat("Fitting models for", pol, "\n")
  results <- fit_subgroup_models(formulas[[pol]], data_subgroup, pol)
  cat("  Total: Got", length(results), "subgroups\n\n")
  subgroup_results <- c(subgroup_results, results)
}

cat("\n========== SUBGROUP ANALYSIS COMPLETE ==========\n")
cat("Total subgroup results:", length(subgroup_results), "\n")

  # Create dataframe for plotting
  subgroup_df <- do.call(rbind, lapply(names(subgroup_results), function(name) {
    res <- subgroup_results[[name]]
    data.frame(
      Group = res$Group,
      Category = res$Category,
      RR = as.numeric(res$RR),
      RR_L = as.numeric(res$L),
      RR_U = as.numeric(res$U),
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(subgroup_df) <- NULL
  
  cat("\nSubgroup dataframe dimensions:", nrow(subgroup_df), "rows x", ncol(subgroup_df), "columns\n")
  print(head(subgroup_df, 10))
  
  # Extract pollutant name from Group (e.g., "NO2 - Male" -> "NO2")
  subgroup_df$Pollutant <- sub(" - .*", "", subgroup_df$Group)
  
  # Reverse sex and age group labels
  subgroup_df$Group <- gsub("Female", "TEMP_FEMALE", subgroup_df$Group)
  subgroup_df$Group <- gsub("Male", "Female", subgroup_df$Group)
  subgroup_df$Group <- gsub("TEMP_FEMALE", "Male", subgroup_df$Group)
  
  subgroup_df$Group <- gsub("Age 50\\+", "TEMP_AGE50", subgroup_df$Group)
  subgroup_df$Group <- gsub("Age <50", "Age 50+", subgroup_df$Group)
  subgroup_df$Group <- gsub("TEMP_AGE50", "Age <50", subgroup_df$Group)
  
  # Create forest plots: one for sex, one for age groups
  
  # Plot 1: Sex stratification (Male and Female together)
  sex_data <- subgroup_df |>
    filter(Category == "Sex") |>
    mutate(
      Pollutant = factor(Pollutant, levels = c("O3", "CO", "Ox", "SO2", "NO", "NO2", "PM10", "PM25")),
      Sex = sub(".*- ", "", Group)
    )
  
  plot_sex <- ggplot(sex_data, aes(x = Pollutant, y = RR, color = Sex, shape = Sex)) +
    geom_point(size = 3, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, position = position_dodge(width = 0.6)) +
    geom_text(aes(label = paste0(sprintf("%.3f", RR), " (", sprintf("%.3f", RR_L), "-", sprintf("%.3f", RR_U), ")")), 
              position = position_dodge(width = 0.6), vjust = 0.2, hjust = -0.15, size = 2.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
    coord_flip() +
    scale_color_manual(values = c("Female" = "#1f77b4", "Male" = "#ff7f0e")) +
    scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
    theme_minimal(base_size = 13) +
    labs(
      title = "Sex Stratification",
      y = "Risk Ratio (95% CI) for Death",
      x = "Pollutant",
      color = "Sex",
      shape = "Sex"
    ) +
    theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold"),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      legend.position = "bottom"
    )
  
  # Plot 2: Age stratification (Age <50 and Age 50+ together)
  age_data <- subgroup_df |>
    filter(Category == "Age Group") |>
    mutate(
      Pollutant = factor(Pollutant, levels = c("O3", "CO", "Ox", "SO2", "NO", "NO2", "PM10", "PM25")),
      Age = sub(".*- ", "", Group)
    )
  
  plot_age <- ggplot(age_data, aes(x = Pollutant, y = RR, color = Age, shape = Age)) +
    geom_point(size = 3, position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = RR_L, ymax = RR_U), width = 0.2, position = position_dodge(width = 0.6)) +
    geom_text(aes(label = paste0(sprintf("%.3f", RR), " (", sprintf("%.3f", RR_L), "-", sprintf("%.3f", RR_U), ")")), 
              position = position_dodge(width = 0.6), vjust = 0.2, hjust = -0.15, size = 2.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
    coord_flip() +
    scale_color_manual(values = c("Age 50+" = "#1f77b4", "Age <50" = "#ff7f0e")) +
    scale_shape_manual(values = c("Age 50+" = 16, "Age <50" = 17)) +
    theme_minimal(base_size = 13) +
    labs(
      title = "Age Stratification",
      y = "Risk Ratio (95% CI) for Death",
      x = "Pollutant",
      color = "Age Group",
      shape = "Age Group"
    ) +
    theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold"),
      panel.grid.major.x = element_line(color = "gray90"),
      panel.grid.minor.x = element_blank(),
      legend.position = "bottom"
    )
  
  # Display both plots side by side
  cat("\n========== FOREST PLOTS BY STRATIFICATION ==========\n")
  cat("Displaying 2 forest plots: Sex and Age Groups\n\n")
plot_sex
plot_age
  # Combine both plots
  library(patchwork)
  combined_subgroup_plot <- plot_sex + plot_age +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Risk Ratios by Stratification for All Pollutants",
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 12))
    )
  
  combined_subgroup_plot


