### reference model and alternative specifications
library(dplyr)
library(fixest)
library(ggplot2)
library(splines)
library(patchwork)
library(zoo)
library(tidyr)
data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")

names(data)
res <- cbind(data$AggregatedDeath, data$mean_NO2_ensemble, data$mean_NOppb, data$mean_SO2_ensemble, data$mean_Ox, data$mean_PM25_ensemble, data$mean_PM10ug, data$mean_O3_pred_cal, data$CO_pred,  data$year)
summary(res)
####################
### Trend Figures for All 8 Pollutants - Combined in One Figure
####################

# NO2 Trend
no2_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_NO2_ensemble, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

no2_trend[no2_trend$year == 2018, ]
NO22024 <- no2_trend[no2_trend$year == 2024, ]
NO22024$month <- factor(NO22024$month, levels = 1:12, labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
p_no2 <- ggplot(no2_trend, aes(x = month, y = mean_val)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  labs(title = "NO2", x = "Year", y = "NO2 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# NO Trend
no_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_NOppb, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

p_no <- ggplot(no_trend, aes(x = date, y = MA12)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  labs(title = "NO", x = "Year", y = "NO (ppb)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# SO2 Trend
so2_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_SO2_ensemble, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

p_so2 <- ggplot(so2_trend, aes(x = date, y = MA12)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  labs(title = "SO2", x = "Year", y = "SO2 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Ox Trend
ox_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_Ox, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

p_ox <- ggplot(ox_trend, aes(x = date, y = MA12)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  labs(title = "Ox", x = "Year", y = "Ox (ppb)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# PM25 Trend
pm25_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_PM25_ensemble, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

# PM2.5 trend by year (monthly progression within each year)
p_pm25 <- ggplot(pm25_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "PM2.5", x = "Month", y = "PM2.5 (µg/m³)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))


# PM2.5 trend by year - faceted plot
p_pm25_faceted <- ggplot(pm25_trend, aes(x = month, y = mean_val, fill = factor(year))) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  facet_wrap(~year, ncol = 3) +
  labs(title = "PM2.5 Trend by Year (Monthly Progression)", x = "Month", y = "PM2.5 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "none")

p_pm25_faceted

# Aggregate PM2.5 by year (sum of monthly values)
pm25_by_year <- pm25_trend |>
  group_by(year) |>
  summarize(
    total_pm25 = sum(mean_val, na.rm = TRUE),
    mean_pm25 = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

pm25_by_year



# Aggregate PM2.5 by year (sum of monthly values)
pm25_by_year <- pm25_trend |>
  group_by(year) |>
  summarize(
    total_pm25 = sum(mean_val, na.rm = TRUE),
    mean_pm25 = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

pm25_by_year


# PM10 Trend
pm10_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_PM10ug, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

p_pm10 <- ggplot(pm10_trend, aes(x = date, y = MA12)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  labs(title = "PM10", x = "Year", y = "PM10 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# O3 Trend
o3_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_O3_pred_cal, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

p_o3 <- ggplot(o3_trend, aes(x = date, y = MA12)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  labs(title = "O3", x = "Year", y = "O3 (ppb)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# CO Trend
co_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(CO_pred, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_val, k = 12, fill = NA, align = "right")
  )

p_co <- ggplot(co_trend, aes(x = date, y = MA12)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.2, alpha = 0.6) +
  labs(title = "CO", x = "Year", y = "CO (ppm)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Combine all plots in one figure
combined_trend_plot <- (p_no2 + p_no + p_so2 + p_ox) / (p_pm25 + p_pm10 + p_o3 + p_co) +
  plot_annotation(title = "Trend of 12-Month Moving Average Concentrations for All Pollutants",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

combined_trend_plot

####################
### Average Exposure Distribution by Districts and Pollutants
####################

# Calculate average exposure by district for all 8 pollutants
exposure_by_district <- data |>
  group_by(District) |>
  summarize(
    NO2 = mean(mean_NO2_ensemble, na.rm = TRUE),
    NO = mean(mean_NOppb, na.rm = TRUE),
    SO2 = mean(mean_SO2_ensemble, na.rm = TRUE),
    Ox = mean(mean_Ox, na.rm = TRUE),
    PM25 = mean(mean_PM25_ensemble, na.rm = TRUE),
    PM10 = mean(mean_PM10ug, na.rm = TRUE),
    O3 = mean(mean_O3_pred_cal, na.rm = TRUE),
    CO = mean(CO_pred, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(cols = -District, names_to = "Pollutant", values_to = "Avg_Exposure")

# Create boxplot for average exposure distribution by district
p_exposure_no2 <- exposure_by_district |>
  filter(Pollutant == "NO2") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "NO2", x = "District", y = "NO2 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_exposure_no <- exposure_by_district |>
  filter(Pollutant == "NO") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "NO", x = "District", y = "NO (ppb)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_exposure_so2 <- exposure_by_district |>
  filter(Pollutant == "SO2") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "SO2", x = "District", y = "SO2 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_exposure_ox <- exposure_by_district |>
  filter(Pollutant == "Ox") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "Ox", x = "District", y = "Ox (ppb)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_exposure_pm25 <- exposure_by_district |>
  filter(Pollutant == "PM25") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "PM2.5", x = "District", y = "PM2.5 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_exposure_pm10 <- exposure_by_district |>
  filter(Pollutant == "PM10") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "PM10", x = "District", y = "PM10 (µg/m³)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_exposure_o3 <- exposure_by_district |>
  filter(Pollutant == "O3") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "O3", x = "District", y = "O3 (ppb)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_exposure_co <- exposure_by_district |>
  filter(Pollutant == "CO") |>
  ggplot(aes(x = reorder(District, Avg_Exposure, FUN = median), y = Avg_Exposure)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  labs(title = "CO", x = "District", y = "CO (ppm)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Combine exposure distribution plots
combined_exposure_plot <- (p_exposure_no2 + p_exposure_no + p_exposure_so2 + p_exposure_ox) / 
  (p_exposure_pm25 + p_exposure_pm10 + p_exposure_o3 + p_exposure_co) +
  plot_annotation(title = "Average Exposure Distribution by Districts and Pollutants",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

combined_exposure_plot

####################
### Average Monthly Mortality
####################

# Calculate average monthly mortality across all districts
monthly_mortality <- data |>
  group_by(year, month) |>
  summarize(
    mean_deaths = mean(AggregatedDeath, na.rm = TRUE),
    total_deaths = sum(AggregatedDeath, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(year, month) |>
  mutate(
    date = as.Date(paste(year, month, "15", sep = "-")),
    MA12 = zoo::rollmean(mean_deaths, k = 12, fill = NA, align = "right")
  )

# Plot average monthly mortality with 12-month moving average
p_mortality_mean <- ggplot(monthly_mortality, aes(x = date, y = mean_deaths)) +
  geom_line(color = "darkred", linewidth = 0.8, alpha = 0.5) +
  geom_point(color = "darkred", size = 1.2, alpha = 0.4) +
  geom_line(aes(y = MA12), color = "darkred", linewidth = 1.2) +
  labs(
    title = "Average Monthly Mortality (Mean Deaths per District)",
    x = "Year",
    y = "Mean Deaths per District"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot total monthly mortality with 12-month moving average
total_ma12 <- zoo::rollmean(monthly_mortality$total_deaths, k = 12, fill = NA, align = "right")

p_mortality_total <- ggplot(monthly_mortality, aes(x = date, y = total_deaths)) +
  geom_line(color = "darkblue", linewidth = 0.8, alpha = 0.5) +
  geom_point(color = "darkblue", size = 1.2, alpha = 0.4) +
  geom_line(aes(y = total_ma12), color = "darkblue", linewidth = 1.2) +
  labs(
    title = "Total Monthly Mortality (Sum of All Districts)",
    x = "Year",
    y = "Total Deaths"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine mortality plots
combined_mortality_plot <- p_mortality_mean / p_mortality_total +
  plot_annotation(
    title = "Average Monthly Mortality Trends",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

combined_mortality_plot


