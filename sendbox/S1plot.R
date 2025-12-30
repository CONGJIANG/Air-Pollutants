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

# Helper function to convert month numbers to month names
convert_months <- function(df) {
  df |>
    mutate(
      month = factor(month, levels = 1:12, 
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
    )
}

####################
### Trend Figures for All 8 Pollutants - Combined in One Figure
####################

# NO2 Trend
no2_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_NO2_ensemble, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  convert_months()

# NO2 trend by year (monthly progression within each year)
p_no2 <- ggplot(no2_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "NO2", x = "Month", y = "NO2 (ppb)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Aggregate NO2 by year (sum of monthly values)
no2_by_year <- no2_trend |>
  group_by(year) |>
  summarize(
    total_no2 = sum(mean_val, na.rm = TRUE),
    mean_no2 = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

no2_by_year

# NO Trend
no_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_NOppb, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  convert_months()

# NO trend by year (monthly progression within each year)
p_no <- ggplot(no_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "NO", x = "Month", y = "NO (ppb)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Aggregate NO by year (sum of monthly values)
no_by_year <- no_trend |>
  group_by(year) |>
  summarize(
    total_no = sum(mean_val, na.rm = TRUE),
    mean_no = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

no_by_year

# SO2 Trend
so2_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_SO2_ensemble, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  convert_months()

# SO2 trend by year (monthly progression within each year)
p_so2 <- ggplot(so2_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "SO2", x = "Month", y = "SO2 (ppb)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Aggregate SO2 by year (sum of monthly values)
so2_by_year <- so2_trend |>
  group_by(year) |>
  summarize(
    total_so2 = sum(mean_val, na.rm = TRUE),
    mean_so2 = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

so2_by_year

# Ox Trend
ox_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_Ox, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  convert_months()

# Ox trend by year (monthly progression within each year)
p_ox <- ggplot(ox_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "Ox", x = "Month", y = "Ox (ppb)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Aggregate Ox by year (sum of monthly values)
ox_by_year <- ox_trend |>
  group_by(year) |>
  summarize(
    total_ox = sum(mean_val, na.rm = TRUE),
    mean_ox = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

ox_by_year

# PM25 Trend
pm25_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_PM25_ensemble, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  convert_months() 

# PM2.5 trend by year (monthly progression within each year)
p_pm25 <- ggplot(pm25_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "PM2.5", x = "Month", y = "PM2.5 (µg/m³)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))


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
  convert_months()

# PM10 trend by year (monthly progression within each year)
p_pm10 <- ggplot(pm10_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "PM10", x = "Month", y = "PM10 (µg/m³)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Aggregate PM10 by year (sum of monthly values)
pm10_by_year <- pm10_trend |>
  group_by(year) |>
  summarize(
    total_pm10 = sum(mean_val, na.rm = TRUE),
    mean_pm10 = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

pm10_by_year

# O3 Trend
o3_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(mean_O3_pred_cal, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  convert_months()

# O3 trend by year (monthly progression within each year)
p_o3 <- ggplot(o3_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "O3", x = "Month", y = "O3 (ppb)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Aggregate O3 by year (sum of monthly values)
o3_by_year <- o3_trend |>
  group_by(year) |>
  summarize(
    total_o3 = sum(mean_val, na.rm = TRUE),
    mean_o3 = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

o3_by_year

# CO Trend
co_trend <- data |>
  group_by(year, month) |>
  summarize(mean_val = mean(CO_pred, na.rm = TRUE), .groups = "drop") |>
  arrange(year, month) |>
  convert_months()

# CO trend by year (monthly progression within each year)
p_co <- ggplot(co_trend, aes(x = month, y = mean_val, color = factor(year), group = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(title = "CO", x = "Month", y = "CO (ppm)", color = "Year") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Aggregate CO by year (sum of monthly values)
co_by_year <- co_trend |>
  group_by(year) |>
  summarize(
    total_co = sum(mean_val, na.rm = TRUE),
    mean_co = mean(mean_val, na.rm = TRUE),
    .groups = "drop"
  )

co_by_year

# Combine all plots in one figure
combined_trend_plot1 <- (p_so2 + p_ox) / (p_o3 + p_co) +
  plot_annotation(title = "Trend of 12-Month Moving Average Concentrations for All Pollutants",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

combined_trend_plot1



# Combine all plots in one figure
combined_trend_plot2 <- (p_no2 + p_no) / (p_pm25 + p_pm10) +
  plot_annotation(title = "Trend of 12-Month Moving Average Concentrations for All Pollutants",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

combined_trend_plot2

# Combine all plots in one figure
combined_trend_plot <- (p_no2 + p_no + p_so2 + p_ox) / (p_pm25 + p_pm10 + p_o3 + p_co) +
  plot_annotation(title = "Trend of 12-Month Moving Average Concentrations for All Pollutants",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

combined_trend_plot