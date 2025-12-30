library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(purrr)

# Load data
data <- readRDS("/Users/cjiang/Downloads/MortalityPM_28082025.RDS")

# Create a simple district map_sf if not already loaded
if (!exists("map_district")) {
  # Use district centroids - extract first Latitude value (same for all rows per district)
  districts_data <- data |>
    group_by(District) |>
    summarize(
      lon = mean(Longitude, na.rm = TRUE),
      lat = as.numeric(first(Latitude)),
      .groups = "drop"
    ) |>
    filter(!is.na(lon) & !is.na(lat))
  
  map_district <- sf::st_as_sf(districts_data, coords = c("lon", "lat"), crs = 4326)
}

# ---- Define pollutant aggregation function ----
create_district_map <- function(data, map_sf, variable_name, display_name) {
  var_sym <- sym(variable_name)
  
  district_summary <- data |>
    group_by(District) |>
    summarize(
      value = mean(!!var_sym, na.rm = TRUE),
      .groups = "drop"
    )
  
  map_sf |>
    left_join(district_summary, by = "District") |>
    ggplot(aes(fill = value)) +
    geom_sf(color = "white", size = 0.2) +
    scale_fill_viridis_c(option = "plasma", name = display_name) +
    theme_void(base_size = 11) +
    labs(title = display_name) +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10)
    )
}

# ---- Aggregate all 8 pollutants and mortality by district ----
pollutant_data <- data.frame(
  name = c("NO₂", "NO", "SO₂", "Ox", "PM₂.₅", "PM₁₀", "CO", "O₃", "Mortality"),
  variable = c("mean_NO2_ensemble", "mean_NOppb", "mean_SO2_ensemble", "mean_Ox", 
               "mean_PM25_ensemble", "mean_PM10ug", "CO_pred", "mean_O3_pred_cal", "AggregatedDeath"),
  stringsAsFactors = FALSE
)

# Create all maps using pmap
maps <- pmap(
  pollutant_data,
  ~ create_district_map(data, map_district, ..2, ..1)  # ..2 = variable, ..1 = name
)
names(maps) <- pollutant_data$name



# ---- Display maps in 3x3 grid ----
combined_map <- (maps[[1]] + maps[[2]] + maps[[3]]) / 
                (maps[[4]] + maps[[5]] + maps[[6]]) / 
                (maps[[7]] + maps[[8]] + maps[[9]]) +
  plot_annotation(
    title = "District-Wise Spatial Distribution: 8 Air Pollutants & Mortality",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10))
    )
  )

combined_map

# Save combined map
ggsave("Spatial_Distribution_All_Pollutants_Mortality.png",
       plot = combined_map, width = 18, height = 14, dpi = 300)

# ---- Individual high-quality maps ----
# Save each pollutant individually with clean filenames
filenames <- c("NO2", "NO", "SO2", "Ox", "PM25", "PM10", "CO", "O3", "Mortality")

walk2(
  filenames,
  maps,
  ~ ggsave(
    paste0("Spatial_Distribution_", .x, ".png"),
    plot = .y,
    width = 8, height = 7, dpi = 300
  )
)
