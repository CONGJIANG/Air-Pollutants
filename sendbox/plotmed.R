library(tidyverse)

# Set the folder path
folder_path <- "/n/home09/c55jiang/res_full 2"

# Get all CSV files
files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# Function to extract metadata from filename
extract_metadata <- function(filename) {
  parts <- str_match(basename(filename), "(.*)_n_(\\d+)rate(0\\.\\d+)\\.csv")
  tibble(
    method_name = parts[,2],
    sample_size = as.integer(parts[,3]),
    censoring_rate = as.numeric(parts[,4]),
    file = filename
  )
}

# Extract metadata for all files
metadata <- map_dfr(files, extract_metadata)
metadata$file

metadata <- metadata %>%
  mutate(method_name = case_when(
    method_name == "Traditional" ~ "LOD/2 imputation",
    method_name == "EM-HAL-tetero" ~ "EM (4): heterosced. CDE, HAL", 
    method_name == "EM-HAL" ~ "EM (5): homosced. CDE, HAL",
    method_name == "EM-correctGLM" ~ "EM (1): true nuisances",
    method_name == "EM-misspecGLM" ~ "EM (2): oracle dens., misspec. GLMs",
    method_name == "EM-misspecDens" ~ "EM (3): misspec. dens., misspec. GLMs",
    TRUE ~ method_name  # Keep other values unchanged
  ))



# Function to read a file and add metadata
read_and_label_data <- function(file, method_name, sample_size, censoring_rate) {
  read_csv(file) %>%
    mutate(
      effect_type = c("NDE", "NIE"),  # Assign effect type to rows
      method_name = method_name,
      sample_size = sample_size,
      censoring_rate = censoring_rate
    )
}

# Apply function to all metadata entries
full_data <- metadata %>%
  pmap_dfr(~ read_and_label_data(..4, ..1, ..2, ..3))




library(tidyverse)
library(ggplot2)
library(latex2exp)
library(MetBrewer)

summary_plot_fun <- function(full_data, effect_type, type = "b") {
  # Transform data to long format
  comb_tbl <- full_data %>%
    pivot_longer(
      cols = c(Mean_Bias, Mean_Var, MSE, Cov_stand, Cov_moutn),
      names_to = "measure",
      values_to = "value"
    ) %>%
    mutate(
      measure = factor(measure, levels = c("Mean_Bias", "Mean_Var", "MSE", "Cov_stand", "Cov_moutn")),
      sample_size = factor(sample_size),
      censoring_rate = factor(censoring_rate),
      y_intercept = case_when(
        measure == "Mean_Bias" ~ 0,
        measure == "Mean_Var" ~ 0,
        measure == "MSE" ~ 0,
        measure == "Cov_stand" ~ 0.95,
        measure == "Cov_moutn" ~ 0.95
      )
    )
  
  # Define labels
  levels(comb_tbl$measure) <- c(
    latex2exp::TeX("Bias"),
    latex2exp::TeX("Variance"),
    latex2exp::TeX("MSE"),
    latex2exp::TeX("CI\\ Coverage"),
    latex2exp::TeX("M-out-n \\ Cov.")
  )
  
  # Define custom colors and point characters
  # Using MetBrewer colorblind-friendly palette: VanGogh1
  # Colorblind-friendly palette inspired by Van Gogh's work
  custom_colors <- met.brewer(name = "VanGogh1", n = 6, type = "discrete")
  #custom_pch <- c(17, 17, 20, 20, 8, 5)
  custom_pch <- c(8, 20, 20, 17, 17, 5)  # Different point characters

  # Create the plot combining geom_jitter and geom_point
  plt <- ggplot(comb_tbl, aes(x = sample_size, y = value, colour = method_name, shape = method_name, group = method_name)) +
    facet_grid(rows = vars(measure), cols = vars(censoring_rate), scales = "free_y", labeller = label_parsed) +
    geom_line(alpha = 0.8) +
    geom_jitter(alpha = 0.3, size = ifelse(type == "o", 2, 3), width = 0.1, height = 0) +
    geom_point(alpha = 0.6, size = ifelse(type == "o", 2, 3)) +
    geom_hline(aes(yintercept = y_intercept), linetype = 2, alpha = 0.5) +
    xlab("Sample Size") +
    ylab("Value (300 Replicates)") +
    scale_color_manual(name = "Method", values = custom_colors) +
    scale_shape_manual(name = "Method", values = custom_pch) +
    ggtitle(paste("", effect_type)) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(plt)
}


# Example usage:
NDE <- full_data[full_data$effect_type == "NDE",]
summary_plot_fun(NDE, "G-computation NDE estimates under different imputation methods, true RD = 0.413")

NIE <- full_data[full_data$effect_type == "NIE",]
summary_plot_fun(NIE, "G-computation NIE estimates under different imputation methods, true RD = 0.396")



