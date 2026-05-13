Data Collection and Analysis Scripts for Climate-Malaria Relationships 
1. DATA COLLECTION
NDVI DATA
/************************************************
 * Annual Mean NDVI (2010–2023)
 * Malaria-endemic countries
 * CLEAN CSV: country | year | ndvi
 ************************************************/
// ---------------------------------------------
// 1. Load MODIS NDVI
// ---------------------------------------------
var ndviIC = ee.ImageCollection('MODIS/061/MOD13Q1')
  .select('NDVI');
// ---------------------------------------------
// 2. Years
// ---------------------------------------------
var years = ee.List.sequence(2010, 2023);
// ---------------------------------------------
// 3. Country boundaries
// ---------------------------------------------
var countries = ee.FeatureCollection('FAO/GAUL/2015/level0');
// ---------------------------------------------
// 4. Malaria-endemic countries (GAUL names)
// ---------------------------------------------
var malariaCountries = ee.List([
  'Angola','Benin','Botswana','Burkina Faso','Burundi','Cameroon',
  'Central African Republic','Chad','Comoros','Congo',
  "Côte d'Ivoire",'Democratic Republic of the Congo','Djibouti','Egypt',
  'Equatorial Guinea','Eritrea','Eswatini','Ethiopia','Gabon','Gambia',
  'Ghana','Guinea','Guinea-Bissau','Kenya','Liberia','Madagascar',
  'Malawi','Mali','Mauritania','Mozambique','Namibia','Niger','Nigeria',
  'Rwanda','Sao Tome and Principe','Senegal','Sierra Leone',
  'South Africa','South Sudan','Sudan','Togo','Uganda',
  'United Republic of Tanzania','Zambia','Zimbabwe',
  'Bangladesh','Bhutan','India','Indonesia','Myanmar','Nepal',
  'Thailand','Timor-Leste',
  'Cambodia',"Lao People's Democratic Republic",'Malaysia',
  'Papua New Guinea','Philippines','Solomon Islands','Vanuatu',
  'Viet Nam',"Democratic People's Republic of Korea",
  'Afghanistan','Pakistan','Saudi Arabia','Somalia','Yemen',
  'Bolivia','Brazil','Colombia','Ecuador','Guyana','Peru',
  'Suriname','Venezuela','Costa Rica','Dominican Republic',
  'Guatemala','Haiti','Honduras','Mexico','Nicaragua','Panama'
]);
// ---------------------------------------------
// 5. Filter countries (SAFE METHOD)
// ---------------------------------------------
var malariaFC = countries.filter(
  ee.Filter.inList('ADM0_NAME', malariaCountries)
);
print('Countries selected:', malariaFC.size()); // EXPECT ~82
// ---------------------------------------------
// 6. Annual NDVI extraction
// ---------------------------------------------
var ndviTable = ee.FeatureCollection(
  years.map(function(year) {
    var img = ndviIC
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .mean()
      .multiply(0.0001);

    var stats = img.reduceRegions({
      collection: malariaFC,
      reducer: ee.Reducer.mean(),
      scale: 250
    });
    return stats.map(function(f) {
      return ee.Feature(null, {
        country: f.get('ADM0_NAME'),
        year: year,
        ndvi: f.get('mean')
      });
    });
  })
).flatten();

print('Final rows:', ndviTable.size()); 
print('Preview:', ndviTable.limit(5));
// ---------------------------------------------
// 7. Export
// ---------------------------------------------
Export.table.toDrive({
  collection: ndviTable,
  description: 'Annual_NDVI_Malaria_Endemic_2010_2023_CLEAN',
  fileFormat: 'CSV'
});
CLIMATE DATA
/********************************************
 * ANNUAL RAINFALL & TEMPERATURE (2010–2023)
 * With validation against CHIRPS (more reliable for precipitation)
 ********************************************/
// 1. Country boundaries
var countries = ee.FeatureCollection('FAO/GAUL/2015/level0');
var countryNames = [
  'Angola','Benin','Botswana','Burkina Faso','Burundi','Cameroon',
  'Central African Republic','Chad','Comoros','Congo',"Côte d'Ivoire",
  'Democratic Republic of the Congo','Djibouti','Egypt','Equatorial Guinea',
  'Eritrea','Eswatini','Ethiopia','Gabon','Gambia','Ghana','Guinea',
  'Guinea-Bissau','Kenya','Liberia','Madagascar','Malawi','Mali',
  'Mauritania','Mozambique','Namibia','Niger','Nigeria','Rwanda',
  'Sao Tome and Principe','Senegal','Sierra Leone','South Africa',
  'South Sudan','Sudan','Togo','Uganda','United Republic of Tanzania',
  'Zambia','Zimbabwe','Bangladesh','Bhutan','India','Indonesia',
  'Myanmar','Nepal','Thailand','Timor-Leste','Cambodia','Lao PDR',
  'Malaysia','Papua New Guinea','Philippines','Solomon Islands',
  'Vanuatu','Viet Nam','DPR Korea','Afghanistan','Pakistan',
  'Saudi Arabia','Somalia','Yemen','Bolivia','Brazil','Colombia',
  'Ecuador','Guyana','Peru','Suriname','Venezuela','Costa Rica',
  'Dominican Republic','Guatemala','Haiti','Honduras','Mexico',
  'Nicaragua','Panama'
];
var malariaCountries = countries.filter(
  ee.Filter.inList('ADM0_NAME', countryNames)
);
// 2. DATA SOURCES
// -----------------
// A. CHIRPS Daily Precipitation (More reliable for rainfall)
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .select('precipitation'); // Units: mm/day
// B. ERA5-Land for Temperature (More reliable for temp)
var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY')
  .select('temperature_2m');
// 3. Year list
var years = ee.List.sequence(2010, 2023);
// 4. Annual extraction 
var annualData = ee.FeatureCollection(
  years.map(function(year) {
    year = ee.Number(year);
        // Define date range
    var startDate = ee.Date.fromYMD(year, 1, 1);
    var endDate = ee.Date.fromYMD(year.add(1), 1, 1);
        // ---- RAINFALL: Using CHIRPS (mm/day) -> Annual total (mm)
    var annualRain = chirps
      .filterDate(startDate, endDate)
      .sum(); // Sum of daily mm -> Annual total in mm
      // NO multiplication needed - CHIRPS is already in mm/day
        // ---- TEMPERATURE: Using ERA5-Land (°C)
    var annualTemp = era5
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .mean()
      .subtract(273.15); // K → °C
        return malariaCountries.map(function(country) {
      var geometry = country.geometry();
            // Extract rainfall
      var rainStats = annualRain.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geometry,
        scale: 5500, // CHIRPS native resolution ~5.5km
        maxPixels: 1e13
      });
            // Extract temperature
      var tempStats = annualTemp.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geometry,
        scale: 10000, // ERA5-Land resolution
        maxPixels: 1e13
      });
            // Get values with null checks
      var rainfall_mm = ee.Number(rainStats.get('precipitation'));
      var temperature_c = ee.Number(tempStats.get('temperature_2m'));
            // Quality checks
      rainfall_mm = ee.Algorithms.If(
        rainfall_mm.gt(0).and(rainfall_mm.lt(15000)), // Realistic range: 0-15,000 mm/year
        rainfall_mm,
        null
      );
            temperature_c = ee.Algorithms.If(
        temperature_c.gte(-50).and(temperature_c.lte(50)), // Realistic range
        temperature_c,
        null
      );
            return ee.Feature(null, {
        country: country.get('ADM0_NAME'),
        year: year,
        rainfall_mm: rainfall_mm,
        temperature_c: temperature_c,
        geometry: geometry
      });
    });
  })
).flatten();
// 5. VALIDATION: Compare with known averages for key countries
// -------------------------------------------------------------
// Add validation statistics for well-known countries
var validationCountries = ['Bangladesh', 'Egypt', 'Brazil', 'India', 'Indonesia'];
var validationData = malariaCountries.filter(
  ee.Filter.inList('ADM0_NAME', validationCountries)
).map(function(country) {
  var name = country.get('ADM0_NAME');
    // Calculate long-term average (2010-2023)
  var longTermRain = annualData
    .filter(ee.Filter.eq('country', name))
    .filter(ee.Filter.lt('year', 2023)) // Exclude 2023 for baseline
    .aggregate_mean('rainfall_mm');
      var longTermTemp = annualData
    .filter(ee.Filter.eq('country', name))
    .filter(ee.Filter.lt('year', 2023))
    .aggregate_mean('temperature_c');
    return ee.Feature(null, {
    country: name,
    avg_rainfall_mm_2010_2022: longTermRain,
    avg_temperature_c_2010_2022: longTermTemp,
    // Known reference values (approximate)
    known_rainfall_mm: ee.Dictionary({
      'Bangladesh': 2200,
      'Egypt': 51,
      'Brazil': 1760,
      'India': 1083,
      'Indonesia': 2702
    }).get(name),
    known_temperature_c: ee.Dictionary({
      'Bangladesh': 25,
      'Egypt': 22,
      'Brazil': 25,
      'India': 25,
      'Indonesia': 26
    }).get(name)
  });
});

// Print validation results
print('VALIDATION: Long-term averages vs known values');
print(validationData);
// 6. Export main data
Export.table.toDrive({
  collection: annualData,
  description: 'Annual_Rainfall_Temp_CHIRPS_ERA5_2010_2023',
  fileFormat: 'CSV',
  selectors: ['country', 'year', 'rainfall_mm', 'temperature_c']
});
// 7. Export validation data
Export.table.toDrive({
  collection: validationData,
  description: 'Validation_Against_Known_Averages',
  fileFormat: 'CSV'
});
// 8. Quick sample check for problematic 2023
print('2023 SAMPLE VALUES (check for anomalies):');
var sample2023 = annualData
  .filter(ee.Filter.eq('year', 2023))
  .limit(10);
print(sample2023);


2. DATA ANALYSIS
########################################################################
# 1. SETUP
########################################################################
library(tidyverse)
library(plm)
library(Hmisc)
library(car)
library(lmtest)
library(sandwich)
library(mediation)
library(ggplot2)
library(corrplot)
library(knitr)
library(margins)
library(scales)
library(interactions)

########################################################################
# Load required packages
library(dplyr)
library(plm)   # for pdata.frame

# 1. Read data
malaria <- read.csv("Data")

# 2. Check column names (IMPORTANT: compare with your CSV)
print(names(malaria))

# 3. Convert specified columns to numeric (only if they exist)
numeric_cols <- c("ind_malaria_inc", "temp_c", "rain_mm", "ndvi", "water_safe",
                  "gdp", "health_exp", "urban_pct", "pop_density")
existing_numeric <- intersect(numeric_cols, names(malaria))
malaria <- malaria %>%
  mutate(across(all_of(existing_numeric), as.numeric))

# 4. Remove rows with any NA (only in those numeric columns, keep index columns intact)
data_clean <- malaria[complete.cases(malaria[, existing_numeric]), ]

# 5. Ensure country and year columns exist and have no NA
if (!all(c("country", "year") %in% names(data_clean))) {
  stop("Columns 'country' and/or 'year' not found in the data.")
}
data_clean <- data_clean %>% filter(!is.na(country), !is.na(year))

# 6. Create panel data frame
malaria_panel <- pdata.frame(data_clean, index = c("country", "year"))

# Confirm success
print(head(malaria_panel, 3))




########################################################################
# 2. DESCRIPTIVE STATISTICS WITH EXPORT TO EXCEL
########################################################################

library(dplyr)
library(knitr)
library(tidyr)

# Ensure writexl is available
if (!require(writexl)) {
  install.packages("writexl")
  library(writexl)
}

# Generate descriptive statistics (wide format)
desc_wide <- data_clean %>%
  summarise(
    across(
      c(ind_malaria_inc, temp_c, rain_mm, ndvi, water_safe,
        gdp, health_exp, urban_pct, pop_density),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE),
        min  = ~min(.x, na.rm = TRUE),
        max  = ~max(.x, na.rm = TRUE)
      )
    )
  )

# Convert to long format
desc_long <- desc_wide %>%
  pivot_longer(everything(), 
               names_to = c("variable", "statistic"), 
               names_pattern = "(.*)_(mean|sd|min|max)") %>%
  pivot_wider(names_from = statistic, values_from = value)

# Display in console
kable(desc_wide, caption = "Table 4.1: Descriptive Statistics (Wide Format)")
cat("\n")
kable(desc_long, caption = "Table 4.1: Descriptive Statistics (Long Format)")

# Export options
export_wide <- TRUE
export_long <- TRUE

if (export_wide) {
  tryCatch({
    write_xlsx(desc_wide, "Descriptive_Statistics_Wide.xlsx")
    cat("Exported wide format to Excel.\n")
  }, error = function(e) {
    write.csv(desc_wide, "Descriptive_Statistics_Wide.csv", row.names = FALSE)
    cat("Excel export failed. Saved as CSV instead: Descriptive_Statistics_Wide.csv\n")
  })
}

if (export_long) {
  tryCatch({
    write_xlsx(desc_long, "Descriptive_Statistics_Long.xlsx")
    cat("Exported long format to Excel.\n")
  }, error = function(e) {
    write.csv(desc_long, "Descriptive_Statistics_Long.csv", row.names = FALSE)
    cat("Excel export failed. Saved as CSV instead: Descriptive_Statistics_Long.csv\n")
  })
}
getwd()




########################################################################
# 3. DESCRIPTIVE PLOTS 
########################################################################

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# ---- Prepare regional summary ----
region_summary <- data_clean %>%
  group_by(region) %>%
  summarise(mean_malaria = mean(ind_malaria_inc, na.rm = TRUE),
            .groups = "drop")

# ---- Clean, publication-ready theme ----
pub_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- Panel A: Histogram of malaria incidence ----
p1 <- ggplot(data_clean, aes(x = ind_malaria_inc)) +
  geom_histogram(bins = 30, fill = "#4D7EA8", color = "white", alpha = 0.8) +
  labs(title = "A", x = "Malaria incidence", y = "Frequency") +
  pub_theme

# ---- Panel B: Regional mean bars (horizontal) ----
p2 <- ggplot(region_summary, aes(x = reorder(region, mean_malaria), y = mean_malaria)) +
  geom_col(fill = "#D95F02", color = "black", width = 0.7) +
  labs(title = "B", x = "Region", y = "Mean malaria incidence") +
  coord_flip() +
  pub_theme +
  theme(axis.text.y = element_text(size = 10))

# ---- Arrange side‑by‑side using gridExtra ----
combined <- grid.arrange(
  p1, p2,
  ncol = 2,
  widths = c(1, 0.9),
  top = textGrob("Distribution of malaria incidence", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# ---- Export ----
ggsave("Figure1_Descriptive.tiff", plot = combined,
       width = 6.5, height = 3.5, units = "in", dpi = 300,
       compression = "lzw")

ggsave("Figure1_Descriptive.pdf", plot = combined,
       width = 6.5, height = 3.5, units = "in", device = pdf)

# Display on screen
print(combined)


########################################################################
# CLIMATE VS MALARIA SCATTER PLOTS 
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# ---- Clean publication theme ----
pub_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 9),
    legend.position = "none",                # no legends needed
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- Panel A: Rainfall vs Malaria (linear) ----
p1 <- ggplot(data_clean, aes(x = rain_mm, y = ind_malaria_inc)) +
  geom_point(alpha = 0.3, size = 1.2, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#D95F02", fill = "#D95F02", alpha = 0.2) +
  labs(title = "A", x = "Rainfall (mm)", y = "Malaria incidence") +
  pub_theme

# ---- Panel B: Rainfall vs Malaria (non-linear, loess) ----
p2 <- ggplot(data_clean, aes(x = rain_mm, y = ind_malaria_inc)) +
  geom_point(alpha = 0.3, size = 1.2, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "#1B9E77", fill = "#1B9E77", alpha = 0.2) +
  labs(title = "B", x = "Rainfall (mm)", y = "Malaria incidence") +
  pub_theme

# ---- Panel C: Temperature vs Malaria (linear) ----
p3 <- ggplot(data_clean, aes(x = temp_c, y = ind_malaria_inc)) +
  geom_point(alpha = 0.3, size = 1.2, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#D95F02", fill = "#D95F02", alpha = 0.2) +
  labs(title = "C", x = "Temperature (°C)", y = "Malaria incidence") +
  pub_theme

# ---- Panel D: Temperature vs Malaria (non-linear, loess) ----
p4 <- ggplot(data_clean, aes(x = temp_c, y = ind_malaria_inc)) +
  geom_point(alpha = 0.3, size = 1.2, color = "grey40") +
  geom_smooth(method = "loess", se = TRUE, color = "#1B9E77", fill = "#1B9E77", alpha = 0.2) +
  labs(title = "D", x = "Temperature (°C)", y = "Malaria incidence") +
  pub_theme

# ---- Combine into 2x2 grid ----
combined <- grid.arrange(
  p1, p2, p3, p4,
  ncol = 2,
  top = textGrob("Relationship between climate variables and malaria incidence",
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# ---- Export ----
# Option 1: Single column width (6.5 inches square)
ggsave("Figure2_Climate_Malaria.tiff", plot = combined,
       width = 6.5, height = 6.5, units = "in", dpi = 300, compression = "lzw")
ggsave("Figure2_Climate_Malaria.pdf", plot = combined,
       width = 6.5, height = 6.5, units = "in")

# Option 2: Full page width (9 inches square – adjust as needed)
# ggsave("Figure2_Climate_Malaria_full.tiff", plot = combined,
#        width = 9, height = 9, units = "in", dpi = 300, compression = "lzw")

# Display on screen
print(combined)






################################################################################
# CLIMATE-MALARIA TRENDS – JOURNAL STANDARD
################################################################################
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)

# -------------------- 1. READ DATA --------------------
data_raw <- read.csv("Data")

# Check column names (diagnostic)
cat("Column names in CSV:\n")
print(names(data_raw))

# Verify required columns exist
required_cols <- c("country", "region", "year", "ind_malaria_inc", "rain_mm", "temp_c")
missing_cols <- required_cols[!required_cols %in% names(data_raw)]
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# -------------------- 2. CLEAN DATA --------------------
data_clean <- data_raw %>%
  dplyr::select(country, region, year, ind_malaria_inc, rain_mm, temp_c) %>%
  filter(region %in% c("Africa", "Asia", "North_America", "South_America")) %>%
  mutate(
    year = as.numeric(year),
    region = factor(region, levels = c("Africa", "Asia", "North_America", "South_America"))
  ) %>%
  filter(!is.na(temp_c), !is.na(rain_mm), !is.na(ind_malaria_inc))

if (nrow(data_clean) == 0) {
  stop("No data after cleaning. Check region names: ", 
       paste(unique(data_raw$region), collapse = ", "))
} else {
  message("Data loaded successfully: ", nrow(data_clean), " rows.")
}

# -------------------- 3. AGGREGATED DATA --------------------
temp_region <- data_clean %>% 
  group_by(year, region) %>% 
  summarise(mean_temp = mean(temp_c), .groups = "drop")
rain_region <- data_clean %>% 
  group_by(year, region) %>% 
  summarise(mean_rain = mean(rain_mm), .groups = "drop")
malaria_region <- data_clean %>% 
  group_by(year, region) %>% 
  summarise(mean_malaria = mean(ind_malaria_inc), .groups = "drop")

temp_overall <- temp_region %>% group_by(year) %>% summarise(mean_temp = mean(mean_temp))
rain_overall <- rain_region %>% group_by(year) %>% summarise(mean_rain = mean(mean_rain))
malaria_overall <- malaria_region %>% group_by(year) %>% summarise(mean_malaria = mean(mean_malaria))

# -------------------- 4. COLOURS --------------------
region_colours <- c(
  "Africa"        = "#1b9e77",
  "Asia"          = "#d95f02",
  "North_America" = "#7570b3",
  "South_America" = "#e7298a",
  "Average"       = "black"
)

# -------------------- 5. BASE THEME --------------------
my_theme <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 8)),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# -------------------- 6. TEMPERATURE PLOT --------------------
p1 <- ggplot(temp_region, aes(x = year, y = mean_temp, colour = region)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_line(data = temp_overall, aes(x = year, y = mean_temp, colour = "Average"),
            linewidth = 1.0, linetype = "dashed", inherit.aes = FALSE) +
  scale_colour_manual(values = region_colours, 
                      breaks = c(levels(temp_region$region), "Average")) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(title = "Temperature Trends", x = "", y = "Mean Temp (°C)", colour = "Key") +
  my_theme +
  guides(colour = guide_legend(override.aes = list(linewidth = c(rep(1.2, 4), 1.0),
                                                   linetype = c(rep("solid", 4), "dashed"))))

# -------------------- 7. RAINFALL PLOT --------------------
p2 <- ggplot(rain_region, aes(x = year, y = mean_rain, colour = region)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_line(data = rain_overall, aes(x = year, y = mean_rain, colour = "Average"),
            linewidth = 1.0, linetype = "dashed", inherit.aes = FALSE) +
  scale_colour_manual(values = region_colours, 
                      breaks = c(levels(rain_region$region), "Average")) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(title = "Rainfall Trends", x = "", y = "Mean Rainfall (mm)", colour = "Key") +
  my_theme +
  guides(colour = guide_legend(override.aes = list(linewidth = c(rep(1.2, 4), 1.0),
                                                   linetype = c(rep("solid", 4), "dashed"))))

# -------------------- 8. MALARIA PLOT (facets, no legend) --------------------
p3 <- ggplot(malaria_region, aes(x = year, y = mean_malaria, colour = region)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_line(data = malaria_overall, aes(x = year, y = mean_malaria, colour = "Average"),
            linewidth = 1.0, linetype = "dashed", inherit.aes = FALSE) +
  scale_colour_manual(values = region_colours) +
  scale_x_continuous(breaks = 2015:2023) +
  labs(title = "Malaria Trends (separate y-axis per region)", 
       x = "Year", y = "Malaria Incidence") +
  facet_wrap(~ region, scales = "free_y", ncol = 2) +
  my_theme +
  theme(strip.text = element_text(face = "bold", size = 11),
        legend.position = "none")

# -------------------- 9. EXTRACT LEGEND SAFELY --------------------
extract_legend <- function(p) {
  tmp <- ggplotGrob(p)
  leg <- tmp$grobs[[which(sapply(tmp$grobs, function(x) x$name) == "guide-box")]]
  if (is.null(leg)) stop("Legend not found in plot. Check that p1 has a legend.")
  return(leg)
}

shared_legend <- extract_legend(p1)

# Remove legends from p1 and p2
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

# -------------------- 10. COMBINE PLOTS --------------------
stack_plots <- plot_grid(p1, p2, p3, ncol = 1, align = "v", 
                         rel_heights = c(1, 1, 1.2),
                         labels = c("A", "B", "C"), label_size = 12)

final_plot <- plot_grid(stack_plots, shared_legend, ncol = 2, rel_widths = c(3.5, 0.8))

title <- ggdraw() + draw_label("Climate-Malaria Trends by Region", 
                               fontface = "bold", size = 16, hjust = 0.5)

final_with_title <- plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.08, 1))

# -------------------- 11. DISPLAY AND SAVE --------------------
# Display in RStudio
print(final_with_title)

# Save as journal-standard TIFF and PDF
ggsave("Climate_Malaria_Trends.tiff", plot = final_with_title,
       width = 7, height = 9, units = "in", dpi = 300, compression = "lzw")

ggsave("Climate_Malaria_Trends.pdf", plot = final_with_title,
       width = 7, height = 9, units = "in", device = pdf)

message("Plot saved as TIFF and PDF in: ", getwd())





########################################################################
# CHECK DISTRIBUTIONS FOR VARIABLE TRANSFORMATION

########################################################################
# 1. LOAD ORIGINAL DATASET
########################################################################

original_data <- read.csv(
  "D:/Documents/SEMESTER 3/Thesis/Data collection/Malaria_Panel_Data_csv.txt",
  stringsAsFactors = FALSE
)

########################################################################
# 2. RESTORE MISSING VARIABLES INTO data_clean
########################################################################

# If data_clean already exists, restore missing variables from original_data
restore_vars <- c(
  "ndvi", "water_safe", "gdp", "health_exp",
  "pop_density", "urban_pct"
)

for (v in restore_vars) {
  if (v %in% names(original_data) && !v %in% names(data_clean)) {
    data_clean[[v]] <- original_data[[v]]
  }
}

########################################################################
# 3. CHECK MISSINGNESS
########################################################################

missing_count <- colSums(is.na(data_clean))
missing_percent <- round(missing_count / nrow(data_clean) * 100, 2)

missing_summary <- data.frame(
  variable = names(missing_count),
  missing_count = as.integer(missing_count),
  missing_percent = as.numeric(missing_percent)
)

print(missing_summary[missing_summary$missing_count > 0, ])

########################################################################
# 4. IMPUTE NUMERIC VARIABLES IF NEEDED
########################################################################

to_numeric_safe <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- gsub("%", "", x)
  x <- gsub("\\$", "", x)
  x <- trimws(x)
  as.numeric(x)
}

vars_to_impute <- c(
  "ind_malaria_inc", "rain_mm", "temp_c", "ndvi",
  "water_safe", "gdp", "health_exp", "pop_density", "urban_pct"
)

for (v in vars_to_impute) {
  if (v %in% names(data_clean)) {
    data_clean[[v]] <- to_numeric_safe(data_clean[[v]])
    if (any(is.na(data_clean[[v]]))) {
      data_clean[[v]][is.na(data_clean[[v]])] <- mean(data_clean[[v]], na.rm = TRUE)
    }
  }
}



########################################################################
# CHECK DISTRIBUTIONS WITH HISTOGRAMS + DENSITY
# Journal-standard export: TIFF (600 dpi) + PDF
########################################################################

# Create output directory if needed
if(!dir.exists("outputs/plots")) dir.create("outputs/plots", recursive = TRUE)

# -----------------------------
# TIFF EXPORT
# -----------------------------
tiff("outputs/plots/distribution_assessment.tiff",
     width = 12, height = 12, units = "in",
     res = 600, compression = "lzw")

par(mfrow = c(3, 3),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 3, 0),
    cex.main = 1,
    cex.lab = 0.9,
    cex.axis = 0.8)

hist(data_clean$ind_malaria_inc, freq = FALSE,
     main = "Malaria",
     xlab = "Indigenous malaria incidence",
     col = "steelblue",
     border = "white")
lines(density(data_clean$ind_malaria_inc, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$rain_mm, freq = FALSE,
     main = "Rainfall",
     xlab = "Rainfall (mm)",
     col = "purple",
     border = "white")
lines(density(data_clean$rain_mm, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$temp_c, freq = FALSE,
     main = "Temperature",
     xlab = "Temperature (°C)",
     col = "goldenrod",
     border = "white")
lines(density(data_clean$temp_c, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$ndvi, freq = FALSE,
     main = "NDVI",
     xlab = "NDVI",
     col = "darkgreen",
     border = "white")
lines(density(data_clean$ndvi, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$water_safe, freq = FALSE,
     main = "Water Access",
     xlab = "Water access",
     col = "cyan4",
     border = "white")
lines(density(data_clean$water_safe, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$gdp, freq = FALSE,
     main = "GDP",
     xlab = "GDP",
     col = "tomato",
     border = "white")
lines(density(data_clean$gdp, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$health_exp, freq = FALSE,
     main = "Health Expenditure",
     xlab = "Health expenditure",
     col = "orange",
     border = "white")
lines(density(data_clean$health_exp, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$pop_density, freq = FALSE,
     main = "Population Density",
     xlab = "Population density",
     col = "seagreen",
     border = "white")
lines(density(data_clean$pop_density, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$urban_pct, freq = FALSE,
     main = "Urban Percentage",
     xlab = "Urban percentage",
     col = "darkred",
     border = "white")
lines(density(data_clean$urban_pct, na.rm = TRUE),
      col = "red", lwd = 2)

mtext("Distribution of Variables for Transformation Assessment",
      outer = TRUE,
      cex = 1.5,
      font = 2)

dev.off()

# -----------------------------
# PDF EXPORT
# -----------------------------
pdf("outputs/plots/distribution_assessment.pdf",
    width = 12, height = 12)

par(mfrow = c(3, 3),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 3, 0),
    cex.main = 1,
    cex.lab = 0.9,
    cex.axis = 0.8)

hist(data_clean$ind_malaria_inc, freq = FALSE,
     main = "Malaria",
     xlab = "Indigenous malaria incidence",
     col = "steelblue",
     border = "white")
lines(density(data_clean$ind_malaria_inc, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$rain_mm, freq = FALSE,
     main = "Rainfall",
     xlab = "Rainfall (mm)",
     col = "purple",
     border = "white")
lines(density(data_clean$rain_mm, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$temp_c, freq = FALSE,
     main = "Temperature",
     xlab = "Temperature (°C)",
     col = "goldenrod",
     border = "white")
lines(density(data_clean$temp_c, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$ndvi, freq = FALSE,
     main = "NDVI",
     xlab = "NDVI",
     col = "darkgreen",
     border = "white")
lines(density(data_clean$ndvi, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$water_safe, freq = FALSE,
     main = "Water Access",
     xlab = "Water access",
     col = "cyan4",
     border = "white")
lines(density(data_clean$water_safe, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$gdp, freq = FALSE,
     main = "GDP",
     xlab = "GDP",
     col = "tomato",
     border = "white")
lines(density(data_clean$gdp, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$health_exp, freq = FALSE,
     main = "Health Expenditure",
     xlab = "Health expenditure",
     col = "orange",
     border = "white")
lines(density(data_clean$health_exp, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$pop_density, freq = FALSE,
     main = "Population Density",
     xlab = "Population density",
     col = "seagreen",
     border = "white")
lines(density(data_clean$pop_density, na.rm = TRUE),
      col = "red", lwd = 2)

hist(data_clean$urban_pct, freq = FALSE,
     main = "Urban Percentage",
     xlab = "Urban percentage",
     col = "darkred",
     border = "white")
lines(density(data_clean$urban_pct, na.rm = TRUE),
      col = "red", lwd = 2)

mtext("Distribution of Variables for Transformation Assessment",
      outer = TRUE,
      cex = 1.5,
      font = 2)

dev.off()

getwd()


################################################################################
# TRANSFORMATION SCRIPT – NATURAL LOG FOR RIGHT-SKEWED VARIABLES
# Added journal-standard exports:
# - TIFF + PDF plots
# - Excel tables
################################################################################

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(openxlsx)

# Create output folders
if(!dir.exists("outputs")) dir.create("outputs")
if(!dir.exists("outputs/plots")) dir.create("outputs/plots", recursive = TRUE)
if(!dir.exists("outputs/tables")) dir.create("outputs/tables", recursive = TRUE)

# -------------------- 1. DEFINE VARIABLE GROUPS --------------------
log_vars <- c("ind_malaria_inc", "gdp", "health_exp", "pop_density")
keep_vars <- c("temp_c", "rain_mm", "ndvi", "water_safe", "urban_pct")

# -------------------- 2. CHECK FOR NON-POSITIVE VALUES --------------------
for (var in log_vars) {
  if (var %in% names(data_clean)) {
    min_val <- min(data_clean[[var]], na.rm = TRUE)
    if (min_val <= 0) {
      warning(sprintf("Variable '%s' has non-positive values (min = %g). Adding +1 to avoid undefined log.", var, min_val))
      data_clean[[var]] <- data_clean[[var]] + 1
    }
  } else {
    warning(sprintf("Variable '%s' not found in data_clean. Skipping.", var))
  }
}

# -------------------- 3. APPLY NATURAL LOG TRANSFORMATION --------------------
data_transformed <- data_clean %>%
  mutate(
    across(all_of(log_vars), ~ log(.), .names = "log_{.col}")
  )

data_transformed <- data_transformed %>%
  dplyr::select(
    country, region, year,
    all_of(log_vars),
    starts_with("log_"),
    all_of(keep_vars)
  )

# -------------------- 4. CHECK RESULTS --------------------
print("Transformed data first rows:")
print(head(data_transformed))

print("Column names after transformation:")
print(names(data_transformed))

# -------------------- 5. SUMMARY STATISTICS --------------------
if (require(e1071, quietly = TRUE)) {
  
  orig_stats <- data_clean %>%
    summarise(across(all_of(log_vars), 
                     list(mean = ~mean(., na.rm = TRUE),
                          sd = ~sd(., na.rm = TRUE),
                          skew = ~skewness(., na.rm = TRUE))))
  
  log_stats <- data_transformed %>%
    summarise(across(starts_with("log_"),
                     list(mean = ~mean(., na.rm = TRUE),
                          sd = ~sd(., na.rm = TRUE),
                          skew = ~skewness(., na.rm = TRUE))))
  
  print("Original variable statistics:")
  print(orig_stats)
  
  print("Log-transformed variable statistics:")
  print(log_stats)
  
  # Save statistics to Excel workbook
  wb <- createWorkbook()
  
  addWorksheet(wb, "Original_Stats")
  writeData(wb, "Original_Stats", orig_stats)
  
  addWorksheet(wb, "Log_Transformed_Stats")
  writeData(wb, "Log_Transformed_Stats", log_stats)
  
  saveWorkbook(wb,
               "outputs/tables/transformation_summary_tables.xlsx",
               overwrite = TRUE)
  
} else {
  print("Install e1071 for skewness (optional): install.packages('e1071')")
}

# -------------------- 6. HISTOGRAM COMPARISON PLOTS --------------------
p_orig <- ggplot(data_clean, aes(x = ind_malaria_inc)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(title = "Original Malaria Incidence",
       x = "ind_malaria_inc",
       y = "Frequency") +
  theme_minimal(base_size = 12)

p_log <- ggplot(data_transformed, aes(x = log_ind_malaria_inc)) +
  geom_histogram(bins = 30, fill = "tomato", alpha = 0.7) +
  labs(title = "Log-transformed Malaria Incidence",
       x = "log(ind_malaria_inc)",
       y = "Frequency") +
  theme_minimal(base_size = 12)

# TIFF export (journal standard)
tiff("outputs/plots/log_transformation_comparison.tiff",
     width = 10, height = 5, units = "in",
     res = 600, compression = "lzw")

grid.arrange(p_orig, p_log, ncol = 2)

dev.off()

# PDF export
pdf("outputs/plots/log_transformation_comparison.pdf",
    width = 10, height = 5)

grid.arrange(p_orig, p_log, ncol = 2)

dev.off()

# Display in session
grid.arrange(p_orig, p_log, ncol = 2)

# -------------------- 7. SAVE TRANSFORMED DATA --------------------
write.csv(data_transformed,
          "outputs/Malaria_Data_Transformed.csv",
          row.names = FALSE)

saveRDS(data_transformed,
        "outputs/Malaria_Data_Transformed.rds")

write.xlsx(data_transformed,
           "outputs/tables/Malaria_Data_Transformed.xlsx",
           overwrite = TRUE)

# -------------------- 8. CREATE PANEL DATA --------------------
if (require(plm, quietly = TRUE)) {
  malaria_panel_transformed <- pdata.frame(data_transformed,
                                           index = c("country", "year"))
  print("Panel data frame created.")
} else {
  message("plm package not installed. To create panel data, run: install.packages('plm')")
}

# -------------------- 9. FINAL DATASET --------------------
data <- na.omit(data_transformed)

message("Transformation complete. Files saved in: ", getwd())


#############################################################################
# 1. DESCRIPTIVE STATISTICS AND CORRELATION ANALYSIS (Objective 1)
#############################################################################

# ==== CLEAN SESSION ====
rm(list = ls())
graphics.off()

# ==== LOAD PACKAGES ====
packages <- c("psych", "corrplot")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# ==== LOAD TRANSFORMED DATA ====
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Transformed.rds"

if (!file.exists(file_path)) {
  stop("File not found. Check the file path.")
}

data <- readRDS(file_path)

cat("✓ Data loaded successfully\n")

# ==== INSPECT DATA STRUCTURE ====
cat("\nAvailable variables:\n")
print(names(data))

# ==== DEFINE VARIABLES ====
cor_vars <- c("log_ind_malaria_inc", "temp_c", "rain_mm", "ndvi", "water_safe",
              "log_gdp", "log_health_exp", "log_pop_density", "urban_pct")

# ==== LOAD PACKAGE FOR EXCEL EXPORT ====
library(openxlsx)

# ==== VALIDATE VARIABLES ====
missing_vars <- cor_vars[!cor_vars %in% names(data)]

if (length(missing_vars) > 0) {
  stop(paste("Missing variables:", paste(missing_vars, collapse = ", ")))
}

cat("✓ All required variables found\n")

# ==== CORRELATION ANALYSIS ====
cor_result <- psych::corr.test(
  data[, cor_vars],
  use = "pairwise.complete.obs"
)

cor_matrix <- cor_result$r
cor_pvals  <- cor_result$p

cat("\n✓ Correlation analysis complete\n")
print(round(cor_matrix, 3))

# ==== SAVE SUMMARY TABLES TO EXCEL ====

# Convert matrices to data frames for export
cor_matrix_df <- as.data.frame(round(cor_matrix, 4))
cor_matrix_df <- cbind(Variable = rownames(cor_matrix_df), cor_matrix_df)

cor_pvals_df <- as.data.frame(round(cor_pvals, 4))
cor_pvals_df <- cbind(Variable = rownames(cor_pvals_df), cor_pvals_df)

# Create workbook
wb <- createWorkbook()

# Correlation coefficients sheet
addWorksheet(wb, "Correlation_Coefficients")
writeData(wb, "Correlation_Coefficients", cor_matrix_df)

# P-values sheet
addWorksheet(wb, "Correlation_P_Values")
writeData(wb, "Correlation_P_Values", cor_pvals_df)

# Save workbook
saveWorkbook(
  wb,
  file = "Table_Correlation_Summary.xlsx",
  overwrite = TRUE
)

cat("✓ Summary tables saved: Table_Correlation_Summary.xlsx\n")

# ==== SAVE PLOT OUTPUTS ====

# TIFF (publication quality)
tiff("Figure1_Correlation_Matrix.tiff",
     width = 8, height = 7, units = "in",
     res = 300, compression = "lzw")

corrplot::corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  p.mat = cor_pvals,
  sig.level = 0.05,
  insig = "blank",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  number.cex = 0.7,
  title = "Correlation Matrix\n(*p < 0.05)",
  mar = c(0, 0, 2, 0)
)

dev.off()

# PDF version
pdf("Figure1_Correlation_Matrix.pdf", width = 8, height = 7)

corrplot::corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  p.mat = cor_pvals,
  sig.level = 0.05,
  insig = "blank",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  number.cex = 0.7,
  title = "Correlation Matrix\n(*p < 0.05)",
  mar = c(0, 0, 2, 0)
)

dev.off()

cat("\n✓ Files saved:\n")
cat("  - Figure1_Correlation_Matrix.tiff\n")
cat("  - Figure1_Correlation_Matrix.pdf\n")
cat("  - Table_Correlation_Summary.xlsx\n")
cat("✓ Non-significant correlations removed from visual plot (p > 0.05)\n")

# ==== SHOW WORKING DIRECTORY ====
getwd()


################################################################################
# FULL PANEL MODEL ESTIMATION, VIF, HAUSMAN TEST, AND DIAGNOSTICS
# Includes all explanatory variables
################################################################################

# ===============================
# LOAD REQUIRED PACKAGES
# ===============================
library(plm)
library(car)
library(lmtest)
library(openxlsx)
library(dplyr)

# ===============================
# LOAD TRANSFORMED DATASET
# ===============================
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Transformed.rds"

data_transformed <- readRDS(file_path)

# Convert to dataframe
data <- as.data.frame(data_transformed)

cat("✓ Data successfully loaded\n")
cat("✓ Observations:", nrow(data), "\n")
cat("✓ Variables:", ncol(data), "\n")

# ===============================
# REMOVE MISSING VALUES
# ===============================
data <- na.omit(data)

# ===============================
# CREATE PANEL DATA FRAME
# ===============================
pdata <- pdata.frame(
  data,
  index = c("country", "year")
)

cat("✓ Panel data frame created\n")

# ===============================
# DEFINE FULL MODEL FORMULA
# ===============================
model_formula <- log_ind_malaria_inc ~ 
  temp_c +
  rain_mm +
  ndvi +
  water_safe +
  log_gdp +
  log_health_exp +
  log_pop_density +
  urban_pct

# ===============================
# VIF ASSESSMENT
# ===============================
vif_model <- lm(
  model_formula,
  data = data
)

vif_values <- car::vif(vif_model)

vif_table <- data.frame(
  Variable = names(vif_values),
  VIF = round(as.numeric(vif_values), 4)
)

cat("✓ Full VIF assessment complete\n")
print(vif_table)

# ===============================
# FIXED EFFECTS MODEL
# ===============================
fe_model <- plm(
  model_formula,
  data = pdata,
  model = "within"
)

# ===============================
# RANDOM EFFECTS MODEL
# ===============================
re_model <- plm(
  model_formula,
  data = pdata,
  model = "random"
)

cat("✓ FE and RE models estimated\n")

# ===============================
# HAUSMAN TEST
# ===============================
hausman <- phtest(fe_model, re_model)

hausman_table <- data.frame(
  Test = "Hausman Test (FE vs RE)",
  Chi_Square = round(as.numeric(hausman$statistic), 4),
  DF = as.numeric(hausman$parameter),
  P_Value = round(as.numeric(hausman$p.value), 4),
  Decision = ifelse(
    hausman$p.value < 0.05,
    "Fixed Effects Preferred",
    "Random Effects Preferred"
  )
)

cat("✓ Hausman test complete\n")
print(hausman_table)

# ===============================
# SERIAL CORRELATION TEST
# ===============================
serial_test <- pbgtest(fe_model)

serial_table <- data.frame(
  Test = "Breusch-Godfrey Serial Correlation",
  Statistic = round(as.numeric(serial_test$statistic), 4),
  DF = as.numeric(serial_test$parameter),
  P_Value = round(as.numeric(serial_test$p.value), 4),
  Decision = ifelse(
    serial_test$p.value < 0.05,
    "Serial Correlation Present",
    "No Serial Correlation"
  )
)

cat("✓ Serial correlation test complete\n")
print(serial_table)

# ===============================
# HETEROSKEDASTICITY TEST
# ===============================
hetero_test <- bptest(
  fe_model,
  studentize = FALSE
)

hetero_table <- data.frame(
  Test = "Breusch-Pagan Heteroskedasticity",
  Statistic = round(as.numeric(hetero_test$statistic), 4),
  DF = as.numeric(hetero_test$parameter),
  P_Value = round(as.numeric(hetero_test$p.value), 4),
  Decision = ifelse(
    hetero_test$p.value < 0.05,
    "Heteroskedasticity Present",
    "Homoskedasticity"
  )
)

cat("✓ Heteroskedasticity test complete\n")
print(hetero_table)

# ===============================
# MODEL COEFFICIENT TABLES
# ===============================
fe_summary <- summary(fe_model)
re_summary <- summary(re_model)

fe_table <- data.frame(
  Variable = rownames(fe_summary$coefficients),
  Coefficient = round(fe_summary$coefficients[, 1], 4),
  Std_Error = round(fe_summary$coefficients[, 2], 4),
  T_Value = round(fe_summary$coefficients[, 3], 4),
  P_Value = round(fe_summary$coefficients[, 4], 4)
)

re_table <- data.frame(
  Variable = rownames(re_summary$coefficients),
  Coefficient = round(re_summary$coefficients[, 1], 4),
  Std_Error = round(re_summary$coefficients[, 2], 4),
  T_Value = round(re_summary$coefficients[, 3], 4),
  P_Value = round(re_summary$coefficients[, 4], 4)
)

# ===============================
# SAVE RESULTS TO EXCEL
# ===============================
output_excel <- "D:/Documents/SEMESTER 3/Thesis/Results/Final_Outputs/Outputs/Tables/Full_Panel_Model_Diagnostics.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "VIF")
writeData(wb, "VIF", vif_table)

addWorksheet(wb, "Fixed_Effects_Model")
writeData(wb, "Fixed_Effects_Model", fe_table)

addWorksheet(wb, "Random_Effects_Model")
writeData(wb, "Random_Effects_Model", re_table)

addWorksheet(wb, "Hausman_Test")
writeData(wb, "Hausman_Test", hausman_table)

addWorksheet(wb, "Serial_Correlation")
writeData(wb, "Serial_Correlation", serial_table)

addWorksheet(wb, "Heteroskedasticity")
writeData(wb, "Heteroskedasticity", hetero_table)

saveWorkbook(
  wb,
  file = output_excel,
  overwrite = TRUE
)

cat("\n✓ Excel file saved:\n", output_excel, "\n")



################################################################################
# H1: CLIMATE → MALARIA
# REVISED STEP 2:
#   A. Adjusted climate-only pooled model
#   B. Adjusted full pathway pooled model
#
# Robust Clustered Standard Errors (group clustered)
#
# Objective:
# Compare adjusted associations:
# - Climate variables only
# - Climate + environmental + socioeconomic pathway
#
# Outputs:
# - Robust pooled models
# - Combined tables
# - Panel coefficient plots
# - Excel exports
# - TIFF/PDF journal-quality figures
################################################################################

# ===============================
# LOAD REQUIRED PACKAGES
# ===============================
library(plm)
library(lmtest)
library(sandwich)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

# ===============================
# LOAD DATA
# ===============================
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Transformed.rds"

data_transformed <- readRDS(file_path)
data <- na.omit(as.data.frame(data_transformed))

cat("✓ Data loaded successfully\n")

# ===============================
# CREATE PANEL DATA
# ===============================
pdata <- pdata.frame(
  data,
  index = c("country", "year")
)

cat("✓ Panel data created\n")

# ===============================
# DEFINE MODELS
# ===============================
climate_formula <- log_ind_malaria_inc ~
  temp_c +
  rain_mm

full_formula <- log_ind_malaria_inc ~
  temp_c +
  rain_mm +
  ndvi +
  water_safe +
  log_gdp +
  log_health_exp +
  log_pop_density +
  urban_pct

# ===============================
# MODEL A: ADJUSTED CLIMATE-ONLY
# ===============================
model_a <- plm(
  formula = climate_formula,
  data = pdata,
  model = "pooling"
)

robust_a <- vcovHC(
  model_a,
  method = "arellano",
  type = "HC1",
  cluster = "group"
)

adjusted_a <- coeftest(
  model_a,
  vcov = robust_a
)

table_a <- data.frame(
  Variable = rownames(adjusted_a),
  Coefficient = round(adjusted_a[,1], 4),
  Robust_SE = round(adjusted_a[,2], 4),
  T_Value = round(adjusted_a[,3], 4),
  P_Value = round(adjusted_a[,4], 4),
  Model = "Adjusted Climate Only"
)

cat("✓ Adjusted climate-only model complete\n")

# ===============================
# MODEL B: ADJUSTED FULL PATHWAY
# ===============================
model_b <- plm(
  formula = full_formula,
  data = pdata,
  model = "pooling"
)

robust_b <- vcovHC(
  model_b,
  method = "arellano",
  type = "HC1",
  cluster = "group"
)

adjusted_b <- coeftest(
  model_b,
  vcov = robust_b
)

table_b <- data.frame(
  Variable = rownames(adjusted_b),
  Coefficient = round(adjusted_b[,1], 4),
  Robust_SE = round(adjusted_b[,2], 4),
  T_Value = round(adjusted_b[,3], 4),
  P_Value = round(adjusted_b[,4], 4),
  Model = "Adjusted Full Pathway"
)

cat("✓ Adjusted full pathway model complete\n")

# ===============================
# COMBINED TABLE
# ===============================
combined_table <- bind_rows(
  table_a,
  table_b
)

print(combined_table)

# ===============================
# SAVE TABLES TO EXCEL
# ===============================
output_excel <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Tables/H1_Adjusted_Pooled_Results.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "Adjusted_Climate_Only")
writeData(wb, "Adjusted_Climate_Only", table_a)

addWorksheet(wb, "Adjusted_Full_Pathway")
writeData(wb, "Adjusted_Full_Pathway", table_b)

addWorksheet(wb, "Combined_Results")
writeData(wb, "Combined_Results", combined_table)

saveWorkbook(
  wb,
  file = output_excel,
  overwrite = TRUE
)

cat("✓ Excel tables saved\n")

# ===============================
# PREPARE PLOT DATA (CORRECTED)
# ===============================

plot_a <- data.frame(
  Variable = rownames(adjusted_a),
  Coefficient = adjusted_a[,1],
  Robust_SE = adjusted_a[,2],
  P_Value = adjusted_a[,4]
) %>%
  mutate(
    Lower_CI = Coefficient - 1.96 * Robust_SE,
    Upper_CI = Coefficient + 1.96 * Robust_SE,
    
    # CORRECT SIGNIFICANCE CLASSIFICATION
    Significance = ifelse(
      P_Value < 0.05,
      "Significant",
      "Not Significant"
    ),
    
    # FIX LEGEND ORDER
    Significance = factor(
      Significance,
      levels = c("Not Significant", "Significant")
    )
  )

plot_b <- data.frame(
  Variable = rownames(adjusted_b),
  Coefficient = adjusted_b[,1],
  Robust_SE = adjusted_b[,2],
  P_Value = adjusted_b[,4]
) %>%
  mutate(
    Lower_CI = Coefficient - 1.96 * Robust_SE,
    Upper_CI = Coefficient + 1.96 * Robust_SE,
    
    # CORRECT SIGNIFICANCE CLASSIFICATION
    Significance = ifelse(
      P_Value < 0.05,
      "Significant",
      "Not Significant"
    ),
    
    # FIX LEGEND ORDER
    Significance = factor(
      Significance,
      levels = c("Not Significant", "Significant")
    )
  )

# ===============================
# CLEAN VARIABLE LABELS
# ===============================

plot_a$Variable <- recode(
  plot_a$Variable,
  "(Intercept)" = "Intercept",
  "temp_c" = "Temp (°C)",
  "rain_mm" = "Rainfall (mm)"
)

plot_b$Variable <- recode(
  plot_b$Variable,
  "(Intercept)" = "Intercept",
  "temp_c" = "Temp (°C)",
  "rain_mm" = "Rainfall (mm)",
  "ndvi" = "NDVI",
  "water_safe" = "Safe Water",
  "log_gdp" = "Log GDP",
  "log_health_exp" = "Log Health Exp",
  "log_pop_density" = "Log Pop Density",
  "urban_pct" = "Urban (%)"
)

# ===============================
# SHAPE DEFINITIONS
# ===============================

shape_values <- c(
  "Not Significant" = 16,
  "Significant" = 17
)

cat("✓ Plot datasets prepared correctly\n")

# ===============================
# PANEL PLOT A
# ===============================
p1 <- ggplot(
  plot_a,
  aes(
    x = Variable,
    y = Coefficient
  )
) +
  geom_point(
    aes(shape = Significance),
    size = 4.5,
    color = "black"
  ) +
  scale_shape_manual(
    values = shape_values
  ) +
  geom_errorbar(
    aes(
      ymin = Lower_CI,
      ymax = Upper_CI
    ),
    width = 0.25,
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red",
    linewidth = 0.9
  ) +
  labs(
    title = "Adjusted Climate-Only Model",
    x = "Variables",
    y = "Coefficient Estimate"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16,
      margin = margin(b = 12)
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12,
      margin = margin(t = 10)
    ),
    axis.text.y = element_text(
      size = 12
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(
      t = 25,
      r = 20,
      b = 45,
      l = 25
    )
  ) +
  coord_cartesian(clip = "off")

# ===============================
# PANEL PLOT B
# ===============================
p2 <- ggplot(
  plot_b,
  aes(
    x = Variable,
    y = Coefficient
  )
) +
  geom_point(
    aes(shape = Significance),
    size = 4.5,
    color = "black"
  ) +
  scale_shape_manual(
    values = shape_values
  ) +
  geom_errorbar(
    aes(
      ymin = Lower_CI,
      ymax = Upper_CI
    ),
    width = 0.25,
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red",
    linewidth = 0.9
  ) +
  labs(
    title = "Adjusted Full Pathway Model",
    x = "Variables",
    y = "Coefficient Estimate"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16,
      margin = margin(b = 12)
    ),
    axis.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 11,
      margin = margin(t = 10)
    ),
    axis.text.y = element_text(
      size = 12
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(
      t = 25,
      r = 25,
      b = 55,
      l = 30
    )
  ) +
  coord_cartesian(clip = "off")

# ===============================
# DISPLAY PANEL
# ===============================
combined_plot <- grid.arrange(
  p1,
  p2,
  ncol = 2,
  widths = c(1, 1.35),
  top = textGrob(
    "Adjusted Pooled Models: Climate and Full Pathway Predictors of Malaria Incidence",
    gp = gpar(
      fontsize = 20,
      fontface = "bold"
    )
  ),
  padding = unit(3, "lines")
)

# ===============================
# SAVE TIFF
# ===============================
tiff(
  "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/H1_Adjusted_Pooled_Panel.tiff",
  width = 22,
  height = 11,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.arrange(
  p1,
  p2,
  ncol = 2,
  widths = c(1, 1.35),
  top = textGrob(
    "Adjusted Pooled Models: Climate and Full Pathway Predictors of Malaria Incidence",
    gp = gpar(
      fontsize = 20,
      fontface = "bold"
    )
  ),
  padding = unit(3, "lines")
)

dev.off()

# ===============================
# SAVE PDF
# ===============================
pdf(
  "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/H1_Adjusted_Pooled_Panel.pdf",
  width = 22,
  height = 11
)

grid.arrange(
  p1,
  p2,
  ncol = 2,
  widths = c(1, 1.35),
  top = textGrob(
    "Adjusted Pooled Models: Climate and Full Pathway Predictors of Malaria Incidence",
    gp = gpar(
      fontsize = 20,
      fontface = "bold"
    )
  ),
  padding = unit(3, "lines")
)

dev.off()

cat("✓ Journal-quality panel plots saved (TIFF + PDF)\n")














################################################################################
################################################################################
# MULTIPLE PREDICTOR MULTIPLE MEDIATION ANALYSIS
# H2/H3:
#   Temperature → NDVI / Safe Water → Malaria
#   Rainfall    → NDVI / Safe Water → Malaria
#
# Predictors:
#   - temp_c
#   - rain_mm
#
# Mediators:
#   - ndvi
#   - water_safe
#
# Outcome:
#   - log_ind_malaria_inc
#
# Outputs:
# - Separate indirect effects for each predictor-mediator pathway
# - Direct effects
# - Total effects
# - Assumption diagnostics
# - Combined summary tables
# - Bootstrapped mediation significance
# - Journal-quality plots
################################################################################

# ===============================
# LOAD REQUIRED PACKAGES
# ===============================
library(mediation)
library(boot)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(gridExtra)
library(lmtest)
library(car)

# ===============================
# LOAD DATA
# ===============================
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Transformed.rds"

data_transformed <- readRDS(file_path)
data <- na.omit(as.data.frame(data_transformed))

# ===============================
# DEFINE VARIABLES
# ===============================
predictors <- c("temp_c", "rain_mm")
mediators  <- c("ndvi", "water_safe")
outcome    <- "log_ind_malaria_inc"

# ===============================
# OUTPUT PATHS
# ===============================
base_table_path <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Tables/"
base_plot_path  <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/"

dir.create(base_table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(base_plot_path, recursive = TRUE, showWarnings = FALSE)

# ===============================
# EXCEL WORKBOOK
# ===============================
wb <- createWorkbook()

# ===============================
# STORE ALL RESULTS
# ===============================
all_results <- list()

# ===============================
# LOOP THROUGH PREDICTORS + MEDIATORS
# ===============================
for (pred in predictors) {
  
  for (med in mediators) {
    
    cat("\n=====================================\n")
    cat("Predictor:", pred, "| Mediator:", med, "\n")
    cat("=====================================\n")
    
    # ---------------------------
    # MODEL FORMULAS
    # ---------------------------
    mediator_formula <- as.formula(
      paste(med, "~", pred)
    )
    
    outcome_formula <- as.formula(
      paste(outcome, "~", pred, "+", med)
    )
    
    # ---------------------------
    # FIT MODELS
    # ---------------------------
    mediator_model <- lm(
      mediator_formula,
      data = data
    )
    
    outcome_model <- lm(
      outcome_formula,
      data = data
    )
    
    # ---------------------------
    # ASSUMPTION TESTS
    # ---------------------------
    shapiro_med <- shapiro.test(residuals(mediator_model))
    shapiro_out <- shapiro.test(residuals(outcome_model))
    
    bp_med <- bptest(mediator_model)
    bp_out <- bptest(outcome_model)
    
    dw_med <- dwtest(mediator_model)
    dw_out <- dwtest(outcome_model)
    
    assumption_table <- data.frame(
      Predictor = pred,
      Mediator = med,
      Assumption = c(
        "Normality (Mediator)",
        "Normality (Outcome)",
        "Homoscedasticity (Mediator)",
        "Homoscedasticity (Outcome)",
        "Independence (Mediator)",
        "Independence (Outcome)"
      ),
      Statistic = round(c(
        shapiro_med$statistic,
        shapiro_out$statistic,
        bp_med$statistic,
        bp_out$statistic,
        dw_med$statistic,
        dw_out$statistic
      ), 4),
      P_Value = round(c(
        shapiro_med$p.value,
        shapiro_out$p.value,
        bp_med$p.value,
        bp_out$p.value,
        dw_med$p.value,
        dw_out$p.value
      ), 4),
      Decision = c(
        ifelse(shapiro_med$p.value > 0.05,
               "Accepted",
               "Potential Non-Normality"),
        ifelse(shapiro_out$p.value > 0.05,
               "Accepted",
               "Potential Non-Normality"),
        ifelse(bp_med$p.value > 0.05,
               "Accepted",
               "Violated"),
        ifelse(bp_out$p.value > 0.05,
               "Accepted",
               "Violated"),
        ifelse(dw_med$p.value > 0.05,
               "Accepted",
               "Violated"),
        ifelse(dw_out$p.value > 0.05,
               "Accepted",
               "Violated")
      )
    )
    
    # ---------------------------
    # MEDIATION ANALYSIS
    # ---------------------------
    set.seed(123)
    
    med_result <- mediation::mediate(
      mediator_model,
      outcome_model,
      treat = pred,
      mediator = med,
      boot = TRUE,
      sims = 5000
    )
    
    med_summary <- summary(med_result)
    
    # ---------------------------
    # MEDIATION RESULTS TABLE
    # ---------------------------
    mediation_table <- data.frame(
      Predictor = pred,
      Mediator = med,
      Effect = c(
        "ACME (Indirect Effect)",
        "ADE (Direct Effect)",
        "Total Effect",
        "Proportion Mediated"
      ),
      Estimate = round(c(
        med_summary$d0,
        med_summary$z0,
        med_summary$tau.coef,
        med_summary$n0
      ), 4),
      CI_Lower = round(c(
        med_summary$d0.ci[1],
        med_summary$z0.ci[1],
        med_summary$tau.ci[1],
        med_summary$n0.ci[1]
      ), 4),
      CI_Upper = round(c(
        med_summary$d0.ci[2],
        med_summary$z0.ci[2],
        med_summary$tau.ci[2],
        med_summary$n0.ci[2]
      ), 4),
      P_Value = round(c(
        med_summary$d0.p,
        med_summary$z0.p,
        med_summary$tau.p,
        med_summary$n0.p
      ), 4),
      Decision = c(
        ifelse(med_summary$d0.p < 0.05,
               "Significant Indirect Effect",
               "No Significant Indirect Effect"),
        ifelse(med_summary$z0.p < 0.05,
               "Significant Direct Effect",
               "No Significant Direct Effect"),
        ifelse(med_summary$tau.p < 0.05,
               "Significant Total Effect",
               "No Significant Total Effect"),
        ifelse(med_summary$n0.p < 0.05,
               "Significant Mediation",
               "No Significant Mediation")
      )
    )
    
    # ---------------------------
    # SAVE TO WORKBOOK
    # ---------------------------
    sheet_base <- paste(pred, med, sep = "_")
    
    addWorksheet(wb, paste0(sheet_base, "_Assumptions"))
    writeData(wb, paste0(sheet_base, "_Assumptions"), assumption_table)
    
    addWorksheet(wb, paste0(sheet_base, "_Mediation"))
    writeData(wb, paste0(sheet_base, "_Mediation"), mediation_table)
    
    # ---------------------------
    # STORE RESULTS
    # ---------------------------
    all_results[[sheet_base]] <- mediation_table
    
    # ---------------------------
    # BOOTSTRAP PLOT
    # ---------------------------
    boot_data <- data.frame(
      Indirect_Effect = med_result$d.avg.sims
    )
    
    boot_plot <- ggplot(
      boot_data,
      aes(x = Indirect_Effect)
    ) +
      geom_histogram(
        bins = 40,
        fill = "steelblue",
        alpha = 0.8
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "red",
        linewidth = 1
      ) +
      labs(
        title = paste("Bootstrapped Indirect Effect:", pred, "→", med),
        subtitle = paste(pred, "→", med, "→ Malaria"),
        x = "Indirect Effect Estimate",
        y = "Frequency"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin = margin(20, 20, 20, 20)
      )
    
    # TIFF
    tiff(
      paste0(base_plot_path, sheet_base, "_Bootstrap_Mediation.tiff"),
      width = 12,
      height = 8,
      units = "in",
      res = 600,
      compression = "lzw"
    )
    print(boot_plot)
    dev.off()
    
    # PDF
    pdf(
      paste0(base_plot_path, sheet_base, "_Bootstrap_Mediation.pdf"),
      width = 12,
      height = 8
    )
    print(boot_plot)
    dev.off()
    
    # CSV
    write.csv(
      mediation_table,
      paste0(base_table_path, sheet_base, "_Mediation_Results.csv"),
      row.names = FALSE
    )
    
    write.csv(
      assumption_table,
      paste0(base_table_path, sheet_base, "_Assumptions.csv"),
      row.names = FALSE
    )
  }
}

# ===============================
# COMBINED SUMMARY TABLE
# ===============================
combined_results <- bind_rows(all_results)

addWorksheet(wb, "Combined_Mediation_Summary")
writeData(wb, "Combined_Mediation_Summary", combined_results)

# ===============================
# SAVE WORKBOOK
# ===============================
saveWorkbook(
  wb,
  paste0(base_table_path, "Multiple_Predictor_Mediation_Analysis.xlsx"),
  overwrite = TRUE
)

# ===============================
# SAVE COMBINED CSV
# ===============================
write.csv(
  combined_results,
  paste0(base_table_path, "Combined_Mediation_Summary.csv"),
  row.names = FALSE
)

# ===============================
# FINAL MESSAGE
# ===============================
message("Complete multiple predictor mediation analysis finished. All outputs saved.")



################################################################################
# MULTIPLE PREDICTOR MULTIPLE MEDIATION ANALYSIS
# Climate Predictors:
#   - Temperature
#   - Rainfall
#
# Mediators:
#   - NDVI
#   - Safe Water Access
#
# Outcome:
#   - Malaria Incidence
#
# Includes:
# - Separate mediation pathways
# - Assumption diagnostics
# - Bootstrapped mediation tests
# - Combined Excel tables
# - Journal-quality panel plots
################################################################################

# ===============================
# LOAD REQUIRED PACKAGES
# ===============================
library(mediation)
library(boot)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(gridExtra)
library(lmtest)
library(car)

# ===============================
# LOAD DATA
# ===============================
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Transformed.rds"

data_transformed <- readRDS(file_path)
data <- na.omit(as.data.frame(data_transformed))

# ===============================
# DEFINE VARIABLES
# ===============================
predictors <- c("temp_c", "rain_mm")
mediators  <- c("ndvi", "water_safe")
outcome    <- "log_ind_malaria_inc"

# ===============================
# OUTPUT PATHS
# ===============================
base_table_path <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Tables/"
base_plot_path  <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/"

dir.create(base_table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(base_plot_path, recursive = TRUE, showWarnings = FALSE)

# ===============================
# CREATE WORKBOOK
# ===============================
wb <- createWorkbook()

# ===============================
# STORAGE
# ===============================
all_results <- list()
panel_plots <- list()

# ===============================
# LOOP THROUGH ALL PATHWAYS
# ===============================
plot_index <- 1

for (pred in predictors) {
  
  for (med in mediators) {
    
    cat("\n=====================================\n")
    cat("Predictor:", pred, "| Mediator:", med, "\n")
    cat("=====================================\n")
    
    # ---------------------------
    # FORMULAS
    # ---------------------------
    mediator_formula <- as.formula(
      paste(med, "~", pred)
    )
    
    outcome_formula <- as.formula(
      paste(outcome, "~", pred, "+", med)
    )
    
    # ---------------------------
    # FIT MODELS
    # ---------------------------
    mediator_model <- lm(
      mediator_formula,
      data = data
    )
    
    outcome_model <- lm(
      outcome_formula,
      data = data
    )
    
    # ---------------------------
    # ASSUMPTION TESTS
    # ---------------------------
    shapiro_med <- shapiro.test(residuals(mediator_model))
    shapiro_out <- shapiro.test(residuals(outcome_model))
    
    bp_med <- bptest(mediator_model)
    bp_out <- bptest(outcome_model)
    
    dw_med <- dwtest(mediator_model)
    dw_out <- dwtest(outcome_model)
    
    assumption_table <- data.frame(
      Predictor = pred,
      Mediator = med,
      Assumption = c(
        "Normality (Mediator)",
        "Normality (Outcome)",
        "Homoscedasticity (Mediator)",
        "Homoscedasticity (Outcome)",
        "Independence (Mediator)",
        "Independence (Outcome)"
      ),
      Statistic = round(c(
        shapiro_med$statistic,
        shapiro_out$statistic,
        bp_med$statistic,
        bp_out$statistic,
        dw_med$statistic,
        dw_out$statistic
      ), 4),
      P_Value = round(c(
        shapiro_med$p.value,
        shapiro_out$p.value,
        bp_med$p.value,
        bp_out$p.value,
        dw_med$p.value,
        dw_out$p.value
      ), 4),
      Decision = c(
        ifelse(shapiro_med$p.value > 0.05,
               "Accepted",
               "Potential Non-Normality"),
        ifelse(shapiro_out$p.value > 0.05,
               "Accepted",
               "Potential Non-Normality"),
        ifelse(bp_med$p.value > 0.05,
               "Accepted",
               "Violated"),
        ifelse(bp_out$p.value > 0.05,
               "Accepted",
               "Violated"),
        ifelse(dw_med$p.value > 0.05,
               "Accepted",
               "Violated"),
        ifelse(dw_out$p.value > 0.05,
               "Accepted",
               "Violated")
      )
    )
    
    # ---------------------------
    # MEDIATION ANALYSIS
    # ---------------------------
    set.seed(123)
    
    med_result <- mediation::mediate(
      mediator_model,
      outcome_model,
      treat = pred,
      mediator = med,
      boot = TRUE,
      sims = 5000
    )
    
    med_summary <- summary(med_result)
    
    # ---------------------------
    # RESULTS TABLE
    # ---------------------------
    mediation_table <- data.frame(
      Predictor = pred,
      Mediator = med,
      Effect = c(
        "ACME (Indirect Effect)",
        "ADE (Direct Effect)",
        "Total Effect",
        "Proportion Mediated"
      ),
      Estimate = round(c(
        med_summary$d0,
        med_summary$z0,
        med_summary$tau.coef,
        med_summary$n0
      ), 4),
      CI_Lower = round(c(
        med_summary$d0.ci[1],
        med_summary$z0.ci[1],
        med_summary$tau.ci[1],
        med_summary$n0.ci[1]
      ), 4),
      CI_Upper = round(c(
        med_summary$d0.ci[2],
        med_summary$z0.ci[2],
        med_summary$tau.ci[2],
        med_summary$n0.ci[2]
      ), 4),
      P_Value = round(c(
        med_summary$d0.p,
        med_summary$z0.p,
        med_summary$tau.p,
        med_summary$n0.p
      ), 4),
      Decision = c(
        ifelse(med_summary$d0.p < 0.05,
               "Significant Indirect Effect",
               "No Significant Indirect Effect"),
        ifelse(med_summary$z0.p < 0.05,
               "Significant Direct Effect",
               "No Significant Direct Effect"),
        ifelse(med_summary$tau.p < 0.05,
               "Significant Total Effect",
               "No Significant Total Effect"),
        ifelse(med_summary$n0.p < 0.05,
               "Significant Mediation",
               "No Significant Mediation")
      )
    )
    
    # ---------------------------
    # SAVE TABLES
    # ---------------------------
    sheet_base <- paste(pred, med, sep = "_")
    
    addWorksheet(wb, paste0(sheet_base, "_Assumptions"))
    writeData(wb, paste0(sheet_base, "_Assumptions"), assumption_table)
    
    addWorksheet(wb, paste0(sheet_base, "_Mediation"))
    writeData(wb, paste0(sheet_base, "_Mediation"), mediation_table)
    
    all_results[[sheet_base]] <- mediation_table
    
    # ---------------------------
    # BOOTSTRAP PLOT
    # ---------------------------
    boot_data <- data.frame(
      Indirect_Effect = med_result$d.avg.sims
    )
    
    panel_plots[[plot_index]] <- ggplot(
      boot_data,
      aes(x = Indirect_Effect)
    ) +
      geom_histogram(
        bins = 35,
        fill = "steelblue",
        alpha = 0.8
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "red",
        linewidth = 0.8
      ) +
      labs(
        title = paste(pred, "→", med),
        subtitle = "Bootstrapped Indirect Effect",
        x = "Indirect Effect",
        y = "Frequency"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 13
        ),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 10
        ),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
      )
    
    plot_index <- plot_index + 1
    
    # CSV
    write.csv(
      mediation_table,
      paste0(base_table_path, sheet_base, "_Mediation_Results.csv"),
      row.names = FALSE
    )
    
    write.csv(
      assumption_table,
      paste0(base_table_path, sheet_base, "_Assumptions.csv"),
      row.names = FALSE
    )
  }
}

# ===============================
# COMBINED SUMMARY
# ===============================
combined_results <- bind_rows(all_results)

addWorksheet(wb, "Combined_Mediation_Summary")
writeData(wb, "Combined_Mediation_Summary", combined_results)

# ===============================
# SAVE EXCEL
# ===============================
saveWorkbook(
  wb,
  paste0(base_table_path, "Multiple_Predictor_Mediation_Analysis.xlsx"),
  overwrite = TRUE
)

# ===============================
# SAVE COMBINED CSV
# ===============================
write.csv(
  combined_results,
  paste0(base_table_path, "Combined_Mediation_Summary.csv"),
  row.names = FALSE
)

# ===============================
# JOURNAL-QUALITY PANEL PLOT
# ===============================

# TIFF
tiff(
  paste0(base_plot_path, "Combined_Mediation_Panel_Plots.tiff"),
  width = 18,
  height = 14,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.arrange(
  grobs = panel_plots,
  ncol = 2,
  top = textGrob(
    "Bootstrapped Indirect Effects for Climate-Malaria Mediation Pathways",
    gp = gpar(fontsize = 18, fontface = "bold")
  )
)

dev.off()

# PDF
pdf(
  paste0(base_plot_path, "Combined_Mediation_Panel_Plots.pdf"),
  width = 18,
  height = 14
)

grid.arrange(
  grobs = panel_plots,
  ncol = 2,
  top = textGrob(
    "Bootstrapped Indirect Effects for Climate-Malaria Mediation Pathways",
    gp = gpar(fontsize = 18, fontface = "bold")
  )
)

dev.off()

cat("✓ Journal-quality panel plots saved\n")

# ===============================
# FINAL MESSAGE
# ===============================
message("Multiple predictor mediation analysis complete. All outputs saved.")


################################################################################
################################################################################
# H3 & H4: MODERATION ANALYSIS
#
# H3: Socioeconomic Adaptive Capacity Moderation
#   Moderators:
#   - log_gdp
#   - log_health_exp
#
# H4: Demographic Moderation
#   Moderators:
#   - log_pop_density
#   - urban_pct
#
# Predictors:
#   - temp_c
#   - rain_mm
#
# Outcome:
#   - log_ind_malaria_inc
#
# Model:
#   Panel pooled regression with robust clustered SE
#
# Outputs:
# - Individual moderation models
# - Combined summary tables
# - VIF diagnostics
# - Serial correlation
# - Heteroskedasticity
# - Journal-quality forest plots
# - Panel plot outputs
# - Excel + CSV exports
################################################################################

# ===============================
# LOAD REQUIRED PACKAGES
# ===============================
library(plm)
library(lmtest)
library(sandwich)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(car)
library(gridExtra)
library(grid)

# ===============================
# LOAD DATA
# ===============================
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Transformed.rds"

data_transformed <- readRDS(file_path)
data <- na.omit(as.data.frame(data_transformed))

pdata <- pdata.frame(
  data,
  index = c("country", "year")
)

cat("✓ Data loaded successfully\n")

# ===============================
# DEFINE MODERATION MODELS
# ===============================
moderation_models <- list(
  
  # H3: GDP
  "Temp_x_GDP" =
    log_ind_malaria_inc ~ temp_c + rain_mm + log_gdp +
    temp_c:log_gdp + rain_mm:log_gdp,
  
  # H3: Health expenditure
  "Temp_x_HealthExp" =
    log_ind_malaria_inc ~ temp_c + rain_mm + log_health_exp +
    temp_c:log_health_exp + rain_mm:log_health_exp,
  
  # H4: Population density
  "Temp_x_PopDensity" =
    log_ind_malaria_inc ~ temp_c + rain_mm + log_pop_density +
    temp_c:log_pop_density + rain_mm:log_pop_density,
  
  # H4: Urbanization
  "Temp_x_UrbanPct" =
    log_ind_malaria_inc ~ temp_c + rain_mm + urban_pct +
    temp_c:urban_pct + rain_mm:urban_pct
)

# ===============================
# OUTPUT PATHS
# ===============================
table_path <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Tables/"
plot_path  <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/"

dir.create(table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

# ===============================
# STORAGE OBJECTS
# ===============================
all_results <- list()
all_vif <- list()
plot_list <- list()

# ===============================
# WORKBOOK
# ===============================
wb <- createWorkbook()

# ===============================
# RUN MODELS
# ===============================
plot_index <- 1

for (model_name in names(moderation_models)) {
  
  cat("\n=============================\n")
  cat("Running:", model_name, "\n")
  cat("=============================\n")
  
  formula_current <- moderation_models[[model_name]]
  
  # ---------------------------
  # MODEL
  # ---------------------------
  model <- plm(
    formula = formula_current,
    data = pdata,
    model = "pooling"
  )
  
  # ---------------------------
  # ROBUST CLUSTERED SE
  # ---------------------------
  robust_cov <- vcovHC(
    model,
    method = "arellano",
    type = "HC1",
    cluster = "group"
  )
  
  robust_results <- coeftest(
    model,
    vcov = robust_cov
  )
  
  # ---------------------------
  # RESULTS TABLE
  # ---------------------------
  results_table <- data.frame(
    Variable = rownames(robust_results),
    Coefficient = round(robust_results[,1], 4),
    Robust_SE = round(robust_results[,2], 4),
    T_Value = round(robust_results[,3], 4),
    P_Value = round(robust_results[,4], 4),
    Significance = ifelse(
      robust_results[,4] < 0.05,
      "Significant",
      "Not Significant"
    ),
    Model = model_name
  )
  
  all_results[[model_name]] <- results_table
  
  # ---------------------------
  # VIF
  # ---------------------------
  vif_model <- lm(
    formula_current,
    data = data
  )
  
  vif_vals <- vif(vif_model)
  
  vif_table <- data.frame(
    Variable = names(vif_vals),
    VIF = round(vif_vals, 4),
    Model = model_name
  )
  
  all_vif[[model_name]] <- vif_table
  
  # ---------------------------
  # SERIAL CORRELATION
  # ---------------------------
  serial_test <- pbgtest(model)
  
  serial_table <- data.frame(
    Test = "Serial Correlation",
    Statistic = round(as.numeric(serial_test$statistic), 4),
    P_Value = round(as.numeric(serial_test$p.value), 4),
    Decision = ifelse(
      serial_test$p.value < 0.05,
      "Present",
      "Absent"
    )
  )
  
  # ---------------------------
  # HETEROSKEDASTICITY
  # ---------------------------
  hetero_test <- bptest(
    model,
    studentize = FALSE
  )
  
  hetero_table <- data.frame(
    Test = "Heteroskedasticity",
    Statistic = round(as.numeric(hetero_test$statistic), 4),
    P_Value = round(as.numeric(hetero_test$p.value), 4),
    Decision = ifelse(
      hetero_test$p.value < 0.05,
      "Present",
      "Absent"
    )
  )
  
  # ---------------------------
  # SAVE TO EXCEL
  # ---------------------------
  addWorksheet(wb, paste0(model_name, "_Results"))
  writeData(wb, paste0(model_name, "_Results"), results_table)
  
  addWorksheet(wb, paste0(model_name, "_VIF"))
  writeData(wb, paste0(model_name, "_VIF"), vif_table)
  
  addWorksheet(wb, paste0(model_name, "_Serial"))
  writeData(wb, paste0(model_name, "_Serial"), serial_table)
  
  addWorksheet(wb, paste0(model_name, "_Hetero"))
  writeData(wb, paste0(model_name, "_Hetero"), hetero_table)
  
  # ---------------------------
  # PLOT DATA
  # ---------------------------
  plot_data <- results_table %>%
    filter(grepl(":", Variable)) %>%
    mutate(
      Lower_CI = Coefficient - 1.96 * Robust_SE,
      Upper_CI = Coefficient + 1.96 * Robust_SE
    )
  
  # ---------------------------
  # FOREST PLOT
  # ---------------------------
  plot_list[[plot_index]] <- ggplot(
    plot_data,
    aes(
      x = Variable,
      y = Coefficient
    )
  ) +
    geom_point(
      aes(shape = Significance),
      size = 4.5,
      color = "black"
    ) +
    geom_errorbar(
      aes(
        ymin = Lower_CI,
        ymax = Upper_CI
      ),
      width = 0.2,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "red",
      linewidth = 0.9
    ) +
    labs(
      title = model_name,
      x = "Interaction Terms",
      y = "Coefficient Estimate"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 15
      ),
      axis.title = element_text(
        face = "bold",
        size = 13
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 11,
        margin = margin(t = 10)
      ),
      axis.text.y = element_text(size = 11),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(
        t = 20,
        r = 20,
        b = 45,
        l = 25
      )
    ) +
    coord_cartesian(clip = "off")
  
  plot_index <- plot_index + 1
}

# ===============================
# COMBINE RESULTS
# ===============================
combined_results <- bind_rows(all_results)
combined_vif <- bind_rows(all_vif)

addWorksheet(wb, "Combined_Moderation_Results")
writeData(wb, "Combined_Moderation_Results", combined_results)

addWorksheet(wb, "Combined_VIF")
writeData(wb, "Combined_VIF", combined_vif)

# ===============================
# SAVE WORKBOOK
# ===============================
saveWorkbook(
  wb,
  file = paste0(table_path, "H3_H4_Moderation_Analysis.xlsx"),
  overwrite = TRUE
)

cat("✓ Excel workbook saved\n")

# ===============================
# SAVE CSV
# ===============================
write.csv(
  combined_results,
  paste0(table_path, "H3_H4_Combined_Moderation_Results.csv"),
  row.names = FALSE
)

# ===============================
# PANEL FOREST PLOTS
# ===============================

# TIFF
tiff(
  paste0(plot_path, "H3_H4_Moderation_Panel.tiff"),
  width = 24,
  height = 16,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.arrange(
  grobs = plot_list,
  ncol = 2,
  top = textGrob(
    "Moderation Effects of Socioeconomic and Demographic Factors on Climate-Malaria Relationships",
    gp = gpar(
      fontsize = 22,
      fontface = "bold"
    )
  ),
  padding = unit(3, "lines")
)

dev.off()

# PDF
pdf(
  paste0(plot_path, "H3_H4_Moderation_Panel.pdf"),
  width = 24,
  height = 16
)

grid.arrange(
  grobs = plot_list,
  ncol = 2,
  top = textGrob(
    "Moderation Effects of Socioeconomic and Demographic Factors on Climate-Malaria Relationships",
    gp = gpar(
      fontsize = 22,
      fontface = "bold"
    )
  ),
  padding = unit(3, "lines")
)

dev.off()

cat("✓ Journal-quality moderation panel plots saved\n")

# ===============================
# FINAL MESSAGE
# ===============================
message("H3 and H4 moderation analyses complete.")
getwd()






################################################################################
# H5: REGIONAL HETEROGENEITY ANALYSIS 
# INCLUDED REGIONS:
# - Africa
# - Asia
# - North America
# - South America
#
# ADDED IMPROVEMENTS:
# 1. Saves fully cleaned regional dataset as:
#    - Malaria_Data_Regional_Cleaned.rds
#    - Malaria_Data_Regional_Cleaned.xlsx
# 2. Ensures all four continents remain visible in:
#    - Combined regional results
#    - Combined VIF
#    - Summary tables
#    - Panel plots
# 3. Saves continent distribution summary
################################################################################

# ===============================
# LOAD REQUIRED PACKAGES
# ===============================
library(plm)
library(lmtest)
library(sandwich)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(car)
library(gridExtra)
library(grid)
library(stringr)
library(forcats)
library(tidyr)

# ===============================
# LOAD DATA
# ===============================
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Transformed.rds"

data_transformed <- readRDS(file_path)

# Preserve full dataset
data <- as.data.frame(data_transformed)

cat("✓ Data loaded successfully\n")

# ===============================
# CLEAN COUNTRY NAMES
# ===============================
data$country <- as.character(data$country)
data$country <- trimws(data$country)
data$country <- str_squish(data$country)
data$country <- str_to_title(data$country)

# ===============================
# CLEAN REGION NAMES
# ===============================
data$region <- as.character(data$region)
data$region <- trimws(data$region)
data$region <- str_squish(data$region)
data$region <- str_replace_all(data$region, "_", " ")
data$region <- str_to_title(data$region)

# ===============================
# STANDARDIZE AFRICA + ASIA
# ===============================
data$region[data$region %in% c(
  "Africa",
  "Sub-Saharan Africa",
  "Eastern Africa",
  "Western Africa",
  "Southern Africa",
  "Central Africa",
  "Northern Africa"
)] <- "Africa"

data$region[data$region %in% c(
  "Asia",
  "South Asia",
  "East Asia",
  "Southeast Asia",
  "Western Asia",
  "Central Asia"
)] <- "Asia"

# ===============================
# NORTH AMERICA COUNTRY ASSIGNMENT
# ===============================
north_america_countries <- c(
  "Costa Rica",
  "Dominican Republic",
  "Guatemala",
  "Honduras",
  "Mexico"
)

data$region[data$country %in% north_america_countries] <- "North America"

# ===============================
# SOUTH AMERICA COUNTRY ASSIGNMENT
# ===============================
south_america_countries <- c(
  "Brazil",
  "Colombia",
  "Ecuador",
  "Guyana",
  "Peru",
  "Suriname"
)

data$region[data$country %in% south_america_countries] <- "South America"

# ===============================
# TARGET REGIONS
# ===============================
target_regions <- c(
  "Africa",
  "Asia",
  "North America",
  "South America"
)

# ===============================
# FILTER TARGET REGIONS ONLY
# ===============================
data <- data %>%
  filter(region %in% target_regions)

# ===============================
# FACTOR ORDER
# ===============================
data$region <- factor(
  data$region,
  levels = target_regions
)

# ===============================
# REGIONAL SUMMARY SECTION
# ===============================
regional_summary <- data %>%
  group_by(region) %>%
  summarise(
    Countries = n_distinct(country),
    Observations = n(),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    region = factor(target_regions, levels = target_regions),
    fill = list(
      Countries = 0,
      Observations = 0
    )
  ) %>%
  arrange(region)

cat("\n✓ Final regional distribution:\n")
print(regional_summary)
print(regional_summary)

# ===============================
# SAVE CLEANED REGIONAL DATASET
# ===============================
saveRDS(
  data,
  "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Regional_Cleaned.rds"
)

# ===============================
# MODEL FORMULA
# ===============================
regional_formula <- log_ind_malaria_inc ~
  temp_c +
  rain_mm +
  ndvi +
  water_safe +
  log_gdp +
  log_health_exp +
  log_pop_density +
  urban_pct

# ===============================
# OUTPUT PATHS
# ===============================
table_path <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Tables/"
plot_path  <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/"

dir.create(table_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)

# ===============================
# STORAGE
# ===============================
all_results <- list()
all_vif <- list()
plot_list <- list()

# ===============================
# CREATE EXCEL WORKBOOK
# ===============================
wb <- createWorkbook()

# ===============================
# SAVE CLEANED DATA SHEETS
# ===============================
addWorksheet(wb, "Regional_Cleaned_Data")
writeData(wb, "Regional_Cleaned_Data", data)

addWorksheet(wb, "Regional_Summary")
writeData(wb, "Regional_Summary", regional_summary)

# ===============================
# LOOP THROUGH REGIONS
# ===============================
for(region_name in target_regions){
  
  cat("\nAnalyzing:", region_name, "\n")
  
  region_data <- data %>%
    filter(region == region_name)
  
  if(nrow(region_data) == 0){
    
    empty_table <- data.frame(
      Variable = NA,
      Coefficient = NA,
      Robust_SE = NA,
      T_Value = NA,
      P_Value = NA,
      Significance = NA,
      Region = region_name
    )
    
    all_results[[region_name]] <- empty_table
    next
  }
  
  pdata_region <- pdata.frame(
    region_data,
    index = c("country", "year")
  )
  
  # ===========================
  # MODEL
  # ===========================
  model <- plm(
    formula = regional_formula,
    data = pdata_region,
    model = "pooling"
  )
  
  robust_cov <- vcovHC(
    model,
    method = "arellano",
    type = "HC1",
    cluster = "group"
  )
  
  robust_results <- coeftest(
    model,
    vcov = robust_cov
  )
  
  # ===========================
  # RESULTS TABLE
  # ===========================
  results_table <- data.frame(
    Variable = rownames(robust_results),
    Coefficient = round(robust_results[,1], 4),
    Robust_SE = round(robust_results[,2], 4),
    T_Value = round(robust_results[,3], 4),
    P_Value = round(robust_results[,4], 4),
    Significance = ifelse(
      robust_results[,4] < 0.05,
      "Significant",
      "Not Significant"
    ),
    Region = region_name
  )
  
  all_results[[region_name]] <- results_table
  
  # ===========================
  # VIF
  # ===========================
  vif_model <- lm(
    regional_formula,
    data = region_data
  )
  
  vif_vals <- vif(vif_model)
  
  vif_table <- data.frame(
    Variable = names(vif_vals),
    VIF = round(vif_vals, 4),
    Decision = ifelse(
      vif_vals < 5,
      "Acceptable",
      ifelse(vif_vals < 10, "Moderate", "High")
    ),
    Region = region_name
  )
  
  all_vif[[region_name]] <- vif_table
  
  # ===========================
  # DIAGNOSTICS
  # ===========================
  serial_test <- pbgtest(model)
  hetero_test <- bptest(model, studentize = FALSE)
  
  serial_table <- data.frame(
    Test = "Breusch-Godfrey Serial Correlation",
    Statistic = round(as.numeric(serial_test$statistic), 4),
    P_Value = round(as.numeric(serial_test$p.value), 4),
    Decision = ifelse(serial_test$p.value < 0.05, "Present", "Absent"),
    Region = region_name
  )
  
  hetero_table <- data.frame(
    Test = "Breusch-Pagan Heteroskedasticity",
    Statistic = round(as.numeric(hetero_test$statistic), 4),
    P_Value = round(as.numeric(hetero_test$p.value), 4),
    Decision = ifelse(hetero_test$p.value < 0.05, "Present", "Absent"),
    Region = region_name
  )
  
  # ===========================
  # EXPORT REGION SHEETS
  # ===========================
  addWorksheet(wb, paste0(region_name, "_Results"))
  writeData(wb, paste0(region_name, "_Results"), results_table)
  
  addWorksheet(wb, paste0(region_name, "_VIF"))
  writeData(wb, paste0(region_name, "_VIF"), vif_table)
  
  addWorksheet(wb, paste0(region_name, "_Serial"))
  writeData(wb, paste0(region_name, "_Serial"), serial_table)
  
  addWorksheet(wb, paste0(region_name, "_Hetero"))
  writeData(wb, paste0(region_name, "_Hetero"), hetero_table)
  
  # ===========================
  # PLOT DATA
  # ===========================
  plot_data <- results_table %>%
    filter(Variable != "(Intercept)") %>%
    mutate(
      Lower_CI = Coefficient - 1.96 * Robust_SE,
      Upper_CI = Coefficient + 1.96 * Robust_SE
    )
  
  plot_list[[region_name]] <- ggplot(
    plot_data,
    aes(x = Variable, y = Coefficient)
  ) +
    geom_point(
      aes(shape = Significance),
      size = 4.5,
      color = "black"
    ) +
    geom_errorbar(
      aes(ymin = Lower_CI, ymax = Upper_CI),
      width = 0.2,
      linewidth = 1
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "red"
    ) +
    labs(
      title = region_name,
      x = "Predictors",
      y = "Coefficient Estimate"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 16
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 11,
        margin = margin(t = 12)
      ),
      plot.margin = margin(
        t = 25,
        r = 25,
        b = 55,
        l = 30
      ),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(clip = "off")
}

# ===============================
# COMBINED OUTPUTS
# ===============================
combined_results <- bind_rows(all_results) %>%
  mutate(
    Region = factor(Region, levels = target_regions)
  ) %>%
  arrange(Region)

combined_vif <- bind_rows(all_vif) %>%
  mutate(
    Region = factor(Region, levels = target_regions)
  ) %>%
  arrange(Region)

# ===============================
# SAVE COMBINED SHEETS
# ===============================
addWorksheet(wb, "Combined_Regional_Results")
writeData(wb, "Combined_Regional_Results", combined_results)

addWorksheet(wb, "Combined_VIF")
writeData(wb, "Combined_VIF", combined_vif)

# ===============================
# SAVE FINAL EXCEL
# ===============================
saveWorkbook(
  wb,
  paste0(table_path, "H5_Regional_Heterogeneity_Results.xlsx"),
  overwrite = TRUE
)

cat("✓ Combined Excel workbook saved successfully\n")

# ===============================
# PANEL PLOTS
# ===============================
ordered_plots <- plot_list[target_regions]
ordered_plots <- ordered_plots[!sapply(ordered_plots, is.null)]
ordered_grobs <- lapply(
  ordered_plots,
  ggplotGrob
)

# TIFF
tiff(
  paste0(plot_path, "H5_Regional_Heterogeneity_Panel.tiff"),
  width = 26,
  height = 20,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.arrange(
  grobs = ordered_grobs,
  ncol = 2,
  top = textGrob(
    "Regional Differences in Climate, Environmental, Adaptive Capacity, and Demographic Effects on Malaria Incidence",
    gp = gpar(
      fontsize = 24,
      fontface = "bold"
    )
  ),
  padding = unit(4, "lines")
)

dev.off()

# PDF
pdf(
  paste0(plot_path, "H5_Regional_Heterogeneity_Panel.pdf"),
  width = 26,
  height = 20
)

grid.arrange(
  grobs = ordered_grobs,
  ncol = 2,
  top = textGrob(
    "Regional Differences in Climate, Environmental, Adaptive Capacity, and Demographic Effects on Malaria Incidence",
    gp = gpar(
      fontsize = 24,
      fontface = "bold"
    )
  ),
  padding = unit(4, "lines")
)

dev.off()

cat("✓ Regional panel plots saved successfully\n")

# ===============================
# FINAL MESSAGE
# ===============================
message("H5 regional heterogeneity analysis complete.")
getwd()



################################################################################
################################################################################
# MALARIA VULNERABILITY INDEX (MVI) CONSTRUCTION SCRIPT
#
# IPCC Framework:
# Vulnerability = Exposure + Sensitivity + Adaptive Capacity
#
# INCLUDED:
# Exposure:
#   - temp_c
#   - rain_mm
#
# Sensitivity:
#   - ndvi
#   - water_safe (reverse coded)
#
# Adaptive Capacity:
#   - log_gdp (reverse coded)
#   - log_health_exp (reverse coded)
#
# EXCLUDED:
#   - log_pop_density
#   - urban_pct
#
# OUTPUTS:
# - Reverse coding
# - Normalization
# - Sub-indices
# - Final MVI
# - Validation diagnostics
# - PCA
# - Cronbach alpha
# - Summary tables
# - Journal-quality plots
# - Excel outputs
################################################################################

# ===============================
# LOAD REQUIRED PACKAGES
# ===============================
library(dplyr)
library(ggplot2)
library(openxlsx)
library(psych)
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(grid)
library(scales)

# ===============================
# LOAD DATA
# ===============================
file_path <- "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_Regional_Cleaned.rds"
data <- readRDS(file_path)

cat("✓ Data loaded successfully\n")

# ===============================
# DEFINE MIN-MAX FUNCTION
# ===============================
min_max <- function(x){
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# ===============================
# STEP 1: REVERSE CODING
# ===============================
data <- data %>%
  mutate(
    rev_water_safe      = max(water_safe, na.rm = TRUE) - water_safe,
    rev_log_gdp         = max(log_gdp, na.rm = TRUE) - log_gdp,
    rev_log_health_exp  = max(log_health_exp, na.rm = TRUE) - log_health_exp
  )

cat("✓ Reverse coding complete\n")

# ===============================
# STEP 2: NORMALIZATION
# ===============================
data <- data %>%
  mutate(
    norm_temp          = min_max(temp_c),
    norm_rain          = min_max(rain_mm),
    norm_ndvi          = min_max(ndvi),
    norm_water         = min_max(rev_water_safe),
    norm_gdp           = min_max(rev_log_gdp),
    norm_health        = min_max(rev_log_health_exp)
  )

cat("✓ Normalization complete\n")

# ===============================
# STEP 3: SUB-INDICES
# ===============================
data <- data %>%
  mutate(
    Exposure_Index = (norm_temp + norm_rain) / 2,
    Sensitivity_Index = (norm_ndvi + norm_water) / 2,
    Adaptive_Capacity_Index = (norm_gdp + norm_health) / 2
  )

cat("✓ Sub-indices constructed\n")

# ===============================
# STEP 4: FINAL MVI
# ===============================
data <- data %>%
  mutate(
    MVI = (
      Exposure_Index +
        Sensitivity_Index +
        Adaptive_Capacity_Index
    ) / 3
  )

cat("✓ Final MVI calculated\n")

# ===============================
# STEP 5: MVI CATEGORIZATION
# ===============================
data <- data %>%
  mutate(
    MVI_Category = case_when(
      MVI <= quantile(MVI, 0.33, na.rm = TRUE) ~ "Low",
      MVI <= quantile(MVI, 0.66, na.rm = TRUE) ~ "Moderate",
      TRUE ~ "High"
    )
  )

# ===============================
# STEP 6: VALIDATION
# ===============================

# Variables used in index
mvi_vars <- data %>%
  select(
    norm_temp,
    norm_rain,
    norm_ndvi,
    norm_water,
    norm_gdp,
    norm_health
  )

# Cronbach Alpha
alpha_result <- psych::alpha(mvi_vars)

alpha_table <- data.frame(
  Cronbach_Alpha = round(alpha_result$total$raw_alpha, 4)
)

cat("✓ Cronbach alpha complete\n")

# PCA
pca_result <- PCA(
  mvi_vars,
  graph = FALSE
)

pca_variance <- data.frame(
  Component = paste0("PC", 1:length(pca_result$eig[,1])),
  Eigenvalue = round(pca_result$eig[,1], 4),
  Variance_Percent = round(pca_result$eig[,2], 4),
  Cumulative_Percent = round(pca_result$eig[,3], 4)
)

cat("✓ PCA complete\n")

# ===============================
# STEP 7: CORRELATION WITH MALARIA
# ===============================
mvi_validation <- cor.test(
  data$MVI,
  data$log_ind_malaria_inc,
  use = "complete.obs"
)

validation_table <- data.frame(
  Correlation = round(as.numeric(mvi_validation$estimate), 4),
  P_Value = round(mvi_validation$p.value, 4),
  Decision = ifelse(
    mvi_validation$p.value < 0.05,
    "Significant",
    "Not Significant"
  )
)

# ===============================
# STEP 8: SUMMARY TABLES
# ===============================
summary_table <- data %>%
  group_by(region) %>%
  summarise(
    Mean_MVI = round(mean(MVI, na.rm = TRUE), 4),
    SD_MVI = round(sd(MVI, na.rm = TRUE), 4),
    Min_MVI = round(min(MVI, na.rm = TRUE), 4),
    Max_MVI = round(max(MVI, na.rm = TRUE), 4),
    Countries = n_distinct(country),
    Observations = n()
  )

# ===============================
# STEP 9: EXPORT TO EXCEL
# ===============================
output_excel <- "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Tables/MVI_Results.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "MVI_Data")
writeData(wb, "MVI_Data", data)

addWorksheet(wb, "Regional_Summary")
writeData(wb, "Regional_Summary", summary_table)

addWorksheet(wb, "Cronbach_Alpha")
writeData(wb, "Cronbach_Alpha", alpha_table)

addWorksheet(wb, "PCA")
writeData(wb, "PCA", pca_variance)

addWorksheet(wb, "Validation")
writeData(wb, "Validation", validation_table)

saveWorkbook(
  wb,
  file = output_excel,
  overwrite = TRUE
)

cat("✓ Excel outputs saved\n")

# ===============================
# STEP 10: PLOTS
# ===============================

# A. MVI Distribution
p1 <- ggplot(data, aes(x = MVI)) +
  geom_histogram(
    bins = 35,
    fill = "steelblue",
    color = "black"
  ) +
  geom_density(
    aes(y = after_stat(count)),
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Distribution of Malaria Vulnerability Index",
    x = "MVI",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.margin = margin(25,25,25,25)
  )

# B. Regional Boxplot
p2 <- ggplot(data, aes(x = region, y = MVI)) +
  geom_boxplot() +
  labs(
    title = "Regional Distribution of MVI",
    x = "Region",
    y = "MVI"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.margin = margin(25,25,40,25)
  )

# C. MVI vs Malaria
p3 <- ggplot(data, aes(x = MVI, y = log_ind_malaria_inc)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "MVI and Malaria Incidence",
    x = "MVI",
    y = "Log Malaria Incidence"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.margin = margin(25,25,25,25)
  )

# D. PCA Scree Plot
p4 <- fviz_eig(
  pca_result,
  addlabels = TRUE
)

# ===============================
# PANEL PLOT
# ===============================
tiff(
  "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/MVI_Panel.tiff",
  width = 22,
  height = 18,
  units = "in",
  res = 600,
  compression = "lzw"
)

grid.arrange(
  p1, p2, p3, p4,
  ncol = 2,
  top = textGrob(
    "Malaria Vulnerability Index Construction and Validation",
    gp = gpar(
      fontsize = 24,
      fontface = "bold"
    )
  ),
  padding = unit(3, "lines")
)

dev.off()

pdf(
  "D:/Documents/SEMESTER 3/Thesis/Final_Outputs/Outputs/Plots/MVI_Panel.pdf",
  width = 22,
  height = 18
)

grid.arrange(
  p1, p2, p3, p4,
  ncol = 2,
  top = textGrob(
    "Malaria Vulnerability Index Construction and Validation",
    gp = gpar(
      fontsize = 24,
      fontface = "bold"
    )
  ),
  padding = unit(3, "lines")
)

dev.off()

cat("✓ Journal-quality MVI plots saved\n")

# ===============================
# SAVE FINAL DATASET
# ===============================
saveRDS(
  data,
  "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_With_MVI.rds"
)

write.csv(
  data,
  "D:/Documents/SEMESTER 3/Thesis/Results/Malaria_Data_With_MVI.csv",
  row.names = FALSE
)

cat("✓ Final MVI dataset saved\n")

# ===============================
# FINAL MESSAGE
# ===============================
message("MVI construction complete.")
getwd()




 MIT License

Copyright (c) 2026 Tadius Chengetai Chitambwe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
[17:48, 13/05/2026] Taddy: # Thesis Scripts: Data Collection and Analysis Pipeline

This repository contains the complete codebase for my thesis, divided into data collection and data analysis phases.

## 💻 System Prerequisites
* *Language Environment:* R (Version 4.3 or higher recommended)
* *Required Packages:* tidyverse, rvest, httr, ggplot2

---

## 🛰️ Phase 1: Data Collection
These scripts gather the raw data used in this research project.

### File Inventory
* collect_data.R - The primary data extraction and API/web scraping script.
* api_config.R - Configuration file containing data endpoints and parameters.

### Execution Instructions
1. Open your R environment (RStudio or terminal).
2. Run the collection script using the following command in your R console:
   source("collect_data.R")
3. The script will fetch the records and save a raw data file named raw_dataset.csv into your project working directory.

---

## 📊 Phase 2: Data Analysis
These scripts process the raw data to generate the final statistics, models, and figures shown in the thesis.

### File Inventory
* clean_data.R - Handles data wrangling, removes missing values, and prepares the dataset.
* generate_plots.R - Generates the charts and data visualizations used in Chapter 4.
* statistical_tests.R - Runs the ANOVA models, regressions, and hypothesis testing.

### Execution Instructions
1. Ensure the raw_dataset.csv file from Phase 1 is present in your working directory.
2. Run the data cleaning script first to format your variables:
   source("clean_data.R")
3. Run the plotting and statistical analysis scripts to recreate the exact outputs used in the thesis:
   source("generate_plots.R")
   source("statistical_tests.R")
