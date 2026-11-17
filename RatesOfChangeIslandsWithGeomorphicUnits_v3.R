# ============================================================
# Wadden Sea elevation, accretion, and geomorphic-unit analysis
# ============================================================
# Required packages
library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(ggbreak)
library(patchwork)
library(knitr)

# ------------------------------------------------------------
# 1. Input paths and data loading
# ------------------------------------------------------------
dem_paths <- list(
  "2016" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2016.tif",
  "2017" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2017.tif",
  "2018" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2018.tif",
  "2019" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2019.tif",
  "2020" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2020.tif",
  "2021" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2021.tif",
  "2022" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2022.tif"
)

polygon_paths <- list(
  "Texel"           = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Texel2.shp",
  "Vlieland"        = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Vlieland2.shp",
  "Terschelling"    = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Terschelling2.shp",
  "Ameland"         = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Ameland2.shp",
  "Schiermonnikoog" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Schiermonnikoog2.shp",
  "Borkum"          = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Borkum2.shp",
  "Juist"           = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Juist2.shp",
  "Norderney"       = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Nordeney2.shp",
  "Baltrum"         = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Baltrum2.shp",
  "Langeoog"        = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Langeoog2.shp",
  "Spiekeroog"      = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Spiekeroog2.shp",
  "Wangerooge"      = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Wangeroog2.shp",
  "Amrum"           = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Amrum2.shp",
  "Sylt"            = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Sylt2.shp"
)

polygons <- map(polygon_paths, vect)
dems <- map(dem_paths, rast)
polygons <- map(polygons, ~ project(.x, crs(dems[[1]])))
years <- names(dems)

island_order <- c(
  "Texel", "Vlieland", "Terschelling", "Ameland", "Schiermonnikoog",
  "Borkum", "Juist", "Norderney", "Baltrum", "Langeoog",
  "Spiekeroog", "Wangerooge", "Amrum", "Sylt"
)
selected_islands <- island_order
island_levels <- c(island_order, "Pooled")

# ------------------------------------------------------------
# 2. Extract annual raster values by geomorphic unit
# ------------------------------------------------------------
extract_pixel_data <- function(polygons, dems, years, geo_filter) {
  map2_dfr(polygons, names(polygons), function(poly, island_name) {
    dem_stack <- rast(dems)
    vals <- terra::extract(dem_stack, poly)
    names(vals) <- c("ID", years)
    
    poly_df <- as.data.frame(poly) %>%
      mutate(ID = seq_len(n())) %>%
      select(ID, GEO)
    
    vals %>%
      left_join(poly_df, by = "ID") %>%
      mutate(
        GEO = factor(
          GEO,
          levels = c("1", "2", "3", "4"),
          labels = c("Beach", "Dune Complex", "Polder", "Back-barrier Marshes")
        ),
        island = island_name,
        pixel_id = paste(island_name, ID, sep = "_")
      ) %>%
      filter(GEO %in% geo_filter) %>%
      pivot_longer(
        cols = all_of(years),
        names_to = "year",
        values_to = "elev_m"
      ) %>%
      mutate(year_num = as.numeric(year)) %>%
      filter(!is.na(elev_m))
  })
}

beach_pixel_data <- extract_pixel_data(polygons, dems, years, "Beach")
dune_pixel_data <- extract_pixel_data(polygons, dems, years, "Dune Complex")
island_mean_pixel_data <- extract_pixel_data(
  polygons, dems, years,
  c("Beach", "Dune Complex", "Polder", "Back-barrier Marshes")
)
backbarrier_pixel_data <- extract_pixel_data(
  polygons, dems, years,
  c("Polder", "Back-barrier Marshes")
)

# ------------------------------------------------------------
# 3. OLS trends in elevation (2016–2022)
# ------------------------------------------------------------
run_island_ols <- function(pixel_data, geo_unit_label) {
  pixel_data %>%
    group_by(island) %>%
    nest() %>%
    mutate(
      n_pixels = map_int(data, ~ n_distinct(.x$pixel_id)),
      model = map(data, ~ lm(elev_m ~ year_num, data = .x)),
      tidied = map(model, tidy)
    ) %>%
    unnest(tidied) %>%
    filter(term == "year_num") %>%
    mutate(
      geo_unit = geo_unit_label,
      accretion_mm_yr = estimate * 1000,
      se_mm_yr = std.error * 1000,
      ci_lower_mm_yr = accretion_mm_yr - 1.96 * se_mm_yr,
      ci_upper_mm_yr = accretion_mm_yr + 1.96 * se_mm_yr,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

run_pooled_ols <- function(pixel_data, geo_unit_label) {
  pixel_data %>%
    lm(elev_m ~ year_num, data = .) %>%
    tidy() %>%
    filter(term == "year_num") %>%
    mutate(
      geo_unit = geo_unit_label,
      island = "Pooled",
      accretion_mm_yr = estimate * 1000,
      se_mm_yr = std.error * 1000,
      ci_lower_mm_yr = accretion_mm_yr - 1.96 * se_mm_yr,
      ci_upper_mm_yr = accretion_mm_yr + 1.96 * se_mm_yr,
      sig = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

backbarrier_label <- "Back-barrier\n(Polder + Marshes)"

beach_ols <- run_island_ols(beach_pixel_data, "Beach")
dune_ols <- run_island_ols(dune_pixel_data, "Dune Complex")
mean_ols <- run_island_ols(island_mean_pixel_data, "Island Mean")
back_ols <- run_island_ols(backbarrier_pixel_data, backbarrier_label)

beach_pool <- run_pooled_ols(beach_pixel_data, "Beach")
dune_pool <- run_pooled_ols(dune_pixel_data, "Dune Complex")
mean_pool <- run_pooled_ols(island_mean_pixel_data, "Island Mean")
back_pool <- run_pooled_ols(backbarrier_pixel_data, backbarrier_label)

plot_data <- bind_rows(
  bind_rows(beach_ols, dune_ols, mean_ols, back_ols) %>%
    filter(island %in% island_order),
  bind_rows(beach_pool, dune_pool, mean_pool, back_pool)
) %>%
  mutate(
    island = factor(island, levels = island_levels),
    is_pooled = island == "Pooled",
    geo_unit = factor(
      geo_unit,
      levels = c("Beach", "Dune Complex", "Island Mean", backbarrier_label)
    )
  )

plot_data_mean <- plot_data %>% filter(geo_unit == "Island Mean")
plot_data_units <- plot_data %>% filter(geo_unit != "Island Mean")

# Accretion-rate plot: island means are shown as background bars;
# beach, dune-complex, and back-barrier rates are shown as points.
ggplot() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 2.9, linetype = "dashed", color = "darkblue", linewidth = 0.7) +
  annotate("text", x = "Sylt", y = 2.9, label = "Current SLR (2.9 mm/yr)",
           hjust = 1, vjust = -0.5, size = 3, color = "darkblue", fontface = "italic") +
  geom_hline(yintercept = 15, linetype = "dashed", color = "darkred", linewidth = 0.7) +
  annotate("text", x = "Sylt", y = 15, label = "RCP 8.5 SLR in 2100 (15 mm/yr)",
           hjust = 1, vjust = -0.5, size = 3, color = "darkred", fontface = "italic") +
  geom_col(data = plot_data_mean, aes(x = island, y = accretion_mm_yr),
           fill = "gray85", alpha = 0.5, width = 0.7) +
  geom_errorbar(data = plot_data_mean,
                aes(x = island, y = accretion_mm_yr,
                    ymin = ci_lower_mm_yr, ymax = ci_upper_mm_yr),
                color = "gray50", width = 0.3, alpha = 0.7) +
  geom_errorbar(data = plot_data_units,
                aes(x = island, y = accretion_mm_yr,
                    ymin = ci_lower_mm_yr, ymax = ci_upper_mm_yr,
                    color = geo_unit),
                position = position_dodge(0.6), width = 0.25) +
  geom_point(data = plot_data_units,
             aes(x = island, y = accretion_mm_yr,
                 color = geo_unit, shape = geo_unit, size = is_pooled),
             position = position_dodge(0.6)) +
  scale_color_manual(values = c(
    "Beach" = "#FFCC33",
    "Dune Complex" = "#FFBB03",
    "Back-barrier\n(Polder + Marshes)" = "#66C23F"
  )) +
  scale_shape_manual(values = c(18, 16, 15)) +
  scale_size_manual(values = c(3, 4.5), guide = "none") +
  theme_minimal() +
  scale_y_break(c(60, 90), scales = 0.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top") +
  labs(title = "Accretion Rates per Island (2016–2022)",
       y = "Accretion Rate (mm/year)", x = NULL)

# ------------------------------------------------------------
# 4. Elevation change from 2016 to 2022 by geomorphic unit
# ------------------------------------------------------------
calc_polygon_change <- function(poly, dem2016, dem2022) {
  dem_stack <- c(dem2016, dem2022)
  vals <- terra::extract(dem_stack, poly)
  names(vals)[2:3] <- c("value_2016", "value_2022")
  vals$GEO <- poly$GEO[vals$ID]
  
  vals %>%
    mutate(change = value_2022 - value_2016) %>%
    group_by(GEO) %>%
    summarise(
      deposition = sum(change[change > 0], na.rm = TRUE),
      erosion = sum(change[change < 0], na.rm = TRUE),
      net_change = sum(change, na.rm = TRUE),
      mean_2016 = mean(value_2016, na.rm = TRUE),
      mean_2022 = mean(value_2022, na.rm = TRUE),
      mean_change = mean(change, na.rm = TRUE),
      .groups = "drop"
    )
}

results <- map2_dfr(polygons, names(polygons), function(poly, island_name) {
  res <- calc_polygon_change(poly, dems[["2016"]], dems[["2022"]])
  res$island <- island_name
  res
})
print(results)

# ------------------------------------------------------------
# 5. Area-weighted elevation contributions by geomorphic unit
# ------------------------------------------------------------
pixel_area <- prod(res(dems[["2016"]]))

contribution_data <- map2_dfr(polygons, names(polygons), function(poly, island_name) {
  dem_stack <- c(dems[["2016"]], dems[["2022"]])
  vals <- terra::extract(dem_stack, poly)
  value_cols <- setdiff(names(vals), "ID")
  names(vals)[match(value_cols, names(vals))] <- c("z2016", "z2022")
  
  poly_df <- as.data.frame(poly) %>%
    mutate(ID = seq_len(n())) %>%
    select(ID, GEO)
  
  vals <- vals %>%
    left_join(poly_df, by = "ID") %>%
    filter(!is.na(z2016) | !is.na(z2022)) %>%
    mutate(GEO = factor(
      GEO,
      levels = c("1", "2", "3", "4"),
      labels = c("Beach", "Dune Complex", "Polder", "Back-barrier Marshes")
    ))
  
  vals %>%
    pivot_longer(cols = c(z2016, z2022), names_to = "year", values_to = "elev") %>%
    mutate(year = sub("z", "", year)) %>%
    filter(!is.na(elev)) %>%
    group_by(year, GEO) %>%
    summarise(
      mean_elev_geo = mean(elev, na.rm = TRUE),
      n_pixels = n(),
      area_geo = n() * pixel_area,
      .groups = "drop"
    ) %>%
    group_by(year) %>%
    mutate(
      total_area = sum(area_geo),
      area_fraction = area_geo / total_area,
      contribution_m = mean_elev_geo * area_fraction,
      island = island_name
    ) %>%
    ungroup()
})

# Contributions within each island and year should sum to its mean elevation.
verification <- contribution_data %>%
  group_by(island, year) %>%
  summarise(total_contribution = sum(contribution_m), .groups = "drop")
print("Total contributions (should equal mean island elevation):")
print(verification)

geo_colors <- c(
  "Beach" = "#FFCC33",
  "Dune Complex" = "#FFBB03",
  "Polder" = "#39b54a",
  "Back-barrier Marshes" = "#99CC33"
)

# Alternative views of the 2016 and 2022 contribution data.
stacked_bar_plot <- ggplot(contribution_data,
                           aes(x = island, y = contribution_m, fill = GEO)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25, position = "stack") +
  facet_wrap(~ year, ncol = 2) +
  scale_fill_manual(values = geo_colors, name = "Geomorphic Unit") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Island Mean Elevation Decomposed by Geomorphic Unit",
       subtitle = "Bar height = area-weighted mean island elevation; stacks = unit contributions",
       x = "Island", y = "Contribution to Mean Elevation (m)", fill = "Geomorphic Unit")
print(stacked_bar_plot)

stacked_bar_plot_alternative <- ggplot(
  contribution_data, aes(x = year, y = contribution_m, fill = GEO)
) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25, position = "stack") +
  facet_wrap(~ island, ncol = 2) +
  scale_fill_manual(values = geo_colors, name = "Geomorphic Unit") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(title = "Island Mean Elevation Decomposed by Geomorphic Unit (2016 vs 2022)",
       subtitle = "Bar height = area-weighted mean island elevation; stacks = unit contributions",
       x = "Year", y = "Contribution to Mean Elevation (m)", fill = "Geomorphic Unit")
print(stacked_bar_plot_alternative)

contribution_data_combined <- contribution_data %>%
  mutate(island_year = paste(island, year, sep = "\n"))

stacked_bar_plot_single <- ggplot(
  contribution_data_combined, aes(x = island_year, y = contribution_m, fill = GEO)
) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25, position = "stack") +
  scale_fill_manual(values = geo_colors, name = "Geomorphic Unit") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Island Mean Elevation Decomposed by Geomorphic Unit",
       subtitle = "Bar height = area-weighted mean island elevation; stacks = unit contributions",
       x = "Island & Year", y = "Contribution to Mean Elevation (m)", fill = "Geomorphic Unit")
print(stacked_bar_plot_single)

summary_table <- contribution_data %>%
  select(island, year, GEO, mean_elev_geo, area_geo, contribution_m) %>%
  arrange(island, year, desc(contribution_m))
print("Summary of contributions:")
print(summary_table)

change_data <- contribution_data %>%
  select(island, GEO, year, contribution_m) %>%
  pivot_wider(names_from = year, values_from = contribution_m,
              names_prefix = "contrib_") %>%
  mutate(change_contribution = contrib_2022 - contrib_2016)

change_plot <- ggplot(change_data,
                      aes(x = island, y = change_contribution, fill = GEO)) +
  geom_col(position = "dodge", color = "grey30", linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = geo_colors, name = "Geomorphic Unit") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(title = "Change in Geomorphic Unit Contribution (2022 - 2016)",
       x = "Island", y = "Change in Contribution (m)", fill = "Geomorphic Unit")
print(change_plot)

contribution_2022 <- contribution_data %>%
  filter(year == "2022", island %in% selected_islands) %>%
  mutate(island = factor(island, levels = selected_islands))

plot_2022_contribution <- ggplot(
  contribution_2022, aes(x = island, y = contribution_m, fill = GEO)
) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25) +
  scale_fill_manual(values = geo_colors, name = "Geomorphic Unit") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Mean Island Elevation", x = "Island",
       y = "Contribution to Mean Elevation (m)")
print(plot_2022_contribution)

# ------------------------------------------------------------
# 6. Convert contributions to volume and summarize across islands
# ------------------------------------------------------------
contrib_volume <- contribution_data %>%
  mutate(
    volume_m3 = contribution_m * area_geo,
    volume_km3 = volume_m3 / 1e9
  )

joint_summary <- contrib_volume %>%
  group_by(GEO, year) %>%
  summarise(total_volume_km3 = sum(volume_km3), .groups = "drop") %>%
  pivot_wider(names_from = year, values_from = total_volume_km3,
              names_prefix = "year_") %>%
  mutate(change_volume_km3 = year_2022 - year_2016) %>%
  arrange(desc(year_2022))
print("Joint summary of geomorphic unit volumes across all islands (km³):")
print(joint_summary)

joint_plot <- joint_summary %>%
  pivot_longer(cols = c(year_2016, year_2022), names_to = "year",
               values_to = "volume_km3") %>%
  ggplot(aes(x = year, y = volume_km3, fill = GEO)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25, position = "stack") +
  scale_fill_manual(values = c(
    "Beach" = "#FFCC33", "Dune Complex" = "#FFBB03",
    "Polder" = "#39b54a", "Back-barrier Marshes" = "#99CC33"
  ), name = "Geomorphic Unit") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  labs(title = "Total Volume of Geomorphic Units Across All Islands",
       subtitle = "Comparison of 2016 and 2022 (km³)", x = "Year",
       y = "Total Volume (km³)")
print(joint_plot)

# ------------------------------------------------------------
# 7. Combined accretion-rate and 2022 elevation-contribution figure
# ------------------------------------------------------------
# Merge polder and back-barrier marsh contributions into one category.
contribution_data_merged <- contribution_data %>%
  mutate(GEO = case_when(
    GEO %in% c("Polder", "Back-barrier Marshes") ~ "Back-barrier\n(Polder + Marshes)",
    TRUE ~ as.character(GEO)
  )) %>%
  group_by(island, year, GEO) %>%
  summarise(
    contribution_m = sum(contribution_m, na.rm = TRUE),
    area_geo = sum(area_geo, na.rm = TRUE),
    n_pixels = sum(n_pixels, na.rm = TRUE),
    .groups = "drop"
  )

geo_colors_merged <- c(
  "Beach" = "#FFCC33",
  "Dune Complex" = "#FFBB03",
  "Back-barrier\n(Polder + Marshes)" = "#66C23F"
)

# The pooled contribution here is the mean of island-level contributions.
pooled_contribution <- contribution_data_merged %>%
  filter(year == "2022", island %in% selected_islands) %>%
  group_by(GEO) %>%
  summarise(contribution_m = mean(contribution_m, na.rm = TRUE), .groups = "drop") %>%
  mutate(island = "Pooled")

contribution_2022_full <- contribution_data_merged %>%
  filter(year == "2022", island %in% selected_islands) %>%
  mutate(island = as.character(island)) %>%
  bind_rows(pooled_contribution) %>%
  mutate(
    island = factor(island, levels = island_levels),
    GEO = factor(GEO, levels = c(
      "Beach", "Dune Complex", "Back-barrier\n(Polder + Marshes)"
    ))
  )

plot_data <- plot_data %>% mutate(island = factor(island, levels = island_levels))
plot_data_mean <- plot_data %>% filter(geo_unit == "Island Mean")
plot_data_units <- plot_data %>% filter(geo_unit != "Island Mean")

# Panel A: accretion rates, including reference SLR rates.
p_top <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = 3.2, linetype = "dashed", color = "darkblue", linewidth = 0.7) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "darkred", linewidth = 0.7) +
  annotate("text", x = "Pooled", y = 3.2,
           label = "3.2 mm/yr (Current SLR)", hjust = 1, vjust = -0.5,
           size = 3, color = "darkblue", fontface = "italic") +
  annotate("text", x = "Pooled", y = 10, label = "10 mm/yr",
           hjust = 1, vjust = -0.5, size = 3, color = "darkred", fontface = "italic") +
  geom_col(data = plot_data_mean, aes(x = island, y = accretion_mm_yr),
           fill = "gray85", alpha = 0.5, width = 0.7) +
  geom_errorbar(data = plot_data_mean,
                aes(x = island, y = accretion_mm_yr,
                    ymin = ci_lower_mm_yr, ymax = ci_upper_mm_yr),
                color = "gray50", width = 0.3) +
  geom_errorbar(data = plot_data_units,
                aes(x = island, y = accretion_mm_yr,
                    ymin = ci_lower_mm_yr, ymax = ci_upper_mm_yr, color = geo_unit),
                position = position_dodge(0.6), width = 0.25) +
  geom_point(data = plot_data_units,
             aes(x = island, y = accretion_mm_yr, color = geo_unit,
                 shape = geo_unit, size = is_pooled),
             position = position_dodge(0.6)) +
  scale_color_manual(
    name = "Geomorphic Unit",
    values = c("Beach" = "#FFCC33", "Dune Complex" = "#FFBB03",
               "Back-barrier\n(Polder + Marshes)" = "#66C23F"),
    breaks = c("Beach", "Dune Complex", "Back-barrier\n(Polder + Marshes)")
  ) +
  scale_shape_manual(
    name = "Geomorphic Unit",
    values = c("Beach" = 18, "Dune Complex" = 16,
               "Back-barrier\n(Polder + Marshes)" = 15),
    breaks = c("Beach", "Dune Complex", "Back-barrier\n(Polder + Marshes)")
  ) +
  scale_size_manual(values = c(3, 4.5), guide = "none") +
  scale_y_break(c(60, 90), scales = 0.5) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        legend.position = "bottom", legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  labs(title = "A) Accretion Rates per Island (2016–2022)",
       y = "Accretion Rate (mm/year)", x = NULL)

# Panel B: 2022 contribution to mean island elevation.
p_bottom <- ggplot(contribution_2022_full,
                   aes(x = island, y = contribution_m, fill = GEO)) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25, na.rm = TRUE) +
  scale_fill_manual(
    values = geo_colors_merged,
    name = "Geomorphic Unit",
    na.translate = FALSE,
    breaks = c("Beach", "Dune Complex", "Back-barrier\n(Polder + Marshes)"),
    guide = guide_legend(reverse = FALSE)
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.x = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "B) Mean Island Elevation by Geomorphic Unit (2022)",
       y = "Contribution to Mean Elevation (m)")

final_plot <- p_top / p_bottom +
  plot_layout(heights = c(1.4, 1), guides = "collect") &
  theme(legend.position = "bottom")
final_plot

# Tables corresponding to the combined figure.
table_top <- plot_data %>%
  filter(island %in% c(selected_islands, "Pooled")) %>%
  mutate(island = factor(island, levels = island_levels)) %>%
  arrange(island, geo_unit) %>%
  select(
    Island = island,
    `Geomorphic Unit` = geo_unit,
    `Accretion (mm/yr)` = accretion_mm_yr,
    `CI Lower (mm/yr)` = ci_lower_mm_yr,
    `CI Upper (mm/yr)` = ci_upper_mm_yr
  ) %>%
  mutate(
    `Uncertainty (± mm/yr)` = round((`CI Upper (mm/yr)` - `CI Lower (mm/yr)`) / 2, 2),
    across(where(is.numeric), ~ round(.x, 2))
  )

kable(table_top,
      caption = "Table A: Accretion rates per island and geomorphic unit (2016–2022)",
      align = c("l", "l", "r", "r", "r", "r"))

table_bottom <- contribution_2022_full %>%
  filter(!is.na(GEO)) %>%
  mutate(island = factor(island, levels = island_levels)) %>%
  arrange(island, GEO) %>%
  select(Island = island, `Geomorphic Unit` = GEO,
         `Contribution (m)` = contribution_m) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

kable(table_bottom,
      caption = "Table B: Contribution to mean island elevation by geomorphic unit (2022)",
      align = c("l", "l", "r"))

# ------------------------------------------------------------
# 8. Beach elevation summaries
# ------------------------------------------------------------
# Overall annual means across all beach pixels and islands.
beach_elevation_overall <- beach_pixel_data %>%
  group_by(year) %>%
  summarise(
    n_pixels = n_distinct(pixel_id),
    mean_elev_m = mean(elev_m, na.rm = TRUE),
    sd_elev_m = sd(elev_m, na.rm = TRUE),
    se_elev_m = sd_elev_m / sqrt(n_pixels),
    ci_lower_m = mean_elev_m - 1.96 * se_elev_m,
    ci_upper_m = mean_elev_m + 1.96 * se_elev_m,
    .groups = "drop"
  )
print(beach_elevation_overall)

# Annual summaries by island.
beach_elevation_by_island <- beach_pixel_data %>%
  group_by(island, year) %>%
  summarise(
    n_pixels = n_distinct(pixel_id),
    mean_elev_m = mean(elev_m, na.rm = TRUE),
    sd_elev_m = sd(elev_m, na.rm = TRUE),
    se_elev_m = sd_elev_m / sqrt(n_pixels),
    ci_lower_m = mean_elev_m - 1.96 * se_elev_m,
    ci_upper_m = mean_elev_m + 1.96 * se_elev_m,
    .groups = "drop"
  ) %>%
  mutate(island = factor(island, levels = island_levels))
print(beach_elevation_by_island)

# Single summary collapsed across years and islands.
beach_elevation_single <- beach_pixel_data %>%
  summarise(
    n_pixels = n_distinct(pixel_id),
    n_obs = n(),
    mean_elev_m = mean(elev_m, na.rm = TRUE),
    sd_elev_m = sd(elev_m, na.rm = TRUE),
    min_elev_m = min(elev_m, na.rm = TRUE),
    max_elev_m = max(elev_m, na.rm = TRUE)
  )
print(beach_elevation_single)

# Pixel-level temporal averages, to avoid treating repeated years as separate pixels.
beach_elevation_perpixel <- beach_pixel_data %>%
  group_by(island, pixel_id) %>%
  summarise(mean_elev_m = mean(elev_m, na.rm = TRUE), .groups = "drop") %>%
  group_by(island) %>%
  summarise(
    n_pixels = n(),
    mean_elev_m = mean(mean_elev_m, na.rm = TRUE),
    sd_elev_m = sd(mean_elev_m, na.rm = TRUE),
    se_elev_m = sd_elev_m / sqrt(n_pixels),
    ci_lower_m = mean_elev_m - 1.96 * se_elev_m,
    ci_upper_m = mean_elev_m + 1.96 * se_elev_m,
    .groups = "drop"
  ) %>%
  mutate(island = factor(island, levels = island_levels)) %>%
  arrange(island)
print(beach_elevation_perpixel)

# Wadden Sea-wide average of pixel-level temporal means.
beach_elevation_wadden <- beach_pixel_data %>%
  group_by(pixel_id) %>%
  summarise(mean_elev_m = mean(elev_m, na.rm = TRUE), .groups = "drop") %>%
  summarise(
    n_pixels = n(),
    mean_elev_m = mean(mean_elev_m),
    sd_elev_m = sd(mean_elev_m),
    se_elev_m = sd(mean_elev_m) / sqrt(n_pixels),
    ci_lower_m = mean_elev_m - 1.96 * se_elev_m,
    ci_upper_m = mean_elev_m + 1.96 * se_elev_m
  )
print(beach_elevation_wadden)

# ------------------------------------------------------------
# 9. Elevation summaries by island, year, and geomorphic unit
# ------------------------------------------------------------
elevation_pixels <- bind_rows(
  beach_pixel_data %>% mutate(geo_unit = "Beach"),
  dune_pixel_data %>% mutate(geo_unit = "Dune Complex"),
  backbarrier_pixel_data %>% mutate(geo_unit = "Back-barrier (Polder + Marshes)")
) %>%
  filter(!is.na(elev_m))

island_elevation_summary <- elevation_pixels %>%
  group_by(island, year, geo_unit) %>%
  summarise(
    mean_elev_m = mean(elev_m, na.rm = TRUE),
    sd_within_island_m = sd(elev_m, na.rm = TRUE),
    n_pixels = n(),
    .groups = "drop"
  ) %>%
  mutate(
    mean_of_island_means_m = NA_real_,
    sd_between_island_means_m = NA_real_,
    n_islands = NA_integer_
  )

# Pooled mean is pixel-weighted; pooled variability is reported as the
# standard deviation between island-level means.
pooled_elevation_summary <- elevation_pixels %>%
  group_by(year, geo_unit) %>%
  summarise(mean_elev_m = mean(elev_m, na.rm = TRUE), n_pixels = n(),
            .groups = "drop") %>%
  left_join(
    island_elevation_summary %>%
      group_by(year, geo_unit) %>%
      summarise(
        mean_of_island_means_m = mean(mean_elev_m, na.rm = TRUE),
        sd_between_island_means_m = sd(mean_elev_m, na.rm = TRUE),
        n_islands = n_distinct(island),
        .groups = "drop"
      ),
    by = c("year", "geo_unit")
  ) %>%
  mutate(island = "Pooled", sd_within_island_m = NA_real_) %>%
  select(island, year, geo_unit, mean_elev_m, sd_within_island_m,
         mean_of_island_means_m, sd_between_island_means_m, n_pixels, n_islands)

elevation_summary_table <- bind_rows(
  island_elevation_summary %>%
    select(island, year, geo_unit, mean_elev_m, sd_within_island_m,
           mean_of_island_means_m, sd_between_island_means_m, n_pixels, n_islands),
  pooled_elevation_summary
) %>%
  arrange(year, geo_unit, island)
print(elevation_summary_table)

# ------------------------------------------------------------
# 10. 2022 elevations by geomorphic unit
# ------------------------------------------------------------
geo_levels <- c("Beach", "Dune Complex", "Polder", "Back-barrier Marshes")

elevation_2022 <- island_mean_pixel_data %>%
  filter(year == "2022", !is.na(elev_m)) %>%
  mutate(geo_unit = factor(as.character(GEO), levels = geo_levels)) %>%
  filter(!is.na(geo_unit))

island_summary_2022 <- elevation_2022 %>%
  group_by(island, geo_unit) %>%
  summarise(
    mean_elev_m = mean(elev_m, na.rm = TRUE),
    sd_within_island_m = sd(elev_m, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    mean_of_island_means_m = NA_real_,
    sd_between_island_means_m = NA_real_,
    n_islands = NA_integer_
  )

pooled_summary_2022 <- elevation_2022 %>%
  group_by(geo_unit) %>%
  summarise(mean_elev_m = mean(elev_m, na.rm = TRUE), n_cells = n(),
            .groups = "drop") %>%
  left_join(
    island_summary_2022 %>%
      group_by(geo_unit) %>%
      summarise(
        mean_of_island_means_m = mean(mean_elev_m, na.rm = TRUE),
        sd_between_island_means_m = sd(mean_elev_m, na.rm = TRUE),
        n_islands = n_distinct(island),
        .groups = "drop"
      ),
    by = "geo_unit"
  ) %>%
  mutate(island = "Pooled", sd_within_island_m = NA_real_) %>%
  select(island, geo_unit, mean_elev_m, sd_within_island_m,
         mean_of_island_means_m, sd_between_island_means_m, n_cells, n_islands)

elevation_summary_2022 <- bind_rows(
  island_summary_2022 %>%
    select(island, geo_unit, mean_elev_m, sd_within_island_m,
           mean_of_island_means_m, sd_between_island_means_m, n_cells, n_islands),
  pooled_summary_2022
) %>%
  arrange(geo_unit, island)
print(elevation_summary_2022)

write.csv(elevation_summary_2022,
          "mean_elevation_by_island_and_unit_2022.csv", row.names = FALSE)

# Distribution plot: violins show cell values; diamonds mark pooled means;
# vertical ranges show pooled minima and maxima.
pooled_cell_stats_2022 <- elevation_2022 %>%
  group_by(geo_unit) %>%
  summarise(
    mean_elev_m = mean(elev_m, na.rm = TRUE),
    min_elev_m = min(elev_m, na.rm = TRUE),
    max_elev_m = max(elev_m, na.rm = TRUE),
    .groups = "drop"
  )

elevation_plot_2022 <- ggplot(
  elevation_2022, aes(x = geo_unit, y = elev_m, fill = geo_unit)
) +
  geom_violin(outlier.shape = NA, alpha = 0.75, width = 0.65) +
  geom_linerange(data = pooled_cell_stats_2022,
                 aes(x = geo_unit, ymin = min_elev_m, ymax = max_elev_m),
                 inherit.aes = FALSE, color = "grey25", linewidth = 0.7) +
  geom_point(data = pooled_cell_stats_2022,
             aes(x = geo_unit, y = mean_elev_m), inherit.aes = FALSE,
             shape = 23, size = 3, fill = "White") +
  scale_fill_manual(values = c(
    "Beach" = "#FFCC33", "Dune Complex" = "#FFBB03",
    "Polder" = "#39b54a", "Back-barrier Marshes" = "#99CC33"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid.minor = element_blank()) +
  labs(title = "2022 Elevation Distribution by Geomorphic Unit",
       subtitle = "Cell distributions; white diamonds show pooled means",
       x = "Geomorphic unit", y = "Elevation (m)")
print(elevation_plot_2022)

