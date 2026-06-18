# ============================================================
# iNaturalist Observer Body Size Analysis
# Within-user first-species body size trends
# Everything is saved in Outputs and subfolders within the Outputs
# ============================================================

library(tidyverse)
library(sf)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggeffects)
set.seed(123)

# ------------------------------------------------------------
# 1. Load data (Data from previous repo)
# ------------------------------------------------------------

birds_raw <- readRDS("Data/body_size_birds.RDS")
butterflies_raw <- readRDS("Data/body_size_butterflies.RDS")

# ------------------------------------------------------------
# 2. Settings (Some basic settings for users, set decline thresholds for figures and then output pathways)
# ------------------------------------------------------------

min_obs <- 20
research_grade_only <- TRUE
mean_line_color <- "maroon"
decline_thresholds <- c(-5, -10, -20)

output_dir <- "Outputs/iNat_First_Species_Results"
fig_dir <- file.path(output_dir, "Figures")
table_dir <- file.path(output_dir, "Tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Clean / prepare function (Function to help prepare body size data)
# ------------------------------------------------------------

prep_body_size_data <- function(df, taxon_name, min_obs = 20, research_grade_only = TRUE) {
  
  df_clean <- df %>%
    st_drop_geometry() %>%
    select(`user.id`, observed_on, body_size, quality_grade, taxon.name) %>%
    rename(user_id = `user.id`, species_name = taxon.name) %>%
    mutate(obs_date = as.Date(observed_on), taxon = taxon_name)
  
  if (research_grade_only) {
    df_clean <- df_clean %>% filter(quality_grade == "research")
  }
  
  df_clean <- df_clean %>%
    filter(
      !is.na(user_id),
      !is.na(obs_date),
      !is.na(body_size),
      !is.na(species_name),
      body_size > 0
    ) %>%
    mutate(
      user_id = as.character(user_id),
      species_name = as.character(species_name),
      body_size = as.numeric(body_size),
      log_body_size = log10(body_size)
    ) %>%
    arrange(user_id, obs_date) %>%
    group_by(user_id) %>%
    mutate(
      observation_number = row_number(),
      n_obs_user = n()
    ) %>%
    filter(n_obs_user >= min_obs) %>%
    mutate(
      log_n_obs_user = log10(n_obs_user),
      relative_time = ifelse(
        n_obs_user == 1,
        0,
        (observation_number - 1) / (n_obs_user - 1)
      ),
      taxon = factor(taxon, levels = c("Birds", "Butterflies"))
    ) %>%
    ungroup()
  
  return(df_clean)
}

# ------------------------------------------------------------
# 4. Prepare all observations (Prepare the observation data read in step 1)
# ------------------------------------------------------------

birds_clean <- prep_body_size_data(birds_raw, "Birds", min_obs, research_grade_only)
butterflies_clean <- prep_body_size_data(butterflies_raw, "Butterflies", min_obs, research_grade_only)

all_obs <- bind_rows(birds_clean, butterflies_clean)

all_obs_summary <- all_obs %>%
  group_by(taxon) %>%
  summarise(
    n_observations = n(),
    n_users = n_distinct(user_id),
    n_species = n_distinct(species_name),
    median_obs_per_user = median(n_obs_user),
    min_obs_per_user = min(n_obs_user),
    max_obs_per_user = max(n_obs_user),
    .groups = "drop"
  )

cat("\n================ ALL OBSERVATIONS SUMMARY ================\n")
print(all_obs_summary, n = Inf)

# ------------------------------------------------------------
# 5. Create first-species dataset (Generate the 1st time a species is observed by a user dataset)
# This is the big thing we should be setting up.
# ------------------------------------------------------------

first_species_dataset <- all_obs %>%
  arrange(user_id, species_name, obs_date) %>%
  group_by(taxon, user_id, species_name) %>%
  slice_min(obs_date, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(taxon, user_id, obs_date) %>%
  group_by(taxon, user_id) %>%
  mutate(
    first_species_number = row_number(),
    n_new_species_user = n(),
    log_n_new_species_user = log10(n_new_species_user),
    first_species_start_date = min(obs_date, na.rm = TRUE),
    days_since_first_species = as.numeric(obs_date - first_species_start_date),
    years_since_first_species = days_since_first_species / 365.25,
    max_years_first_species_user = max(years_since_first_species, na.rm = TRUE)) %>%
  filter(n_new_species_user >= 20) %>%
  ungroup()

first_species_summary <- first_species_dataset %>%
  group_by(taxon) %>%
  summarise(
    n_first_species_records = n(),
    n_users = n_distinct(user_id),
    n_species = n_distinct(species_name),
    median_new_species_per_user = median(n_new_species_user),
    min_new_species_per_user = min(n_new_species_user),
    max_new_species_per_user = max(n_new_species_user),
    median_years_on_app = median(max_years_first_species_user, na.rm = TRUE),
    max_years_on_app = max(max_years_first_species_user, na.rm = TRUE),
    mean_body_size = mean(body_size, na.rm = TRUE),
    median_body_size = median(body_size, na.rm = TRUE),
    mean_log_body_size = mean(log_body_size, na.rm = TRUE),
    median_log_body_size = median(log_body_size, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================ FIRST-SPECIES DATASET SUMMARY ================\n")
print(first_species_summary, n = Inf)

# ------------------------------------------------------------
# 6. Raw within-user first-species slopes
# ------------------------------------------------------------

raw_first_species_slopes <- first_species_dataset %>%
  group_by(taxon, user_id) %>%
  summarise(
    raw_slope = coef(lm(log_body_size ~ first_species_number))[2],
    raw_intercept = coef(lm(log_body_size ~ first_species_number))[1],
    n_new_species = n(),
    max_years_first_species_user = max(years_since_first_species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    predicted_start_log = raw_intercept,
    predicted_end_log = raw_intercept + raw_slope,
    predicted_start_body_size = 10^predicted_start_log,
    predicted_end_body_size = 10^predicted_end_log,
    absolute_change_body_size = predicted_end_body_size - predicted_start_body_size,
    percent_change_body_size = ((predicted_end_body_size / predicted_start_body_size) - 1) * 100
  ) %>%
  group_by(taxon) %>%
  arrange(raw_slope, .by_group = TRUE) %>%
  mutate(user_rank = row_number()) %>%
  ungroup()

raw_first_species_slope_summary <- raw_first_species_slopes %>%
  group_by(taxon) %>%
  summarise(
    n_users = n(),
    mean_raw_slope = mean(raw_slope, na.rm = TRUE),
    median_raw_slope = median(raw_slope, na.rm = TRUE),
    sd_raw_slope = sd(raw_slope, na.rm = TRUE),
    min_raw_slope = min(raw_slope, na.rm = TRUE),
    max_raw_slope = max(raw_slope, na.rm = TRUE),
    n_negative = sum(raw_slope < 0, na.rm = TRUE),
    percent_negative = mean(raw_slope < 0, na.rm = TRUE) * 100,
    median_percent_change_body_size = median(percent_change_body_size, na.rm = TRUE),
    mean_percent_change_body_size = mean(percent_change_body_size, na.rm = TRUE),
    median_start_body_size = median(predicted_start_body_size, na.rm = TRUE),
    median_end_body_size = median(predicted_end_body_size, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================ RAW FIRST-SPECIES SLOPE SUMMARY ================\n")
print(raw_first_species_slope_summary, n = Inf)

# ------------------------------------------------------------
# 7. Robust statistics for raw slopes
# ------------------------------------------------------------

wilcox_summary <- raw_first_species_slopes %>%
  group_by(taxon) %>%
  summarise(
    test = "Wilcoxon signed-rank test",
    n_users = n(),
    statistic = wilcox.test(raw_slope, mu = 0, exact = FALSE)$statistic,
    p_value = wilcox.test(raw_slope, mu = 0, exact = FALSE)$p.value,
    .groups = "drop"
  )

binomial_summary <- raw_first_species_slopes %>%
  group_by(taxon) %>%
  summarise(
    test = "Exact binomial test",
    n_users = n(),
    n_negative = sum(raw_slope < 0, na.rm = TRUE),
    proportion_negative = mean(raw_slope < 0, na.rm = TRUE),
    p_value = binom.test(
      x = sum(raw_slope < 0, na.rm = TRUE),
      n = sum(!is.na(raw_slope)),
      p = 0.5,
      alternative = "greater"
    )$p.value,
    .groups = "drop"
  )

cat("\n================ WILCOXON TESTS: RAW SLOPES VS ZERO ================\n")
print(wilcox_summary, n = Inf)

cat("\n================ BINOMIAL TESTS: PROPORTION NEGATIVE SLOPES ================\n")
print(binomial_summary, n = Inf)

# ------------------------------------------------------------
# 8. Figure: Raw user-level first-species slopes
# ------------------------------------------------------------

fig_raw_first_species_slopes_ranked <- ggplot(
  raw_first_species_slopes,
  aes(x = user_rank, y = raw_slope)
) +
  geom_point(alpha = 0.65, size = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  facet_wrap(~ taxon, scales = "free_x") +
  theme_classic(base_size = 15) +
  labs(
    x = "Users ranked by within-user slope",
    y = expression("Raw within-user slope in log"[10]*"(body size)"),
    title = "Within-user body size trends in newly detected species",
    subtitle = "Each point is one user; negative slopes indicate smaller-bodied new species over time"
  )

print(fig_raw_first_species_slopes_ranked)

# ------------------------------------------------------------
# 9. Figure: Individual fitted trajectories
# ------------------------------------------------------------

first_species_user_lines <- raw_first_species_slopes %>%
  select(taxon, user_id, raw_slope, raw_intercept, n_new_species) %>%
  group_by(taxon, user_id) %>%
  mutate(
    first_species_number = list(seq(1, n_new_species, length.out = n_new_species))
  ) %>%
  unnest(first_species_number) %>%
  mutate(
    fitted_log_body_size = raw_intercept + raw_slope * first_species_number,
    fitted_body_size = 10^fitted_log_body_size
  ) %>%
  ungroup()

mean_time_line <- first_species_user_lines %>%
  group_by(taxon, first_species_number) %>%
  summarise(
    mean_fitted_log_body_size = mean(fitted_log_body_size, na.rm = TRUE),
    mean_fitted_body_size = mean(fitted_body_size, na.rm = TRUE),
    .groups = "drop"
  )


fig_first_species_user_lines_backtransformed <- ggplot(
  first_species_user_lines,
  aes(x = first_species_number, y = fitted_body_size, group = user_id)
) +
  geom_line(alpha = 0.035, linewidth = 0.25) +
  geom_line(
    data = mean_time_line,
    aes(x = first_species_number, y = mean_fitted_body_size, group = taxon),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_bw(base_size = 15) +
  labs(
    x = "Within First-species Observation History",
    y = "Body Size",
  ) +
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.spacing = unit(1.2, "lines"),
    panel.grid = element_blank()
  )

print(fig_first_species_user_lines_backtransformed)



# ------------------------------------------------------------
# 12. Mixed-effects model for first-species dataset (Statistics)
# ------------------------------------------------------------

fit_first_species_model <- function(df) {
  
  df <- df %>%
    filter(
      is.finite(log_body_size),
      is.finite(first_species_number),
      is.finite(log_n_new_species_user)
    )
  
  lmer(
    log_body_size ~ first_species_number + log_n_new_species_user +
      (1 + first_species_number | user_id),
    data = df,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )
}

first_species_models <- list()

for (tx in c("Birds", "Butterflies")) {
  
  message("\nFitting first-species mixed model for: ", tx)
  
  first_species_models[[tx]] <- first_species_dataset %>%
    filter(taxon == tx) %>%
    fit_first_species_model()
}

first_species_singular_summary <- tibble(
  taxon = names(first_species_models),
  singular_fit = map_lgl(first_species_models, ~ isSingular(.x, tol = 1e-4))
)

cat("\n================ FIRST-SPECIES MODEL SINGULAR FIT CHECK ================\n")
print(first_species_singular_summary, n = Inf)

first_species_fixed_effects <- map_dfr(names(first_species_models), function(tx) {
  
  broom.mixed::tidy(
    first_species_models[[tx]],
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(taxon = tx)
}) %>%
  select(taxon, term, estimate, std.error, conf.low, conf.high, statistic, p.value)

cat("\n================ FIRST-SPECIES MIXED MODEL FIXED EFFECTS ================\n")
print(first_species_fixed_effects, n = Inf)

first_species_number_effects <- first_species_fixed_effects %>%
  filter(term == "first_species_number")

cat("\n================ FIRST-SPECIES MIXED MODEL RELATIVE TIME EFFECTS ================\n")
print(first_species_number_effects, n = Inf)

clean_first_species_model_summary <- first_species_number_effects %>%
  mutate(
    estimate = round(estimate, 4),
    std.error = round(std.error, 4),
    conf.low = round(conf.low, 4),
    conf.high = round(conf.high, 4),
    statistic = round(statistic, 3),
    p.value = signif(p.value, 3)
  )

cat("\n================ CLEAN FIRST-SPECIES MODEL SUMMARY ================\n")
print(clean_first_species_model_summary, n = Inf)

# ------------------------------------------------------------
# 9. Figure: Individual fitted trajectories
# ------------------------------------------------------------

first_species_user_lines <- raw_first_species_slopes %>%
  select(taxon, user_id, raw_slope, raw_intercept, n_new_species) %>%
  group_by(taxon, user_id) %>%
  mutate(
    first_species_number = list(seq(1, n_new_species, length.out = n_new_species))
  ) %>%
  unnest(first_species_number) %>%
  mutate(
    fitted_log_body_size = raw_intercept + raw_slope * first_species_number,
    fitted_body_size = 10^fitted_log_body_size
  ) %>%
  ungroup()

mixed_model_line <- map_dfr(names(first_species_models), function(tx) {
  
  model <- first_species_models[[tx]]
  
  taxon_data <- first_species_dataset %>%
    filter(taxon == tx)
  
  mean_log_n_new_species <- mean(
    taxon_data$log_n_new_species_user,
    na.rm = TRUE
  )
  
  newdata <- tibble(
    first_species_number = seq(
      min(taxon_data$first_species_number, na.rm = TRUE),
      max(taxon_data$first_species_number, na.rm = TRUE),
      length.out = 200
    ),
    log_n_new_species_user = mean_log_n_new_species
  )
  
  newdata$pred_log_body_size <- predict(
    model,
    newdata = newdata,
    re.form = NA
  )
  
  # Back-transform
  newdata$pred_body_size <- 10^(newdata$pred_log_body_size)
  
  newdata$taxon <- tx
  
  newdata
})

# crop the x-axis for better readability
first_species_user_lines_plot <- first_species_user_lines %>%
  filter(
    fitted_body_size < 500,
    (taxon == "Birds" & first_species_number <= 500) |
      (taxon == "Butterflies" & first_species_number <= 300)
  )

mixed_model_line_plot <- mixed_model_line %>%
  filter(
    (taxon == "Birds" & first_species_number <= 500) |
      (taxon == "Butterflies" & first_species_number <= 300)
  )

fig_first_species_user_lines_backtransformed <- ggplot(
  first_species_user_lines_plot,
  aes(x = first_species_number, y = fitted_body_size, group = user_id)
) +
  geom_line(alpha = 0.035, linewidth = 0.25) +
  geom_line(
    data = mixed_model_line_plot,
    aes(
      x = first_species_number,
      y = pred_body_size
    ),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free") +
  theme_bw(base_size = 15) +
  labs(
    x = "Within First-species Observation History",
    y = "Body Size",
  ) +
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.spacing = unit(1.2, "lines"),
    panel.grid = element_blank()
  )

print(fig_first_species_user_lines_backtransformed)

ggsave(file.path(fig_dir, "Fig_03_individual_first_species_trajectories_log_obs_full_data.png"),
       fig_first_species_user_lines_backtransformed, width = 12, height = 7, dpi = 300)

# ------------------------------------------------------------
# Percent decrease in body size
# ------------------------------------------------------------

percent_decrease_body_size_bird <- mixed_model_line %>%
  filter(taxon=="Birds") %>%
  mutate(
    start_body_size = first(pred_body_size),
    percent_change_from_start =
      (pred_body_size / start_body_size - 1) * 100
  )

percent_decrease_body_size_butterflies <- mixed_model_line %>%
  filter(taxon == "Butterflies") %>%
  mutate(
    start_body_size = first(pred_body_size),
    percent_change_from_start =
      (pred_body_size / start_body_size - 1) * 100
  )

percent_decrease_body_size <- rbind(percent_decrease_body_size_bird, percent_decrease_body_size_butterflies)

fig_percent_decline_by_species_number <- ggplot(
  percent_decrease_body_size,
  aes(x = first_species_number, y = percent_change_from_start)
) +
  geom_line(linewidth = 1.8, color = mean_line_color) +
  geom_hline(yintercept = decline_thresholds, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ taxon, scales = "free_x") +
  theme_bw(base_size = 20) +
  labs(
    x = "Unique species detected",
    y = "Change in Body Size (%)",
  ) +   
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 12, 10, 10),
    panel.grid = element_blank()
  )

print(fig_percent_decline_by_species_number)

ggsave(file.path(fig_dir, "Fig_06_percent_decline_by_unique_species_all_data.png"),
       fig_percent_decline_by_species_number, width = 12, height = 7, dpi = 300)


# ------------------------------------------------------------
# Slope histogram figure
# ------------------------------------------------------------

mean_slopes <- raw_first_species_slopes %>%
  group_by(taxon) %>%
  summarise(
    mean_raw_slope = mean(raw_slope, na.rm = TRUE),
    .groups = "drop"
  )

fig_raw_first_species_slopes_density <- ggplot(
  raw_first_species_slopes,
  aes(x = raw_slope)
) +
  geom_histogram(color = "white", fill = mean_line_color) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, colour = "grey35") +
  geom_text(
    data = mean_slopes,
    aes(
      x = Inf,
      y = Inf,
      label = paste0("Mean slope = ", round(mean_raw_slope, 3))
    ),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.2,
    size = 5
  ) +
  facet_wrap(~ taxon, scales = "free") +
  theme_classic(base_size = 15) +
  labs(
    x = "Within-user Body Size Slope",
    y = "Density of Users"
  ) +
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 12, 10, 10)
  )

print(fig_raw_first_species_slopes_density)

ggsave(file.path(fig_dir, "Fig_11_raw_slope_density_all_data.png"),
       fig_raw_first_species_slopes_density, width = 12, height = 7, dpi = 300)
