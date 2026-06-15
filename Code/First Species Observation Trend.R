# ============================================================
# iNaturalist Observer Body Size Analysis
# Within-user first-species body size trends
# Saves figures and tables for sharing
# ============================================================

library(tidyverse)
library(sf)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggeffects)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

birds_raw <- readRDS("Data/body_size_birds_clustered.RDS")
butterflies_raw <- readRDS("Data/body_size_butterflies_clustered.RDS")

# ------------------------------------------------------------
# 2. Settings
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
# 3. Clean / prepare function
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
# 4. Prepare all observations
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
# 5. Create first-species dataset
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
    max_years_first_species_user = max(years_since_first_species, na.rm = TRUE),
    relative_time_first_species = ifelse(
      n_new_species_user == 1,
      0,
      (first_species_number - 1) / (n_new_species_user - 1)
    )
  ) %>%
  filter(n_new_species_user >= min_obs) %>%
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
    raw_slope = coef(lm(log_body_size ~ relative_time_first_species))[2],
    raw_intercept = coef(lm(log_body_size ~ relative_time_first_species))[1],
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
  select(taxon, user_id, raw_slope, raw_intercept) %>%
  tidyr::crossing(relative_time_first_species = seq(0, 1, length.out = 50)) %>%
  mutate(
    fitted_log_body_size = raw_intercept + raw_slope * relative_time_first_species,
    fitted_body_size = 10^fitted_log_body_size
  )

mean_time_line <- first_species_user_lines %>%
  group_by(taxon, relative_time_first_species) %>%
  summarise(
    mean_fitted_log_body_size = mean(fitted_log_body_size, na.rm = TRUE),
    mean_fitted_body_size = mean(fitted_body_size, na.rm = TRUE),
    .groups = "drop"
  )

fig_first_species_user_lines <- ggplot(
  first_species_user_lines,
  aes(x = relative_time_first_species, y = fitted_log_body_size, group = user_id)
) +
  geom_line(alpha = 0.035, linewidth = 0.25) +
  geom_line(
    data = mean_time_line,
    aes(x = relative_time_first_species, y = mean_fitted_log_body_size, group = taxon),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Relative time within first-species observation history",
    y = expression("Fitted log"[10]*"(body size)"),
    title = "Individual first-species body size trajectories",
    subtitle = "Thin lines are users; maroon line is the average trajectory"
  )

print(fig_first_species_user_lines)




fig_first_species_user_lines_backtransformed <- ggplot(
  first_species_user_lines %>%
    group_by(taxon) %>%
    mutate(body_size_y_limit_99 = quantile(fitted_body_size, 0.99, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(fitted_body_size <= body_size_y_limit_99),
  aes(x = relative_time_first_species, y = fitted_body_size, group = user_id)
) +
  geom_line(alpha = 0.035, linewidth = 0.25) +
  geom_line(
    data = mean_time_line,
    aes(x = relative_time_first_species, y = mean_fitted_body_size, group = taxon),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
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
    panel.spacing = unit(1.2, "lines")
  )

print(fig_first_species_user_lines_backtransformed)










# ------------------------------------------------------------
# 10. Figure: Unique species number on x-axis
#     Plot-capped version for both x- and y-axes
# ------------------------------------------------------------

first_species_count_lines <- raw_first_species_slopes %>%
  select(taxon, user_id, raw_slope, raw_intercept, n_new_species) %>%
  group_by(taxon, user_id) %>%
  summarise(
    raw_slope = first(raw_slope),
    raw_intercept = first(raw_intercept),
    n_new_species = first(n_new_species),
    first_species_number = list(seq(1, n_new_species)),
    .groups = "drop"
  ) %>%
  unnest(first_species_number) %>%
  mutate(
    relative_time_first_species = ifelse(
      n_new_species == 1,
      0,
      (first_species_number - 1) / (n_new_species - 1)
    ),
    fitted_log_body_size = raw_intercept + raw_slope * relative_time_first_species,
    fitted_body_size = 10^fitted_log_body_size
  )

species_plot_limits <- first_species_count_lines %>%
  group_by(taxon) %>%
  summarise(
    species_x_limit_95 = as.numeric(quantile(n_new_species, 0.95, na.rm = TRUE)),
    body_size_y_limit_99 = as.numeric(quantile(fitted_body_size, 0.99, na.rm = TRUE)),
    max_n_new_species = max(n_new_species, na.rm = TRUE),
    max_fitted_body_size = max(fitted_body_size, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================ SPECIES PLOT X- AND Y-AXIS LIMITS ================\n")
print(species_plot_limits, n = Inf)

first_species_count_lines_plot <- first_species_count_lines %>%
  left_join(species_plot_limits, by = "taxon") %>%
  filter(
    first_species_number <= species_x_limit_95,
    fitted_body_size <= body_size_y_limit_99
  )

mean_count_line <- first_species_count_lines_plot %>%
  mutate(
    species_bin = case_when(
      first_species_number <= 100 ~ floor(first_species_number / 5) * 5,
      first_species_number <= 500 ~ floor(first_species_number / 25) * 25,
      TRUE ~ floor(first_species_number / 100) * 100
    )
  ) %>%
  group_by(taxon, species_bin) %>%
  summarise(
    mean_species_number = mean(first_species_number, na.rm = TRUE),
    mean_fitted_log_body_size = mean(fitted_log_body_size, na.rm = TRUE),
    mean_fitted_body_size = mean(fitted_body_size, na.rm = TRUE),
    median_fitted_body_size = median(fitted_body_size, na.rm = TRUE),
    n_users_at_species_number = n_distinct(user_id),
    .groups = "drop"
  ) %>%
  group_by(taxon) %>%
  arrange(mean_species_number, .by_group = TRUE) %>%
  mutate(
    start_body_size = first(mean_fitted_body_size),
    percent_change_from_start = ((mean_fitted_body_size / start_body_size) - 1) * 100
  ) %>%
  ungroup()

fig_first_species_count_lines_log <- ggplot(
  first_species_count_lines_plot,
  aes(x = first_species_number, y = fitted_log_body_size, group = user_id)
) +
  geom_line(alpha = 0.018, linewidth = 0.2) +
  geom_line(
    data = mean_count_line,
    aes(x = mean_species_number, y = mean_fitted_log_body_size, group = taxon),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free") +
  theme_classic(base_size = 15) +
  labs(
    x = "Number of unique species detected",
    y = expression("Fitted log"[10]*"(body size)"),
    title = "First-species trajectories by number of unique species",
    subtitle = "Thin lines are users; maroon line is the binned average trajectory; x-axis capped at the 95th percentile"
  )

print(fig_first_species_count_lines_log)

fig_first_species_count_lines_backtransformed <- ggplot(
  first_species_count_lines_plot,
  aes(x = first_species_number, y = fitted_body_size, group = user_id)
) +
  geom_line(alpha = 0.018, linewidth = 0.2) +
  geom_line(
    data = mean_count_line,
    aes(x = mean_species_number, y = mean_fitted_body_size, group = taxon),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free") +
  theme_classic(base_size = 15) +
  labs(
    x = "Number of unique species detected",
    y = "Fitted body size",
  )

print(fig_first_species_count_lines_backtransformed)

# ------------------------------------------------------------
# 11. Total size decrease by number of unique species
#     Based on capped visualization data
# ------------------------------------------------------------

species_decline_threshold_table <- map_dfr(decline_thresholds, function(thresh) {
  
  mean_count_line %>%
    filter(percent_change_from_start <= thresh) %>%
    group_by(taxon) %>%
    summarise(
      threshold_percent_decline = thresh,
      unique_species_needed = min(mean_species_number, na.rm = TRUE),
      body_size_at_threshold = mean_fitted_body_size[which.min(mean_species_number)],
      users_contributing_at_threshold = n_users_at_species_number[which.min(mean_species_number)],
      .groups = "drop"
    )
})

cat("\n================ SPECIES NUMBER NEEDED FOR BODY-SIZE DECLINE ================\n")
print(species_decline_threshold_table, n = Inf)

fig_percent_decline_by_species_number <- ggplot(
  mean_count_line,
  aes(x = mean_species_number, y = percent_change_from_start)
) +
  geom_line(linewidth = 1.8, color = mean_line_color) +
  geom_hline(yintercept = decline_thresholds, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ taxon, scales = "free_x") +
  theme_classic(base_size = 15) +
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
    plot.margin = margin(10, 12, 10, 10)
  )

print(fig_percent_decline_by_species_number)

# ------------------------------------------------------------
# 11B. Relative-time position of body-size decline thresholds
#     Table only
# ------------------------------------------------------------

species_time_curve <- first_species_dataset %>%
  group_by(taxon, user_id, first_species_number) %>%
  summarise(
    relative_time_first_species = first(relative_time_first_species),
    years_since_first_species = first(years_since_first_species),
    .groups = "drop"
  ) %>%
  group_by(taxon, first_species_number) %>%
  summarise(
    median_relative_time = median(relative_time_first_species, na.rm = TRUE),
    q25_relative_time = quantile(relative_time_first_species, 0.25, na.rm = TRUE),
    q75_relative_time = quantile(relative_time_first_species, 0.75, na.rm = TRUE),
    median_years_since_first_species = median(years_since_first_species, na.rm = TRUE),
    q25_years_since_first_species = quantile(years_since_first_species, 0.25, na.rm = TRUE),
    q75_years_since_first_species = quantile(years_since_first_species, 0.75, na.rm = TRUE),
    n_users_contributing = n_distinct(user_id),
    .groups = "drop"
  )

decline_threshold_time_table <- species_decline_threshold_table %>%
  mutate(first_species_number_join = round(unique_species_needed)) %>%
  left_join(
    species_time_curve,
    by = c("taxon", "first_species_number_join" = "first_species_number")
  ) %>%
  select(
    taxon,
    threshold_percent_decline,
    unique_species_needed,
    body_size_at_threshold,
    users_contributing_at_threshold,
    median_relative_time,
    q25_relative_time,
    q75_relative_time,
    median_years_since_first_species,
    q25_years_since_first_species,
    q75_years_since_first_species,
    n_users_contributing
  )

cat("\n================ DECLINE THRESHOLDS BY SPECIES AND RELATIVE TIME ================\n")
print(decline_threshold_time_table, n = Inf)

# ------------------------------------------------------------
# 12. Mixed-effects model for first-species dataset
# ------------------------------------------------------------

fit_first_species_model <- function(df) {
  
  df <- df %>%
    filter(
      is.finite(log_body_size),
      is.finite(relative_time_first_species),
      is.finite(log_n_new_species_user)
    )
  
  lmer(
    body_size ~ relative_time_first_species + log_n_new_species_user +
      (1 + relative_time_first_species | user_id),
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

first_species_relative_time_effects <- first_species_fixed_effects %>%
  filter(term == "relative_time_first_species")

cat("\n================ FIRST-SPECIES MIXED MODEL RELATIVE TIME EFFECTS ================\n")
print(first_species_relative_time_effects, n = Inf)

clean_first_species_model_summary <- first_species_relative_time_effects %>%
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
# 13. Model effect size in multiplicative / percent terms
# ------------------------------------------------------------

model_effect_size_summary <- first_species_relative_time_effects %>%
  mutate(
    multiplicative_change = 10^estimate,
    percent_change = (multiplicative_change - 1) * 100,
    multiplicative_change_low = 10^conf.low,
    multiplicative_change_high = 10^conf.high,
    percent_change_low = (multiplicative_change_low - 1) * 100,
    percent_change_high = (multiplicative_change_high - 1) * 100
  ) %>%
  select(
    taxon,
    estimate,
    conf.low,
    conf.high,
    p.value,
    multiplicative_change,
    percent_change,
    percent_change_low,
    percent_change_high
  )

cat("\n================ MODEL EFFECT SIZE ON ORIGINAL BODY-SIZE SCALE ================\n")
print(model_effect_size_summary, n = Inf)

# ------------------------------------------------------------
# 14. Figure: First-species model predictions
# ------------------------------------------------------------

first_species_prediction_data <- map_dfr(names(first_species_models), function(tx) {
  
  ggpredict(
    first_species_models[[tx]],
    terms = "relative_time_first_species [0:1 by=0.01]"
  ) %>%
    as.data.frame() %>%
    mutate(
      taxon = tx,
      predicted_body_size = 10^predicted,
      conf.low_body_size = 10^conf.low,
      conf.high_body_size = 10^conf.high
    )
})

fig_first_species_model_predictions <- ggplot(
  first_species_prediction_data,
  aes(x = x, y = predicted)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.20) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Relative time within first-species observation history",
    y = expression("Predicted log"[10]*"(body size)"),
    title = "Mixed-model predictions for newly detected species",
    subtitle = "Predictions show the average within-user body size trajectory"
  )

print(fig_first_species_model_predictions)

fig_first_species_model_predictions_backtransformed <- ggplot(
  first_species_prediction_data,
  aes(x = x, y = predicted_body_size)
) +
  geom_ribbon(aes(ymin = conf.low_body_size, ymax = conf.high_body_size), alpha = 0.20) +
  geom_line(linewidth = 1.5, color = mean_line_color) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Relative time within first-species observation history",
    y = "Predicted body size",
  )

print(fig_first_species_model_predictions_backtransformed)

# ------------------------------------------------------------
# 15. Supplementary comparison: All observations vs first species
# ------------------------------------------------------------

all_obs_dataset <- all_obs %>%
  mutate(dataset_compare = "All observations") %>%
  select(
    taxon,
    user_id,
    log_body_size,
    relative_time_all = relative_time,
    log_n_obs_user,
    dataset_compare
  ) %>%
  rename(
    relative_time_compare = relative_time_all,
    log_n_compare = log_n_obs_user
  )

first_species_compare_dataset <- first_species_dataset %>%
  mutate(dataset_compare = "First species") %>%
  select(
    taxon,
    user_id,
    log_body_size,
    relative_time_first_species,
    log_n_new_species_user,
    dataset_compare
  ) %>%
  rename(
    relative_time_compare = relative_time_first_species,
    log_n_compare = log_n_new_species_user
  )

comparison_data <- bind_rows(all_obs_dataset, first_species_compare_dataset) %>%
  mutate(
    dataset_compare = factor(
      dataset_compare,
      levels = c("All observations", "First species")
    )
  ) %>%
  filter(
    is.finite(log_body_size),
    is.finite(relative_time_compare),
    is.finite(log_n_compare)
  )

interaction_models <- list()

for (tx in c("Birds", "Butterflies")) {
  
  message("\nFitting comparison interaction model for: ", tx)
  
  interaction_models[[tx]] <- comparison_data %>%
    filter(taxon == tx) %>%
    lmer(
      log_body_size ~ relative_time_compare * dataset_compare + log_n_compare +
        (1 + relative_time_compare | user_id),
      data = .,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 100000)
      )
    )
}

interaction_fixed_effects <- map_dfr(names(interaction_models), function(tx) {
  
  broom.mixed::tidy(
    interaction_models[[tx]],
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(taxon = tx)
}) %>%
  select(taxon, term, estimate, std.error, conf.low, conf.high, statistic, p.value)

cat("\n================ COMPARISON MODEL FIXED EFFECTS ================\n")
print(interaction_fixed_effects, n = Inf)

comparison_key_interactions <- interaction_fixed_effects %>%
  filter(str_detect(term, "relative_time_compare:dataset_compare"))

cat("\n================ KEY INTERACTION TERMS ================\n")
print(comparison_key_interactions, n = Inf)

clean_comparison_model_summary <- comparison_key_interactions %>%
  mutate(
    estimate = round(estimate, 4),
    std.error = round(std.error, 4),
    conf.low = round(conf.low, 4),
    conf.high = round(conf.high, 4),
    statistic = round(statistic, 3),
    p.value = signif(p.value, 3)
  )

cat("\n================ CLEAN COMPARISON MODEL SUMMARY ================\n")
print(clean_comparison_model_summary, n = Inf)

comparison_prediction_data <- map_dfr(names(interaction_models), function(tx) {
  
  ggpredict(
    interaction_models[[tx]],
    terms = c("relative_time_compare [0:1 by=0.01]", "dataset_compare")
  ) %>%
    as.data.frame() %>%
    mutate(
      taxon = tx,
      predicted_body_size = 10^predicted,
      conf.low_body_size = 10^conf.low,
      conf.high_body_size = 10^conf.high
    )
})

fig_comparison_predictions <- ggplot(
  comparison_prediction_data,
  aes(x = x, y = predicted, linetype = group)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Relative time within user observation history",
    y = expression("Predicted log"[10]*"(body size)"),
    linetype = "Dataset",
    title = "All observations vs newly detected species",
    subtitle = "Comparison model tests whether first-species trends differ from all observations"
  )

print(fig_comparison_predictions)

fig_comparison_predictions_backtransformed <- ggplot(
  comparison_prediction_data,
  aes(x = x, y = predicted_body_size, linetype = group)
) +
  geom_ribbon(aes(ymin = conf.low_body_size, ymax = conf.high_body_size), alpha = 0.15) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Relative time within user observation history",
    y = "Predicted body size",
    linetype = "Dataset",
  )

print(fig_comparison_predictions_backtransformed)

# ------------------------------------------------------------
# 16. Supplementary slope figures
# ------------------------------------------------------------

fig_raw_first_species_slopes_density <- ggplot(
  raw_first_species_slopes,
  aes(x = raw_slope)
) +
  geom_density(linewidth = 1.2, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, colour = "grey35") +
  facet_wrap(~ taxon, scales = "free_y") +
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







fig_raw_first_species_slopes_box <- ggplot(
  raw_first_species_slopes,
  aes(x = taxon, y = raw_slope)
) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  theme_classic(base_size = 15) +
  labs(
    x = NULL,
    y = expression("Raw within-user slope in log"[10]*"(body size)"),
    title = "Raw first-species slopes across users"
  )

print(fig_raw_first_species_slopes_box)

# ------------------------------------------------------------
# 17. Save figures
# ------------------------------------------------------------

ggsave(file.path(fig_dir, "Fig_01_raw_first_species_slopes_ranked.png"),
       fig_raw_first_species_slopes_ranked, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_02_individual_first_species_trajectories_log.png"),
       fig_first_species_user_lines, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_03_individual_first_species_trajectories_backtransformed.png"),
       fig_first_species_user_lines_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_04_first_species_by_unique_species_log.png"),
       fig_first_species_count_lines_log, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_05_first_species_by_unique_species_backtransformed.png"),
       fig_first_species_count_lines_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_06_percent_decline_by_unique_species.png"),
       fig_percent_decline_by_species_number, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_07_first_species_model_predictions_log.png"),
       fig_first_species_model_predictions, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_08_first_species_model_predictions_backtransformed.png"),
       fig_first_species_model_predictions_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_09_all_observations_vs_first_species_log.png"),
       fig_comparison_predictions, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_10_all_observations_vs_first_species_backtransformed.png"),
       fig_comparison_predictions_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_11_raw_slope_density.png"),
       fig_raw_first_species_slopes_density, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_12_raw_slope_boxplot.png"),
       fig_raw_first_species_slopes_box, width = 10, height = 7, dpi = 300)

# ------------------------------------------------------------
# 18. Save tables
# ------------------------------------------------------------

write_csv(all_obs_summary,
          file.path(table_dir, "all_observations_summary.csv"))

write_csv(first_species_summary,
          file.path(table_dir, "first_species_dataset_summary.csv"))

write_csv(raw_first_species_slopes,
          file.path(table_dir, "raw_first_species_slopes_by_user.csv"))

write_csv(raw_first_species_slope_summary,
          file.path(table_dir, "raw_first_species_slope_summary.csv"))

write_csv(wilcox_summary,
          file.path(table_dir, "wilcoxon_raw_slope_tests.csv"))

write_csv(binomial_summary,
          file.path(table_dir, "binomial_negative_slope_tests.csv"))

write_csv(mean_count_line,
          file.path(table_dir, "mean_body_size_by_unique_species_bin.csv"))

write_csv(species_decline_threshold_table,
          file.path(table_dir, "species_decline_threshold_table.csv"))

write_csv(species_plot_limits,
          file.path(table_dir, "species_plot_x_y_axis_limits.csv"))

write_csv(first_species_count_lines_plot,
          file.path(table_dir, "first_species_count_lines_plot_capped.csv"))

write_csv(first_species_singular_summary,
          file.path(table_dir, "first_species_singular_fit_summary.csv"))

write_csv(first_species_fixed_effects,
          file.path(table_dir, "first_species_mixed_model_fixed_effects.csv"))

write_csv(clean_first_species_model_summary,
          file.path(table_dir, "clean_first_species_model_summary.csv"))

write_csv(model_effect_size_summary,
          file.path(table_dir, "model_effect_size_summary.csv"))

write_csv(interaction_fixed_effects,
          file.path(table_dir, "comparison_model_fixed_effects.csv"))

write_csv(clean_comparison_model_summary,
          file.path(table_dir, "clean_comparison_model_summary.csv"))

write_csv(species_time_curve,
          file.path(table_dir, "species_number_to_relative_time_curve.csv"))

write_csv(decline_threshold_time_table,
          file.path(table_dir, "decline_threshold_time_table.csv"))

cat("\n============================================================\n")
cat("Analysis complete.\n")
cat("Figures saved to: ", fig_dir, "\n", sep = "")
cat("Tables saved to: ", table_dir, "\n", sep = "")
cat("============================================================\n")








# ============================================================
# Sensitivity check: random subset of 500 users per taxon
# ============================================================

set.seed(123)

subset_users_per_taxon <- 500

sub_users <- first_species_dataset %>%
  distinct(taxon, user_id) %>%
  group_by(taxon) %>%
  group_modify(~ slice_sample(.x, n = min(subset_users_per_taxon, nrow(.x)))) %>%
  ungroup()

sub_first_species_dataset <- first_species_dataset %>%
  semi_join(sub_users, by = c("taxon", "user_id"))

cat("\n================ SUBSET FIRST-SPECIES DATASET SUMMARY ================\n")
sub_first_species_dataset %>%
  group_by(taxon) %>%
  summarise(
    n_records = n(),
    n_users = n_distinct(user_id),
    n_species = n_distinct(species_name),
    median_new_species_per_user = median(n_new_species_user),
    min_new_species_per_user = min(n_new_species_user),
    max_new_species_per_user = max(n_new_species_user),
    .groups = "drop"
  ) %>%
  print(n = Inf)

# ------------------------------------------------------------
# Subset raw within-user slopes on ACTUAL body-size scale
# ------------------------------------------------------------

sub_raw_first_species_slopes <- sub_first_species_dataset %>%
  group_by(taxon, user_id) %>%
  summarise(
    raw_slope = coef(lm(body_size ~ relative_time_first_species))[2],
    raw_intercept = coef(lm(body_size ~ relative_time_first_species))[1],
    n_new_species = n(),
    max_years_first_species_user = max(years_since_first_species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    predicted_start_body_size = raw_intercept,
    predicted_end_body_size = raw_intercept + raw_slope,
    absolute_change_body_size = predicted_end_body_size - predicted_start_body_size,
    percent_change_body_size = ((predicted_end_body_size / predicted_start_body_size) - 1) * 100
  ) %>%
  group_by(taxon) %>%
  arrange(raw_slope, .by_group = TRUE) %>%
  mutate(user_rank = row_number()) %>%
  ungroup()

sub_raw_first_species_slope_summary <- sub_raw_first_species_slopes %>%
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

cat("\n================ SUBSET RAW FIRST-SPECIES SLOPE SUMMARY: ACTUAL BODY SIZE ================\n")
print(sub_raw_first_species_slope_summary, n = Inf)

# ------------------------------------------------------------
# Subset robust statistics
# ------------------------------------------------------------

sub_wilcox_summary <- sub_raw_first_species_slopes %>%
  group_by(taxon) %>%
  summarise(
    test = "Wilcoxon signed-rank test",
    n_users = n(),
    statistic = wilcox.test(raw_slope, mu = 0, exact = FALSE)$statistic,
    p_value = wilcox.test(raw_slope, mu = 0, exact = FALSE)$p.value,
    .groups = "drop"
  )

sub_binomial_summary <- sub_raw_first_species_slopes %>%
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

cat("\n================ SUBSET WILCOXON TESTS: ACTUAL BODY-SIZE SLOPES VS ZERO ================\n")
print(sub_wilcox_summary, n = Inf)

cat("\n================ SUBSET BINOMIAL TESTS: PROPORTION NEGATIVE ACTUAL BODY-SIZE SLOPES ================\n")
print(sub_binomial_summary, n = Inf)

# ------------------------------------------------------------
# Subset user trajectories on actual body-size scale
# ------------------------------------------------------------

sub_first_species_user_lines <- sub_raw_first_species_slopes %>%
  select(taxon, user_id, raw_slope, raw_intercept) %>%
  tidyr::crossing(relative_time_first_species = seq(0, 1, length.out = 50)) %>%
  mutate(
    fitted_body_size = raw_intercept + raw_slope * relative_time_first_species
  ) %>%
  filter(is.finite(fitted_body_size), fitted_body_size > 0)

sub_mean_time_line <- sub_first_species_user_lines %>%
  group_by(taxon, relative_time_first_species) %>%
  summarise(
    mean_fitted_body_size = mean(fitted_body_size, na.rm = TRUE),
    median_fitted_body_size = median(fitted_body_size, na.rm = TRUE),
    .groups = "drop"
  )

sub_fig_first_species_user_lines <- ggplot(
  sub_first_species_user_lines %>%
    group_by(taxon) %>%
    mutate(body_size_y_limit_99 = quantile(fitted_body_size, 0.99, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(fitted_body_size <= body_size_y_limit_99),
  aes(x = relative_time_first_species, y = fitted_body_size, group = user_id)
) +
  geom_line(alpha = 0.035, linewidth = 0.25) +
  geom_line(
    data = sub_mean_time_line,
    aes(x = relative_time_first_species, y = mean_fitted_body_size, group = taxon),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Within First-species Observation History",
    y = "Body Size",
  ) +   theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.spacing = unit(1.2, "lines")
  )

print(sub_fig_first_species_user_lines)

# ------------------------------------------------------------
# Subset unique-species trajectories on actual body-size scale
# ------------------------------------------------------------

sub_first_species_count_lines <- sub_raw_first_species_slopes %>%
  select(taxon, user_id, raw_slope, raw_intercept, n_new_species) %>%
  group_by(taxon, user_id) %>%
  summarise(
    raw_slope = first(raw_slope),
    raw_intercept = first(raw_intercept),
    n_new_species = first(n_new_species),
    first_species_number = list(seq(1, n_new_species)),
    .groups = "drop"
  ) %>%
  unnest(first_species_number) %>%
  mutate(
    relative_time_first_species = ifelse(
      n_new_species == 1,
      0,
      (first_species_number - 1) / (n_new_species - 1)
    ),
    fitted_body_size = raw_intercept + raw_slope * relative_time_first_species
  ) %>%
  filter(is.finite(fitted_body_size), fitted_body_size > 0)

sub_species_plot_limits <- sub_first_species_count_lines %>%
  group_by(taxon) %>%
  summarise(
    species_x_limit_95 = as.numeric(quantile(n_new_species, 0.95, na.rm = TRUE)),
    body_size_y_limit_99 = as.numeric(quantile(fitted_body_size, 0.99, na.rm = TRUE)),
    max_n_new_species = max(n_new_species, na.rm = TRUE),
    max_fitted_body_size = max(fitted_body_size, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================ SUBSET SPECIES PLOT X- AND Y-AXIS LIMITS: ACTUAL BODY SIZE ================\n")
print(sub_species_plot_limits, n = Inf)

sub_first_species_count_lines_plot <- sub_first_species_count_lines %>%
  left_join(sub_species_plot_limits, by = "taxon") %>%
  filter(
    first_species_number <= species_x_limit_95,
    fitted_body_size <= body_size_y_limit_99
  )

sub_mean_count_line <- sub_first_species_count_lines_plot %>%
  mutate(
    species_bin = case_when(
      first_species_number <= 100 ~ floor(first_species_number / 5) * 5,
      first_species_number <= 500 ~ floor(first_species_number / 25) * 25,
      TRUE ~ floor(first_species_number / 100) * 100
    )
  ) %>%
  group_by(taxon, species_bin) %>%
  summarise(
    mean_species_number = mean(first_species_number, na.rm = TRUE),
    mean_fitted_body_size = mean(fitted_body_size, na.rm = TRUE),
    median_fitted_body_size = median(fitted_body_size, na.rm = TRUE),
    n_users_at_species_number = n_distinct(user_id),
    .groups = "drop"
  ) %>%
  group_by(taxon) %>%
  arrange(mean_species_number, .by_group = TRUE) %>%
  mutate(
    start_body_size = first(mean_fitted_body_size),
    percent_change_from_start = ((mean_fitted_body_size / start_body_size) - 1) * 100
  ) %>%
  ungroup()

sub_fig_first_species_count_lines <- ggplot(
  sub_first_species_count_lines_plot,
  aes(x = first_species_number, y = fitted_body_size, group = user_id)
) +
  geom_line(alpha = 0.018, linewidth = 0.2) +
  geom_line(
    data = sub_mean_count_line,
    aes(x = mean_species_number, y = mean_fitted_body_size, group = taxon),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free") +
  theme_classic(base_size = 15) +
  labs(
    x = "Number of unique species detected",
    y = "Fitted body size",
  )

print(sub_fig_first_species_count_lines)

sub_fig_percent_decline_by_species_number <- ggplot(
  sub_mean_count_line,
  aes(x = mean_species_number, y = percent_change_from_start)
) +
  geom_line(linewidth = 1.8, color = mean_line_color) +
  geom_hline(yintercept = decline_thresholds, linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ taxon, scales = "free_x") +
  theme_classic(base_size = 15) +
  labs(
    x = "Unique Species Detected",
    y = "Change in Body Size (%)",
  ) +  
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.spacing = unit(1.2, "lines")
  )

print(sub_fig_percent_decline_by_species_number)

# ------------------------------------------------------------
# Subset decline threshold table
# ------------------------------------------------------------

sub_species_decline_threshold_table <- map_dfr(decline_thresholds, function(thresh) {
  
  sub_mean_count_line %>%
    filter(percent_change_from_start <= thresh) %>%
    group_by(taxon) %>%
    summarise(
      threshold_percent_decline = thresh,
      unique_species_needed = min(mean_species_number, na.rm = TRUE),
      body_size_at_threshold = mean_fitted_body_size[which.min(mean_species_number)],
      users_contributing_at_threshold = n_users_at_species_number[which.min(mean_species_number)],
      .groups = "drop"
    )
})

cat("\n================ SUBSET SPECIES NUMBER NEEDED FOR BODY-SIZE DECLINE: ACTUAL BODY SIZE ================\n")
print(sub_species_decline_threshold_table, n = Inf)

# ------------------------------------------------------------
# Subset mixed model on actual body-size scale
# ------------------------------------------------------------

fit_subset_actual_size_model <- function(df) {
  
  df <- df %>%
    filter(
      is.finite(body_size),
      is.finite(relative_time_first_species),
      is.finite(log_n_new_species_user),
      body_size > 0
    )
  
  lmer(
    body_size ~ relative_time_first_species + log_n_new_species_user +
      (1 + relative_time_first_species | user_id),
    data = df,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )
}

sub_first_species_models <- list()

for (tx in c("Birds", "Butterflies")) {
  
  message("\nFitting subset actual-size mixed model for: ", tx)
  
  sub_first_species_models[[tx]] <- sub_first_species_dataset %>%
    filter(taxon == tx) %>%
    fit_subset_actual_size_model()
}

sub_first_species_singular_summary <- tibble(
  taxon = names(sub_first_species_models),
  singular_fit = map_lgl(sub_first_species_models, ~ isSingular(.x, tol = 1e-4))
)

cat("\n================ SUBSET ACTUAL-SIZE MODEL SINGULAR FIT CHECK ================\n")
print(sub_first_species_singular_summary, n = Inf)

sub_first_species_fixed_effects <- map_dfr(names(sub_first_species_models), function(tx) {
  
  broom.mixed::tidy(
    sub_first_species_models[[tx]],
    effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(taxon = tx)
}) %>%
  select(taxon, term, estimate, std.error, conf.low, conf.high, statistic, p.value)

cat("\n================ SUBSET ACTUAL-SIZE MIXED MODEL FIXED EFFECTS ================\n")
print(sub_first_species_fixed_effects, n = Inf)

sub_model_effect_size_summary <- map_dfr(names(sub_first_species_models), function(tx) {
  
  pred <- ggpredict(
    sub_first_species_models[[tx]],
    terms = "relative_time_first_species [0,1]"
  ) %>%
    as.data.frame()
  
  tibble(
    taxon = tx,
    predicted_start_body_size = pred$predicted[pred$x == 0],
    predicted_end_body_size = pred$predicted[pred$x == 1],
    absolute_change_body_size = predicted_end_body_size - predicted_start_body_size,
    percent_change = ((predicted_end_body_size / predicted_start_body_size) - 1) * 100
  )
})

cat("\n================ SUBSET ACTUAL-SIZE MODEL EFFECT SIZE ================\n")
print(sub_model_effect_size_summary, n = Inf)

# ------------------------------------------------------------
# Save subset sensitivity-check outputs to Supplement folder
# ------------------------------------------------------------

supplement_dir <- file.path(output_dir, "Supplement")
supplement_fig_dir <- file.path(supplement_dir, "Figures")
supplement_table_dir <- file.path(supplement_dir, "Tables")

dir.create(supplement_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplement_table_dir, recursive = TRUE, showWarnings = FALSE)

# Save supplement figures
ggsave(
  file.path(supplement_fig_dir, "Supplement_subset_individual_first_species_trajectories_actual_size.png"),
  sub_fig_first_species_user_lines,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(supplement_fig_dir, "Supplement_subset_first_species_by_unique_species_actual_size.png"),
  sub_fig_first_species_count_lines,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(supplement_fig_dir, "Supplement_subset_percent_decline_by_unique_species_actual_size.png"),
  sub_fig_percent_decline_by_species_number,
  width = 12,
  height = 7,
  dpi = 300
)

# Save supplement tables
write_csv(
  sub_raw_first_species_slope_summary,
  file.path(supplement_table_dir, "Supplement_subset_raw_first_species_slope_summary_actual_size.csv")
)

write_csv(
  sub_wilcox_summary,
  file.path(supplement_table_dir, "Supplement_subset_wilcoxon_actual_size.csv")
)

write_csv(
  sub_binomial_summary,
  file.path(supplement_table_dir, "Supplement_subset_binomial_actual_size.csv")
)

write_csv(
  sub_species_plot_limits,
  file.path(supplement_table_dir, "Supplement_subset_species_plot_x_y_axis_limits_actual_size.csv")
)

write_csv(
  sub_species_decline_threshold_table,
  file.path(supplement_table_dir, "Supplement_subset_species_decline_threshold_table_actual_size.csv")
)

write_csv(
  sub_first_species_singular_summary,
  file.path(supplement_table_dir, "Supplement_subset_actual_size_singular_fit_summary.csv")
)

write_csv(
  sub_first_species_fixed_effects,
  file.path(supplement_table_dir, "Supplement_subset_actual_size_mixed_model_fixed_effects.csv")
)

write_csv(
  sub_model_effect_size_summary,
  file.path(supplement_table_dir, "Supplement_subset_actual_size_model_effect_size_summary.csv")
)

cat("\n================ SUBSET SUPPLEMENT OUTPUTS SAVED ================\n")
cat("Supplement figures saved to: ", supplement_fig_dir, "\n", sep = "")
cat("Supplement tables saved to: ", supplement_table_dir, "\n", sep = "")
cat("============================================================\n")














# ============================================================
# Sensitivity check: slope stability across random sample sizes
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

n_reps <- 100

bird_max_users <- raw_first_species_slopes %>%
  filter(taxon == "Birds") %>%
  summarise(n = n_distinct(user_id)) %>%
  pull(n)

butterfly_max_users <- raw_first_species_slopes %>%
  filter(taxon == "Butterflies") %>%
  summarise(n = n_distinct(user_id)) %>%
  pull(n)

# Use the same subsampling sizes for both taxa so panels are directly comparable
common_sample_sizes <- c(50, 75, 100, 150, 200, 300, 400, 500, 750)
common_sample_sizes <- common_sample_sizes[
  common_sample_sizes <= min(bird_max_users, butterfly_max_users)
]

# ------------------------------------------------------------
# Helper: calculate mean-line slope from a sampled set of users
# ------------------------------------------------------------

calc_mean_line_slope <- function(slope_df, sampled_users) {
  
  sampled_lines <- slope_df %>%
    filter(user_id %in% sampled_users) %>%
    select(taxon, user_id, raw_slope, raw_intercept) %>%
    tidyr::crossing(relative_time_first_species = seq(0, 1, length.out = 50)) %>%
    mutate(
      fitted_log_body_size = raw_intercept + raw_slope * relative_time_first_species,
      fitted_body_size = 10^fitted_log_body_size
    )
  
  mean_line <- sampled_lines %>%
    group_by(taxon, relative_time_first_species) %>%
    summarise(
      mean_fitted_log_body_size = mean(fitted_log_body_size, na.rm = TRUE),
      mean_fitted_body_size = mean(fitted_body_size, na.rm = TRUE),
      .groups = "drop"
    )
  
  mean_line_slope_log <- coef(
    lm(mean_fitted_log_body_size ~ relative_time_first_species, data = mean_line)
  )[2]
  
  mean_line_slope_body <- coef(
    lm(mean_fitted_body_size ~ relative_time_first_species, data = mean_line)
  )[2]
  
  start_body_size <- mean_line %>%
    filter(relative_time_first_species == min(relative_time_first_species)) %>%
    pull(mean_fitted_body_size)
  
  end_body_size <- mean_line %>%
    filter(relative_time_first_species == max(relative_time_first_species)) %>%
    pull(mean_fitted_body_size)
  
  percent_change_body <- ((end_body_size / start_body_size) - 1) * 100
  
  tibble(
    mean_line_slope_log = as.numeric(mean_line_slope_log),
    mean_line_slope_body = as.numeric(mean_line_slope_body),
    start_body_size = as.numeric(start_body_size),
    end_body_size = as.numeric(end_body_size),
    percent_change_body = as.numeric(percent_change_body)
  )
}

# ------------------------------------------------------------
# Birds
# ------------------------------------------------------------

bird_users_all <- raw_first_species_slopes %>%
  filter(taxon == "Birds") %>%
  distinct(user_id) %>%
  pull(user_id)

bird_slope_stability <- purrr::map_dfr(common_sample_sizes, function(n_users_sampled) {
  
  purrr::map_dfr(seq_len(n_reps), function(rep_id) {
    
    sampled_users <- sample(
      bird_users_all,
      size = n_users_sampled,
      replace = FALSE
    )
    
    calc_mean_line_slope(
      slope_df = raw_first_species_slopes %>% filter(taxon == "Birds"),
      sampled_users = sampled_users
    ) %>%
      mutate(
        taxon = "Birds",
        n_users_sampled = n_users_sampled,
        replicate = rep_id
      )
  })
})

# ------------------------------------------------------------
# Butterflies
# ------------------------------------------------------------

butterfly_users_all <- raw_first_species_slopes %>%
  filter(taxon == "Butterflies") %>%
  distinct(user_id) %>%
  pull(user_id)

butterfly_slope_stability <- purrr::map_dfr(common_sample_sizes, function(n_users_sampled) {
  
  purrr::map_dfr(seq_len(n_reps), function(rep_id) {
    
    sampled_users <- sample(
      butterfly_users_all,
      size = n_users_sampled,
      replace = FALSE
    )
    
    calc_mean_line_slope(
      slope_df = raw_first_species_slopes %>% filter(taxon == "Butterflies"),
      sampled_users = sampled_users
    ) %>%
      mutate(
        taxon = "Butterflies",
        n_users_sampled = n_users_sampled,
        replicate = rep_id
      )
  })
})

slope_stability_results <- bind_rows(
  bird_slope_stability,
  butterfly_slope_stability
)

cat("\n================ SLOPE STABILITY RESULTS ================\n")
print(slope_stability_results, n = 20)

# ------------------------------------------------------------
# Full-data reference slopes
# ------------------------------------------------------------

full_data_mean_line_slope_reference <- raw_first_species_slopes %>%
  group_split(taxon) %>%
  purrr::map_dfr(function(slope_df) {
    
    tx <- unique(slope_df$taxon)
    sampled_users <- unique(slope_df$user_id)
    
    calc_mean_line_slope(
      slope_df = slope_df,
      sampled_users = sampled_users
    ) %>%
      mutate(taxon = tx)
  })

cat("\n================ FULL-DATA REFERENCE MEAN-LINE SLOPES ================\n")
print(full_data_mean_line_slope_reference, n = Inf)

# ------------------------------------------------------------
# Summary table by sample size
# ------------------------------------------------------------

slope_stability_summary <- slope_stability_results %>%
  group_by(taxon, n_users_sampled) %>%
  summarise(
    mean_mean_line_slope_log = mean(mean_line_slope_log, na.rm = TRUE),
    sd_mean_line_slope_log = sd(mean_line_slope_log, na.rm = TRUE),
    mean_mean_line_slope_body = mean(mean_line_slope_body, na.rm = TRUE),
    sd_mean_line_slope_body = sd(mean_line_slope_body, na.rm = TRUE),
    mean_percent_change_body = mean(percent_change_body, na.rm = TRUE),
    sd_percent_change_body = sd(percent_change_body, na.rm = TRUE),
    q025_slope_log = quantile(mean_line_slope_log, 0.025, na.rm = TRUE),
    q975_slope_log = quantile(mean_line_slope_log, 0.975, na.rm = TRUE),
    q025_slope_body = quantile(mean_line_slope_body, 0.025, na.rm = TRUE),
    q975_slope_body = quantile(mean_line_slope_body, 0.975, na.rm = TRUE),
    q025_percent_change = quantile(percent_change_body, 0.025, na.rm = TRUE),
    q975_percent_change = quantile(percent_change_body, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================ SLOPE STABILITY SUMMARY ================\n")
print(slope_stability_summary, n = Inf)

# ------------------------------------------------------------
# Formal stability analysis:
# Does deviation from the full-data estimate decline with sample size?
# ------------------------------------------------------------

bird_full_pct <- full_data_mean_line_slope_reference %>%
  filter(taxon == "Birds") %>%
  pull(percent_change_body)

butterfly_full_pct <- full_data_mean_line_slope_reference %>%
  filter(taxon == "Butterflies") %>%
  pull(percent_change_body)

stability_results_for_test <- slope_stability_results %>%
  mutate(
    full_data_percent_change = case_when(
      taxon == "Birds" ~ bird_full_pct,
      taxon == "Butterflies" ~ butterfly_full_pct
    ),
    abs_dev_percent_change = abs(percent_change_body - full_data_percent_change)
  )

bird_stability_lm <- lm(
  abs_dev_percent_change ~ log10(n_users_sampled),
  data = stability_results_for_test %>% filter(taxon == "Birds")
)

butterfly_stability_lm <- lm(
  abs_dev_percent_change ~ log10(n_users_sampled),
  data = stability_results_for_test %>% filter(taxon == "Butterflies")
)

cat("\n================ FORMAL STABILITY TESTS ================\n")
cat("\nBirds:\n")
print(summary(bird_stability_lm))

cat("\nButterflies:\n")
print(summary(butterfly_stability_lm))

bird_stability_coef <- as.data.frame(summary(bird_stability_lm)$coefficients)
bird_stability_coef$term <- rownames(bird_stability_coef)
rownames(bird_stability_coef) <- NULL
bird_stability_coef$taxon <- "Birds"

butterfly_stability_coef <- as.data.frame(summary(butterfly_stability_lm)$coefficients)
butterfly_stability_coef$term <- rownames(butterfly_stability_coef)
rownames(butterfly_stability_coef) <- NULL
butterfly_stability_coef$taxon <- "Butterflies"

stability_test_coefficients <- bind_rows(
  bird_stability_coef,
  butterfly_stability_coef
) %>%
  select(taxon, term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)

stability_test_model_fit <- tibble(
  taxon = c("Birds", "Butterflies"),
  response = "abs_dev_percent_change",
  r_squared = c(
    summary(bird_stability_lm)$r.squared,
    summary(butterfly_stability_lm)$r.squared
  ),
  adj_r_squared = c(
    summary(bird_stability_lm)$adj.r.squared,
    summary(butterfly_stability_lm)$adj.r.squared
  ),
  model_p_value = c(
    pf(
      summary(bird_stability_lm)$fstatistic[1],
      summary(bird_stability_lm)$fstatistic[2],
      summary(bird_stability_lm)$fstatistic[3],
      lower.tail = FALSE
    ),
    pf(
      summary(butterfly_stability_lm)$fstatistic[1],
      summary(butterfly_stability_lm)$fstatistic[2],
      summary(butterfly_stability_lm)$fstatistic[3],
      lower.tail = FALSE
    )
  )
)

cat("\n================ STABILITY TEST COEFFICIENTS ================\n")
print(stability_test_coefficients)

cat("\n================ STABILITY TEST MODEL FIT ================\n")
print(stability_test_model_fit)

# ------------------------------------------------------------
# Prep discrete x-axis positions so actual sample sizes are shown
# but spacing is cleaner and less clustered
# ------------------------------------------------------------

plot_results_discrete <- slope_stability_results %>%
  mutate(
    n_users_sampled_f = factor(n_users_sampled, levels = common_sample_sizes)
  )

plot_summary_discrete <- slope_stability_summary %>%
  mutate(
    n_users_sampled_f = factor(n_users_sampled, levels = common_sample_sizes)
  )

# ------------------------------------------------------------
# Figure: slope stability
# ------------------------------------------------------------

fig_slope_stability_body <- ggplot(
  plot_results_discrete,
  aes(x = n_users_sampled_f, y = mean_line_slope_body)
) +
  geom_point(
    alpha = 0.30,
    size = 1.6,
    position = position_jitter(width = 0.12, height = 0)
  ) +
  geom_line(
    data = plot_summary_discrete,
    aes(x = n_users_sampled_f, y = mean_mean_line_slope_body, group = taxon),
    linewidth = 1.2,
    colour = mean_line_color
  ) +
  geom_hline(
    data = full_data_mean_line_slope_reference,
    aes(yintercept = mean_line_slope_body),
    linetype = "dashed",
    linewidth = 0.9,
    colour = "grey35"
  ) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Number of Users Subsampled",
    y = "Mean Body Size Slope"
  ) +
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 12, 10, 10)
  )

print(fig_slope_stability_body)

# ------------------------------------------------------------
# Figure: percent change stability
# ------------------------------------------------------------

full_data_percent_change_reference <- full_data_mean_line_slope_reference %>%
  select(taxon, percent_change_body)

fig_percent_change_stability <- ggplot(
  plot_results_discrete,
  aes(x = n_users_sampled_f, y = percent_change_body)
) +
  geom_point(
    alpha = 0.30,
    size = 1.6,
    position = position_jitter(width = 0.12, height = 0)
  ) +
  geom_line(
    data = plot_summary_discrete,
    aes(x = n_users_sampled_f, y = mean_percent_change_body, group = taxon),
    linewidth = 1.2,
    colour = mean_line_color
  ) +
  geom_hline(
    data = full_data_percent_change_reference,
    aes(yintercept = percent_change_body),
    linetype = "dashed",
    linewidth = 0.9,
    colour = "grey35"
  ) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_classic(base_size = 15) +
  labs(
    x = "Number of Users Subsampled",
    y = "Start-to-End Percent Change in Mean Body Size"
  ) +
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 12, 10, 10)
  )

print(fig_percent_change_stability)

# ------------------------------------------------------------
# Save slope-stability sensitivity outputs
# ------------------------------------------------------------

ggsave(
  file.path(supplement_fig_dir, "Supplement_slope_stability_body_size.png"),
  fig_slope_stability_body,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(supplement_fig_dir, "Supplement_percent_change_stability.png"),
  fig_percent_change_stability,
  width = 12,
  height = 7,
  dpi = 300
)

write_csv(
  slope_stability_results,
  file.path(supplement_table_dir, "Supplement_slope_stability_results.csv")
)

write_csv(
  slope_stability_summary,
  file.path(supplement_table_dir, "Supplement_slope_stability_summary.csv")
)

write_csv(
  full_data_mean_line_slope_reference,
  file.path(supplement_table_dir, "Supplement_full_data_mean_line_slope_reference.csv")
)

write_csv(
  stability_test_coefficients,
  file.path(supplement_table_dir, "Supplement_stability_test_coefficients.csv")
)

write_csv(
  stability_test_model_fit,
  file.path(supplement_table_dir, "Supplement_stability_test_model_fit.csv")
)

cat("\n================ SAMPLE-SIZE STABILITY OUTPUTS SAVED ================\n")
cat("Supplement figures saved to: ", supplement_fig_dir, "\n", sep = "")
cat("Supplement tables saved to: ", supplement_table_dir, "\n", sep = "")
cat("============================================================\n")