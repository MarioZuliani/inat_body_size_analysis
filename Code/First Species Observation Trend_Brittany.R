
# iNaturalist Observer Body Size Analysis
# Within-user first-species body size trends


library(tidyverse)
library(sf)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggeffects)
library(patchwork)
library(plotly)
set.seed(123)


# Load data ---------------------------------------------------------------

birds_raw <- readRDS("Data/body_size_birds.RDS")
butterflies_raw <- readRDS("Data/body_size_butterflies.RDS")


# Settings ----------------------------------------------------------------

min_obs <- 20
research_grade_only <- TRUE
mean_line_color <- "maroon"
decline_thresholds <- c(-5, -10, -20)

output_dir <- "Outputs/iNat_First_Species_Results_Brittany"
fig_dir <- file.path(output_dir, "Figures")
table_dir <- file.path(output_dir, "Tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# Clean data function -----------------------------------------------------

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


# Prepare all observations ------------------------------------------------

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


# Overall body mass given user experience ---------------------------------

# Get count of observations per user and define observation groups
obs_levels <- c(
  "1001-5000",
  "501-1000",
  "251-500",
  "101-250",
  "50-100",
  "Less Than 50"
)

assign_obs_group <- function(dat) {
  user_groups <- dat %>%
    count(user_id, name = "count") %>%
    mutate(
      Group = case_when(
        count > 1000 ~ "1001-5000",
        count >= 501 ~ "501-1000",
        count >= 251 ~ "251-500",
        count >= 101 ~ "101-250",
        count >= 50  ~ "50-100",
        TRUE         ~ "Less Than 50"
      ),
      Group = factor(Group, levels = obs_levels)
    )
  
  left_join(dat, user_groups %>% select(user_id, Group), by = "user_id")
}

birds_clean <- assign_obs_group(birds_clean)
butterflies_clean <- assign_obs_group(butterflies_clean)

comb_clean <- bind_rows(birds_clean, butterflies_clean)

comb_clean_plot <- comb_clean %>%
  group_by(user_id, taxon) %>%
  summarise(body_size=mean(body_size),
            Group=first(Group))

comb_clean_plot %>% 
  group_by(taxon, Group) %>%
  summarise(mean_body_size=mean(body_size))

user_group_trends <- ggplot(comb_clean_plot, aes(x = Group, y = body_size)) +
  geom_violin(fill = mean_line_color) +
  geom_boxplot(width=0.1) +
  labs(x = "Number of Total Observations", y = "Body Size") +
  theme_bw(base_size=20) +
  facet_wrap(~taxon, scales="free_x") +
  theme(text = element_text(size = 20),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8)) +
  scale_y_log10() +
  coord_flip()

print(user_group_trends)



# Create first-species dataset --------------------------------------------

# Create first-species dataset (Generate the 1st time a species is observed by a user dataset)
# This is the big thing we should be setting up.

first_species_dataset <- all_obs %>%
  arrange(user_id, species_name, obs_date) %>%
  group_by(taxon, user_id, species_name) %>%
  slice_min(obs_date, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(taxon, user_id, obs_date) %>%
  group_by(taxon, user_id) %>%
  slice_head(n=100) %>%
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
  filter(n_new_species_user >= 100) %>%
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


# Raw within-user first-species slopes ------------------------------------

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


# Robust statistics for raw slopes ----------------------------------------

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


# Mixed-effects model for first-species dataset ---------------------------

fit_first_species_model <- function(df) {
  
  df <- df %>%
    filter(
      is.finite(log_body_size),
      is.finite(relative_time_first_species),
      is.finite(log_n_new_species_user)
    )
  
  lmer(
    log_body_size ~ relative_time_first_species + 
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


# Figure: Individual fitted trajectories ----------------------------------

first_species_user_lines <- raw_first_species_slopes %>%
  select(taxon, user_id, raw_slope, raw_intercept) %>%
  tidyr::crossing(relative_time_first_species = seq(0, 1, length.out = 50)) %>%
  mutate(
    fitted_log_body_size = raw_intercept + raw_slope * relative_time_first_species,
    fitted_body_size = 10^fitted_log_body_size
  )

mixed_model_line <- map_dfr(names(first_species_models), function(tx) {
  
  model <- first_species_models[[tx]]
  
  mean_log_n_new_species <- first_species_dataset %>%
    filter(taxon == tx) %>%
    summarise(mean_val = mean(log_n_new_species_user, na.rm = TRUE)) %>%
    pull(mean_val)
  
  newdata <- tibble(
    relative_time_first_species = seq(0, 1, length.out = 100),
    log_n_new_species_user = mean_log_n_new_species
  )
  
  newdata$pred_log_body_size <- predict(
    model,
    newdata = newdata,
    re.form = NA
  )
  
  newdata$pred_body_size <- 10^(newdata$pred_log_body_size)
  
  newdata$taxon <- tx
  
  newdata
})


fig_first_species_user_lines_backtransformed <- ggplot(
  first_species_user_lines %>%
    group_by(taxon) %>%
    mutate(body_size_y_limit_99 = quantile(fitted_body_size, 0.99, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(fitted_body_size <= body_size_y_limit_99),
  aes(x = relative_time_first_species*100, y = fitted_body_size, group = user_id)
) +
  geom_line(alpha = 0.035, linewidth = 0.5) +
  geom_line(
    data = mixed_model_line,
    aes(
      x = relative_time_first_species*100,
      y = pred_body_size
    ),
    linewidth = 2,
    color = mean_line_color,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ taxon, scales = "free_y") +
  theme_bw(base_size = 20) +
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



# Relative-time position of body-size decline thresholds ------------------

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



# Percent decrease in body size -------------------------------------------

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

# how much is body size expected to decrease after 100 unique observations?
percent_decrease_body_size %>%
  filter(relative_time_first_species==1) %>%
  select(relative_time_first_species, taxon, percent_change_from_start)

fig_percent_decline_by_species_number <- ggplot(
  percent_decrease_body_size,
  aes(x = relative_time_first_species*100, y = percent_change_from_start)
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


# Example plots -----------------------------------------------------------

# let's choose 3 users from birds and 3 users from butterflies with negative, neutral, and positive slopes
users_birds <- c(6067690, 524378)
users_butterflies <- c(3038372, 330735)

first_species_dataset_birds <- first_species_dataset %>%
  filter(taxon=="Birds")

first_species_dataset_butterflies <- first_species_dataset %>%
  filter(taxon=="Butterflies")

example_birds <- first_species_dataset_birds %>%
  filter(user_id %in% users_birds) %>%
  left_join(raw_first_species_slopes %>% 
              filter(taxon=="Birds") %>%
              select(user_id, raw_slope), by="user_id") %>%
  mutate(
    user_label = factor(
      user_id,
      labels = c("User A", "User B")
    )
  )

user_lines <- example_birds %>%
  group_by(user_id) %>%
  summarise(
    raw_slope = first(raw_slope),
    intercept =
      mean(log_body_size, na.rm = TRUE) -
      first(raw_slope) * mean(relative_time_first_species, na.rm = TRUE),
    max_obs = max(relative_time_first_species),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    relative_time_first_species = list(seq(0, 1, length.out = 100))
  ) %>%
  tidyr::unnest(relative_time_first_species) %>%
  mutate(
    fitted_log_body_size =
      intercept + raw_slope * relative_time_first_species,
    fitted_body_size =
      10^fitted_log_body_size
  ) %>%
  mutate(
    user_label = factor(
      user_id,
      labels = c("User A", "User B")
    )
  )

birds_example <- ggplot(
  example_birds,
  aes(
    x = relative_time_first_species*100,
    y = body_size
  )
) +
  geom_point(
    alpha = 0.4,
    size = 1
  ) +
  scale_y_log10() +
  geom_line(
    data = user_lines,
    aes(
      x = relative_time_first_species*100,
      y = fitted_body_size
    ),
    color = mean_line_color,
    linewidth = 1
  ) +
  facet_wrap(~user_label, scales = "free_x") +
  theme_bw(base_size = 15) +
  labs(
    x = "Within First-species Observation History",
    y = "Body Size",
    title = "A) Birds"
  ) + 
  theme(
    panel.grid = element_blank()
  )

print(birds_example)


example_butterflies <- first_species_dataset_butterflies %>%
  filter(user_id %in% users_butterflies) %>%
  left_join(raw_first_species_slopes %>% 
              filter(taxon=="Butterflies") %>%
              select(user_id, raw_slope), by="user_id") %>%
  mutate(
    user_label = factor(
      user_id,
      labels = c("User A", "User B")
    )
  )

user_lines <- example_butterflies %>%
  group_by(user_id) %>%
  summarise(
    raw_slope = first(raw_slope),
    intercept =
      mean(log_body_size, na.rm = TRUE) -
      first(raw_slope) * mean(relative_time_first_species, na.rm = TRUE),
    max_obs = max(relative_time_first_species),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    relative_time_first_species = list(seq(0, 1, length.out = 100))
  ) %>%
  tidyr::unnest(relative_time_first_species) %>%
  mutate(
    fitted_log_body_size =
      intercept + raw_slope * relative_time_first_species,
    fitted_body_size =
      10^fitted_log_body_size
  ) %>%
  mutate(
    user_label = factor(
      user_id,
      labels = c("User A", "User B")
    )
  )

butterfly_example <- ggplot(
  example_butterflies,
  aes(
    x = relative_time_first_species*100,
    y = body_size
  )
) +
  geom_point(
    alpha = 0.4,
    size = 1
  ) +
  geom_line(
    data = user_lines,
    aes(
      x = relative_time_first_species*100,
      y = fitted_body_size
    ),
    color = mean_line_color,
    linewidth = 1
  ) +
  scale_y_log10() +
  facet_wrap(~user_label, scales = "free_x") +
  theme_bw(base_size=15) +
  labs(
    x = "Within First-species Observation History",
    y = "Body Size",
    title = "B) Butterflies"
  ) +
  theme(
    panel.grid = element_blank()
  )

print(butterfly_example)

comb_example <- birds_example/butterfly_example
comb_example

p <- ggplot(
  example_butterflies,
  aes(
    x = relative_time_first_species,
    y = body_size,
    text = paste(
      "Species:", species_name,
      "<br>Body size:", round(body_size, 2),
      "<br>Obs #:", observation_number
    )
  )
) +
  geom_point(
    alpha = 0.4,
    size = 1
  ) +
  geom_line(
    data = user_lines,
    aes(
      x = relative_time_first_species,
      y = fitted_body_size
    ),
    color = "red",
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  scale_y_log10() +
  facet_wrap(~user_id, scales = "free_x") +
  theme_bw(base_size=15) +
  labs(
    x = "Observation Number",
    y = "Body Size"
  )

ggplotly(p, tooltip = "text")


# Model effect size in multiplicative / percent terms ---------------------

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


# Figure: First-species model predictions ---------------------------------

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

# Supplementary comparison: All observations vs first species -------------

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


# Slope figures ---------------------------------------------

# percentage of slopes that are negative
raw_first_species_slopes %>% 
  group_by(taxon) %>%
  summarise(
    negative_slopes = sum(raw_slope < 0, na.rm = TRUE),
    proportion_negative = mean(raw_slope < 0, na.rm = TRUE)*100
  )

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
    size = 6
  ) +
  facet_wrap(~ taxon, scales = "free") +
  theme_bw(base_size = 20) +
  labs(
    x = "Within-user Body Size Slope",
    y = "Density of Users"
  ) +
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 20),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 12, 10, 10),
    panel.grid = element_blank()
  )

print(fig_raw_first_species_slopes_density)



# Supplemental comparison of user group tends -----------------------------

# for this instead of selecting the first 100 unique species observations,
# we will randomly select 100 observations from each group and see how the 
# mixed effects results vary
user_groups <- comb_clean %>%
  group_by(taxon, user_id) %>%
  summarise(
    n_new_species_user = n_distinct(species_name),
    .groups = "drop"
  ) %>%
  mutate(
    Group = case_when(
      n_new_species_user >= 1000 ~ "1000+",
      n_new_species_user >= 501  ~ "501-1000",
      n_new_species_user >= 251  ~ "251-500",
      n_new_species_user >= 100  ~ "100-250",
      TRUE ~ NA_character_
    ),
    Group = factor(
      Group,
      levels = c("100-250", "251-500", "501-1000", "1000+")
    )
  )

data_with_groups <- comb_clean %>%
  select(-Group) %>%
  left_join(user_groups, by = c("taxon", "user_id")) %>%
  filter(!is.na(Group))

data_ordered <- data_with_groups %>%
  arrange(user_id, obs_date) %>%
  group_by(taxon, user_id, species_name) %>%
  slice_min(obs_date, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(taxon, user_id, obs_date) %>%
  group_by(taxon, user_id) %>%
  mutate(
    first_species_number = row_number(),
    n_new_species_user = n()
  ) %>%
  mutate(
    relative_time_first_species =
      ifelse(n_new_species_user == 1,
             0,
             (first_species_number - 1) / (n_new_species_user - 1))
  ) %>%
  ungroup()

first_species_dataset_user_comp <- data_ordered %>%
  group_by(taxon, user_id) %>%
  slice_sample(n = 100) %>%
  ungroup()

first_species_dataset_user_comp %>%
  group_by(Group, taxon) %>%
  summarise(count=n(),
            users=n_distinct(user_id))

# remove groups without enough data
first_species_dataset_user_comp <- first_species_dataset_user_comp %>%
  filter(!Group %in% c("1000+", "501-1000"))

# create grouping grid
model_grid <- expand_grid(
  taxon = c("Birds", "Butterflies"),
  Group = levels(first_species_dataset_user_comp$Group)
)

first_species_models_grouped <- pmap(
  model_grid,
  ~ {
    
    tx <- ..1
    grp <- ..2
    
    message("Fitting: ", tx, " | ", grp)
    
    df <- first_species_dataset_user_comp %>%
      filter(
        taxon == tx,
        Group == grp
      ) %>%
      filter(
        is.finite(log_body_size),
        is.finite(relative_time_first_species)
      )
    
    if (nrow(df) < 30) return(NULL)
    
    lmer(
      log_body_size ~ relative_time_first_species +
        (1 + relative_time_first_species | user_id),
      data = df,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 100000)
      )
    )
  }
)

names(first_species_models_grouped) <- paste(
  model_grid$taxon,
  model_grid$Group,
  sep = "__"
)

first_species_models_grouped <- compact(first_species_models_grouped)

first_species_group_effects <- imap_dfr(
  first_species_models_grouped,
  ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = .y)
)

relative_time_effects <- first_species_group_effects %>%
  filter(term == "relative_time_first_species")

first_species_prediction_data <- imap_dfr(
  first_species_models_grouped,
  function(model, name) {
    
    parts <- strsplit(name, "__")[[1]]
    tx <- parts[1]
    grp <- parts[2]
    
    ggpredict(
      model,
      terms = "relative_time_first_species [0:1 by=0.01]"
    ) %>%
      as.data.frame() %>%
      mutate(
        taxon = tx,
        Group = grp,
        predicted_body_size = 10^predicted,
        conf.low_body_size = 10^conf.low,
        conf.high_body_size = 10^conf.high
      )
  }
)

fig_first_species <- ggplot(
  first_species_prediction_data,
  aes(x = x*100, y = predicted_body_size)
) +
  geom_ribbon(
    aes(ymin = conf.low_body_size, ymax = conf.high_body_size),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  facet_grid(taxon ~ Group, scales = "free_y") +
  theme_bw(base_size = 15) +
  labs(
    x = "Relative time within first-species observation history",
    y = "Predicted body size"
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8))

print(fig_first_species)


# Save figures ------------------------------------------------------------

ggsave(file.path(fig_dir, "Fig_01_raw_first_species_slopes_ranked.png"),
       fig_raw_first_species_slopes_ranked, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_03_individual_first_species_trajectories_backtransformed_equal_obs_first_100_sp.png"),
       fig_first_species_user_lines_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_04_first_species_by_unique_species_log.png"),
       fig_first_species_count_lines_log, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_05_first_species_by_unique_species_backtransformed.png"),
       fig_first_species_count_lines_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_06_percent_decline_by_unique_species_100_sp.png"),
       fig_percent_decline_by_species_number, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_07_first_species_model_predictions_log.png"),
       fig_first_species_model_predictions, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_08_first_species_model_predictions_backtransformed.png"),
       fig_first_species_model_predictions_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_09_all_observations_vs_first_species_log.png"),
       fig_comparison_predictions, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_10_all_observations_vs_first_species_backtransformed.png"),
       fig_comparison_predictions_backtransformed, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_11_raw_slope_density_first_100.png"),
       fig_raw_first_species_slopes_density, width = 12, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "Fig_12_raw_slope_boxplot.png"),
       fig_raw_first_species_slopes_box, width = 10, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "User_Example.png"),
       comb_example, width = 7, height = 7, dpi = 300)

ggsave(file.path(fig_dir, "User_Group_Trends.png"),
       user_group_trends, width = 10, height = 8, dpi = 300)

ggsave(file.path(fig_dir, "User_Group_Mixed_Effects_Result.png"),
       fig_first_species, width = 7, height = 7, dpi = 300)



# Save tables -------------------------------------------------------------

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



# Sensitivity check: sample size comparision ------------------------------

# run mixed effects model

sample_sizes <- c(50, 75, 100, 200, 300, 400, 500, 600, 700)
n_runs <- 50

sample_users <- function(df, n_users) {
  
  users <- df %>%
    distinct(user_id) %>%
    slice_sample(n = n_users)
  
  df %>%
    semi_join(users, by = "user_id")
}

fit_model_from_data <- function(df) {
  
  df <- df %>%
    filter(
      is.finite(log_body_size),
      is.finite(relative_time_first_species)
    )
  
  lmer(
    log_body_size ~ relative_time_first_species +
      (1 + relative_time_first_species | user_id),
    data = df,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )
}

results <- list()

set.seed(123)

for (tx in c("Birds", "Butterflies")) {
  
  df_tx <- first_species_dataset %>%
    filter(taxon == tx)
  
  for (sample_n in sample_sizes) {
    
    for (run in 1:n_runs) {
      
      message("Fitting: ", tx,
              " | n = ", sample_n,
              " | run = ", run)
      
      df_sample <- sample_users(df_tx, sample_n)
      
      model <- tryCatch(
        fit_model_from_data(df_sample),
        error = function(e) return(NULL)
      )
      
      if (is.null(model)) next
      
      tidy_model <- broom.mixed::tidy(
        model,
        effects = "fixed",
        conf.int = TRUE
      ) %>%
        mutate(
          taxon = tx,
          sample_size = sample_n,
          run = run
        )
      
      results[[paste(tx, sample_n, run, sep = "_")]] <- tidy_model
    }
  }
}

final_results_table <- bind_rows(results) %>%
  select(
    taxon,
    sample_size,
    run,
    term,
    estimate,
    std.error,
    conf.low,
    conf.high,
    statistic,
    p.value
  )

first_species_dataset %>% group_by(taxon) %>% summarise(users=n_distinct(user_id))

# filter the data according to our sample size
final_results_table <- final_results_table %>%
  filter(
    !(taxon == "Birds" & sample_size >= 600),
    !(taxon == "Butterflies" & sample_size >= 200)
  )

sensitivity_analysis <- ggplot(final_results_table %>% filter(term == "relative_time_first_species"),
       aes(x = factor(sample_size), y = estimate)) +
  geom_boxplot(fill = mean_line_color, alpha = 0.6, outlier.alpha = 0.3) +
  facet_wrap(~ taxon, scale = "free") +
  theme_bw(base_size = 15) +
  labs(
    x = "Number of users subsampled",
    y = "Estimated slope"
  ) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8)
  )

print(sensitivity_analysis)


ggsave(
  file.path(output_dir, "Supplement_sensitivity_analysis.png"),
  sensitivity_analysis,
  width = 8,
  height = 5,
  dpi = 300
)

# save table
write_csv(final_results_table,
          file.path(supplement_table_dir, "sensitivity_analysis.csv"))

