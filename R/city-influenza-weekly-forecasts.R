library(tidyverse)
library(lubridate)
library(cowplot)
library(gam)
library(MMWRweek)
theme_set(theme_cowplot())
source('R/city-copycat-function.R')

# Read-in recent metrocast data --------------------------------------------
curr_resp_season <- 2025
## Switch to the date you are making the forecast
forecast_date <- Sys.Date() 
# forecast_date <- ymd('2025-12-24')

quantiles_needed <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
weeks_to_drop <- 0

## Pull data in from github and format properly
flu_data <- read_csv('https://raw.githubusercontent.com/reichlab/flu-metrocast/refs/heads/main/target-data/latest-data.csv')
location_info <- read_csv('https://raw.githubusercontent.com/reichlab/flu-metrocast/refs/heads/main/auxiliary-data/locations.csv')

flu_data |> 
  mutate(week = MMWRweek(target_end_date)$MMWRweek,
         year = MMWRweek(target_end_date)$MMWRyear)  |> 
  mutate(resp_season = ifelse(week>=40, year, year-1)) |> 
  filter(!resp_season %in% c(2008, 2009, 2020)) |>
  group_by(location, target, resp_season) |> 
  arrange(target_end_date) |> 
  mutate(resp_season_week = seq_along(week)) |> 
  ungroup() |> 
  left_join(location_info |> 
              select(location, state, state_abb, location_name, population), by = 'location') -> flu_data

## Separate into recent and historic
flu_data |> 
  filter(resp_season == curr_resp_season) |> 
  arrange(resp_season_week) -> recent

flu_data |> 
  filter(resp_season != curr_resp_season) |> 
  arrange(resp_season_week) -> historic


## Set the week of the forecast
forecast_week <- recent |> pull(resp_season_week) |> max()



# Create seasonal trajectory splines -------------------------------------
get_seasonal_spline_vals <- function(season_weeks, value){
  # add_before <- 
  # season_weeks
  # browser()
  padding <- 5
  new_value <- c(rep(head(value, 1),padding), value, rep(tail(value, 1),padding))
  new_season_weeks <- c(rev(min(season_weeks) - 1:padding), season_weeks, max(season_weeks)+1:padding)
  weekly_change <- lead(new_value+1)/(new_value+1)
  weekly_change <- ifelse(is.na(weekly_change), 1, weekly_change)
  df <- tibble(new_season_weeks, weekly_change) 
  
  mod <- gam(log(weekly_change) ~ s(new_season_weeks, 12), data = df)
  
  tibble(weeks=new_season_weeks,
         pred = mod$fitted.values,
         pred_se = as.numeric(predict(mod, se = T)$se.fit)) |> 
    filter(weeks %in% season_weeks)
}
library(gam)

historic |> 
  # filter(metric == 'fluview', location == 'wy', resp_season == 2017) |> 
  group_by(target, location, resp_season) |> 
  filter(n() >= 30) |> ## Makes sure every season is full
  arrange(resp_season_week) |> 
  mutate(get_seasonal_spline_vals(resp_season_week, observation)) |> 
  # filter(!any(weekly_change>10)) |> 
  ungroup() |> 
  select(target, location, resp_season, resp_season_week, pred, pred_se) -> traj_db



traj_db |> 
  ggplot(aes(resp_season_week, exp(pred), group = interaction(target, location, resp_season), color = as.factor(resp_season))) +
  geom_line(alpha = .8) 

save(traj_db, historic, file = 'processed-data/city-historic-flu-spline.rda')


# Make forecasts ----------------------------------------------------------
load('processed-data/city-historic-flu-spline.rda')
# traj_db <- traj_db |> filter(target != 'Flu ED visits pct')

locations <- recent |> 
  pull(location) |> unique()


location_forecasts <- vector('list', length = length(locations))
for(curr_location in locations){
  ## curr_location = 'NYC'
  ## Subset data to correct time period and location
  recent |>
    ungroup() |>
    filter(location == curr_location,
           resp_season_week <= forecast_week - weeks_to_drop) |> ## Removes most recent data point
    mutate(value = observation+1) |> ## Makes sure no zeroes, need to correct later
    mutate(curr_weekly_change = log(lead(value)/value)) |> 
    select(resp_season_week, value, curr_weekly_change) |>
    city_copycat(db = traj_db,
                 recent_weeks_touse = 100,
                 top_matches = 20,
                 resp_week_range = 1,
                 target_type = 'percentage',
                 nbinom_disp = .25, ## When it's not a count this is standard deviation
                 forecast_horizon = 5 + weeks_to_drop) |>
    mutate(forecast = forecast-1) |>
    mutate(forecast = ifelse(forecast < 0, 0, forecast)) -> forecast_trajectories
  
  
  forecast_trajectories |>
    group_by(resp_season_week) |>
    summarize(qs = list(value = quantile(forecast, probs = quantiles_needed))) |>
    mutate(horizon = seq_along(resp_season_week)-weeks_to_drop-2) |>
    unnest_wider(qs) |>
    gather(quantile, value, -resp_season_week, -horizon) |>
    ungroup() |>
    mutate(quantile = as.numeric(gsub("[\\%,]", "", quantile))/100) |>
    mutate(location = curr_location,
           target = ifelse(curr_location == 'nyc', paste0("ILI ED visits pct"), paste0("Flu ED visits pct")),
           reference_date = forecast_date + 3,
           target_end_date = forecast_date + 3 + horizon*7,
           output_type_id = quantile,
           output_type = 'quantile',
           value = value) %>%
    arrange(location, horizon, quantile) |>
    dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value) -> cleaned_forecasts_quantiles
  
  
  location_forecasts[[match(curr_location, locations)]] <- cleaned_forecasts_quantiles |>
    mutate(output_type_id = as.character(output_type_id)) 
}



location_forecasts |>
  bind_rows() |> 
  filter(horizon >=0) -> location_forecasts

location_forecasts |> 
  write_csv(paste0("processed-data/city-rt-forecasts/", forecast_date + 3, "-NAU-Copycat.csv"))


# Plot forecasts ----------------------------------------------------------
make_individual_plot <- function(curr_location,
                                 curr_season_data, 
                                 prev_season_data,
                                 forecast_df){
  # browser()
  curr_df <- curr_season_data |> 
    filter(location == curr_location)
  
  # curr_df |> 
  #   filter(resp_season_week == 1) |> pull(date) -> week1_date
  prev_df <- prev_season_data |> 
    filter(location == curr_location)
  
  forecast_df |> 
    filter(location == curr_location) -> forecast_df
  
  forecast_df |> 
    mutate(week = horizon + rec_resp_season_week+1) |> 
    ggplot(aes(week, `0.5`)) +
    geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), alpha = .2) +
    geom_ribbon(aes(ymin = `0.25`, ymax = `0.75`), alpha = .2) +
    geom_line() +
    geom_point(data = curr_df, aes(resp_season_week, observation)) +
    geom_point(data = prev_df, aes(resp_season_week, observation), color = 'red',alpha = .6) +
    labs(title = curr_location, x = NULL, y ='Admits') +
    background_grid(major = 'xy', minor = 'y') +
    coord_cartesian(ylim = c(0, max(c(curr_df$observation, forecast_df$`0.75`))))
}


curr_season_data <- recent |>
  ungroup() |>
  filter(resp_season_week <= forecast_week, resp_season == curr_resp_season)


rec_resp_season_week <- forecast_week

prev_season_data <- historic |>
  ungroup() |>
  filter(resp_season == curr_resp_season-1, 
         resp_season_week >= 1,
         resp_season_week < rec_resp_season_week+7)

location_forecasts |>
  filter(output_type == 'quantile') |> 
  filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |> 
  spread(output_type_id, value) -> forecast_df


unique(forecast_df$location) |> 
  map(make_individual_plot, 
      curr_season_data = curr_season_data, 
      prev_season_data = prev_season_data,
      forecast_df = forecast_df) -> plots

plot_grid(plotlist = plots) |> 
  save_plot(filename = paste0('figures/city-rt/', forecast_date + 3, '_rt-forecast.png'), base_height = 12, base_asp = 1.6, bg = 'white')


# Double check files ------------------------------------------------------
library(hubValidations)

file.copy(from=paste0("processed-data/city-rt-forecasts/", forecast_date + 3, "-NAU-Copycat.csv"), 
          to=paste0("../flu-metrocast/model-output/NAU-Copycat/", 
                    forecast_date + 3, 
                    "-NAU-Copycat.csv"), copy.mode = TRUE, overwrite = T)

hubValidations::validate_submission(hub_path = '~/projects/flu-metrocast',
                                    file_path = paste0('NAU-Copycat/', forecast_date + 3, '-NAU-Copycat.csv')) -> sub_validation

# hubValidations::validate_submission(hub_path = '~/projects/flu-metrocast',
#                                     file_path = paste0('NAU-INFLAenza/', forecast_date + 3, '-NAU-INFLAenza.csv')) -> sub_validation

# Want all \green checkmarks
sub_validation

## Want to make sure there are no missing required values
sub_validation$req_vals$missing


# 
# read_csv(paste0('../flu-metrocast/model-output/epiENGAGE-INFLAenza/', forecast_date + 3, '-epiENGAGE-INFLAenza.csv')) ->temp
# 
# location_info
# load('raw-data/locations-data.rda')
# 
# temp |> 
#   left_join(locations |> 
#               select(location, location_name)) |> 
#   mutate(location = ifelse(!is.na(location_name), tolower(location_name), location)) |> 
#   select(-location_name) |> 
#   write_csv(paste0('../flu-metrocast/model-output/epiENGAGE-INFLAenza/', forecast_date + 3, '-epiENGAGE-INFLAenza.csv'))
# 
# 
# temp |> 
#   mutate(target = ifelse(target == 'wk inc flu prop ed visits', 'Flu ED visits pct', target)) |> 
#   filter(target %in% c('Flu ED visits pct', 'ILI ED visits pct'),
#          location != 'new york') |> 
#   write_csv(paste0('../flu-metrocast/model-output/epiENGAGE-INFLAenza/', forecast_date + 3, '-epiENGAGE-INFLAenza.csv'))
# 
# temp |> 
#   filter(output_type_id %in% quantiles_needed) |> 
#   write_csv(paste0('../flu-metrocast/model-output/epiENGAGE-INFLAenza/', forecast_date + 3, '-epiENGAGE-INFLAenza.csv'))


# temp |> 
#   filter(location %in% location_info$location) |> 
#   write_csv(paste0('../flu-metrocast/model-output/epiENGAGE-INFLAenza/', forecast_date + 3, '-epiENGAGE-INFLAenza.csv'))
# 
# 
# temp |> 
#   mutate(value = value*100) |> 
#   write_csv(paste0('../flu-metrocast/model-output/epiENGAGE-INFLAenza/', forecast_date + 3, '-epiENGAGE-INFLAenza.csv'))
# 
# list.files('../flu-metrocast/model-output/NAU-Copycat', full.names = T) -> file_names
# 
# change_name_resave <- function(file_path){
#   # browser()
#   read_csv(file_path) |> 
#     write_csv(str_replace(file_path, pattern = 'epiENGAGE', replacement = 'NAU'))
# }
# file_names |> map(change_name_resave)
