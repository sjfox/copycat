
library(tidyverse)
library(lubridate)
library(MMWRweek)
library(cowplot)
theme_set(theme_cowplot())

source('R/get-traj-forecast.R')
load('raw-data/locations-data.rda')
load('processed-data/clean-historic-flu-spline.rda')


# Grab NHSN influenza data ----------------------------------------------
flu_data <- read_csv('https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/target-data/target-hospital-admissions.csv')

curr_resp_season <- 2024 ## Set the current year
first_week_of_season <- '2024-10-05' ## This is used to start the season of the recent data
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

flu_data |> 
  mutate(year = MMWRweek(date)$MMWRyear,
         week = MMWRweek(date)$MMWRweek) |> 
  mutate(resp_season = ifelse(week>=40, year, year-1))  |> 
  filter(date >= first_week_of_season) |> 
  group_by(location, location_name, resp_season) |>
  arrange(date) |>
  mutate(resp_season_week = seq_along(week)) -> recent_flu_cleaned

## Switch to the date of the most recent HHS weekly data point
forecast_date <- recent_flu_cleaned |> pull(date) |> max()
forecast_date ## Should be most recent date of data - can be useful to check

# Run and save forecasts --------------------------------------------------
weeks_to_use <- 100
week_range <- 2
fcast_horizon <- 4

loc_forecasts <- vector('list', length = length(locations$location_name))
for(loc in locations$location_name){
  # i = 6; loc = 'Illinois'
  ## Subset data to correct time period and location
  recent_flu_cleaned |>
    ungroup() |>
    filter(location_name == loc,
           resp_season == curr_resp_season,
           date <= forecast_date) |>
    mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
    mutate(curr_weekly_change = log(lead(value)/value)) |> 
    select(resp_season_week, value, curr_weekly_change) |>
    get_traj_forecast(db = traj_db,
                      recent_weeks_touse = weeks_to_use,
                      resp_week_range = week_range,
                      forecast_horizon = fcast_horizon,
                      nsamps = 1000,
                      min_allowed_weight = 0.02,
                      top_matches = 500,
                      error_exponentiation = 2,
                      nbinom_disp = 100) |>
    mutate(forecast = forecast-1) |>
    mutate(forecast = ifelse(forecast < 0, 0, forecast)) -> forecast_trajectories
  
  
  forecast_trajectories |>
    group_by(resp_season_week) |>
    summarize(qs = list(value = quantile(forecast, probs = quantiles_needed))) |>
    mutate(horizon = seq_along(resp_season_week)-1) |>
    unnest_wider(qs) |>
    gather(quantile, value, -resp_season_week, -horizon) |>
    ungroup() |>
    mutate(quantile = as.numeric(gsub("[\\%,]", "", quantile))/100) |>
    mutate(location = locations$location[match(loc, locations$location_name)],
           target = paste0("wk inc flu hosp"),
           reference_date = forecast_date + 7,
           target_end_date = forecast_date + 7 + horizon*7,
           output_type_id = quantile,
           output_type = 'quantile',
           value = round(value)) %>%
    arrange(location, horizon, quantile) |>
    dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value) -> cleaned_forecasts_quantiles
  
  locations |> 
    filter(location_name == loc) |> pull(population) -> loc_population
  
  forecast_trajectories |>
    group_by(id) |> 
    arrange(resp_season_week) |> 
    mutate(horizon = seq_along(resp_season_week)-1) |> 
    ungroup() |> 
    mutate(reference_val = recent_flu_cleaned |>
             ungroup() |>
             filter(location_name == loc,
                    resp_season == curr_resp_season,
                    date == forecast_date) |> pull(value)) |>
    mutate(forecast = round(forecast)) |> 
    mutate(change = forecast - reference_val,
           pop_change = (forecast - reference_val)/loc_population*100000) |>
    mutate(classification = case_when(horizon == -1 & (abs(change) < 10 | abs(pop_change) < 1) ~ 'stable',
                                      horizon == -1 & (pop_change >= 1 & pop_change < 2) ~ 'increase',
                                      horizon == -1 & (pop_change >= 2) ~ 'large_increase',
                                      horizon == -1 & (pop_change <= -1 & pop_change > -2) ~ 'decrease',
                                      horizon == -1 & (pop_change <= -2) ~ 'large_decrease',
                                      horizon == 0 & (abs(change) < 10 | abs(pop_change) < 1) ~ 'stable',
                                      horizon == 0 & (pop_change >= 1 & pop_change < 3) ~ 'increase',
                                      horizon == 0 & (pop_change >= 3) ~ 'large_increase',
                                      horizon == 0 & (pop_change <= -1 & pop_change > -3) ~ 'decrease',
                                      horizon == 0 & (pop_change <= -3) ~ 'large_decrease',
                                      horizon == 1 & (abs(change) < 10 | abs(pop_change) < 2) ~ 'stable',
                                      horizon == 1 & (pop_change >= 2 & pop_change < 4) ~ 'increase',
                                      horizon == 1 & (pop_change >= 4) ~ 'large_increase',
                                      horizon == 1 & (pop_change <= -2 & pop_change > -4) ~ 'decrease',
                                      horizon == 1 & (pop_change <= -4) ~ 'large_decrease',
                                      horizon == 2 & (abs(change) < 10 | abs(pop_change) < 2.5) ~ 'stable',
                                      horizon == 2 & (pop_change >= 2.5 & pop_change < 5) ~ 'increase',
                                      horizon == 2 & (pop_change >= 5) ~ 'large_increase',
                                      horizon == 2 & (pop_change <= -2.5 & pop_change > -5) ~ 'decrease',
                                      horizon == 2 & (pop_change <= -5) ~ 'large_decrease',
                                      horizon == 3 & (abs(change) < 10 | abs(pop_change) < 2.5) ~ 'stable',
                                      horizon == 3 & (pop_change >= 2.5 & pop_change < 5) ~ 'increase',
                                      horizon == 3 & (pop_change >= 5) ~ 'large_increase',
                                      horizon == 3 & (pop_change <= -2.5 & pop_change > -5) ~ 'decrease',
                                      horizon == 3 & (pop_change <= -5) ~ 'large_decrease',
                                      T ~ NA
    )) |>
    count(horizon, classification) |>
    group_by(horizon) |>
    mutate(value = n/sum(n)) |>
    mutate(location = locations$location[match(loc, locations$location_name)],
           target = 'wk flu hosp rate change',
           reference_date = forecast_date + 7,
           target_end_date = forecast_date + 7 + horizon*7,
           output_type_id = classification,
           output_type = 'pmf',
           value = value) %>%
    dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value) -> cleaned_forecasts_ratechange
  cleaned_forecasts_ratechange |>
    distinct(reference_date, target, horizon, target_end_date, location, output_type) |>
    mutate(output_type_id= list(c('large_decrease', 'decrease', 'stable', 'increase', 'large_increase'))) |>
    unnest(output_type_id) |>
    left_join(cleaned_forecasts_ratechange) |>
    mutate(value = ifelse(is.na(value), 0, value)) -> cleaned_forecasts_ratechange
  
  
  loc_forecasts[[match(loc, locations$location_name)]] <- cleaned_forecasts_quantiles |>
    mutate(output_type_id = as.character(output_type_id)) |>
    bind_rows(cleaned_forecasts_ratechange)
}
loc_forecasts |>
  bind_rows() |>
  write_csv(paste0("processed-data/us-rt-forecasts/", forecast_date + 7, "-UGA_flucast-Copycat.csv"))


# Visualize the forecasts -------------------------------------------------
make_individual_plot <- function(curr_location_name,
                                 curr_season_data, 
                                 prev_season_data,
                                 forecast_df){
  # browser()
  curr_df <- curr_season_data |> 
    filter(location_name == curr_location_name)
  
  curr_df |> 
    filter(resp_season_week == 1) |> pull(date) -> week1_date
  prev_df <- prev_season_data |> 
    filter(location_name == curr_location_name) |> 
    mutate(date = week1_date + 7*(resp_season_week)-7)
  
  forecast_df |> 
    filter(location_name == curr_location_name) -> forecast_df
  
  forecast_df |> 
    ggplot(aes(target_end_date, `0.5`)) +
    geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), alpha = .2) +
    geom_ribbon(aes(ymin = `0.25`, ymax = `0.75`), alpha = .2) +
    geom_line() +
    geom_point(data = curr_df, aes(date, value)) +
    geom_point(data = prev_df, aes(date, value), color = 'red',alpha = .6) +
    scale_x_date(date_breaks = '1 month', date_labels = '%b') +
    labs(title = curr_location_name, x = NULL, y ='Admits') +
    background_grid(major = 'xy', minor = 'y') +
    coord_cartesian(ylim = c(0, max(c(curr_df$value, forecast_df$`0.75`))))
}


curr_season_data <- recent_flu_cleaned |>
  ungroup() |>
  filter(date<=forecast_date, resp_season == curr_resp_season)

# curr_season_data |> 
#   filter(location_name == 'Washington')
rec_resp_season_week <- recent_flu_cleaned |>
  ungroup() |>
  filter(date<forecast_date) |> 
  arrange(date) |> tail(1) |> pull(resp_season_week)
prev_season_data <- recent_flu_cleaned |>
  ungroup() |>
  filter(resp_season == curr_resp_season-1, 
         resp_season_week > rec_resp_season_week-10,
         resp_season_week < rec_resp_season_week+7)

loc_forecasts |>
  bind_rows() |> 
  filter(output_type == 'quantile') |> 
  filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |> 
  left_join(locations, by = 'location') |> 
  spread(output_type_id, value) -> forecast_df



locations$location_name |> 
  map(make_individual_plot, 
      curr_season_data = curr_season_data, 
      prev_season_data = prev_season_data,
      forecast_df = forecast_df) -> plots


plot_grid(plotlist = plots) |> 
  save_plot(filename = paste0('figures/us-rt/', forecast_date + 7, '_rt-forecast.png'), base_height = 14, base_asp = 1.8, bg = 'white')


# Move file to correct location -------------------------------------------------
library(hubValidations)

file.copy(from=paste0("processed-data/us-rt-forecasts/", forecast_date + 7, "-UGA_flucast-Copycat.csv"), 
          to=paste0("../FluSight-forecast-hub/model-output/UGA_flucast-Copycat/", 
                    forecast_date + 7, 
                    "-UGA_flucast-Copycat.csv"), copy.mode = TRUE, overwrite = T)


hubValidations::validate_submission(hub_path = '~/projects/FluSight-forecast-hub/',
                                    file_path = paste0('UGA_flucast-Copycat/', forecast_date + 7, '-UGA_flucast-Copycat.csv')) -> sub_validation

# Want all \green checkmarks
sub_validation

## Want to make sure there are no missing required values
sub_validation$req_vals$missing
