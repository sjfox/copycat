## Generate weekly trajectory influenza forecasts

library(tidyverse)
library(lubridate)
library(MMWRweek)

# Clean HHS influenza data ------------------------------------------------
refresh_data <- T
flu_hhs_data_loc <- 'raw-data/flu-hhs-data.rda'
if(!file.exists(flu_hhs_data_loc)  | refresh_data){
  source('R/fetch-flu-hhs-data.R')
  flu_data <- fetch_flu()
  save(flu_data, file = flu_hhs_data_loc)
} else{
  load(flu_hhs_data_loc)
}

flu_data |> 
  mutate(year = MMWRweek(date)$MMWRyear,
         week = MMWRweek(date)$MMWRweek) |> 
  mutate(resp_season = ifelse(week>=40, year, year-1))  |> 
  group_by(location, location_name, resp_season) |> 
  arrange(date) |> 
  mutate(resp_season_week = seq_along(week)) -> recent_flu_cleaned

## Should be the latest expected date, and that date should match
## forecast_date variable value below
recent_flu_cleaned |> pull(date) |> max()


# Setup paramaeters ---------------------------------------------------
source('R/get-traj-forecast.R')
load('raw-data/locations-data.rda')
load('processed-data/clean-historic-flu-spline.rda')

## Switch to the date of the most recent HHS weekly data point
forecast_date <- ymd('2024-01-27')
curr_resp_season <- 2023
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)




# Run and save forecasts --------------------------------------------------
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
                      recent_weeks_touse = 100,
                      resp_week_range = 5,
                      forecast_horizon = 4) |>
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
  write_csv(paste0("processed-data/rt-forecasts/", forecast_date + 7,"-UGA_flucast-Copycat.csv"))


# Double check files ------------------------------------------------------
## Make sure to copy and put the csv in the right folder before this step
library(hubValidations)
hubValidations::validate_submission(hub_path = '~/projects/FluSight-forecast-hub/',
                                    file_path = 'UGA_flucast-Copycat/2024-02-03-UGA_flucast-Copycat.csv') -> sub_validation

hubValidations::validate_submission(hub_path = '~/projects/FluSight-forecast-hub/',
                  file_path = 'UGA_flucast-INFLAenza/2024-02-03-UGA_flucast-INFLAenza.csv') -> sub_validation

# temp <- read_csv("../FluSight-forecast-hub/model-output/UGA_flucast-INFLAenza/2023-11-25-UGA_flucast-INFLAenza.csv")
# temp |>
#   filter(output_type == 'pmf') |>
#   distinct(reference_date, target, horizon, target_end_date, location, output_type) |>
#   mutate(output_type_id= list(c('large_decrease', 'decrease', 'stable', 'increase', 'large_increase'))) |>
#   unnest(output_type_id) |>
#   left_join(temp) |>
#   mutate(value = ifelse(is.na(value), 0, value)) -> temp2
# 
# temp |>
#   filter(output_type == 'quantile') |>
#   bind_rows(temp2) |>
#   write_csv("../FluSight-forecast-hub/model-output/UGA_flucast-INFLAenza/2023-11-25-UGA_flucast-INFLAenza.csv")

# Want all \green checkmarks
sub_validation

## Want to make sure there are no missing required values
sub_validation$req_vals$missing


# Plot forecasts ----------------------------------------------------------
library(cowplot)
theme_set(theme_cowplot())
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
  save_plot(filename = 'figures/rt-forecast.png', base_height = 12, base_asp = 1.6, bg = 'white')


# Plot trajectories -------------------------------------------------------

traj_db |> 
  filter(resp_season_week<=30) |> 
  ggplot(aes(resp_season_week, pred, group = interaction(location,resp_season))) +
  geom_line(alpha = .1) +
  facet_wrap(~metric) +
  labs(x = 'Week of respiratory season', y = 'Growth rate') +
  background_grid(major = 'xy', minor = 'x') -> traj_plot
traj_plot  
save_plot('figures/trajectories.png', traj_plot, base_height = 4.5, base_asp = 1.8, bg = 'white')



traj_db |> 
  filter(resp_season_week<=30) |> 
  ggplot(aes(resp_season_week, pred, group = interaction(location,resp_season))) +
  geom_line(alpha = .1) +
  facet_wrap(~metric) +
  labs(x = 'Week of respiratory season', y = 'Growth rate') +
  background_grid(major = 'xy', minor = 'x') +
  geom_point(data = recent_flu_cleaned |>
               ungroup() |>
               filter(location_name == 'US',
                      resp_season == 2023,
                      date <= forecast_date) |>
               mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
               mutate(curr_weekly_change = log(lead(value)/value)) |> 
               select(resp_season_week, value, curr_weekly_change), 
             aes(x = resp_season_week, curr_weekly_change),color = 'red', inherit.aes=F) +
  geom_line(data = recent_flu_cleaned |>
               ungroup() |>
               filter(location_name == 'US',
                      resp_season == 2023,
                      date <= forecast_date) |>
               mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
               mutate(curr_weekly_change = log(lead(value)/value)) |> 
               select(resp_season_week, value, curr_weekly_change), 
             aes(x = resp_season_week, curr_weekly_change),color = 'red', inherit.aes=F) -> traj_plot2
traj_plot2
save_plot('figures/trajectories2.png', traj_plot2, base_height = 4.5, base_asp = 1.8, bg = 'white')



# Matching trajectories ---------------------------------------------------
curr_data <- recent_flu_cleaned |>
  ungroup() |>
  filter(location_name == 'US',
         resp_season == 2023,
         date <= forecast_date) |>
  mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
  mutate(curr_weekly_change = log(lead(value)/value)) |> 
  select(resp_season_week, value, curr_weekly_change)
most_recent_week <- max(curr_data$resp_season_week)
most_recent_value <- tail(curr_data$value,1)

recent_weeks_touse = 100
resp_week_range = 5
forecast_horizon = 4


cleaned_data <- curr_data |> 
  select(resp_season_week, curr_weekly_change) |> 
  filter(!is.na(curr_weekly_change)) |> 
  tail(recent_weeks_touse)  

if(resp_week_range != 0){  
  matching_data <- cleaned_data |> 
    mutate(week_change = list(c(-(1:resp_week_range), 0, (1:resp_week_range)))) |> 
    unnest(week_change) |> 
    mutate(resp_season_week = resp_season_week + week_change) |> 
    filter(resp_season_week>0)
} else {
  matching_data <- cleaned_data |> 
    mutate(week_change = 0) |> 
    filter(resp_season_week>0)
}

traj_db |> 
  inner_join(matching_data, 
             by = 'resp_season_week', relationship = "many-to-many") |> 
  group_by(week_change, metric, location, resp_season) |> 
  filter(n() == nrow(cleaned_data) | n() >= 4) |> ## Makes sure you've matched as many as the cleaned data or at least a full month
  # filter(metric == 'flusurv') |>  ## If you want to limit to specific database
  summarize(weight = sum((pred - curr_weekly_change)^2)/n(), .groups = "drop") |> 
  ungroup() |> 
  filter(!is.na(weight)) -> traj_temp

min_allowed_weight <- 0.02

traj_temp |>   
  mutate(weight = ifelse(weight < min_allowed_weight, min_allowed_weight, weight)) |>
  arrange(weight) |>
  slice(1:200) |> 
  sample_n(size = 1000, replace = T, weight = 1/weight^2) |> 
  mutate(id = seq_along(metric)) |> 
  select(id, metric, location, resp_season, week_change) -> trajectories

trajectories |> 
  left_join(traj_db |> 
              nest(data = c('resp_season_week', 'pred', 'pred_se')),
            by = c('metric', 'location', 'resp_season')) |> 
  unnest(data) |> 
  mutate(resp_season_week = resp_season_week - week_change) |> 
  filter(resp_season_week<30) |> 
  ggplot(aes(resp_season_week, pred, group = id)) +
  geom_line(alpha = .05) +
  labs(x = 'Week of respiratory season', y = 'Growth rate') +
  background_grid(major = 'xy', minor = 'x') +
  geom_point(data = recent_flu_cleaned |>
               ungroup() |>
               filter(location_name == 'US',
                      resp_season == 2023,
                      date <= forecast_date) |>
               mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
               mutate(curr_weekly_change = log(lead(value)/value)) |> 
               select(resp_season_week, value, curr_weekly_change), 
             aes(x = resp_season_week, curr_weekly_change),color = 'red', inherit.aes=F) +
  geom_line(data = recent_flu_cleaned |>
              ungroup() |>
              filter(location_name == 'US',
                     resp_season == 2023,
                     date <= forecast_date) |>
              mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
              mutate(curr_weekly_change = log(lead(value)/value)) |> 
              select(resp_season_week, value, curr_weekly_change), 
            aes(x = resp_season_week, curr_weekly_change),color = 'red', inherit.aes=F) -> traj_plot3
traj_plot3
save_plot('figures/trajectories3.png', traj_plot3, base_height = 4.5, base_asp = 1.8, bg = 'white')



  