## Retrospective forecasts for 2022 to test model
library(tidyverse)
library(lubridate)
library(MMWRweek)

# Clean HHS influenza data ------------------------------------------------
refresh_data = F
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


# Setup paramaeters ---------------------------------------------------
load('raw-data/locations-data.rda')
load('processed-data/clean-historic-flu-spline.rda')

forecast_weeks <- ymd('2022-10-17')+days((0:24)*7)
curr_resp_season <- 2022
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

source('R/get-traj-forecast.R')


# Run and save forecasts --------------------------------------------------
for(i in 1:length(forecast_weeks)){
  loc_forecasts <- vector('list', length = length(locations$location_name))
  for(loc in locations$location_name){
    # i = 7; loc = 'California'
    ## Subset data to correct time period and location
    recent_flu_cleaned |> 
      ungroup() |> 
      filter(location_name == loc, 
             resp_season == curr_resp_season, 
             date<forecast_weeks[i]) |> 
      mutate(value = value + 1) |> ## Makes sure no zeroes, need to correct later
      mutate(curr_weekly_change = log(lead(value)/value)) |> 
      # filter(!is.na(curr_weekly_change)) |> 
      # mutate(resp_season_week  = resp_season_week) |> 
      select(resp_season_week, value, curr_weekly_change) |> 
      get_traj_forecast(db = traj_db, 
                        recent_weeks_touse = 100, 
                        resp_week_range = 5,
                        forecast_horizon = 5) |> 
      mutate(forecast = forecast - 1) |> 
      mutate(forecast = ifelse(forecast < 0, 0, forecast)) -> forecast_trajectories
    
    
    forecast_trajectories |> 
      group_by(resp_season_week) |> 
      summarize(qs = list(value = quantile(forecast, probs = quantiles_needed))) |>
      mutate(forecast_horizon = seq_along(resp_season_week)) |> 
      unnest_wider(qs) |>
      gather(quantile, value, -resp_season_week, -forecast_horizon) |> 
      ungroup() |>
      mutate(quantile = as.numeric(gsub("[\\%,]", "", quantile))/100) |> 
      mutate(location = locations$location[match(loc, locations$location_name)],
             target = paste0(forecast_horizon, " wk inc flu hosp"),
             forecast_date = forecast_weeks[i],
             target_end_date = forecast_weeks[i] + 5 + 7*(forecast_horizon-1),
             type = 'quantile',
             value = round(value)) %>%
      arrange(location, target, quantile) |>
      dplyr::select(forecast_date, target, target_end_date, location, type, quantile, value) -> cleaned_forecasts_quantiles
    
    loc_forecasts[[match(loc, locations$location_name)]] <- cleaned_forecasts_quantiles
  }
  loc_forecasts |> 
    bind_rows() |> 
    write_csv(file = paste0('processed-data/2016forecasts-spline/', forecast_weeks[i], '-UGA', 'trajrate.csv'))
}



# recent_flu_cleaned |>
#   ungroup() |>
#   filter(location_name == loc,
#          resp_season == curr_resp_season,
#          date<forecast_weeks[i+5]) |>
#   mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
#   mutate(curr_weekly_change = lead(value)/value) |>
#   ggplot(aes(resp_season_week, value)) +
#   geom_point() +
#   geom_line(data = forecast_trajectories,
#             aes(resp_season_week, forecast, group = id), alpha = .1) +
#   coord_cartesian(ylim = c(0,10000))


# Score and compare forecasts ---------------------------------------------


library(covidHubUtils)
get_forecast_score <- function(file_path, flu_truth){
  forecasts <- read_csv(file_path, show_col_types = F)
  score_forecasts(forecasts |> 
                    mutate(model = 'UGAtrajrate',
                           target_variable = 'inc flu hosp',
                           horizon = str_replace(target, 
                                                 pattern = ' wk inc flu hosp', 
                                                 replacement = ''),
                           temporal_resolution = 'wk') |> 
                    select(-target), 
                  truth = flu_truth, 
                  use_median_as_point = T) |> 
    select(model, location, horizon, forecast_date, target_end_date, wis, quantile_coverage_0.95)
}

flu_truth <- load_truth(truth_source = 'HealthData',
                        target_variable = 'inc flu hosp', 
                        hub = 'FluSight')

list.files('processed-data/2016forecasts-spline', full.names = T) |> 
  map(get_forecast_score, flu_truth = flu_truth) |> 
  bind_rows() -> scored_forecasts_all

scored_forecasts_all |> 
  write_csv('processed-data/scored-2016forecasts-spline.csv')




# Analyze scores ----------------------------------------------------------
comp_file_path <- 'processed-data/2022-comparison-scores.rda'
if(!file.exists(comp_file_path) | refresh_data){
  list.files('processed-data/Flusight-baseline', full.names = T) |> 
    map(get_forecast_score, flu_truth = flu_truth) |> 
    bind_rows() |> 
    mutate(model = 'baseline') -> scored_baseline_all
  
  list.files('processed-data/Flusight-ensemble', full.names = T) |> 
    map(get_forecast_score, flu_truth = flu_truth) |> 
    bind_rows() |> 
    mutate(model = 'ensemble') -> scored_ensemble_all
  save(scored_baseline_all, scored_ensemble_all, file = comp_file_path) 
}else{
  load(comp_file_path)
}



# Compare scores ----------------------------------------------------------
scored_forecasts_all |> 
  bind_rows(scored_baseline_all,
            scored_ensemble_all) |> 
  filter(location != 'US',
         forecast_date %in% unique(scored_forecasts_all$forecast_date),
         horizon %in% 1:4 | horizon %in% paste0(1:4, ' wk ahead inc flu hosp')) -> combined_scores

combined_scores |> 
  group_by(model) |> 
  summarize(wis = mean(wis))

combined_scores |> 
  group_by(model,forecast_date) |> 
  summarize(wis = mean(wis)) |> 
  ggplot(aes(forecast_date, wis, color = model)) +
  geom_line()

# 
combined_scores |>
  arrange(desc(wis))

# 
forecasts <- read_csv('processed-data/2016forecasts-spline/2022-12-05-UGAtrajrate.csv')
forecasts |>
  filter(location == '06') |>
  filter(quantile %in% c('0.025', '0.5','0.975')) |>
  spread(quantile, value) |>
  ggplot(aes(target_end_date, `0.5`)) +
  geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), alpha = .2) +
  geom_line() +
  geom_point(data = flu_truth |>
               filter(location == '06',
                      target_end_date > '2022-10-01',
                      target_end_date < '2023-03-01'),
             aes(x = target_end_date, y = value), inherit.aes=F)


flu_truth |>
  filter(location == '02',
         target_end_date > '2022-10-01',
         target_end_date < '2023-03-01')

recent_flu_cleaned |> 
  ungroup() |> 
  filter(location_name == 'California', 
         resp_season == 2022, 
         date<'2023-03-01') |> 
  mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
  mutate(curr_weekly_change = lead(value)/value) |> 
  ggplot(aes(resp_season_week, curr_weekly_change)) + geom_line()


recent_flu_cleaned |> 
  ungroup() |> 
  filter(location_name == 'Wyoming', 
         resp_season == 2022, 
         date<'2023-03-01') |> 
  mutate(value = value+1) |> ## Makes sure no zeroes, need to correct later
  mutate(curr_weekly_change = lead(value)/value) |> 
  ggplot(aes(resp_season_week, value)) + geom_line()

traj_db  |> 
  ggplot(aes(resp_season_week, weekly_change, group = interaction(metric, location, resp_season))) +
  geom_line(alpha = .05) 
