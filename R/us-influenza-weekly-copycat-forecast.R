
library(tidyverse)
library(lubridate)
library(MMWRweek)
library(cowplot)
theme_set(theme_cowplot())

source('R/get-traj-forecast.R')
source('R/get-traj-forecast-percentage.R')
load('raw-data/locations-data.rda')
load('processed-data/clean-historic-flu-spline.rda')


# Grab NHSN influenza data ----------------------------------------------
flu_data <- read_csv('https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/target-data/target-hospital-admissions.csv')

curr_resp_season <- 2025 ## Set the current year
first_week_of_season <- '2025-10-04' ## This is used to start the season of the recent data
quantiles_needed <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

flu_data |> 
  mutate(year = MMWRweek(date)$MMWRyear,
         week = MMWRweek(date)$MMWRweek) |> 
  mutate(resp_season = ifelse(week>=40, year, year-1))  |> 
  filter(date >= first_week_of_season) |> 
  group_by(location, location_name, resp_season) |>
  arrange(date) |>
  mutate(resp_season_week = seq_along(week)) -> recent_flu_cleaned


## Get the NHSN Reporting percentage
nhsn <- RSocrata::read.socrata(url = "https://data.cdc.gov/resource/mpgq-jmmr.json") %>% 
  dplyr::filter(weekendingdate >= as.Date("2022-02-01"))
nhsn %>% 
  dplyr::filter(!jurisdiction %in% c("VI", "AS", "GU", "MP")) %>% as_tibble() |> 
  dplyr::select(jurisdiction, weekendingdate, totalconffluhosppatsperc) |> 
  dplyr::rename("pct_reporting" = "totalconffluhosppatsperc", "date"="weekendingdate", "state"="jurisdiction") %>% 
  dplyr::mutate(epiweek = as.Date(date), 
                pct_reporting = as.numeric(pct_reporting),
                location = str_replace(state, "USA", "US")) |>
  select(epiweek,location, pct_reporting) |> 
  group_by(location) |> 
  arrange(epiweek) |> 
  summarize(ratio = max(tail(pct_reporting,15))/tail(pct_reporting,1)) |> 
  mutate(scaling_factor = (1+ratio)/2) -> reporting_scale_factor

## Correct the recent NHSN data point with the scale factor for all locations
recent_flu_cleaned |> 
  left_join(reporting_scale_factor |> 
              left_join(locations |> 
                          select(abbreviation, location_name), by = c('location' = 'abbreviation')) |> 
              select(-location), by = 'location_name') |> 
  mutate(value = ifelse(date == max(date), round(value*scaling_factor), value)) -> recent_flu_cleaned

## Get NSSP data
nssp_data <- read_csv('https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/refs/heads/main/target-data/target-ed-visits-prop.csv')
nssp_data |> 
  mutate(year = MMWRweek(date)$MMWRyear,
         week = MMWRweek(date)$MMWRweek) |> 
  mutate(resp_season = ifelse(week>=40, year, year-1))  |> 
  filter(date >= first_week_of_season) |> 
  group_by(location, location_name, resp_season) |>
  arrange(date) |>
  mutate(resp_season_week = seq_along(week)) -> nssp_recent_cleaned

## Switch to the date of the most recent HHS weekly data point
forecast_date <- recent_flu_cleaned |> pull(date) |> max()
forecast_date ## Should be most recent date of data - can be useful to check


# Run and save forecasts --------------------------------------------------
weeks_to_use <- 100
week_range <- 2
fcast_horizon <- 25


## Get most recent season's peak
prev_season <- flu_data |> 
  mutate(year = MMWRweek(date)$MMWRyear,
         week = MMWRweek(date)$MMWRweek) |> 
  mutate(resp_season = ifelse(week>=40, year, year-1))  |> 
  filter(resp_season == curr_resp_season-1) |> 
  group_by(location, location_name, resp_season) |>
  arrange(date) |>
  mutate(resp_season_week = seq_along(week)) |> 
  ungroup()



loc_forecasts <- vector('list', length = length(locations$location_name))
for(loc in locations$location_name){
  # i = 6; loc = 'Illinois'
  
  ## Get previous season's peak
  prev_season_peak <- prev_season |> 
    ungroup() |>
    filter(location_name == loc) |> 
    pull(value) |> max(na.rm=T)
  
  ## Subset data to correct time period and location
  recent_flu_cleaned |>
    # mutate(value = ifelse(date==max(date), round(value*curr_loc_scale_factor), value)) |> 
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
                      top_matches = 100,
                      error_exponentiation = 2,
                      nbinom_disp = 100) |>
    mutate(forecast = forecast-1) |>
    mutate(forecast = ifelse(forecast < 0, 0, forecast)) -> forecast_trajectories
  
  
  ## Filter out trajectories that go nuts.
  forecast_trajectories |> 
    filter(id %in% (forecast_trajectories |> 
                      group_by(id) |> 
                      summarize(peak_val = max(forecast)) |> 
                      filter(peak_val < 3*prev_season_peak) |> pull(id))) -> forecast_trajectories

  ## Get the observed and forecast trajectory combination for peak week and magnitude predictrions
  expand_grid(id = unique(forecast_trajectories$id),
              recent_flu_cleaned |>
                ungroup() |>
                filter(location_name == loc,
                       resp_season == curr_resp_season,
                       date <= forecast_date) |> 
                select(resp_season_week, value)) |> 
    bind_rows(forecast_trajectories |> 
                rename(value = forecast)) -> full_trajectories
  
  
  
  
  ## Get the standard quantile time-series predictions
  forecast_trajectories |>
    group_by(resp_season_week) |>
    summarize(qs = list(value = quantile(forecast, probs = quantiles_needed))) |>
    mutate(horizon = seq_along(resp_season_week)-1) |>
    filter(horizon<=3) |> 
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
  
  ## Get the rate change forecasts
  forecast_trajectories |>
    group_by(id) |> 
    arrange(resp_season_week) |> 
    mutate(horizon = seq_along(resp_season_week)-1) |> 
    ungroup() |> 
    filter(horizon<=3) |> 
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
  
  
  ## Get the peak week summary (week and magnitude)
  full_trajectories |> 
    group_by(id) %>%
    summarise(
      peak_week = resp_season_week[which.max(value)],
      peak_magnitude = max(value),
      .groups = "drop"
    ) -> peak_trajectory_summary
  
  ## Get the peak incidence quantiles
  peak_trajectory_summary |> 
    summarize(qs = list(value = quantile(peak_magnitude, probs = quantiles_needed))) |>
    unnest_wider(qs) |>
    gather(quantile, value) |>
    ungroup() |>
    mutate(quantile = as.numeric(gsub("[\\%,]", "", quantile))/100) |>
    mutate(location = locations$location[match(loc, locations$location_name)],
           target = paste0("peak inc flu hosp"),
           reference_date = forecast_date + 7,
           target_end_date = NA,
           horizon = NA,
           output_type_id = quantile,
           output_type = 'quantile',
           value = round(value)) %>%
    arrange(location, quantile) |> 
    dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value) -> cleaned_peak_incidence_quantiles
  
  
  tibble(output_type_id = seq.Date(from = ymd('2025-11-22'), to = ymd('2026-05-23'), by = 7)) |> 
    left_join(peak_trajectory_summary |> 
                count(peak_week) |> 
                mutate(value = n/sum(n),
                       output_type_id = ymd(first_week_of_season) + days((peak_week-1)*7)) |> 
                select(output_type_id, value)) |> 
    mutate(value = value/sum(value,na.rm=T)) |> 
    mutate(value = ifelse(is.na(value), 0 , value)) |> 
    mutate(location = locations$location[match(loc, locations$location_name)],
           target = paste0("peak week inc flu hosp"),
           reference_date = forecast_date + 7,
           target_end_date = NA,
           horizon = NA,
           output_type = 'pmf',
           value = value) %>%
    dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value) -> cleaned_peak_date_quantiles
  
  
  loc_forecasts[[match(loc, locations$location_name)]] <- cleaned_forecasts_quantiles |>
    mutate(output_type_id = as.character(output_type_id)) |>
    bind_rows(cleaned_forecasts_ratechange,
              cleaned_peak_incidence_quantiles |> 
                mutate(output_type_id = as.character(output_type_id)),
              cleaned_peak_date_quantiles |> 
                mutate(output_type_id = as.character(output_type_id))) 
}
loc_forecasts |>
  bind_rows() -> nhsn_forecasts



# NSSP Predictions --------------------------------------------------------
loc_forecasts_nssp <- vector('list', length = length(unique(nssp_recent_cleaned$location_name)))
for(loc_nssp in unique(nssp_recent_cleaned$location_name)){
  # loc_nssp = 'Illinois'
  # print(loc_nssp)
  
  ## Subset data to correct time period and location
  nssp_recent_cleaned |>
    ungroup() |>
    filter(location_name == loc_nssp,
           resp_season == curr_resp_season,
           date <= forecast_date) |>
    mutate(value = value + 1) |> ## Makes sure no zeroes, need to correct later - maybe change to smaller number
    mutate(curr_weekly_change = log(lead(value)/value)) |> 
    select(resp_season_week, value, curr_weekly_change) |>
    get_traj_forecast_percentage(db = traj_db |> 
                                   filter(metric == 'nssp'),
                 recent_weeks_touse = weeks_to_use,
                 resp_week_range = week_range,
                 forecast_horizon = fcast_horizon,
                 nsamps = 1000,
                 min_allowed_weight = 0.000000001,
                 # min_allowed_weight = 0.000000001,
                 top_matches = 250,
                 error_exponentiation = 2,
                 norm_sd = .0003) |> 
    mutate(forecast = forecast-1) |>
    mutate(forecast = ifelse(forecast < 0, 0, forecast)) |> 
    mutate(forecast = ifelse(forecast > 1, 1, forecast)) -> nssp_forecast_trajectories
  
  
  ## Get the standard quantile time-series predictions
  nssp_forecast_trajectories |>
    group_by(resp_season_week) |>
    summarize(qs = list(value = quantile(forecast, probs = quantiles_needed))) |>
    mutate(horizon = seq_along(resp_season_week)-1) |>
    filter(horizon<=3) |> 
    unnest_wider(qs) |>
    gather(quantile, value, -resp_season_week, -horizon) |>
    ungroup() |>
    mutate(quantile = as.numeric(gsub("[\\%,]", "", quantile))/100) |>
    mutate(location = locations$location[match(loc_nssp, locations$location_name)],
           target = paste0("wk inc flu prop ed visits"),
           reference_date = forecast_date + 7,
           target_end_date = forecast_date + 7 + horizon*7,
           output_type_id = quantile,
           output_type = 'quantile',
           value = value) %>%
    arrange(location, horizon, quantile) |>
    dplyr::select(reference_date, target, horizon, target_end_date, location, output_type, output_type_id, value) -> cleaned_forecasts_quantiles_nssp
  
  
  loc_forecasts_nssp[[match(loc_nssp, locations$location_name)]] <- cleaned_forecasts_quantiles_nssp |>
    mutate(output_type_id = as.character(output_type_id))
}


nhsn_forecasts |> 
  bind_rows(loc_forecasts_nssp) -> all_forecasts

all_forecasts |> 
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
    filter(!is.na(horizon),
           location_name == curr_location_name) -> forecast_df
  
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

prev_season_data <- prev_season |>
  filter(resp_season_week > rec_resp_season_week-10,
         resp_season_week < rec_resp_season_week+7)

nhsn_forecasts |> 
  filter(target == 'wk inc flu hosp') |> 
  filter(output_type == 'quantile') |> 
  filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |> 
  left_join(locations, by = 'location') |> 
  spread(output_type_id, value) -> forecast_df

## Plot out NHSN forecasts
locations$location_name |> 
  map(make_individual_plot, 
      curr_season_data = curr_season_data, 
      prev_season_data = prev_season_data,
      forecast_df = forecast_df) -> plots


plot_grid(plotlist = plots) |> 
  save_plot(filename = paste0('figures/us-rt/', forecast_date + 7, '_rt-forecast.png'), base_height = 14, base_asp = 1.8, bg = 'white')


## Plot out NSSP forecasts
curr_season_data_nssp <- nssp_recent_cleaned |>
  ungroup() |>
  filter(date<=forecast_date, resp_season == curr_resp_season)

prev_season_data_nssp <- nssp_data |> 
  mutate(year = MMWRweek(date)$MMWRyear,
         week = MMWRweek(date)$MMWRweek) |> 
  mutate(resp_season = ifelse(week>=40, year, year-1))  |> 
  filter(resp_season == curr_resp_season-1) |> 
  group_by(location, location_name, resp_season) |>
  arrange(date) |>
  mutate(resp_season_week = seq_along(week)) |> 
  ungroup() |>
  filter(resp_season_week > rec_resp_season_week-10,
         resp_season_week < rec_resp_season_week+7)

loc_forecasts_nssp |> 
  bind_rows() |> 
  filter(target == 'wk inc flu prop ed visits') |> 
  filter(output_type == 'quantile') |> 
  filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |> 
  left_join(locations, by = 'location') |> 
  spread(output_type_id, value) -> forecast_df_nssp

unique(forecast_df_nssp$location_name) |> 
  map(make_individual_plot, 
      curr_season_data = curr_season_data_nssp, 
      prev_season_data = prev_season_data_nssp,
      forecast_df = forecast_df_nssp) -> plots_nssp

plot_grid(plotlist = plots_nssp) |> 
  save_plot(filename = paste0('figures/us-rt/', forecast_date + 7, '_rt-forecast_nssp.png'), base_height = 14, base_asp = 1.8, bg = 'white')


# Plot peak metrics -------------------------------------------------------
loc_forecasts |>
  bind_rows() |> 
  filter(target == 'peak inc flu hosp') |> 
  filter(output_type_id %in% c(0.025, 0.25, 0.5, 0.75, 0.975)) |> 
  left_join(locations, by = 'location') |> 
  ggplot(aes(T, value, color = as.factor(output_type_id))) +
  geom_point() +
  geom_hline(data = prev_season |> 
               group_by(location_name) |> 
               summarize(peak_val = max(value,na.rm=T)),
             aes(yintercept = peak_val), color = 'red')+
  facet_wrap(~location_name, scales = 'free_y') -> peak_pred_plot

save_plot('figures/us-rt/peak_pred_plot.png',peak_pred_plot, base_height = 8, base_asp = 1.6, bg = 'white')


loc_forecasts |>
  bind_rows() |> 
  filter(target == 'peak week inc flu hosp') |> 
  left_join(locations, by = 'location') |> 
  ggplot(aes(ymd(output_type_id), value)) +
  geom_col() +
  geom_vline(data = prev_season |> 
               group_by(location_name) |>
               filter(value == max(value, na.rm=T)) |> 
               mutate(peak_date = date + days(365)),
             aes(xintercept = peak_date), color = 'red')+
  facet_wrap(~location_name) -> peak_week_plot

save_plot('figures/us-rt/peak_week_plot.png',peak_week_plot, base_height = 8, base_asp = 1.6, bg = 'white')


# nhsn <- RSocrata::read.socrata(url = "https://data.cdc.gov/resource/mpgq-jmmr.json") %>% 
#   dplyr::filter(weekendingdate >= as.Date("2022-02-01"))
# 
# #remove  VI and AS as they are not included for FluSight, keep only necessary vars and add epiweek and epiyear 
# nhsn <- nhsn %>% 
#   dplyr::filter(!jurisdiction %in% c("VI", "AS", "GU", "MP")) %>% as_tibble()
#   dplyr::select(jurisdiction, weekendingdate, totalconfflunewadm) %>% 
#   dplyr::rename("value" = "totalconfflunewadm", "date"="weekendingdate", "state"="jurisdiction") %>% 
#   dplyr::mutate(epiweek = as.Date(date), 
#                 value = as.numeric(value),
#                 location = str_replace(state, "USA", "US")) |> 
#   mutate(metric = 'nhsn') |> 
#   select(metric, location, epiweek, value)
# 
#   nhsn %>% 
#     dplyr::filter(!jurisdiction %in% c("VI", "AS", "GU", "MP")) %>% as_tibble() |> 
#     dplyr::select(jurisdiction, weekendingdate, totalconffluhosppatsperc) |> 
#     dplyr::rename("pct_reporting" = "totalconffluhosppatsperc", "date"="weekendingdate", "state"="jurisdiction") %>% 
#     dplyr::mutate(epiweek = as.Date(date), 
#                   pct_reporting = as.numeric(pct_reporting),
#                   location = str_replace(state, "USA", "US")) |>
#     select(epiweek,location, pct_reporting) |> 
#     group_by(location) |> 
#     arrange(epiweek) |> 
#     summarize(ratio = max(pct_reporting)/tail(pct_reporting,1)) |> 
#     mutate(scaling_factor = (1+ratio)/2) |> 
#     filter(scaling_factor>1.3)
#     ggplot(aes(location, scaling_factor)) + geom_col()
  
  
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
