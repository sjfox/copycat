library(tidyverse)
library(epidatr)
library(MMWRweek)
# set_api_key('00f36c125549')

# HHS data ------------------------------------------------------------
# flusurv_locs <- read_csv('https://raw.githubusercontent.com/cmu-delphi/delphi-epidata/main/labels/flusurv_locations.txt',col_names = F) |> 
#   pull(X1)
# pub_flusurv(locations = flusurv_locs, 
#             epiweeks = epirange(199001, 202320)) -> flusurv
# flusurv |> 
#   mutate(metric = 'flusurv') |> 
#   select(metric, location, epiweek,  value = rate_overall) -> flusurv_clean

# Flusurv data ------------------------------------------------------------
flusurv_locs <- read_csv('https://raw.githubusercontent.com/cmu-delphi/delphi-epidata/main/labels/flusurv_locations.txt',col_names = F) |> 
  pull(X1)
pub_flusurv(locations = flusurv_locs, 
            epiweeks = epirange(199001, 202320)) -> flusurv
flusurv |> 
  mutate(metric = 'flusurv') |> 
  select(metric, location, epiweek,  value = rate_overall) -> flusurv_clean


# ILI data ----------------------------------------------------------------
ili_locs <- c(tolower(state.abb), 'pr', 'dc', paste0('hhs', 1:10), paste0('cen', 1:9), 'nat')
pub_fluview(regions = ili_locs,
            epiweeks = epirange(199001, 202320)) -> ili
## Clinical data only after 2015, so stick to wili for now so we have as many time periods as we can
pub_fluview_clinical(regions = ili_locs,
                     epiweeks = epirange(199001, 202320)) -> ili_clinical_epidat


### Get all of the percent positive data from before 2016
census_join_tbl <- c("New England", "cen1",
                     "Mid-Atlantic", "cen2",
                     "East North Central", "cen3",
                     "West North Central", "cen4",
                     "South Atlantic", "cen5",
                     "East South Central", "cen6",
                     "West South Central", "cen7",
                     "Mountain", "cen8",
                     "Pacific", "cen9") |> matrix(ncol = 2, byrow = T) |> as_tibble()
colnames(census_join_tbl) <- c('REGION', 'region')

## Pull state level locations
load('raw-data/locations-data.rda')

read_csv('raw-data/nrevss-pre2016/hhs.csv', skip = 1) |> 
  mutate(region = paste0('hhs', str_replace(REGION, 'Region ', replacement = '')),
         epiweek = MMWRweek2Date(YEAR, WEEK),
         percent_positive = `PERCENT POSITIVE`) |> 
  select(region, epiweek, percent_positive) |> 
  bind_rows(read_csv('raw-data/nrevss-pre2016/census.csv', skip = 1) |> 
              left_join(census_join_tbl, by = 'REGION') |> 
              mutate(epiweek = MMWRweek2Date(YEAR, WEEK),
                     percent_positive = `PERCENT POSITIVE`) |> 
              select(region, epiweek, percent_positive),
            read_csv('raw-data/nrevss-pre2016/state.csv', skip = 1) |> 
              inner_join(locations |> 
                           select(region = abbreviation,
                                  REGION = location_name), by = 'REGION') |> 
              mutate(region = tolower(region),
                     epiweek = MMWRweek2Date(YEAR, WEEK),
                     percent_positive = as.numeric(`PERCENT POSITIVE`)) |> 
              select(region, epiweek, percent_positive),
            read_csv('raw-data/nrevss-pre2016/national.csv', skip = 1) |> 
              mutate(region = 'nat',
                     epiweek = MMWRweek2Date(YEAR, WEEK),
                     percent_positive = as.numeric(`PERCENT POSITIVE`)) |> 
              select(region, epiweek, percent_positive)) |> 
  mutate(percent_positive = ifelse(is.na(percent_positive), 0, percent_positive)) -> ili_clinical_pre2016





# bind_rows(ili_clinical_epidat,
#           ili_clinical_pre2016) |> 
#   select(region, epiweek, percent_positive)


ili |> 
  inner_join(bind_rows(ili_clinical_epidat,
                       ili_clinical_pre2016) |> 
               select(region, epiweek, percent_positive),
             by = c('region', 'epiweek')) |> 
  mutate(metric = 'fluview') |> 
  mutate(ilip = wili*percent_positive/100) |> 
  select(metric, location = region, epiweek, value = ilip)  -> ili_clean 


# NHSN data set -----------------------------------------------------------
nhsn <- RSocrata::read.socrata(url = "https://data.cdc.gov/resource/mpgq-jmmr.json") %>% 
  dplyr::filter(weekendingdate >= as.Date("2022-02-01"))

#remove  VI and AS as they are not included for FluSight, keep only necessary vars and add epiweek and epiyear 
nhsn <- nhsn %>% 
  dplyr::filter(!jurisdiction %in% c("VI", "AS", "GU", "MP")) %>% 
  dplyr::select(jurisdiction, weekendingdate, totalconfflunewadm) %>% 
  dplyr::rename("value" = "totalconfflunewadm", "date"="weekendingdate", "state"="jurisdiction") %>% 
  dplyr::mutate(epiweek = as.Date(date), 
                value = as.numeric(value),
                location = str_replace(state, "USA", "US")) |> 
  mutate(metric = 'nhsn') |> 
  select(metric, location, epiweek, value)


# Combine datasets --------------------------------------------------------
flusurv_clean |> 
  bind_rows(ili_clean,
            nhsn) |> 
  mutate(year = MMWRweek(epiweek)$MMWRyear,
         week = MMWRweek(epiweek)$MMWRweek) |> 
  mutate(resp_season = ifelse(week>=40, year, year-1)) |> 
  filter(!resp_season %in% c(2020)) |> 
  group_by(metric, location, resp_season) |> 
  arrange(epiweek) |> 
  mutate(resp_season_week = seq_along(week)) |> 
  ungroup() -> all_flu

get_seasonal_spline_vals <- function(season_weeks, value){
  # add_before <- 
  # season_weeks
  # browser()
  padding <- 4
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
  # plot(exp(mod$fitted.values))
  # as.numeric(predict(mod, se = T)$se.fit)
  # sqrt(mod$var)
  # sqrt(mod$nl.df)
  # plot(mod, se = TRUE, col = "blue")
  # plot(predict(mod), type = 'l', col = 'blue')
  # lines(predict(mod) + 2*as.numeric(predict(mod,se=T)$se.fit))
  # lines(predict(mod) - 2*as.numeric(predict(mod,se=T)$se.fit))
  # 
  # plx <- predict(mod, x=season_weeks, se=T)
  # 
  # rerun(500,tibble(val = exp(rnorm(length(season_weeks), mean = plx$fit, sd = plx$se*sqrt(plx$df)))) |> 
  #         mutate(step = seq_along(val))) |>  
  #   bind_rows() |> 
  #   ggplot(aes(step, val)) + geom_point(alpha = .05) + geom_smooth() +
  #   geom_hline(yintercept = 1, color = 'red')
  # 
  # plot(season_weeks, weekly_change)
  # lines(season_weeks,exp(plx$fit))
  # lines(season_weeks, exp(plx$fit - qt(0.975,plx$df)*plx$se), lty=2)
  # lines(season_weeks, exp(plx$fit + qt(0.975,plx$df)*plx$se), lty=2)
  
  # indices_for_spline <- !is.na(weekly_change) & !is.infinite(weekly_change)
  # if(sum(indices_for_spline) < 6){
  #   return(rep(NA,length(season_weeks)))
  # } else {
  #   # smooth.spline(season_weeks[indices_for_spline], 
  #   #               log(weekly_change[indices_for_spline]), 
  #   #               w = value[indices_for_spline], nknots = 5) -> sp
  # 
  #   preds = exp(predict(sp, x=season_weeks)$y)-1
  #   ifelse(preds<0,0,preds)
  # }
}
library(gam)
all_flu |> 
  # filter(metric == 'fluview', location == 'wy', resp_season == 2017) |> 
  group_by(metric, location, resp_season) |> 
  filter(n() >= 30) |> 
  arrange(epiweek) |> 
  mutate(get_seasonal_spline_vals(resp_season_week, value)) |> 
  # filter(!any(weekly_change>10)) |> 
  ungroup() |> 
  select(metric, location, resp_season, resp_season_week, pred, pred_se) -> traj_db


traj_db |> 
  filter(metric == 'nhsn')
  ggplot(aes(resp_season_week, exp(pred), group = interaction(metric, location, resp_season))) +
  geom_line(alpha = .1) 


save(traj_db, all_flu, file = 'processed-data/clean-historic-flu-spline.rda')

# load('processed-data/clean-historic-flu-spline.rda')

