## Growth rate trajectory matching


get_traj_forecast <- function(curr_data,
                              forecast_horizon = 5, ## How many weeks forecast and plotted?
                              recent_weeks_touse = 5, ## 100 means all data from season are used
                              nsamps = 1000,
                              resp_week_range = 0,
                              db = traj_db){
  # browser()
  
  most_recent_week <- max(curr_data$resp_season_week)
  most_recent_value <- tail(curr_data$value,1)
  
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
    
  db |> 
    inner_join(matching_data, 
               by = 'resp_season_week', relationship = "many-to-many") |> 
    group_by(week_change, metric, location, resp_season) |> 
    filter(n() == nrow(cleaned_data) | n() >= 4) |> ## Makes sure you've matched as many as the cleaned data or at least a full month
    # filter(metric == 'flusurv') |>  ## If you want to limit to specific database
    summarize(weight = sum((pred - curr_weekly_change)^2)/n(), .groups = "drop") |> 
    ungroup() |> 
    filter(!is.na(weight)) -> traj_temp
  
  
  # traj_temp |>
  #   arrange(weight) |>
  #   slice(1:1000) |>
  #   left_join(db |>
  #               nest(data = c('resp_season_week', 'pred','pred_se')),
  #             by = c('metric', 'location', 'resp_season')) |>
  #   unnest(data) |>
  #   mutate(resp_season_week = resp_season_week - week_change) |>
  #   ggplot(aes(resp_season_week, pred, group = interaction(metric, location, week_change, resp_season))) +
  #   geom_line(alpha = .05) +
  #   geom_point(data = cleaned_data, aes(resp_season_week, curr_weekly_change), color = 'red', inherit.aes=F)

  
  # traj_temp |> 
  #   arrange(weight) |> 
  #   slice(1:50) |> 
  #   summarize(min_weight = mean(weight)) |> pull(min_weight) -> min_allowed_weight
  # print(min_allowed_weight)
  min_allowed_weight <- 0.02
  
  traj_temp |>   
    mutate(weight = ifelse(weight < min_allowed_weight, min_allowed_weight, weight)) |>
    arrange(weight) |>
    slice(1:200) |> 
    sample_n(size = nsamps, replace = T, weight = 1/weight^2) |> 
    mutate(id = seq_along(metric)) |> 
    select(id, metric, location, resp_season, week_change) -> trajectories
  
  
  # trajectories |>
  #   left_join(db |>
  #               nest(data = c('resp_season_week', 'pred', 'pred_se')),
  #             by = c('metric', 'location', 'resp_season')) |>
  #   unnest(data) |>
  #   mutate(resp_season_week = resp_season_week - week_change) |>
  #   mutate(pred = rnorm(n(), pred, pred_se)) |>
  #   ggplot(aes(resp_season_week, pred, group = interaction(metric, location, week_change, resp_season))) +
  #   geom_point(alpha = .1) +
  #   geom_point(data = cleaned_data, aes(resp_season_week, curr_weekly_change), color = 'red', inherit.aes=F)

  
  
  trajectories |> 
    left_join(db |> 
                nest(data = c('resp_season_week', 'pred', 'pred_se')),
              by = c('metric', 'location', 'resp_season')) |> 
    unnest(data) |> 
    mutate(resp_season_week = resp_season_week - week_change) |> 
    filter(resp_season_week %in% most_recent_week:(most_recent_week+forecast_horizon-1)) |>
    mutate(weekly_change = exp(rnorm(n(), pred, pred_se))) |>
    group_by(id) |> 
    arrange(resp_season_week) |>
    mutate(mult_factor = cumprod(weekly_change)) |> 
    ungroup() |> 
    mutate(forecast = rpois(n(),most_recent_value*mult_factor)) |> 
    mutate(resp_season_week = resp_season_week + 1) |> 
    select(id, resp_season_week,forecast) 
  
}


# 
# 
# ## Plot out the growth rate trajectories alongside the data
# traj_plot <- trajectories |> 
#   left_join(db |> 
#               nest(data = c('resp_season_week', 'weekly_change'))) |> 
#   unnest(data) |> 
#   ggplot(aes(resp_season_week, weekly_change, group = id)) + 
#   geom_line(alpha = .1) +
#   geom_point(data = curr_data |> tail(recent_weeks_touse),
#              aes(resp_season_week, curr_weekly_change), color = 'red', inherit.aes=F)
# if(!is.null(full_data)){
#   forecast_df |> 
#     ggplot(aes(resp_season_week, q50)) +
#     geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .1)   +
#     geom_line() +
#     geom_point(data = full_data,
#                aes(resp_season_week, value), color = 'red', inherit.aes=F) -> forecast_plot
# } else {
#   forecast_df |> 
#     ggplot(aes(resp_season_week, q50)) +
#     geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .1)   +
#     geom_line() +
#     geom_point(data = curr_data,
#                aes(resp_season_week, value), color = 'red', inherit.aes=F)  -> forecast_plot
# }