source("nowcast.R")

rts <- c(rep(2, 20), (2 - 1:15 * 0.1), rep(0.5, 10), (0.5 + 1:7 * 0.1), rep(1.2, 10))
generation_interval <- rowMeans(EpiNow::covid_generation_times)

incubation_period <- list(mean = EpiNow::covid_incubation_period[1, ]$mean,
                          mean_sd = EpiNow::covid_incubation_period[1, ]$mean_sd,
                          sd = EpiNow::covid_incubation_period[1, ]$sd,
                          sd_sd = EpiNow::covid_incubation_period[1, ]$sd_sd,
                          max = 30)

reporting_delay <- list(mean = log(5),
                        mean_sd = log(1.1),
                        sd = log(2),
                        sd_sd = log(1.2),
                        max = 30)

## Sample a report delay as a lognormal - take 10 samples
delay_defs <- EpiNow::lognorm_dist_def(mean = reporting_delay$mean, mean_sd = reporting_delay$mean_sd,
                                       sd = reporting_delay$sd, sd_sd = reporting_delay$sd_sd,
                                       max_value = reporting_delay$max, samples = 1000)


## Sample a incubation period (again using the default for covid) - take 10 samples
incubation_defs <- EpiNow::lognorm_dist_def(mean = incubation_period$mean,
                                            mean_sd = incubation_period$mean_sd,
                                            sd = incubation_period$sd,
                                            sd_sd = incubation_period$sd_sd,
                                            max_value = 30, samples = 1000)

simulated_cases <- EpiNow::simulate_cases(rts, initial_cases = 100 , initial_date = as.Date("2020-03-01"),
                                          generation_interval = generation_interval, delay_def = delay_defs[1, ],
                                          incubation_def = incubation_defs[1, ],
                                          reporting_effect = c(1.4, rep(1, 4), 0.8, 0.8))

simulated_reports <- simulated_cases[reference == "report"][, confirm := cases][, cases := NULL]

# Median covid generation interval - used to smooth prior cases by date of infection
generation_interval <- rowMeans(EpiNow::covid_generation_times)
generation_interval <- sum(!(cumsum(generation_interval) > 0.5)) + 1   

## Reconstruction via backwards sampling
sampling_cases <- nowcast_pipeline(reported_cases = simulated_reports[, import_status := "local"], 
                                   target_date = max(simulated_reports$date),
                                   delay_defs = delay_defs,
                                   incubation_defs = incubation_defs,
                                   nowcast_lag = 0, approx_delay = TRUE)

## Non-parameteric reconstruction
model <- rstan::stan_model("nowcast.stan")
non_parametric_cases <- nowcast(simulated_reports,
                                family = "poisson", incubation_period = incubation_period,
                                reporting_delay = reporting_delay,
                                generation_interval = generation_interval, cores = 4, chains = 4,
                                samples = 1000, return_all = TRUE, model = model, verbose = TRUE)




simulated_cases2 <- simulated_cases[reference %in% c("infection", "report")][, median := cases][,
                                                                                               type := ifelse(reference == "infection", 
                                                                                                              "Simulated infections", 
                                                                                                              "Simulated reported cases")][,
                                                                                                                                           `:=`(cases = NULL, reference = NULL)]

summarise_nowcasting_approaches <- function(sampling_cases, non_parametric_cases) {
  sampling_cases <- sampling_cases[type %in% "infection_upscaled"][,
                                                                   .(median = median(cases), bottom = quantile(cases, 0.025),
                                                                     lower = quantile(cases, 0.25), upper = quantile(cases, 0.75),
                                                                     top = quantile(cases, 0.975)), by = c("date", "type")][,
                                                                                                                            type := "Sampled"]
  
  non_parametric_infections <- non_parametric_cases$infections[,
                                                               .(median = median(value), bottom = quantile(value, 0.025),
                                                                 lower = quantile(value, 0.25), upper = quantile(value, 0.75),
                                                                 top = quantile(value, 0.975)), by = c("date")][,
                                                                                                                type := "Non-parametric"]
  
  out <- data.table::rbindlist(list(sampling_cases, non_parametric_infections), fill = TRUE)
  return(out)
}

summarised_nowcasting_approaches <- summarise_nowcasting_approaches(sampling_cases, non_parametric_cases)

simulated_cases2 <- data.table::rbindlist(list(simulated_cases2,
                                              summarised_nowcasting_approaches),
                                         fill = TRUE)

plot_data <- simulated_cases2[date >= as.Date("2020-03-01")][,
                                                            type := factor(type, levels = c("Simulated infections",
                                                                                            "Simulated reported cases",
                                                                                            "Non-parametric",
                                                                                            "Sampled"))]
plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = date, col = type, fill = type)) +
  ggplot2::geom_col(data = plot_data[type %in% "Simulated infections"],
                    ggplot2::aes(y = median), fill = "grey", col = "white", alpha = 0.8) +
  ggplot2::geom_line(data = plot_data[type %in% "Simulated reported cases"],
                     ggplot2::aes(y = median), size = 1.1) +
  ggplot2::geom_linerange(data = plot_data[!type %in% "Simulated infections"],
                          ggplot2::aes(ymin = bottom, ymax = top), 
                          alpha = 0.4, size = 1.5) +
  ggplot2::geom_linerange(data = plot_data[!type %in% "Simulated infections"],
                          ggplot2::aes(ymin = lower, ymax = upper), 
                          alpha = 0.6, size = 1.5) +
  cowplot::theme_cowplot() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::labs(y = "Cases", x = "Date", col = "Type")

plot
